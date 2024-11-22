// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "agc_decompressor.h"
#include "genome_io.h"
#include <filesystem>
#include <chrono>

using namespace std::filesystem;

// *******************************************************************************************
CAGCDecompressor::CAGCDecompressor(bool _is_app_mode) : CAGCDecompressorLibrary(_is_app_mode)
{
}

// *******************************************************************************************
CAGCDecompressor::~CAGCDecompressor()
{
}

// *******************************************************************************************
void CAGCDecompressor::gzip_contig(contig_t& ctg, contig_t& working_space, refresh::gz_in_memory& gzip_compressor)
{
	auto overhead = gzip_compressor.get_overhead(ctg.size());
	working_space.resize(ctg.size() + overhead);

	auto gzipped_size = gzip_compressor.compress(ctg.data(), ctg.size(), working_space.data(), working_space.size());
	working_space.resize(gzipped_size);

	swap(ctg, working_space);
}
 
// *******************************************************************************************
void CAGCDecompressor::start_decompressing_threads(vector<thread>& v_threads, const uint32_t n_t, uint32_t gzip_level, uint32_t line_len, bool fast)
{
	for (uint32_t i = 0; i < n_t; ++i)
		v_threads.emplace_back([&, i, gzip_level, line_len, fast] {

		auto zstd_ctx = ZSTD_createDCtx();

		contig_t ctg, working_space;
		contig_task_t contig_desc;
		refresh::gz_in_memory gzip_compressor(gzip_level);

		while (!q_contig_tasks->IsCompleted())
		{
			if (!q_contig_tasks->Pop(contig_desc))
				break;

			size_t priority = contig_desc.priority;

			if (!decompress_contig(contig_desc, zstd_ctx, ctg, fast))
				continue;

			if(line_len == 0)
				CNumAlphaConverter::convert_to_alpha(ctg);
			else
				CNumAlphaConverter::convert_and_split_into_lines(ctg, working_space, line_len);

			if (gzip_level)
				gzip_contig(ctg, working_space, gzip_compressor);

			name_range_t contig_name_range = contig_desc.name_range;

			pq_contigs_to_save->Emplace(priority, sample_contig_data_t{ contig_desc.sample_name, contig_name_range.str(), move(ctg) });
			ctg.clear();
		}

		pq_contigs_to_save->MarkCompleted();

		ZSTD_freeDCtx(zstd_ctx);
		});
}

// *******************************************************************************************
// Assign archive from already opened one
bool CAGCDecompressor::AssignArchive(const CAGCBasic& agc_basic)
{
	working_mode = working_mode_t::decompression;
	is_app_mode = false;

	kmer_length = agc_basic.kmer_length;
	pack_cardinality = agc_basic.pack_cardinality;
	min_match_len = agc_basic.min_match_len;

	in_archive_name = "";
	prefetch_archive = false;
	archive_version = agc_basic.archive_version;

	in_archive = agc_basic.in_archive;
	collection_desc = agc_basic.collection_desc;

	m_file_type_info = agc_basic.m_file_type_info;

	compression_params = agc_basic.compression_params;

	verbosity = agc_basic.verbosity;

	return true;
}

// *******************************************************************************************
bool CAGCDecompressor::GetCollectionFiles(const string& _path, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, bool no_ref, bool fast, uint32_t verbosity)
{
	if (working_mode != working_mode_t::decompression)
		return false;

	vector<string> v_samples;
	vector<sample_desc_t> v_sample_desc;

	path p(_path);

	if (!_path.empty() && !is_directory(p))
	{
		cerr << "Path must point to an existing directory\n";
		return false;
	}

	collection_desc->get_samples_list(v_samples, false);

	if (no_ref && !v_samples.empty())
		v_samples.erase(v_samples.begin());

	q_contig_tasks = make_unique<CBoundedQueue<contig_task_t>>(1, 1);
	pq_contigs_to_save = make_unique<CPriorityQueue<sample_contig_data_t>>(no_threads);

	vector<contig_task_t> v_tasks;
	vector<thread> v_threads;

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		sample_contig_data_t ctg;
		uint32_t global_id = 0;
		string prev_sample_name;
		bool is_gio_opened = false;
		size_t file_id = 0;

		string eol = "";

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;

			if (ctg.sample_name != prev_sample_name)
			{
				if (is_gio_opened)
					gio.Close();

				prev_sample_name = ctg.sample_name;
				path cur_path = _path;
				cur_path.append(prev_sample_name + ".fa" + (gzip_level ? ".gz" : ""));
//				cur_path.append(prev_sample_name + (gzip_level ? ".gz" : ""));

				if (_path.empty())
					gio.Open("", true);
				else
				{
					gio.Open(cur_path.string(), true);
					if (verbosity > 0)
					{
						cerr << eol << cur_path.string() << "  (" << ++file_id << " of " << v_samples.size() << ")";
						eol = "\n";
					}
				}

				is_gio_opened = true;
			}

			++global_id;

			gio.SaveContigDirectly(ctg.contig_name, ctg.contig_data, gzip_level);
		}

		if (is_gio_opened)
		{
			gio.Close();
			if (verbosity > 0)
				cerr << endl;
		}
		});

	v_threads.clear();
	v_threads.reserve(no_threads);

	uint32_t global_i = 0;
	uint32_t prev_no_contigs = no_threads;

	q_contig_tasks->Restart(1);

	start_decompressing_threads(v_threads, no_threads, gzip_level, _line_length, fast);

	sample_desc_t sample_desc;

	bool res = true;

	for (const auto& s : v_samples)
	{
		if (!collection_desc->get_sample_desc(s, sample_desc))
		{
			cerr << "There is no sample " << s << endl;

			res = false;
			break;
		}

		v_tasks.clear();
		v_tasks.shrink_to_fit();
		v_tasks.reserve(sample_desc.size());

		for (uint32_t i = 0; i < sample_desc.size(); ++i, ++global_i)
			v_tasks.emplace_back(global_i, s, sample_desc[i].first, sample_desc[i].second);

		sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

		while (pq_contigs_to_save->GetSize() > 3 * prev_no_contigs || q_contig_tasks->GetSize().first > 3 * no_threads)
			this_thread::sleep_for(std::chrono::microseconds(25));

		prev_no_contigs = sample_desc.size();

		for (auto& task : v_tasks)
			q_contig_tasks->Push(task, 0);
	}

	q_contig_tasks->MarkCompleted();

	join_threads(v_threads);

	gio_thread.join();

	q_contig_tasks.release();
	pq_contigs_to_save.release();

	return res;
}

// *******************************************************************************************
bool CAGCDecompressor::GetSampleFile(const string& _file_name, const vector<string>& sample_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity)
{
	if (working_mode != working_mode_t::decompression)
		return false;

	vector<sample_desc_t> v_sample_desc;

	for (const auto &s : sample_names)
	{
		sample_desc_t sample_desc;

		if (!collection_desc->get_sample_desc(s, sample_desc))
		{
			cerr << "There is no sample " << s << endl;

			return false;
		}

		v_sample_desc.emplace_back(sample_desc);
	}

	q_contig_tasks = make_unique<CBoundedQueue<contig_task_t>>(1, 1);
	pq_contigs_to_save = make_unique<CPriorityQueue<sample_contig_data_t>>(no_threads * v_sample_desc.size());

	vector<contig_task_t> v_tasks;
	vector<thread> v_threads;

	string eol = "";

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		sample_contig_data_t ctg;

		gio.Open(_file_name, true);

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;

			gio.SaveContigDirectly(ctg.contig_name, ctg.contig_data, gzip_level);

			if (!_file_name.empty() && verbosity > 0)
			{
				cerr << eol << ctg.contig_name;
				eol = "\n";
			}
		}

		gio.Close();
		});

	v_threads.clear();
	v_threads.reserve(no_threads);

	uint32_t global_i = 0;

	for (auto& sample_desc : v_sample_desc)
	{
		v_tasks.clear();
		v_tasks.shrink_to_fit();
		v_tasks.reserve(sample_desc.size());

		for (uint32_t i = 0; i < sample_desc.size(); ++i, ++global_i)
			v_tasks.emplace_back(global_i, "", sample_desc[i].first, sample_desc[i].second);

		std::sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

		q_contig_tasks->Restart(1);

		start_decompressing_threads(v_threads, no_threads, gzip_level, _line_length);

		for (auto& task : v_tasks)
			q_contig_tasks->Push(task, 0);
		q_contig_tasks->MarkCompleted();

		join_threads(v_threads);

		v_threads.clear();
	}

	gio_thread.join();

	if (!_file_name.empty() && verbosity > 0)
		cerr << eol;

	q_contig_tasks.release();
	pq_contigs_to_save.release();

	return true;
}

// *******************************************************************************************
bool CAGCDecompressor::GetSampleForStreaming(const string& _file_name, const vector<string>& sample_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity)
{
	if (working_mode != working_mode_t::decompression)
		return false;

	contig_t ctg, working_space;
	contig_task_t contig_desc;

	FILE* stream;

	if (_file_name.empty())
	{
		stream = stdout;
#ifdef _WIN32
		_setmode(_fileno(stream), _O_BINARY);
#endif
	}
	else
	{
		stream = fopen(_file_name.c_str(), "wb");
		if (!stream)
		{
//			if(app_mode)
			// !!! TODO: check app mode
			cerr << "Cannot open destination file: " << _file_name << endl;
			return false;
		}
		setvbuf(stream, nullptr, _IOFBF, 1 << 20);
	}

	auto zstd_ctx = ZSTD_createDCtx();

	CStreamWrapper stream_wrapper(stream, _line_length, gzip_level);

	for (const auto &s : sample_names)
	{
		sample_desc_t sample_desc;

		if (!collection_desc->get_sample_desc(s, sample_desc))
		{
			cerr << "There is no sample " << s << endl;

			return false;
		}

		for (auto contig_desc : sample_desc)
		{
			contig_task_t contig_task(0, "", contig_desc.first, contig_desc.second);

			stream_wrapper.start_contig(contig_desc.first);

			if (!decompress_contig_streaming(contig_task, zstd_ctx, stream_wrapper, false))
				continue;
		}
	}

	if (!_file_name.empty())
		fclose(stdout);

	ZSTD_freeDCtx(zstd_ctx);

	return true;
}

// *******************************************************************************************
bool CAGCDecompressor::GetSampleSequences(const string& sample_name, vector<pair<string, vector<uint8_t>>>& v_contig_seq, const uint32_t no_threads)
{
	if (working_mode != working_mode_t::decompression)
		return false;

	sample_desc_t sample_desc;

	if (!collection_desc->get_sample_desc(sample_name, sample_desc))
	{
		cerr << "There is no sample " << sample_name << endl;

		return false;
	}

	q_contig_tasks = make_unique<CBoundedQueue<contig_task_t>>(1, 1);
	pq_contigs_to_save = make_unique<CPriorityQueue<sample_contig_data_t>>(no_threads);

	vector<contig_task_t> v_tasks;
	vector<thread> v_threads;

	v_contig_seq.reserve(sample_desc.size());

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		sample_contig_data_t ctg;

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;

			v_contig_seq.emplace_back(ctg.contig_name, ctg.contig_data);
		}

		gio.Close();
		});

	v_threads.clear();
	v_threads.reserve(no_threads);

	uint32_t global_i = 0;

	v_tasks.clear();
	v_tasks.shrink_to_fit();
	v_tasks.reserve(sample_desc.size());

	for (uint32_t i = 0; i < sample_desc.size(); ++i, ++global_i)
		v_tasks.emplace_back(global_i, "", sample_desc[i].first, sample_desc[i].second);

	sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

	q_contig_tasks->Restart(1);

	start_decompressing_threads(v_threads, no_threads);

	for (auto& task : v_tasks)
		q_contig_tasks->Push(task, 0);
	q_contig_tasks->MarkCompleted();

	join_threads(v_threads);

	v_threads.clear();

	gio_thread.join();

	q_contig_tasks.release();
	pq_contigs_to_save.release();

	return true;
}

// *******************************************************************************************
bool CAGCDecompressor::GetContigFile(const string& _file_name, const vector<string>& contig_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity)
{
	if (working_mode != working_mode_t::decompression)
		return false;

	vector<pair<string, name_range_t>> v_sample_contig;
	name_range_t name_range;
	string sample;

	for (const auto &cn : contig_names)
	{
		if (!analyze_contig_query(cn, sample, name_range))
		{
			cerr << "Wrong contig format: " << cn << endl;
			continue;
		}

		if (!sample.empty())
		{
			if (!collection_desc->is_contig_desc(sample, name_range.name))
			{
				cerr << "There is no sample:contig pair: " << sample << " : " << name_range.name << endl;
				return false;
			}

			v_sample_contig.emplace_back(sample, name_range_t(name_range.name, name_range.from, name_range.to));
		}
		else
		{
			auto v_cand_samples = collection_desc->get_samples_for_contig(name_range.name);
			if (v_cand_samples.size() == 0)
			{
				cerr << "There is no contig: " << name_range.name << endl;
				return false;
			}
			if (v_cand_samples.size() > 1)
			{
				cerr << "There are " << v_cand_samples.size() << " samples with conting: " << name_range.name << endl;
				return false;
			}

			v_sample_contig.emplace_back(v_cand_samples.front(), name_range_t(name_range.name, name_range.from, name_range.to));
		}
	}

	q_contig_tasks = make_unique<CBoundedQueue<contig_task_t>>(1, 1);
	pq_contigs_to_save = make_unique<CPriorityQueue<sample_contig_data_t>>(no_threads);

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		sample_contig_data_t ctg;

		gio.Open(_file_name, true);

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;
			gio.SaveContigDirectly(ctg.contig_name, ctg.contig_data, gzip_level);
		}

		gio.Close();
		});

	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	start_decompressing_threads(v_threads, no_threads, gzip_level, _line_length);

	uint32_t id = 0;
	vector<segment_desc_t> contig_desc;

	for (auto& p_sc : v_sample_contig)
	{
		collection_desc->get_contig_desc(p_sc.first, p_sc.second.name, contig_desc);

		contig_task_t task(id++, "", p_sc.second, contig_desc);
		q_contig_tasks->Push(task, 0);
	}

	q_contig_tasks->MarkCompleted();

	join_threads(v_threads);
	gio_thread.join();

	q_contig_tasks.release();
	pq_contigs_to_save.release();

	return true;
}

// *******************************************************************************************
bool CAGCDecompressor::GetContigForStreaming(const string& _file_name, const vector<string>& contig_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity)
{
	if (working_mode != working_mode_t::decompression)
		return false;

	vector<pair<string, name_range_t>> v_sample_contig;
	name_range_t name_range;
	string sample;

	for (const auto &cn : contig_names)
	{
		if (!analyze_contig_query(cn, sample, name_range))
		{
			cerr << "Wrong contig format: " << cn << endl;
			continue;
		}

		if (!sample.empty())
		{
			if (!collection_desc->is_contig_desc(sample, name_range.name))
			{
				cerr << "There is no sample:contig pair: " << sample << " : " << name_range.name << endl;
				return false;
			}

			v_sample_contig.emplace_back(sample, name_range_t(name_range.name, name_range.from, name_range.to));
		}
		else
		{
			auto v_cand_samples = collection_desc->get_samples_for_contig(name_range.name);
			if (v_cand_samples.size() == 0)
			{
				cerr << "There is no contig: " << name_range.name << endl;
				return false;
			}
			if (v_cand_samples.size() > 1)
			{
				cerr << "There are " << v_cand_samples.size() << " samples with conting: " << name_range.name << endl;
				return false;
			}

			v_sample_contig.emplace_back(v_cand_samples.front(), name_range_t(name_range.name, name_range.from, name_range.to));
		}
	}

	FILE* stream;

	if (_file_name.empty())
	{
		stream = stdout;
#ifdef _WIN32
		_setmode(_fileno(stream), _O_BINARY);
#endif
	}
	else
	{
		stream = fopen(_file_name.c_str(), "wb");
		if (!stream)
		{
			//			if(app_mode)
						// !!! TODO: check app mode
			cerr << "Cannot open destination file: " << _file_name << endl;
			return false;
		}
		setvbuf(stream, nullptr, _IOFBF, 1 << 20);
	}

	auto zstd_ctx = ZSTD_createDCtx();

	CStreamWrapper stream_wrapper(stream, _line_length, gzip_level);

	for (auto& p_sc : v_sample_contig)
	{
		vector<segment_desc_t> contig_desc;
		collection_desc->get_contig_desc(p_sc.first, p_sc.second.name, contig_desc);

//		contig_task_t contig_task(0, p_sc.first, p_sc.second.name, contig_desc);
		contig_task_t contig_task(0, p_sc.first, p_sc.second, contig_desc);

//		stream_wrapper.start_contig(p_sc.second.name);
		stream_wrapper.start_contig(p_sc.second.str());

		if (!decompress_contig_streaming(contig_task, zstd_ctx, stream_wrapper, false))
			continue;
	}

	if (!_file_name.empty())
		fclose(stdout);

	ZSTD_freeDCtx(zstd_ctx);

	return true;
}

// EOF
