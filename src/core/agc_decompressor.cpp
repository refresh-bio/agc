// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include "../core/agc_decompressor.h"
#include "../core/genome_io.h"
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
bool CAGCDecompressor::GetCollectionFiles(const string& _path, const uint32_t _line_length, const uint32_t no_threads)
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

	collection_desc->get_samples_list(v_samples);

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
				cur_path.append(prev_sample_name + ".fa");

				if (_path.empty())
					gio.Open("", true);
				else
					gio.Open(cur_path.string(), true);

				is_gio_opened = true;
			}

			++global_id;

			gio.SaveContig(ctg.contig_name, ctg.contig_data, _line_length);
		}

		if (is_gio_opened)
			gio.Close();

		});

	v_threads.clear();
	v_threads.reserve(no_threads);

	uint32_t global_i = 0;
	uint32_t prev_no_contigs = no_threads;

	q_contig_tasks->Restart(1);

	start_decompressing_threads(v_threads, no_threads, true);

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
bool CAGCDecompressor::GetSampleFile(const string& _file_name, const vector<string>& sample_names, const uint32_t _line_length, const uint32_t no_threads)
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

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		sample_contig_data_t ctg;

		gio.Open(_file_name, true);

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;

			gio.SaveContig(ctg.contig_name, ctg.contig_data, _line_length);
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

		sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

		q_contig_tasks->Restart(1);

		start_decompressing_threads(v_threads, no_threads, true);

		for (auto& task : v_tasks)
			q_contig_tasks->Push(task, 0);
		q_contig_tasks->MarkCompleted();

		join_threads(v_threads);

		v_threads.clear();
	}

	gio_thread.join();

	q_contig_tasks.release();
	pq_contigs_to_save.release();

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

	start_decompressing_threads(v_threads, no_threads, true);

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
bool CAGCDecompressor::GetContigFile(const string& _file_name, const vector<string>& contig_names, const uint32_t _line_length, const uint32_t no_threads)
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
			gio.SaveContig(ctg.contig_name, ctg.contig_data, _line_length);
		}

		gio.Close();
		});

	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	start_decompressing_threads(v_threads, no_threads, true);

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

// EOF
