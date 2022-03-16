// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-02-24
// *******************************************************************************************

#include "../core/agc_decompressor.h"
#include "../core/genome_io.h"

// *******************************************************************************************
CAGCDecompressor::CAGCDecompressor(bool _is_app_mode) : CAGCDecompressorLibrary(_is_app_mode)
{
}

// *******************************************************************************************
CAGCDecompressor::~CAGCDecompressor()
{
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
		if (!collection_desc.get_sample_desc(s, sample_desc))
		{
			cerr << "There is no sample " << s << endl;

			return false;
		}

		v_sample_desc.emplace_back(sample_desc);
	}

	q_contig_tasks = make_unique<CBoundedQueue<tuple<size_t, name_range_t, vector<segment_desc_t>>>>(1, 1);
	pq_contigs_to_save = make_unique<CPriorityQueue<pair<string, contig_t>>>(no_threads * v_sample_desc.size());

	vector<tuple<size_t, name_range_t, vector<segment_desc_t>>> v_tasks;
	vector<thread> v_threads;

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		pair<string, contig_t> ctg;

		gio.Open(_file_name, true);

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;
			gio.SaveContigConverted(ctg.first, ctg.second, _line_length);
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
			v_tasks.emplace_back(global_i, sample_desc[i].first, sample_desc[i].second);

		sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return get<2>(x).size() > get<2>(y).size(); });

		q_contig_tasks->Restart(1);
		//		pq_contigs_to_save->Restart(no_threads);

		start_decompressing_threads(v_threads, no_threads);

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
			if (!collection_desc.is_contig_desc(sample, name_range.name))
			{
				cerr << "There is no sample:contig pair: " << sample << " : " << name_range.name << endl;
				return false;
			}

			v_sample_contig.emplace_back(sample, name_range_t(name_range.name, name_range.from, name_range.to));
		}
		else
		{
			auto v_cand_samples = collection_desc.get_samples_for_contig(name_range.name);
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

	q_contig_tasks = make_unique<CBoundedQueue<tuple<size_t, name_range_t, vector<segment_desc_t>>>>(1, 1);
	pq_contigs_to_save = make_unique<CPriorityQueue<pair<string, contig_t>>>(no_threads);

	// Saving thread
	thread gio_thread([&] {
		CGenomeIO gio;
		pair<string, contig_t> ctg;

		gio.Open(_file_name, true);

		while (!pq_contigs_to_save->IsCompleted())
		{
			if (!pq_contigs_to_save->Pop(ctg))
				break;
			gio.SaveContigConverted(ctg.first, ctg.second, _line_length);
		}

		gio.Close();
		});

	vector<thread> v_threads;
	v_threads.reserve(no_threads);

	start_decompressing_threads(v_threads, no_threads);

	uint32_t id = 0;
	vector<segment_desc_t> contig_desc;

	for (auto& p_sc : v_sample_contig)
	{
		collection_desc.get_contig_desc(p_sc.first, p_sc.second.name, contig_desc);

		tuple<size_t, name_range_t, vector<segment_desc_t>> task(id++, p_sc.second, contig_desc);
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
