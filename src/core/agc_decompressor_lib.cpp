// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-04-05
// *******************************************************************************************

#include "../core/agc_decompressor_lib.h"

// *******************************************************************************************
CAGCDecompressorLibrary::CAGCDecompressorLibrary(bool _is_app_mode) : CAGCBasic()
{
	is_app_mode = _is_app_mode;
}

// *******************************************************************************************
CAGCDecompressorLibrary::~CAGCDecompressorLibrary()
{
	if (working_mode == working_mode_t::decompression)
		close_decompression();
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::ListSamples(vector<string>& v_sample_names)
{
	return collection_desc.get_samples_list(v_sample_names);
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::ListContigs(const string& sample_name, vector<string>& v_contig_names)
{
	return collection_desc.get_contig_list_in_sample(sample_name, v_contig_names);
}

// *******************************************************************************************
void CAGCDecompressorLibrary::GetFileTypeInfo(map<string, string>& _m_file_type_info)
{
	_m_file_type_info = m_file_type_info;
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::IsOpened()
{
	return working_mode != working_mode_t::none;
}

// *******************************************************************************************
int32_t CAGCDecompressorLibrary::GetNoSamples()
{
	return static_cast<int32_t>(collection_desc.get_no_samples());
}

// *******************************************************************************************
int32_t CAGCDecompressorLibrary::GetNoContigs(const string& sample_name)
{
	return collection_desc.get_no_contigs(sample_name);
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::analyze_contig_query(const string& query, string& sample, name_range_t& name_range)
{
	smatch sm;

	name_range.from = -1;
	name_range.to = -1;

	sample.clear();

	if (regex_match(query, sm, re_csr))
	{
		name_range.name = sm[1].str();
		sample = sm[2].str();
		name_range.from = atoll(sm[3].str().c_str());
		name_range.to = atoll(sm[4].str().c_str());
	}
	else if (regex_match(query, sm, re_cs))
	{
		name_range.name = sm[1].str();
		sample = sm[2].str();
	}
	else if (regex_match(query, sm, re_cr))
	{
		name_range.name = sm[1].str();
		name_range.from = atoll(sm[2].str().c_str());
		name_range.to = atoll(sm[3].str().c_str());
	}
	else if (regex_match(query, sm, re_c))
	{
		name_range.name = sm[1].str();
	}
	else
	{
		return false;
	}

	return true;
}

// *******************************************************************************************
int CAGCDecompressorLibrary::GetContigString(const string& sample_name, const string& contig_name, const int start, const int end, string& contig_data)
{
	if (working_mode != working_mode_t::decompression)
		return -1;

	uint32_t id = 0;
	vector<segment_desc_t> contig_desc;
	string det_sample_name = sample_name;

	if (sample_name.empty())
	{
		auto v_cand_samples = collection_desc.get_samples_for_contig(contig_name);
		if (v_cand_samples.size() == 0)
			return -1;
		if (v_cand_samples.size() > 1)
			return -2;

		det_sample_name = v_cand_samples.front();
	}

	string full_contig_name = contig_name;

	if (!collection_desc.get_contig_desc(det_sample_name, full_contig_name, contig_desc))
		return -1;

	tuple<size_t, name_range_t, vector<segment_desc_t>> task(id++, name_range_t(full_contig_name, start, end), contig_desc);
	contig_t ctg;

	decompress_contig(task, nullptr, ctg);

	contig_data.reserve(ctg.size());
	for (auto& c : ctg)
		contig_data.push_back(cnv_num[static_cast<uint8_t>(c)]);

	return 0;
}

// *******************************************************************************************
int64_t CAGCDecompressorLibrary::GetContigLength(const string& sample_name, const string& contig_name)
{
	vector<segment_desc_t> contig_desc;
	string det_sample_name = sample_name;

	if (sample_name.empty())
	{
		auto v_cand_samples = collection_desc.get_samples_for_contig(contig_name);
		if (v_cand_samples.size() == 0)
			return -1;
		if (v_cand_samples.size() > 1)
			return -2;

		det_sample_name = v_cand_samples.front();
	}

	string full_contig_name = contig_name;

	if (!collection_desc.get_contig_desc(det_sample_name, full_contig_name, contig_desc))
		return -1;

	int64_t len = 0;
	for (auto& x : contig_desc)
		len += x.raw_length;

	return len - (contig_desc.size() - 1) * kmer_length;
}

// *******************************************************************************************
void CAGCDecompressorLibrary::start_decompressing_threads(vector<thread>& v_threads, const uint32_t n_t, bool converted_to_alpha)
{
	for (uint32_t i = 0; i < n_t; ++i)
		v_threads.emplace_back([&, converted_to_alpha] {

		auto zstd_ctx = ZSTD_createDCtx();

		contig_t ctg;
		tuple<size_t, name_range_t, vector<segment_desc_t>> contig_desc;

		while (!q_contig_tasks->IsCompleted())
		{
			if (!q_contig_tasks->Pop(contig_desc))
				break;

			size_t priority = get<0>(contig_desc);

			if (!decompress_contig(contig_desc, zstd_ctx, ctg))
				continue;

			if (converted_to_alpha)
				convert_to_alpha(ctg);

			name_range_t contig_name_range = get<1>(contig_desc);

			pq_contigs_to_save->Emplace(priority, make_pair(contig_name_range.str(), move(ctg)));
		}

		pq_contigs_to_save->MarkCompleted();

		ZSTD_freeDCtx(zstd_ctx);
		});
}

// *******************************************************************************************
void CAGCDecompressorLibrary::convert_to_alpha(contig_t& ctg)
{
	size_t size = ctg.size();
	size_t i = size % 8;

	switch (i)
	{
	case 7: ctg[6] = cnv_num[ctg[6]];
	case 6: ctg[5] = cnv_num[ctg[5]];
	case 5: ctg[4] = cnv_num[ctg[4]];
	case 4: ctg[3] = cnv_num[ctg[3]];
	case 3: ctg[2] = cnv_num[ctg[2]];
	case 2: ctg[1] = cnv_num[ctg[1]];
	case 1: ctg[0] = cnv_num[ctg[0]];
	}

	for (; i < size; i += 8)
	{
		ctg[i] = cnv_num[ctg[i]];
		ctg[i+1] = cnv_num[ctg[i+1]];
		ctg[i+2] = cnv_num[ctg[i+2]];
		ctg[i+3] = cnv_num[ctg[i+3]];
		ctg[i+4] = cnv_num[ctg[i+4]];
		ctg[i+5] = cnv_num[ctg[i+5]];
		ctg[i+6] = cnv_num[ctg[i+6]];
		ctg[i+7] = cnv_num[ctg[i+7]];
	}
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::decompress_contig(task_desc_t& contig_desc, ZSTD_DCtx* zstd_ctx, contig_t& ctg)
{
	name_range_t &contig_name_range = get<1>(contig_desc);
	vector<contig_t> v_segments_loc;

	bool need_free_zstd = false;

	if (!zstd_ctx)
	{
		zstd_ctx = ZSTD_createDCtx();
		need_free_zstd = true;
	}

	int64_t from = contig_name_range.from;
	int64_t to = contig_name_range.to;
	int64_t curr_pos = 0;

	if (from < 0 && to < 0)
	{
		from = 0;
		to = 0x7fffffffffffffffu;
	}
	else
	{
		if (from < 0)
		{
			if (is_app_mode)
				cerr << "Warning: Start of range (" + to_string(from) + ") is below 0, so changed to 0\n";
			from = 0;
			contig_name_range.from = 0;
		}
		if (to < 0)
		{
			if (is_app_mode)
				cerr << "Warning: End of range (" + to_string(to) + ") is below 0, so changed to max value\n";
			to = 0x7fffffffffffffffu;
			contig_name_range.to = 0x7fffffffffffffffu;
		}
		if (from > to)
		{
			if (is_app_mode)
				cerr << "Warning: End of range (" + to_string(to) + ") is prior to start of range (" + to_string(from) + ") so changed to whole contig\n";

			from = 0;
			to = 0x7fffffffffffffffu;
			contig_name_range.from = -1;
			contig_name_range.to = -1;
		}
	}

	for (auto seg : get<2>(contig_desc))
	{
		int32_t seg_len = seg.raw_length;

		if (curr_pos + seg_len < from)
		{
			from -= seg_len - kmer_length;
			to -= seg_len - kmer_length;
			continue;
		}
		else if (curr_pos > to)
			break;

		decompress_segment(seg.group_id, seg.in_group_id, ctg, zstd_ctx);
		if (seg.is_rev_comp)
			reverse_complement(ctg);

		v_segments_loc.emplace_back(ctg);

		curr_pos += seg_len - kmer_length;
	}

	if (!v_segments_loc.empty())
	{
		ctg = v_segments_loc.front();

		for (uint32_t j = 1; j < v_segments_loc.size(); ++j)
			if (v_segments_loc[j].size() < compression_params.kmer_length)
			{
				if (is_app_mode)
					cerr << "Corrupted archive!" << endl;
			}
			else
			{
				ctg.insert(ctg.end(), v_segments_loc[j].begin() + compression_params.kmer_length, v_segments_loc[j].end());
				v_segments_loc[j].clear();
				v_segments_loc[j].shrink_to_fit();
			}

		if (ctg.size() > (uint64_t)to + 1)
			ctg.resize((uint64_t)to + 1);
		if (from != 0)
			ctg.erase(ctg.begin(), ctg.begin() + from);
	}
	else
		ctg.clear();

	if (need_free_zstd)
		ZSTD_freeDCtx(zstd_ctx);

	return true;
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::Open(const string& _archive_fn, const bool _prefetch_archive)
{
	if (working_mode != working_mode_t::none)
		return false;

	in_archive_name = _archive_fn;
	prefetch_archive = _prefetch_archive;

	working_mode = working_mode_t::decompression;

	if (!load_file_type_info(in_archive_name))
		return false;

	if (!load_metadata())
		return false;

	return true;
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::close_decompression()
{
	if (working_mode != working_mode_t::decompression)
		return false;

	return true;
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::decompress_segment(const uint32_t group_id, const uint32_t in_group_id, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
	CSegment segment("seg-" + to_string(group_id), in_archive, nullptr, compression_params.pack_cardinality, compression_params.min_match_len, false, archive_version);

	if (group_id < no_raw_groups)
		return segment.get_raw(in_group_id, ctg, zstd_ctx);
	else
		return segment.get(in_group_id, ctg, zstd_ctx);
}

// *******************************************************************************************
void CAGCDecompressorLibrary::GetCmdLines(vector<pair<string, string>>& _cmd_lines)
{
	_cmd_lines.clear();

	if (working_mode != working_mode_t::decompression)
		return;

	collection_desc.get_cmd_lines(_cmd_lines);
}

// *******************************************************************************************
void CAGCDecompressorLibrary::GetParams(uint32_t& _kmer_length, uint32_t& _min_match_len, uint32_t& _pack_cardinality)
{
	if (working_mode != working_mode_t::decompression)
		return;

	_kmer_length = kmer_length;
	_min_match_len = min_match_len;
	_pack_cardinality = pack_cardinality;
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::Close()
{
	bool r = true;

	if (working_mode == working_mode_t::none)
		r = false;
	else if (working_mode == working_mode_t::decompression)
		r = close_decompression();

	working_mode = working_mode_t::none;

	return r;
}

// EOF
