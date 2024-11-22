// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "agc_decompressor_lib.h"
#include <cassert>

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
	return collection_desc->get_samples_list(v_sample_names);
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::ListContigs(const string& sample_name, vector<string>& v_contig_names)
{
	return collection_desc->get_contig_list_in_sample(sample_name, v_contig_names);
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
	return static_cast<int32_t>(collection_desc->get_no_samples());
}

// *******************************************************************************************
int32_t CAGCDecompressorLibrary::GetNoContigs(const string& sample_name)
{
	return collection_desc->get_no_contigs(sample_name);
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
		auto v_cand_samples = collection_desc->get_samples_for_contig(contig_name);
		if (v_cand_samples.size() == 0)
			return -1;
		if (v_cand_samples.size() > 1)
			return -2;

		det_sample_name = v_cand_samples.front();
	}

	string full_contig_name = contig_name;

	if (!collection_desc->get_contig_desc(det_sample_name, full_contig_name, contig_desc))
		return -1;

	contig_task_t task{ id++, "", name_range_t(full_contig_name, start, end), contig_desc };
	contig_t ctg;

	decompress_contig(task, nullptr, ctg);

	contig_data.clear();
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
		auto v_cand_samples = collection_desc->get_samples_for_contig(contig_name);
		if (v_cand_samples.size() == 0)
			return -1;
		if (v_cand_samples.size() > 1)
			return -2;

		det_sample_name = v_cand_samples.front();
	}

	string full_contig_name = contig_name;

	if (!collection_desc->get_contig_desc(det_sample_name, full_contig_name, contig_desc))
		return -1;

	int64_t len = 0;
	for (auto& x : contig_desc)
		len += x.raw_length;

	return len - (contig_desc.size() - 1) * kmer_length;
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::decompress_contig(contig_task_t& contig_desc, ZSTD_DCtx* zstd_ctx, contig_t& ctg, bool fast)
{
	name_range_t &contig_name_range = contig_desc.name_range;
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

	v_segments_loc.reserve(contig_desc.segments.size());

	for (auto seg : contig_desc.segments)
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

		if(!fast)
			decompress_segment(seg.group_id, seg.in_group_id, ctg, zstd_ctx);
		else
			decompress_segment_fast(seg.group_id, seg.in_group_id, ctg, zstd_ctx);

		if (seg.is_rev_comp)
			reverse_complement(ctg);

		v_segments_loc.emplace_back(move(ctg));

		curr_pos += seg_len - kmer_length;
	}

	if (!v_segments_loc.empty())
	{
		size_t req_size = 0;
		for (uint32_t j = 0; j < v_segments_loc.size(); ++j)
			req_size += v_segments_loc[j].size();

		ctg.clear();
		ctg.reserve(req_size);
		ctg.insert(ctg.end(), v_segments_loc.front().begin(), v_segments_loc.front().end());

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
bool CAGCDecompressorLibrary::decompress_contig_streaming(contig_task_t& contig_desc, ZSTD_DCtx* zstd_ctx, CStreamWrapper& stream_wrapper, bool fast)
{
	name_range_t &contig_name_range = contig_desc.name_range;
//	vector<contig_t> v_segments_loc;

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

//	v_segments_loc.reserve(contig_desc.segments.size());
	contig_t ctg;
	bool first_processed = false;

	contig_t local_ctg;

	for (auto seg : contig_desc.segments)
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

		if(!fast)
			decompress_segment(seg.group_id, seg.in_group_id, ctg, zstd_ctx);
		else
			decompress_segment_fast(seg.group_id, seg.in_group_id, ctg, zstd_ctx);

		if (seg.is_rev_comp)
			reverse_complement(ctg);

		if (first_processed)
		{
			if (ctg.size() < compression_params.kmer_length)
			{
				if (is_app_mode)
					cerr << "Corrupted archive!" << endl;
			}
			else
				local_ctg.assign(ctg.begin() + compression_params.kmer_length, ctg.end());
//				stream_wrapper.append(ctg.begin() + compression_params.kmer_length, ctg.end());
		}
		else
//			stream_wrapper.append(ctg.begin() + from, ctg.end());
			local_ctg.assign(ctg.begin(), ctg.end());

		if (local_ctg.size() > (uint64_t)to + 1)
			local_ctg.resize((uint64_t)to + 1);

		stream_wrapper.append(local_ctg.begin() + from, local_ctg.end());

		curr_pos += seg_len - kmer_length;
		first_processed = true;
	}

	if (need_free_zstd)
		ZSTD_freeDCtx(zstd_ctx);

	stream_wrapper.complete_contig();

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

	if (archive_version < 3000)
	{
		if (!load_metadata())
			return false;
	}
	else if (archive_version < 4000)
	{
		if (!load_metadata() ||
			!dynamic_pointer_cast<CCollection_V3>(collection_desc)->set_archives(in_archive, nullptr, 1, pack_cardinality, segment_size, kmer_length))
			return false;
	}

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
	CSegment segment(ss_base(archive_version, group_id), in_archive, nullptr, compression_params.pack_cardinality, compression_params.min_match_len, false, archive_version);

	if (group_id < no_raw_groups)
		return segment.get_raw(in_group_id, ctg, zstd_ctx);
	else
		return segment.get(in_group_id, ctg, zstd_ctx);
}

// *******************************************************************************************
bool CAGCDecompressorLibrary::decompress_segment_fast(const uint32_t group_id, const uint32_t in_group_id, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
	shared_ptr<CSegment> segment;

	{
		shared_lock<shared_mutex> lck(mtx_segment);

		auto p = v_segment.find(group_id);

		if (p != v_segment.end())
			segment = p->second;
	}

	if (segment == nullptr)
	{
		unique_lock<shared_mutex> lck(mtx_segment);

		auto p = v_segment.find(group_id);		// can happen that other thread will add this segment just before this lock

		if (p != v_segment.end())
			segment = p->second;
		else
		{
			segment = make_shared<CSegment>(ss_base(archive_version, group_id), in_archive, nullptr, compression_params.pack_cardinality, compression_params.min_match_len, false, archive_version, true);
			v_segment[group_id] = segment;
		}
	}

	if (group_id < no_raw_groups)
		return segment->get_raw_locked(in_group_id, ctg, zstd_ctx);
	else
		return segment->get_locked(in_group_id, ctg, zstd_ctx);
}

// *******************************************************************************************
void CAGCDecompressorLibrary::GetCmdLines(vector<pair<string, string>>& _cmd_lines)
{
	_cmd_lines.clear();

	if (working_mode != working_mode_t::decompression)
		return;

	collection_desc->get_cmd_lines(_cmd_lines);
}

// *******************************************************************************************
void CAGCDecompressorLibrary::GetParams(uint32_t& _kmer_length, uint32_t& _min_match_len, uint32_t& _pack_cardinality, uint32_t &_segment_size)
{
	if (working_mode != working_mode_t::decompression)
		return;

	_kmer_length = kmer_length;
	_min_match_len = min_match_len;
	_pack_cardinality = pack_cardinality;
	_segment_size = segment_size;
}

// *******************************************************************************************
void CAGCDecompressorLibrary::GetReferenceSample(string& ref_name)
{
	ref_name.clear();

	if (working_mode != working_mode_t::decompression)
		return;

	collection_desc->get_reference_name(ref_name);
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

// *******************************************************************************************
void CAGCDecompressorLibrary::CNumAlphaConverter::convert_to_alpha(contig_t& ctg)
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
		ctg[i + 1] = cnv_num[ctg[i + 1]];
		ctg[i + 2] = cnv_num[ctg[i + 2]];
		ctg[i + 3] = cnv_num[ctg[i + 3]];
		ctg[i + 4] = cnv_num[ctg[i + 4]];
		ctg[i + 5] = cnv_num[ctg[i + 5]];
		ctg[i + 6] = cnv_num[ctg[i + 6]];
		ctg[i + 7] = cnv_num[ctg[i + 7]];
	}
}

// *******************************************************************************************
size_t CAGCDecompressorLibrary::CNumAlphaConverter::convert_and_split_into_lines(contig_t& ctg, contig_t& working_space, uint32_t line_len, uint32_t no_symbols_in_non_complete_line, bool append_eol)
{
	if (ctg.empty())
		return 0;

	size_t dest_size = ctg.size() + (ctg.size() + line_len - 1) / line_len + 2;
	working_space.resize(dest_size);

	auto p = ctg.data();
	auto q = working_space.data();

	size_t to_save = ctg.size();

	if (no_symbols_in_non_complete_line)
	{
		while (to_save && no_symbols_in_non_complete_line++ < line_len)
		{
			*q++ = cnv_num[*p++];
			to_save--;
		}
		*q++ = '\n';
	}

	if (!to_save)
	{
		working_space.resize(q - working_space.data());

		if (!append_eol)
			working_space.pop_back();

		std::swap(ctg, working_space);

		return no_symbols_in_non_complete_line;
	}

	for (; to_save > line_len; to_save -= line_len)
	{
		uint32_t i;

		switch (i = line_len % 8)
		{
		case 7:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 6:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 5:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 4:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 3:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 2:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 1:	*q++ = cnv_num[*p++];
		}

		for (; i < line_len; i += 8)
		{
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
		}

		*q++ = '\n';
	}

	size_t r = to_save;

	if (to_save)
	{
		while (to_save--)
			*q++ = cnv_num[*p++];
		*q++ = '\n';
	}

	assert(q <= working_space.data() + working_space.size());
	working_space.resize(q - working_space.data());

	if (!append_eol)
		working_space.pop_back();

	std::swap(ctg, working_space);

	return r;		// No of symbols in last line
}

// EOF
