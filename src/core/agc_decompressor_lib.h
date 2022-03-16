#ifndef _AGC_DECOMPRESSOR_LIB_H
#define _AGC_DECOMPRESSOR_LIB_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-02-24
// *******************************************************************************************

#include <regex>
#include "../core/agc_basic.h"

// *******************************************************************************************
// Class supporting only decompression of AGC files - library version
class CAGCDecompressorLibrary : public CAGCBasic
{
protected:
	const uint8_t cnv_num[128] = {
		//		0    1	   2    3    4    5    6    7    8    9    10   11   12  13   14    15
				'A', 'C', 'G', 'T', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'U',
				' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
				' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
				' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
				' ',   0,  11,   1,  12,  30,  30,   2,  13,  30,  30,   9,  30,  10,   4,  30,
				 30,  30,   5,   7,   3,  15,  14,   8,  30,   6,  30,  30,  30,  30,  30,  30,
				' ',   0,  11,   1,  12,  30,  30,   2,  13,  30,  30,   9,  30,  10,   4,  30,
				 30,  30,   5,   7,   3,  15,  14,   8,  30,   6,  30,  30,  30,  30,  30,  30
	};

	// *******************************************************************************************
	struct name_range_t
	{
		string name;
		int64_t from;
		int64_t to;

		name_range_t(const string _name = "", const int64_t _from = -1, const int64_t _to = -1) :
			name(_name), from(_from), to(_to)
		{}

		string str() const
		{
			if (from >= 0 && to >= 0)
				return name + ":" + to_string(from) + "-" + to_string(to);
			else
				return name;
		}
	};

	typedef tuple<size_t, name_range_t, vector<segment_desc_t>> task_desc_t;

	const regex re_csr = regex("(.+)@(.+):(.+)-(.+)");
	const regex re_cs = regex("(.+)@(.+)");
	const regex re_cr = regex("(.+):(.+)-(.+)");
	const regex re_c = regex("(.+)");

	unique_ptr<CBoundedQueue<tuple<size_t, name_range_t, vector<segment_desc_t>>>> q_contig_tasks;
	unique_ptr<CPriorityQueue<pair<string, contig_t>>> pq_contigs_to_save;

	bool analyze_contig_query(const string& query, string& sample, name_range_t& name_range);
	void start_decompressing_threads(vector<thread>& v_threads, const uint32_t n_t);
	bool decompress_segment(const uint32_t group_id, const uint32_t in_group_id, contig_t& ctg, ZSTD_DCtx* zstd_ctx);

	bool decompress_contig(task_desc_t& task, ZSTD_DCtx *zstd_ctx, contig_t& ctg);

	bool close_decompression();

public:
	CAGCDecompressorLibrary(bool _is_app_mode);
	~CAGCDecompressorLibrary();

	bool Open(const string& _archive_fn, const bool _prefetch_archive = false);

	void GetCmdLines(vector<pair<string, string>>& _cmd_lines);
	void GetParams(uint32_t& kmer_length, uint32_t& min_match_len, uint32_t& pack_cardinality);

	bool Close();

	int GetContigString(const string& sample_name, const string& contig_name, const int start, const int end, string& contig_data);
	int64_t GetContigLength(const string& sample_name, const string& contig_name);

	bool ListSamples(vector<string>& v_sample_names);
	bool ListContigs(const string& sample_name, vector<string>& v_contig_names);
	int32_t GetNoSamples();
	int32_t GetNoContigs(const string& sample_name);

	void GetFileTypeInfo(map<string, string>& _m_file_type_info);

	bool IsOpened();
};

// EOF
#endif