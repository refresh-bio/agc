#ifndef _AGC_BASIC_H
#define _AGC_BASIC_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <tuple>
#include <thread>
#include <atomic>
#include <shared_mutex>
#include "../core/archive.h"
#include "../core/segment.h"
#include "../core/collection_v1.h"
#include "../core/collection_v2.h"
#include "../core/collection_v3.h"
#include "../core/queue.h"

using namespace std;

// *******************************************************************************************
// Basic compression class
class CAGCBasic
{
	friend class CAGCDecompressor;

protected:
	enum class working_mode_t { none, compression, decompression, appending, pre_appending };

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
	struct compression_params_t
	{
		uint32_t kmer_length;
		uint32_t min_match_len;
		uint32_t pack_cardinality;
		uint32_t segment_size;
	};

	working_mode_t working_mode;
	bool is_app_mode;

	uint32_t kmer_length;
	uint32_t min_match_len;
	uint32_t pack_cardinality;
	uint32_t segment_size;

	string in_archive_name;
	bool prefetch_archive = false;
	uint32_t archive_version;

	shared_ptr<CArchive> in_archive;															// internal mutexes

	shared_ptr<CCollection> collection_desc;

	map<string, string> m_file_type_info;

	compression_params_t compression_params;

	const uint32_t no_raw_groups = 16;

	uint32_t verbosity;

	// *******************************************************************************************
	void read(vector<uint8_t>::iterator& p, uint32_t& num)
	{
		num = 0;

		for (int i = 0; i < 4; ++i)
			num += ((uint32_t)p[i]) << (8 * i);

		p += 4;
	}

	// *******************************************************************************************
	void read64(vector<uint8_t>::iterator& p, uint64_t& num)
	{
		num = 0;

		for (int i = 0; i < 8; ++i)
			num += ((uint64_t)p[i]) << (8 * i);

		p += 8;
	}

	// *******************************************************************************************
	void read(vector<uint8_t>::iterator& p, string& str)
	{
		str.clear();

		for (; *p != 0; ++p)
			str.push_back((char)*p);
		++p;
	}

	// *******************************************************************************************
	void join_threads(vector<thread> &v_threads);
	bool load_metadata_impl_v1();
	bool load_metadata_impl_v2();
	bool load_metadata_impl_v3();

	bool load_metadata();
	bool load_file_type_info(const string& archive_name);

	void reverse_complement(contig_t& contig);
	void reverse_complement_copy(contig_t& src_contig, contig_t& dest_contig);

public:
	CAGCBasic();
	~CAGCBasic();
};

// EOF
#endif