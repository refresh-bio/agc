#ifndef _GENOME_IO_H
#define _GENOME_IO_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include <cstdio>
#include <vector>
#include <string>
#include <cinttypes>
#include "../common/defs.h"
#include <refresh/compression/lib/file_wrapper.h>
#include <refresh/compression/lib/gz_wrapper.h>

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;

class CGenomeIO
{
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

	string file_name;
	bool writing;

	FILE* out;
	bool is_gzipped;
	bool use_stdout;

	refresh::stream_in_file *sif = nullptr;
	refresh::stream_decompression* sdf = nullptr;

	refresh::gz_in_memory gzip_zero_compressor{ 1 };
	vector<uint8_t> gzip_zero_compressor_buffer;

	uint8_t* buffer;
	const size_t write_buffer_size = 32 << 20;
	const size_t read_buffer_size = 4 << 20;
	size_t buffer_filled;
	size_t buffer_pos;

	bool fill_buffer();
	bool eof() { return buffer_pos == buffer_filled; }
	int find_contig_end();

	bool read_contig(string& id, contig_t& contig, const bool converted);
	bool read_contig_raw(string& id, contig_t& contig);

	bool save_contig_directly(const string& id, const contig_t& contig, const uint32_t gzip_level);

#if 0
	bool save_contig(const string& id, const contig_t& contig, const uint32_t line_length, const bool converted);

	void save_contig_imp(const contig_t& contig, const uint32_t line_length);
	void save_contig_imp_cnv(const contig_t& contig, const uint32_t line_length);
#endif

public:
	CGenomeIO();
	~CGenomeIO();

	bool Open(const string &_file_name, const bool _writing);
	bool Close();
	size_t FileSize();

	bool ReadContig(string &id, contig_t&contig);
	bool ReadContigConverted(string& id, contig_t& contig);
	bool ReadContigRaw(string& id, contig_t& contig);

	bool SaveContigDirectly(const string& id, const contig_t& contig, const uint32_t gzip_level);
#if 0
	bool SaveContig(const string& id, const contig_t& contig, const uint32_t line_length);
	bool SaveContigConverted(const string& id, const contig_t& contig, const uint32_t line_length);
#endif
};
#endif

// EOF
