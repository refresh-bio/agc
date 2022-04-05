#ifndef _GENOME_IO_H
#define _GENOME_IO_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-04-05
// *******************************************************************************************

#include <cstdio>
#include <vector>
#include <string>
#include <cinttypes>
#include "defs.h"
#include "../../libs/zlib.h"

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

	FILE* in;
	FILE* out;
	gzFile gz_in;
	gzFile gz_out;
	bool is_gzipped;
	bool use_stdout;

	uint8_t* buffer;
	const size_t buffer_size = 128 << 20;
	const size_t gz_buffer_size = 32 << 20;
	size_t buffer_filled;
	size_t buffer_pos;

	bool fill_buffer();
	bool eof() { return buffer_pos == buffer_filled; }
	int find_contig_end();

	bool read_contig(string& id, contig_t& contig, const bool converted);
	bool read_contig_raw(string& id, contig_t& contig);
	bool save_contig(const string& id, const contig_t& contig, const uint32_t line_length, const bool converted);

	void save_contig_imp(const contig_t& contig, const uint32_t line_length);
	void save_contig_imp_cnv(const contig_t& contig, const uint32_t line_length);

public:
	CGenomeIO();
	~CGenomeIO();

	bool Open(const string &_file_name, const bool _writing);
	bool Close();
	size_t FileSize();

	bool ReadContig(string &id, contig_t&contig);
	bool ReadContigConverted(string& id, contig_t& contig);
	bool ReadContigRaw(string& id, contig_t& contig);

	bool SaveContig(const string& id, const contig_t& contig, const uint32_t line_length);
	bool SaveContigConverted(const string& id, const contig_t& contig, const uint32_t line_length);
};

// EOF
#endif