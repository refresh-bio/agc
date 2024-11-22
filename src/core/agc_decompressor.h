#ifndef _AGC_DECOMPRESSOR_H
#define _AGC_DECOMPRESSOR_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "../common/agc_decompressor_lib.h"
#include <refresh/compression/lib/gz_wrapper.h>

// *******************************************************************************************
// Class supporting only decompression of AGC files - extended version (can store also in gzipped files)
class CAGCDecompressor : public CAGCDecompressorLibrary
{
	void start_decompressing_threads(vector<thread>& v_threads, const uint32_t n_t, uint32_t gzip_level = 0, uint32_t line_len = 0, bool fast = false);

	void gzip_contig(contig_t& ctg, contig_t& working_space, refresh::gz_in_memory& gzip_compressor);

public:
	CAGCDecompressor(bool _is_app_mode);
	~CAGCDecompressor();

	bool GetCollectionFiles(const string& _path, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, bool no_ref, bool fast, uint32_t verbosity);
	bool GetSampleFile(const string& _file_name, const vector<string>& sample_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity);
	bool GetContigFile(const string& _file_name, const vector<string>& contig_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity);
	bool GetSampleForStreaming(const string& _file_name, const vector<string>& sample_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity);
	bool GetContigForStreaming(const string& _file_name, const vector<string>& contig_names, const uint32_t _line_length, const uint32_t no_threads, const uint32_t gzip_level, uint32_t verbosity);

	bool GetSampleSequences(const string& sample_name, vector<pair<string, vector<uint8_t>>> &v_contig_seq, const uint32_t no_threads);

	bool AssignArchive(const CAGCBasic &agc_basic);
};

// EOF
#endif