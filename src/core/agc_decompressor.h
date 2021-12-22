#ifndef _AGC_DECOMPRESSOR_H
#define _AGC_DECOMPRESSOR_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include "../core/agc_decompressor_lib.h"

// *******************************************************************************************
// Class supporting only decompression of AGC files - extended version (can store also in gzipped files)
class CAGCDecompressor : public CAGCDecompressorLibrary
{
public:
	CAGCDecompressor(bool _is_app_mode);
	~CAGCDecompressor();

	bool GetSampleFile(const string& _file_name, const vector<string>& sample_names, const uint32_t _line_length, const uint32_t no_threads);
	bool GetContigFile(const string& _file_name, const vector<string>& contig_names, const uint32_t _line_length, const uint32_t no_threads);
};

// EOF
#endif