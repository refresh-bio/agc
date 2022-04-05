#ifndef _DEFS_H
#define _DEFS_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-04-05
// *******************************************************************************************

#include <string>
#include <vector>
using namespace std;

typedef vector<uint8_t> contig_t;
typedef vector<uint8_t> packed_block_t;

const uint32_t AGC_VER_MAJOR = 2;
const uint32_t AGC_VER_MINOR = 0;
const string AGC_VER_BUILD = "20220405.1"s;
const uint32_t AGC_FILE_MAJOR = 2;
const uint32_t AGC_FILE_MINOR = 0;


const std::string AGC_VERSION = std::string("AGC (Assembled Genomes Compressor) v. ") + 
	to_string(AGC_VER_MAJOR) + "." + to_string(AGC_VER_MINOR) +
	" [build " + AGC_VER_BUILD + "]";

#define IMPROVED_LZ_ENCODING

// EOF
#endif