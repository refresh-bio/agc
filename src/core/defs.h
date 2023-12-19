#ifndef _DEFS_H
#define _DEFS_H

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
#include <cstdint>
using namespace std;

typedef vector<uint8_t> contig_t;
typedef vector<uint8_t> packed_block_t;

const uint32_t AGC_VER_MAJOR = 3;
const uint32_t AGC_VER_MINOR = 0;
const string AGC_VER_BUILD = "20221222.1"s;

const uint32_t AGC_FILE_MAJOR = 3;
const uint32_t AGC_FILE_MINOR = 0;

const std::string AGC_VERSION = std::string("AGC (Assembled Genomes Compressor) v. ") + 
	to_string(AGC_VER_MAJOR) + "." + to_string(AGC_VER_MINOR) +
	" [build " + AGC_VER_BUILD + "]";

#define IMPROVED_LZ_ENCODING

#include <iostream>

#if defined(_MSC_VER)  /* Visual Studio */
#define REFRESH_FORCE_INLINE __forceinline
#define REFRESH_NO_INLINE __declspec(noinline)
#elif defined(__GNUC__)
#define REFRESH_FORCE_INLINE __inline__ __attribute__((always_inline, unused))
#define REFRESH_NO_INLINE __attribute__((noinline))
#else
#define REFRESH_FORCE_INLINE
#define REFRESH_NO_INLINE
#endif

// EOF
#endif