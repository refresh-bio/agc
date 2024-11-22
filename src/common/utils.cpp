// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "utils.h"
#include <algorithm>

// *******************************************************************************************
string int_to_hex(uint32_t n)
{
	const char dig[] = "0123456789ABCDEF";

	string res;

	do
	{
		res.push_back(dig[n & 0xfu]);
		n /= 16;
	} while (n);

//	res.reserve();

	return res;
}

// *******************************************************************************************
string int_to_base64(uint32_t n)
{
	const char dig[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_#";

	string res;

	do
	{
		res.push_back(dig[n & 0x3fu]);
		n /= 64;
	} while (n);

//	res.reserve();

	return res;
}

// *******************************************************************************************
string ss_prefix(uint32_t archive_version)
{
	if (archive_version < 3000)
		return "seg-";
	else
		return "x";
}

// *******************************************************************************************
string ss_base(uint32_t archive_version, uint32_t n)
{
	if (archive_version < 3000)
		return "seg-" + to_string(n);
	else
		return "x" + int_to_base64(n);
}

// *******************************************************************************************
string ss_ref_name(uint32_t archive_version, uint32_t n)
{
	if (archive_version < 3000)
		return "seg-" + to_string(n) + "-ref";
	else
		return "x" + int_to_base64(n) + "r";
}

// *******************************************************************************************
string ss_delta_name(uint32_t archive_version, uint32_t n)
{
	if (archive_version < 3000)
		return "seg-" + to_string(n) + "-delta";
	else
		return "x" + int_to_base64(n) + "d";
}

// *******************************************************************************************
string ss_ref_ext(uint32_t archive_version)
{
	if (archive_version < 3000)
		return "-ref";
	else
		return "r";
}

// *******************************************************************************************
string ss_delta_ext(uint32_t archive_version)
{
	if (archive_version < 3000)
		return "-delta";
	else
		return "d";
}

// EOF
