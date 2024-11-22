#ifndef _COLLECTION_H
#define _COLLECTION_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include <map>
#include <array>
#include <vector>
#include <string>
#include <mutex>
#include <chrono>
#include "../common/utils.h"
#include <zstd/lib/zstd.h>

using namespace std;
using namespace std::chrono;

// *******************************************************************************************
struct segment_desc_t
{
	uint32_t group_id;
	uint32_t in_group_id;
	bool is_rev_comp;
	uint32_t raw_length;

	segment_desc_t() :
		group_id(~0u), in_group_id(~0u), is_rev_comp(false), raw_length(0)
	{}

	segment_desc_t(const uint32_t _group_id, const uint32_t _in_group_id, const bool _is_rev_comp, const uint32_t _raw_length) :
		group_id(_group_id), in_group_id(_in_group_id), is_rev_comp(_is_rev_comp), raw_length(_raw_length)
	{}
};
	
// *******************************************************************************************
struct pair_segment_desc_t
{
	segment_desc_t first;
	segment_desc_t second;
	bool contains_second;

	pair_segment_desc_t(segment_desc_t _first, segment_desc_t _second = segment_desc_t{}, bool _contains_second = false) :
		first(_first), second(_second), contains_second(_contains_second)
	{}
};

// *******************************************************************************************
struct contig_info_t
{
	string sample_name;
	string contig_name;
	uint32_t id;
	uint32_t no_seg;

	contig_info_t(string _sample_name, string _contig_name, uint32_t _id, uint32_t _no_seg) :
		sample_name(_sample_name), contig_name(_contig_name), id(_id), no_seg(_no_seg)
	{};
};

// *******************************************************************************************
struct segments_to_place_t {
	string sample_name;
	string contig_name;
	uint32_t seg_part_no;
	int group_id;
	int in_group_id;
	bool is_rev_comp;
	uint32_t data_size;

	segments_to_place_t(string _sample_name, string _contig_name, uint32_t _seg_part_no, int _group_id, int _in_group_id, bool _is_rev_comp, uint32_t _data_size) :
		sample_name(_sample_name),
		contig_name(_contig_name),
		seg_part_no(_seg_part_no),
		group_id(_group_id),
		in_group_id(_in_group_id),
		is_rev_comp(_is_rev_comp),
		data_size(_data_size) {}

	segments_to_place_t() = default;
	segments_to_place_t(const segments_to_place_t&) = default;
};

// *******************************************************************************************
using sample_desc_t = vector<pair<string, vector<segment_desc_t>>>;

// *******************************************************************************************
class CCollection
{
protected:
	mutex mtx;

	const uint32_t thr_1 = 1u << 7;
	const uint32_t thr_2 = thr_1 + (1u << 14);
	const uint32_t thr_3 = thr_2 + (1u << 21);
	const uint32_t thr_4 = thr_3 + (1u << 28);
	const uint8_t pref_1 = 0;
	const uint8_t pref_2 = 0b10000000u;
	const uint8_t pref_3 = 0b11000000u;
	const uint8_t pref_4 = 0b11100000u;
	const uint8_t pref_5 = 0b11110000u;
	const uint8_t mask_1 = 0b10000000u;
	const uint8_t mask_2 = 0b11000000u;
	const uint8_t mask_3 = 0b11100000u;
	const uint8_t mask_4 = 0b11110000u;

	const uint8_t pref_arr[16] = { 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5};

	ZSTD_DCtx* zstd_dctx = nullptr;

	vector<pair<string, string>> cmd_lines;

	void append(vector<uint8_t>& data, const string& str)
	{
		data.insert(data.end(), str.begin(), str.end());
		data.emplace_back(0);
	}

	void append(vector<uint8_t>& data, uint32_t num)
	{
		if (num < thr_1)
			data.emplace_back(pref_1 + num);
		else if (num < thr_2)
		{
			num -= thr_1;
			data.emplace_back(pref_2 + (num >> 8));
			data.emplace_back(num & 0xffu);
		}
		else if (num < thr_3)
		{
			num -= thr_2;
			data.emplace_back(pref_3 + (num >> 16));
			data.emplace_back((num >> 8) & 0xffu);
			data.emplace_back(num & 0xffu);
		}
		else if (num < thr_4)
		{
			num -= thr_3;
			data.emplace_back(pref_4 + (num >> 24));
			data.emplace_back((num >> 16) & 0xffu);
			data.emplace_back((num >> 8) & 0xffu);
			data.emplace_back(num & 0xffu);
		}
		else
		{
			num -= thr_4;
			data.emplace_back(pref_5);
			data.emplace_back((num >> 24) & 0xffu);
			data.emplace_back((num >> 16) & 0xffu);
			data.emplace_back((num >> 8) & 0xffu);
			data.emplace_back(num & 0xffu);
		}
	}

	void read(uint8_t*& p, string& str)
	{
		auto q = p;
		while (*q)
			++q;

		str.assign((char*)p, q - p);

		p = q + 1;
	}

	void read(uint8_t*& p, uint32_t& num)
	{
		if ((*p & mask_1) == pref_1)
			num = *p++ - pref_1;
		else if ((*p & mask_2) == pref_2)
		{
			num = ((uint32_t)p[0] << 8) + p[1] + thr_1 - (pref_2 << 8);
			p += 2;
		}
		else if ((*p & mask_3) == pref_3)
		{
			num = ((uint32_t) p[0] << 16) + ((uint32_t) p[1] << 8) + p[2] + thr_2 - (pref_3 << 16);
			p += 3;
		}
		else if ((*p & mask_4) == pref_4)
		{
			num = ((uint32_t)p[0] << 24) + ((uint32_t)p[1] << 16) + ((uint32_t)p[2] << 8) + p[3] + thr_3 - (pref_4 << 24);
			p += 4;
		}
		else
		{
			p++;		// skip pref_5
			num = *p++;
			num <<= 8;		num += *p++;
			num <<= 8;		num += *p++;
			num <<= 8;		num += *p++;
			num += thr_4;
		}
	}

	void read_fixed32(vector<uint8_t>::iterator& p, uint32_t& num)
	{
		num = 0;

		for (int i = 0; i < 4; ++i)
			num += ((uint32_t)p[i]) << (8 * i);

		p += 4;
	}

	void skip(uint8_t*& p)
	{
		auto x = pref_arr[*p >> 4];
		p += x;
	}

	string extract_contig_name(const string& s);
	bool is_equal_sample_contig(const pair<string, string>& x, const pair<string, string>& y);

public:
	CCollection() {};
	virtual ~CCollection() { 
		if (zstd_dctx)
			ZSTD_freeDCtx(zstd_dctx);
	};

	virtual bool register_sample_contig(const string& sample_name, const string& contig_name) = 0;

	virtual void add_segment_placed(const string& sample_name, const string& contig_name, const uint32_t place, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length) = 0;
	virtual void add_segments_placed(vector<segments_to_place_t>& segments_to_place) = 0;
	virtual bool get_reference_name(string& reference_name) = 0;
	virtual bool get_samples_list(vector<string>& v_samples, bool sorted = true) = 0;
	virtual bool get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names) = 0;

	virtual bool get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc) = 0;
	virtual bool get_contig_desc(const string& sample_name, string& contig_name, vector<segment_desc_t>& contig_desc) = 0;

	virtual bool is_contig_desc(const string& sample_name, const string& contig_name) = 0;
	virtual vector<string> get_samples_for_contig(const string& contig_name) = 0;

	void add_cmd_line(const string &cmd);
	void get_cmd_lines(vector<pair<string, string>>& _cmd_lines);

	virtual size_t get_no_samples() = 0;
	virtual int32_t get_no_contigs(const string& sample_name) = 0;
};

// EOF
#endif