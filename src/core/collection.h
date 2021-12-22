#ifndef _COLLECTION_H
#define _COLLECTION_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include <map>
#include <vector>
#include <string>
#include <mutex>
#include <chrono>
#include "../core/utils.h"

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
using sample_desc_t = vector<pair<string, vector<segment_desc_t>>>;

// *******************************************************************************************
class CCollection
{
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

	map<string, vector<pair<string, vector<segment_desc_t>>>> col;
	map<pair<string, string>, uint32_t> contig_ids;

	multimap<string, string> mm_contig2sample;

	vector<pair<string, string>> cmd_lines;

	void append(vector<uint8_t>& data, const string& str);
	void append(vector<uint8_t>& data, uint32_t num);

	void read(vector<uint8_t>::iterator& p, string& str);
	void read(vector<uint8_t>::iterator& p, uint32_t& num);

	string extract_contig_name(const string &s);
	bool is_equal_sample_contig(const pair<string, string>& x, const pair<string, string>& y);

	vector<segment_desc_t> & add_segment_basic(const string& sample_name, const string& contig_name, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length);

public:
	CCollection() {};

	bool register_sample_contig(const string& sample_name, const string& contig_name);
	void add_segment_append(const string& sample_name, const string& contig_name, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length);
	void add_segment_placed(const string& sample_name, const string& contig_name, const uint32_t place, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length);

	bool get_samples_list(vector<string>& v_samples);
	bool get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names);
	bool get_samples_info(map<string, vector<string>>& v_samples);
	bool get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc);
	bool get_contig_desc(const string& sample_name, const string& contig_name, vector<segment_desc_t>& contig_desc);
	bool is_contig_desc(const string& sample_name, const string& contig_name);
	vector<string> get_samples_for_contig(const string& contig_name);

	void add_cmd_line(const string &cmd);

	void get_cmd_lines(vector<pair<string, string>>& _cmd_lines);

	size_t get_no_samples();
	int32_t get_no_contigs(const string& sample_name);

	void serialize(vector<uint8_t>& data, bool store_date_time);
	bool deserialize(vector<uint8_t>& data);
};

// EOF
#endif