#ifndef _APPLICATION_H
#define _APPLICATION_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-02-24
// *******************************************************************************************

#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <algorithm>

using namespace std;

// *******************************************************************************************
template<typename T> class b_value
{
	T value;
	T min_value;
	T max_value;

public:
	b_value(const T _value, const T _min_value, const T _max_value) : value(_value), min_value(_min_value), max_value(_max_value) {};

	void assign(const T _value)
	{
		value = clamp(_value, min_value, max_value);
	}

	string info() const
	{
		return "(default: " + to_string(value) + "; min: " + to_string(min_value) + "; max: " + to_string(max_value) + ")";
	}

	T operator()() const
	{
		return value;
	}
};

// *******************************************************************************************
struct CParams
{
	vector<string> input_names;
	string in_archive_name;
	string out_archive_name;
	string kmc_db_name;
	string reference_name;
	string splitters_name;
	string output_name;
	vector<string> sample_names;
	vector<string> contig_names;
	string contig_name;
	string mode;

	b_value<uint32_t> k{ 31, 17, 32 };
	b_value<uint32_t> pack_cardinality{ 50, 1, 1'000'000'000 };
	b_value<uint32_t> segment_size{ 60'000, 100, 1'000'000 };
	b_value<uint32_t> min_match_length{ 20, 15, 32 };
	b_value<uint32_t> no_threads{ max<uint32_t>(1, thread::hardware_concurrency() / 2), 1, max<uint32_t>(16, thread::hardware_concurrency()) };
	b_value<uint32_t> line_length{ 80, 40, 2'000'000'000 };
	b_value<uint32_t> verbosity{ 0, 0, 2 };

	uint32_t no_segments = 0;
	bool concatenated_genomes = false;
	bool reproducibility_mode = true;
	bool use_stdout = true;
	bool store_cmd_line = true;
	bool prefetch = true;
	bool adaptive_compression = false;

	CParams() = default;
};

// *******************************************************************************************
class CApplication
{
	CParams execution_params;
	string cmd_line;

	bool parse_params(const int argc, const char** argv);
	void usage() const;
	void usage_create() const;
	void usage_append() const;
	void usage_getset() const;
	void usage_getctg() const;
	void usage_listset() const;
	void usage_listctg() const;
	void usage_info() const;

	string bool2string(bool x) const;

	bool load_file_names(const string & fn, vector<string>& v_file_names);

	bool parse_params_create(const int argc, const char** argv);
	bool parse_params_append(const int argc, const char** argv);
	bool parse_params_getset(const int argc, const char** argv);
	bool parse_params_getctg(const int argc, const char** argv);
	bool parse_params_listset(const int argc, const char** argv);
	bool parse_params_listctg(const int argc, const char** argv);
	bool parse_params_info(const int argc, const char** argv);

	bool create();
	bool append();
	bool getset();
	bool getctg();
	bool listset();
	bool listctg();
	bool info();

public:
	CApplication() = default;
	~CApplication() = default;

	int Run(const int argc, const char** argv);
};

// EOF
#endif