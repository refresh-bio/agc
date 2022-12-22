// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include <ctime>
#include <iomanip>
#include <algorithm>
#include "../core/collection.h"

#include <iostream>

// *******************************************************************************************
string CCollection::extract_contig_name(const string& s)
{
	string::const_iterator p;

	for (p = s.begin(); p != s.end(); ++p)
		if ((*p < '0') && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t'))
			break;

	return string(s.begin(), p);
}

// *******************************************************************************************
bool CCollection::is_equal_sample_contig(const pair<string, string>& x, const pair<string, string>& y)
{
	return x.first == y.first && extract_contig_name(x.second) == extract_contig_name(y.second);
}


#if 0
// *******************************************************************************************
bool CCollection::get_samples_info(map<string, vector<string>>& v_samples)
{
	lock_guard<mutex> lck(mtx);

	v_samples.clear();

	for (auto& p : col)
	{
		auto& q = v_samples.emplace(p.first, vector<string>()).first->second;

		for (auto& r : p.second)
			q.emplace_back(r.first);
	}

	return true;
}
#endif

// *******************************************************************************************
void CCollection::add_cmd_line(const string &cmd)
{
	lock_guard<mutex> lck(mtx);

	auto tc = time(nullptr);
	char tmp[64];
	string s_time;

	if(strftime(tmp, sizeof(tmp), "%A %c", std::gmtime(&tc)))
		s_time = tmp;

	cmd_lines.emplace_back(cmd, s_time);
}

// *******************************************************************************************
void CCollection::get_cmd_lines(vector<pair<string, string>>& _cmd_lines)
{
	lock_guard<mutex> lck(mtx);

	_cmd_lines = cmd_lines;
}

// EOF
