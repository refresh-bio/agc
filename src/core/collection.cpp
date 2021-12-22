// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include <ctime>
#include <iomanip>
#include "../core/collection.h"

// *******************************************************************************************
void CCollection::append(vector<uint8_t>& data, const string& str)
{
	data.insert(data.end(), str.begin(), str.end());
	data.emplace_back(0);
}

// *******************************************************************************************
void CCollection::append(vector<uint8_t>& data, uint32_t num)
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

// *******************************************************************************************
void CCollection::read(vector<uint8_t>::iterator& p, string& str)
{
	str.clear();

	for (; *p; ++p)
		str.push_back(*p);

	++p;
}

// *******************************************************************************************
void CCollection::read(vector<uint8_t>::iterator& p, uint32_t& num)
{
	num = 0;

	if ((*p & mask_1) == pref_1)
		num = *p++ - pref_1;
	else if ((*p & mask_2) == pref_2)
	{
		num = *p++ - pref_2;
		num <<= 8;		num += *p++;
		num += thr_1;
	}
	else if ((*p & mask_3) == pref_3)
	{
		num = *p++ - pref_3;
		num <<= 8;		num += *p++;
		num <<= 8;		num += *p++;
		num += thr_2;
	}
	else if ((*p & mask_4) == pref_4)
	{
		num = *p++ - pref_4;
		num <<= 8;		num += *p++;
		num <<= 8;		num += *p++;
		num <<= 8;		num += *p++;
		num += thr_3;
	}
	else
	{
		p++;		// skip pref_5
		num <<= 8;		num += *p++;
		num <<= 8;		num += *p++;
		num <<= 8;		num += *p++;
		num <<= 8;		num += *p++;
		num += thr_4;
	}
}

// *******************************************************************************************
string CCollection::extract_contig_name(const string& s)
{
	string::const_iterator p;

	for (p = s.begin(); p != s.end(); ++p)
		if (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')
			break;

	return string(s.begin(), p);
}

// *******************************************************************************************
bool CCollection::is_equal_sample_contig(const pair<string, string>& x, const pair<string, string>& y)
{
	return x.first == y.first && extract_contig_name(x.second) == extract_contig_name(y.second);
}

// *******************************************************************************************
bool CCollection::register_sample_contig(const string& sample_name, const string& contig_name)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);
	string stored_sample_name = sample_name;

	if (sample_name.empty())
		stored_sample_name = short_contig_name;

	auto p = contig_ids.find(make_pair(stored_sample_name, short_contig_name));
	if (p == contig_ids.end())
	{
		uint32_t contig_id = (uint32_t)col[stored_sample_name].size();
		contig_ids[make_pair(stored_sample_name, short_contig_name)] = contig_id;
		col[stored_sample_name].emplace_back(contig_name, vector<segment_desc_t>());

		mm_contig2sample.emplace(short_contig_name, stored_sample_name);

		return true;
	}
	else
		return false;		// The pair sample_name:contig_name is not unique
}

// *******************************************************************************************
vector<segment_desc_t> &CCollection::add_segment_basic(const string& sample_name, const string& contig_name,
	const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length)
{
	uint32_t contig_id;

	string short_contig_name = extract_contig_name(contig_name);
	string stored_sample_name = sample_name;

	if (sample_name.empty())
		stored_sample_name = short_contig_name;

	auto p_contig_ids = contig_ids.find(make_pair(stored_sample_name, short_contig_name));
	if (p_contig_ids != contig_ids.end())
		contig_id = p_contig_ids->second;
	else
	{
		contig_id = (uint32_t)col[stored_sample_name].size();
		contig_ids[make_pair(stored_sample_name, short_contig_name)] = contig_id;
		col[stored_sample_name].emplace_back(contig_name, vector<segment_desc_t>());

		mm_contig2sample.emplace(short_contig_name, stored_sample_name);
	}

	return col[stored_sample_name][contig_id].second;
}

// *******************************************************************************************
void CCollection::add_segment_append(const string& sample_name, const string& contig_name,
	const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length)
{
	lock_guard<mutex> lck(mtx);

	auto &vec = add_segment_basic(sample_name, contig_name, group_id, in_group_id, is_rev_comp, raw_length);

	vec.emplace_back(group_id, in_group_id, is_rev_comp, raw_length);
}

// *******************************************************************************************
void CCollection::add_segment_placed(const string& sample_name, const string& contig_name, const uint32_t place,
	const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length)
{
	lock_guard<mutex> lck(mtx);

	auto& vec = add_segment_basic(sample_name, contig_name, group_id, in_group_id, is_rev_comp, raw_length);

	if (place >= vec.size())
		vec.resize(place + 1, segment_desc_t(555555555, 555555555, true, 555555555));

	vec[place] = segment_desc_t(group_id, in_group_id, is_rev_comp, raw_length);
}

// *******************************************************************************************
bool CCollection::get_samples_list(vector<string>& v_samples)
{
	lock_guard<mutex> lck(mtx);

	v_samples.clear();

	for (auto& p : col)
		v_samples.emplace_back(p.first);

	return true;
}

// *******************************************************************************************
bool CCollection::get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names)
{
	lock_guard<mutex> lck(mtx);

	v_contig_names.clear();

	auto p = col.find(sample_name);
	if (p == col.end())
		return false;

	for (auto q : p->second)
		v_contig_names.emplace_back(q.first);

	return true;
}

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

// *******************************************************************************************
bool CCollection::get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc)
{
	lock_guard<mutex> lck(mtx);

	auto p = col.find(sample_name);
	if (p == col.end())
		return false;

	sample_desc = p->second;

	return true;
}

// *******************************************************************************************
bool CCollection::get_contig_desc(const string& sample_name, const string& contig_name, vector<segment_desc_t>& contig_desc)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);

	auto p_contig_id = contig_ids.find(make_pair(sample_name, short_contig_name));
	if (p_contig_id == contig_ids.end())
		return false;

	uint32_t contig_id = p_contig_id->second;

	contig_desc = col[sample_name][contig_id].second;

	return true;
}

// *******************************************************************************************
bool CCollection::is_contig_desc(const string& sample_name, const string& contig_name)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);

	return contig_ids.find(make_pair(sample_name, short_contig_name)) != contig_ids.end();
}

// *******************************************************************************************
vector<string> CCollection::get_samples_for_contig(const string& contig_name)
{
	auto pq = mm_contig2sample.equal_range(extract_contig_name(contig_name));

	vector<string> vs;

	for (auto p = pq.first; p != pq.second; ++p)
		vs.emplace_back(p->second);

	return vs;
}

// *******************************************************************************************
size_t CCollection::get_no_samples()
{
	lock_guard<mutex> lck(mtx);

	return col.size();
}

// *******************************************************************************************
int32_t CCollection::get_no_contigs(const string& sample_name)
{
	lock_guard<mutex> lck(mtx);

	auto p = col.find(sample_name);
	if (p == col.end())
		return -1;

	return static_cast<int32_t>(p->second.size());
}

// *******************************************************************************************
void CCollection::add_cmd_line(const string &cmd)
{
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
	_cmd_lines = cmd_lines;
}

// *******************************************************************************************
void CCollection::serialize(vector<uint8_t>& data, bool store_date_time)
{
	data.clear();

	append(data, (uint32_t)col.size());

	for (auto& sample : col)
	{
		append(data, sample.first);
		append(data, (uint32_t)sample.second.size());

		for (auto& contig : sample.second)
		{
			append(data, contig.first);
			append(data, (uint32_t)contig.second.size());

			int32_t prev_group_id = 0;
			int32_t prev_in_group_id = 0;
			int32_t prev_raw_length = 0;

			for (auto& seg : contig.second)
			{
				uint32_t e_group_id = (uint32_t)zigzag_encode((int32_t)seg.group_id - prev_group_id);
				uint32_t e_in_group_id = (uint32_t)zigzag_encode((int32_t)seg.in_group_id - prev_in_group_id);
				uint32_t e_raw_length = (uint32_t)zigzag_encode((int32_t)seg.raw_length - prev_raw_length);

				append(data, e_group_id);
				append(data, e_in_group_id);
				append(data, e_raw_length);
				append(data, (uint32_t)seg.is_rev_comp);

				prev_group_id = seg.group_id;
				prev_in_group_id = seg.in_group_id;
				prev_raw_length = seg.raw_length;
			}
		}
	}

	append(data, (uint32_t) cmd_lines.size());

	for (auto& cmd : cmd_lines)
	{
		append(data, cmd.first);
		if (store_date_time)
			append(data, cmd.second);
		else
			append(data, "");
	}
}

// *******************************************************************************************
bool CCollection::deserialize(vector<uint8_t>& data)
{
	auto p = data.begin();

	col.clear();

	uint32_t no_samples;

	read(p, no_samples);

	for (uint32_t i = 0; i < no_samples; ++i)
	{
		string sample_name;
		read(p, sample_name);

		uint32_t no_contigs;
		read(p, no_contigs);

		col[sample_name].resize(no_contigs);

		for (uint32_t j = 0; j < no_contigs; ++j)
		{
			string contig_name;
			read(p, contig_name);

			string short_contig_name = extract_contig_name(contig_name);

			contig_ids[make_pair(sample_name, short_contig_name)] = j;

			mm_contig2sample.emplace(short_contig_name, sample_name);

			uint32_t no_seg;
			read(p, no_seg);

			auto& q = col[sample_name][j];
			q.first = contig_name;

			uint32_t e_group_id;
			uint32_t e_in_group_id;
			uint32_t e_raw_length;
			uint32_t e_orientation;

			int32_t prev_group_id = 0;
			int32_t prev_in_group_id = 0;
			int32_t prev_raw_length = 0;

			for (uint32_t k = 0; k < no_seg; ++k)
			{
				read(p, e_group_id);
				read(p, e_in_group_id);
				read(p, e_raw_length);
				read(p, e_orientation);

				uint32_t c_group_id = (uint32_t)((int32_t)prev_group_id + zigzag_decode(e_group_id));
				uint32_t c_in_group_id = (uint32_t)((int32_t)prev_in_group_id + zigzag_decode(e_in_group_id));
				uint32_t c_raw_length = (uint32_t)((int32_t)prev_raw_length + zigzag_decode(e_raw_length));

				q.second.emplace_back(c_group_id, c_in_group_id, (bool)e_orientation, c_raw_length);

				prev_group_id = c_group_id;
				prev_in_group_id = c_in_group_id;
				prev_raw_length = c_raw_length;
			}
		}
	}

	uint32_t no_cmds;

	read(p, no_cmds);
	cmd_lines.clear();
	cmd_lines.resize(no_cmds);

	for (uint32_t i = 0; i < no_cmds; ++i)
	{
		read(p, cmd_lines[i].first);
		read(p, cmd_lines[i].second);
	}

	return true;
}

// EOF
