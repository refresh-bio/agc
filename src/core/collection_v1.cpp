// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include "../core/collection_v1.h"

// *******************************************************************************************
void CCollection_V1::serialize(vector<uint8_t>& data, bool store_date_time)
{
	data.clear();

	append(data, (uint32_t)col.size());

	auto col_order = get_sample_original_order();

	for (auto& sample_name : col_order)
	{
		auto& sample = col[sample_name];

		append(data, sample_name);
		append(data, (uint32_t)sample.size());

		for (auto& contig : sample)
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

	append(data, (uint32_t)cmd_lines.size());

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
bool CCollection_V1::deserialize(vector<uint8_t>& data)
{
	uint8_t* p = data.data();

	col.clear();

	uint32_t no_samples;
	string sample_name;
	string contig_name;

	read(p, no_samples);

	v_sample_name.reserve(no_samples);

	for (uint32_t i = 0; i < no_samples; ++i)
	{
		read(p, sample_name);

		v_sample_name.emplace_back(sample_name);

		uint32_t no_contigs;
		read(p, no_contigs);

		col[sample_name].resize(no_contigs);

		uint32_t sample_id = (uint32_t)sample_ids.size();
		sample_ids[sample_name] = sample_id;

		for (uint32_t j = 0; j < no_contigs; ++j)
		{
			read(p, contig_name);
			uint32_t no_seg;
			read(p, no_seg);

			string short_contig_name = extract_contig_name(contig_name);

			contig_ids_no_seg[make_pair(sample_name, short_contig_name)] = make_pair(j, no_seg);

			mm_contig2sample.emplace(short_contig_name, sample_name);

			auto& q = col[sample_name][j];
			q.first = contig_name;
			auto& q_second = q.second;

			uint32_t e_group_id;
			uint32_t e_in_group_id;
			uint32_t e_raw_length;
			uint32_t e_orientation;

			int32_t prev_group_id = 0;
			int32_t prev_in_group_id = 0;
			int32_t prev_raw_length = 0;

			q_second.reserve(no_seg);

			for (uint32_t k = 0; k < no_seg; ++k)
			{
				read(p, e_group_id);
				read(p, e_in_group_id);
				read(p, e_raw_length);
				read(p, e_orientation);

				uint32_t c_group_id = (uint32_t)((int32_t)prev_group_id + zigzag_decode(e_group_id));
				uint32_t c_in_group_id = (uint32_t)((int32_t)prev_in_group_id + zigzag_decode(e_in_group_id));
				uint32_t c_raw_length = (uint32_t)((int32_t)prev_raw_length + zigzag_decode(e_raw_length));

				q_second.emplace_back(c_group_id, c_in_group_id, (bool)e_orientation, c_raw_length);

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

	maps_built = true;

	return true;
}

// *******************************************************************************************
vector<string> CCollection_V1::get_sample_original_order()
{
	vector<string> vec(sample_ids.size());

	for (auto& x : sample_ids)
		vec[x.second] = x.first;

	return vec;
}

// *******************************************************************************************
bool CCollection_V1::register_sample_contig(const string& sample_name, const string& contig_name)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);
	string stored_sample_name = sample_name;

	if (sample_name.empty())
		stored_sample_name = short_contig_name;

	auto q = sample_ids.find(stored_sample_name);
	if (q == sample_ids.end())
	{
		uint32_t sample_id = (uint32_t)sample_ids.size();
		sample_ids[stored_sample_name] = sample_id;
	}

	auto p = contig_ids_no_seg.find(make_pair(stored_sample_name, short_contig_name));
	if (p == contig_ids_no_seg.end())
	{
		uint32_t contig_id = (uint32_t)col[stored_sample_name].size();
		contig_ids_no_seg[make_pair(stored_sample_name, short_contig_name)] = make_pair(contig_id, 0);
		col[stored_sample_name].emplace_back(contig_name, vector<segment_desc_t>());

		mm_contig2sample.emplace(short_contig_name, stored_sample_name);

		return true;
	}
	else
		return false;		// The pair sample_name:contig_name is not unique
}

// *******************************************************************************************
vector<segment_desc_t>& CCollection_V1::add_segment_basic(const string& sample_name, const string& contig_name,
	const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length)
{
	uint32_t contig_id;

	string short_contig_name = extract_contig_name(contig_name);
	string stored_sample_name = sample_name;

	if (sample_name.empty())
		stored_sample_name = short_contig_name;

	auto p_contig_ids = contig_ids_no_seg.find(make_pair(stored_sample_name, short_contig_name));
	if (p_contig_ids != contig_ids_no_seg.end())
		contig_id = p_contig_ids->second.first;
	else
	{
		contig_id = (uint32_t)col[stored_sample_name].size();
		contig_ids_no_seg[make_pair(stored_sample_name, short_contig_name)] = make_pair(contig_id, 0);
		col[stored_sample_name].emplace_back(contig_name, vector<segment_desc_t>());

		mm_contig2sample.emplace(short_contig_name, stored_sample_name);
	}

	return col[stored_sample_name][contig_id].second;
}

// *******************************************************************************************
void CCollection_V1::add_segment_placed(const string& sample_name, const string& contig_name, const uint32_t place,
	const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length)
{
	lock_guard<mutex> lck(mtx);

	auto& vec = add_segment_basic(sample_name, contig_name, group_id, in_group_id, is_rev_comp, raw_length);

	if (place >= vec.size())
		vec.resize(place + 1, segment_desc_t(555555555, 555555555, true, 555555555));

	vec[place] = segment_desc_t(group_id, in_group_id, is_rev_comp, raw_length);
}

// *******************************************************************************************
bool CCollection_V1::get_reference_name(string& reference_name)
{
	lock_guard<mutex> lck(mtx);

	if (v_sample_name.empty())
		return false;

	reference_name = v_sample_name.front();

	return true;
}

// *******************************************************************************************
bool CCollection_V1::get_samples_list(vector<string>& v_samples)
{
	lock_guard<mutex> lck(mtx);

	v_samples.clear();

	for (auto& p : col)
		v_samples.emplace_back(p.first);

	return true;
}

// *******************************************************************************************
bool CCollection_V1::get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names)
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
bool CCollection_V1::get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc)
{
	lock_guard<mutex> lck(mtx);

	auto p = col.find(sample_name);
	if (p == col.end())
		return false;

	if (p->second.front().second.empty())
		decompress_sample_details(sample_ids[sample_name]);

	sample_desc = p->second;

	return true;
}

// *******************************************************************************************
bool CCollection_V1::get_contig_desc(const string& sample_name, string& contig_name, vector<segment_desc_t>& contig_desc)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);

	uint32_t contig_id;

	if (maps_built)
	{
		auto p_contig_id = contig_ids_no_seg.find(make_pair(sample_name, short_contig_name));
		if (p_contig_id == contig_ids_no_seg.end())
			return false;

		contig_id = p_contig_id->second.first;
	}
	else
	{
		auto p_contig_info = find_if(v_contig_info.begin(), v_contig_info.end(), [sample_name, short_contig_name](auto& x) {
			return sample_name == x.sample_name && short_contig_name == x.contig_name; });

		if (p_contig_info == v_contig_info.end())
			return false;

		contig_id = p_contig_info->id;
	}

	if (col[sample_name][contig_id].second.empty())
		decompress_sample_details(sample_ids[sample_name]);

	contig_desc = col[sample_name][contig_id].second;

	// Reconstruct full contig name
	contig_name = col[sample_name][contig_id].first;

	return true;
}

// *******************************************************************************************
bool CCollection_V1::is_contig_desc(const string& sample_name, const string& contig_name)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);

	if (maps_built)
		return contig_ids_no_seg.find(make_pair(sample_name, short_contig_name)) != contig_ids_no_seg.end();
	else
	{
		for (auto p = v_contig_info.begin(); p != v_contig_info.end(); ++p)
			if (p->contig_name == contig_name && p->sample_name == sample_name)
				return true;

		return false;
	}
}

// *******************************************************************************************
vector<string> CCollection_V1::get_samples_for_contig(const string& contig_name)
{
	vector<string> vs;

	if (maps_built)
	{
		auto pq = mm_contig2sample.equal_range(extract_contig_name(contig_name));

		for (auto p = pq.first; p != pq.second; ++p)
			vs.emplace_back(p->second);
	}
	else
	{
		for (auto p = v_contig_info.begin(); p != v_contig_info.end(); ++p)
			if (p->contig_name == contig_name)
				vs.emplace_back(p->sample_name);
	}

	return vs;
}

// *******************************************************************************************
size_t CCollection_V1::get_no_samples()
{
	lock_guard<mutex> lck(mtx);

	return col.size();
}

// *******************************************************************************************
int32_t CCollection_V1::get_no_contigs(const string& sample_name)
{
	lock_guard<mutex> lck(mtx);

	auto p = col.find(sample_name);
	if (p == col.end())
		return -1;

	return static_cast<int32_t>(p->second.size());
}

// *******************************************************************************************
void CCollection_V1::decompress_sample_details(uint32_t i_sample)
{
	if (!zstd_dctx)
		zstd_dctx = ZSTD_createDCtx();

	uint32_t i_part = i_sample / details_batch_size;

	vector<uint8_t> v_tmp;

	v_tmp.resize(v_zstd_batches[i_part].second);
	ZSTD_decompressDCtx(zstd_dctx, v_tmp.data(), v_tmp.size(), v_zstd_batches[i_part].first.data(), v_zstd_batches[i_part].first.size());

	v_zstd_batches[i_part].first.clear();
	v_zstd_batches[i_part].first.shrink_to_fit();
	v_zstd_batches[i_part].second = 0;

	uint8_t* p = v_tmp.data();

	vector<col_t::iterator> v_p_sample;
	v_p_sample.reserve(details_batch_size);

	for (uint32_t i = i_part * details_batch_size; i < (i_part + 1) * details_batch_size && i < col.size(); ++i)
		v_p_sample.push_back(col.find(v_sample_name[i]));

	if (maps_built)
	{
		for (auto p_sam = v_p_sample.begin(); p_sam != v_p_sample.end(); ++p_sam)
			for (auto p_ctg = (*p_sam)->second.begin(); p_ctg != (*p_sam)->second.end(); ++p_ctg)
			{
				uint32_t no_seg = contig_ids_no_seg[make_pair((*p_sam)->first, extract_contig_name(p_ctg->first))].second;
				p_ctg->second.resize(no_seg);
			}
	}
	else
	{
		auto p_contig_info = find_if(v_contig_info.begin(), v_contig_info.end(), [&v_p_sample](auto& x) {
			return v_p_sample.front()->first == x.sample_name; });

		for (auto p_sam = v_p_sample.begin(); p_sam != v_p_sample.end(); ++p_sam)
			for (auto p_ctg = (*p_sam)->second.begin(); p_ctg != (*p_sam)->second.end(); ++p_ctg, ++p_contig_info)
			{
				uint32_t no_seg = p_contig_info->no_seg;

				p_ctg->second.resize(no_seg);
			}
	}

	for (auto p_sam = v_p_sample.begin(); p_sam != v_p_sample.end(); ++p_sam)
		for (auto p_ctg = (*p_sam)->second.begin(); p_ctg != (*p_sam)->second.end(); ++p_ctg)
			deserialize_contig_details_group_id(p, p_ctg->second);

	for (auto p_sam = v_p_sample.begin(); p_sam != v_p_sample.end(); ++p_sam)
		for (auto p_ctg = (*p_sam)->second.begin(); p_ctg != (*p_sam)->second.end(); ++p_ctg)
			deserialize_contig_details_in_group_id(p, p_ctg->second);

	for (auto p_sam = v_p_sample.begin(); p_sam != v_p_sample.end(); ++p_sam)
		for (auto p_ctg = (*p_sam)->second.begin(); p_ctg != (*p_sam)->second.end(); ++p_ctg)
			deserialize_contig_details_raw_length(p, p_ctg->second);

	for (auto p_sam = v_p_sample.begin(); p_sam != v_p_sample.end(); ++p_sam)
		for (auto p_ctg = (*p_sam)->second.begin(); p_ctg != (*p_sam)->second.end(); ++p_ctg)
			deserialize_contig_details_orientation(p, p_ctg->second);
}

// *******************************************************************************************
void CCollection_V1::deserialize_contig_details_group_id(uint8_t*& p, vector<segment_desc_t>& contig_segments)
{
	uint32_t e_group_id;
	int32_t prev_group_id = 0;

	for (uint32_t k = 0; k < contig_segments.size(); ++k)
	{
		read(p, e_group_id);
		uint32_t c_group_id = (uint32_t)zigzag_decode(e_group_id, prev_group_id);
		contig_segments[k].group_id = c_group_id;
		prev_group_id = c_group_id;
	}
}

// *******************************************************************************************
void CCollection_V1::deserialize_contig_details_in_group_id(uint8_t*& p, vector<segment_desc_t>& contig_segments)
{
	uint32_t e_in_group_id;
	int32_t prev_in_group_id = 0;

	for (uint32_t k = 0; k < contig_segments.size(); ++k)
	{
		read(p, e_in_group_id);
		uint32_t c_in_group_id = (uint32_t)zigzag_decode(e_in_group_id, prev_in_group_id);
		contig_segments[k].in_group_id = c_in_group_id;
		prev_in_group_id = c_in_group_id;
	}
}

// *******************************************************************************************
void CCollection_V1::deserialize_contig_details_raw_length(uint8_t*& p, vector<segment_desc_t>& contig_segments)
{
	uint32_t e_raw_length;
	int32_t prev_raw_length = 0;

	for (uint32_t k = 0; k < contig_segments.size(); ++k)
	{
		read(p, e_raw_length);
		uint32_t c_raw_length = (uint32_t)zigzag_decode(e_raw_length, prev_raw_length);
		contig_segments[k].raw_length = c_raw_length;
		prev_raw_length = c_raw_length;
	}
}

// *******************************************************************************************
void CCollection_V1::deserialize_contig_details_orientation(uint8_t*& p, vector<segment_desc_t>& contig_segments)
{
	uint32_t e_orientation;

	for (uint32_t k = 0; k < contig_segments.size(); ++k)
	{
		read(p, e_orientation);
		contig_segments[k].is_rev_comp = (bool)e_orientation;
	}
}

// EOF
