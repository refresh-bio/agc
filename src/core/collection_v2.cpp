// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.1
// Date   : 2024-03-12
// *******************************************************************************************

#include "../core/collection_v2.h"

// *******************************************************************************************
void CCollection_V2::serialize(vector<uint8_t>& data_main, vector<vector<uint8_t>>& data_details, bool store_date_time, uint32_t _details_batch_size)
{
	data_main.clear();
	data_details.clear();

	details_batch_size = _details_batch_size;

	auto col_order = get_sample_original_order();

	append(data_main, details_batch_size);
	append(data_main, (uint32_t)col.size());

	data_details.resize((col_order.size() + (details_batch_size - 1)) / details_batch_size);

	array<vector<uint8_t>, 4> v_det;
	size_t no_samples_in_batch = 0;
	uint32_t i_details_part = 0;

	for (auto& sample_name : col_order)
	{
		auto& sample = col[sample_name];

		append(data_main, sample_name);
		append(data_main, (uint32_t)sample.size());

		if (no_samples_in_batch == details_batch_size)
		{
			for (auto& v : v_det)
			{
				data_details[i_details_part].insert(data_details[i_details_part].end(), v.begin(), v.end());
				v.clear();
			}

			no_samples_in_batch = 0;
			++i_details_part;
		}

		for (auto& contig : sample)
		{
			append(data_main, contig.first);
			append(data_main, (uint32_t)contig.second.size());

			int32_t prev_group_id = 0;
			int32_t prev_in_group_id = 0;
			int32_t prev_raw_length = 0;

			for (auto& seg : contig.second)
			{
				uint32_t e_group_id = (uint32_t)zigzag_encode(seg.group_id, prev_group_id);
				uint32_t e_in_group_id = (uint32_t)zigzag_encode(seg.in_group_id, prev_in_group_id);
				uint32_t e_raw_length = (uint32_t)zigzag_encode(seg.raw_length, prev_raw_length);

				append(v_det[0], e_group_id);
				append(v_det[1], e_in_group_id);
				append(v_det[2], e_raw_length);
				append(v_det[3], (uint32_t)seg.is_rev_comp);

				prev_group_id = seg.group_id;
				prev_in_group_id = seg.in_group_id;
				prev_raw_length = seg.raw_length;
			}
		}

		++no_samples_in_batch;
	}

	for (auto& v : v_det)
		data_details[i_details_part].insert(data_details[i_details_part].end(), v.begin(), v.end());

	append(data_main, (uint32_t)cmd_lines.size());

	for (auto& cmd : cmd_lines)
	{
		append(data_main, cmd.first);
		if (store_date_time)
			append(data_main, cmd.second);
		else
			append(data_main, "");
	}
}

// *******************************************************************************************
bool CCollection_V2::deserialize_main(vector<uint8_t>& data_main, bool create_maps)
{
	uint8_t* p = data_main.data();

	col.clear();

	uint32_t no_samples;
	string sample_name;
	string contig_name;

	read(p, details_batch_size);
	read(p, no_samples);

	v_sample_name.reserve(no_samples);

	for (uint32_t i = 0; i < no_samples; ++i)
	{
		read(p, sample_name);

		v_sample_name.emplace_back(sample_name);

		uint32_t no_contigs;
		read(p, no_contigs);

		col[sample_name].resize(no_contigs);
		auto& col_sample = col[sample_name];

		uint32_t sample_id = (uint32_t)sample_ids.size();
		sample_ids[sample_name] = sample_id;

		for (uint32_t j = 0; j < no_contigs; ++j)
		{
			read(p, contig_name);

			uint32_t no_seg;
			read(p, no_seg);

			string short_contig_name = extract_contig_name(contig_name);

			if (create_maps)
			{
				contig_ids_no_seg[make_pair(sample_name, short_contig_name)] = make_pair(j, no_seg);
				mm_contig2sample.emplace(short_contig_name, sample_name);
			}
			else
				v_contig_info.emplace_back(sample_name, short_contig_name, j, no_seg);

			col_sample[j].first = contig_name;
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

	maps_built = create_maps;

	return true;
}

// *******************************************************************************************
bool CCollection_V2::deserialize_details(vector<uint8_t>& zstd_data_details, size_t raw_size, bool deserialize_details)
{
	v_zstd_batches.emplace_back(move(zstd_data_details), raw_size);

	if (deserialize_details)
		decompress_sample_details((v_zstd_batches.size() - 1) * details_batch_size);

	return true;
}

// EOF
