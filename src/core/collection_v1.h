#ifndef _COLLECTION_V1_H
#define _COLLECTION_V1_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.1
// Date   : 2024-03-12
// *******************************************************************************************

#include "collection.h"

class CCollection_V1 : public CCollection
{
protected:
	uint32_t details_batch_size = 1;

	typedef map<string, vector<pair<string, vector<segment_desc_t>>>> col_t;

	vector<pair<vector<uint8_t>, size_t>> v_zstd_batches;

	col_t col;
	map<pair<string, string>, pair<uint32_t, uint32_t>> contig_ids_no_seg;
	multimap<string, string> mm_contig2sample;
	map<string, uint32_t> sample_ids;
	vector<string> v_sample_name;
	vector<contig_info_t> v_contig_info;
	bool maps_built = false;

	vector<string> get_sample_original_order();

	void decompress_sample_details(uint32_t i_sample);
	void deserialize_contig_details_group_id(uint8_t*& p, vector<segment_desc_t>& contig_segments);
	void deserialize_contig_details_in_group_id(uint8_t*& p, vector<segment_desc_t>& contig_segments);
	void deserialize_contig_details_raw_length(uint8_t*& p, vector<segment_desc_t>& contig_segments);
	void deserialize_contig_details_orientation(uint8_t*& p, vector<segment_desc_t>& contig_segments);

	vector<segment_desc_t>& add_segment_basic(const string& sample_name, const string& contig_name, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length);

public:
	CCollection_V1() : CCollection() {}

	virtual ~CCollection_V1() {};


	void serialize(vector<uint8_t>& data, bool store_date_time);
	bool deserialize(vector<uint8_t>& data);

	virtual bool register_sample_contig(const string& sample_name, const string& contig_name);
	virtual void add_segment_placed(const string& sample_name, const string& contig_name, const uint32_t place, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length);
	virtual bool get_reference_name(string& reference_name);
	virtual bool get_samples_list(vector<string>& v_samples);
	virtual bool get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names);
	virtual bool get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc);
	virtual bool get_contig_desc(const string& sample_name, string& contig_name, vector<segment_desc_t>& contig_desc);

	virtual bool is_contig_desc(const string& sample_name, const string& contig_name);
	virtual vector<string> get_samples_for_contig(const string& contig_name);

	virtual size_t get_no_samples();
	virtual int32_t get_no_contigs(const string& sample_name);
};

// EOF
#endif