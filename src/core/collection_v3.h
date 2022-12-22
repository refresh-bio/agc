#ifndef _COLLECTION_V3_H
#define _COLLECTION_V3_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include "collection.h"
#include "archive.h"

class CCollection_V3 : public CCollection
{
	struct contig_desc_t {
		string name;
		vector<segment_desc_t> segments;

		contig_desc_t() : name("") {};
		
		contig_desc_t(const contig_desc_t& x) {
			name = x.name;
			segments = x.segments;
		}

		contig_desc_t(contig_desc_t&& x) noexcept {
			name = move(x.name);
			segments = move(x.segments);
		}

		contig_desc_t(const string& _name) : name(_name) {};

		contig_desc_t& operator=(const contig_desc_t& x)
		{
			name = x.name;
			segments = x.segments;

			return *this;
		}
	};

	struct sample_desc_t {
		string name;
		vector<contig_desc_t> contigs;

		sample_desc_t() : name("") {};
		sample_desc_t(const string& _name, const vector<contig_desc_t>& _contigs) : name(_name), contigs(_contigs) {};
		sample_desc_t(const sample_desc_t&x) {
			name = x.name;
			contigs = x.contigs;
		}

		sample_desc_t(sample_desc_t&&x) noexcept {
			name = move(x.name);
			contigs = move(x.contigs);
		}
		sample_desc_t(const string& _name) : name(_name) {};

		sample_desc_t& operator=(const sample_desc_t& x) {
			name = x.name;
			contigs = x.contigs;

			return *this;
		}
	};

	ZSTD_CCtx* zstd_cctx_samples = nullptr;
	ZSTD_CCtx* zstd_cctx_contigs = nullptr;
	array<ZSTD_CCtx*, 5> zstd_cctx_details = { nullptr, nullptr, nullptr, nullptr, nullptr };
	
	ZSTD_DCtx* zstd_dctx_samples = nullptr;
	ZSTD_DCtx* zstd_dctx_contigs = nullptr;
	array<ZSTD_DCtx*, 5> zstd_dctx_details = { nullptr, nullptr, nullptr, nullptr, nullptr };

	unordered_map<string, uint32_t> sample_ids;
	vector<sample_desc_t> sample_desc;

	int unpacked_contig_data_batch_id = -1;

	uint32_t no_threads;

	size_t batch_size;
	uint32_t segment_size;
	uint32_t kmer_length;
	size_t no_samples_in_last_batch;
	string prev_sample_name;
	
	string placing_sample_name;
	uint32_t placing_sample_id;

	int collection_samples_id;
	int collection_contig_id;
	int collection_details_id;

	shared_ptr<CArchive> in_archive;
	shared_ptr<CArchive> out_archive;
	vector<int> v_in_group_ids;

	void store_batch_sample_names();
	void store_batch_contig_names(uint32_t id_from, uint32_t id_to);
	void store_batch_contig_details(uint32_t id_from, uint32_t id_to);

	void load_batch_sample_names();
	void load_batch_contig_names(size_t id_batch);
	void load_batch_contig_details(size_t id_batch);
	void clear_batch_contig(size_t id_batch);

	void serialize_sample_names(vector<uint8_t> &v_data);
	void serialize_contig_names(vector<uint8_t>& v_data, uint32_t id_from, uint32_t id_to);
	void serialize_contig_details(array<vector<uint8_t>, 5>& v_data, uint32_t id_from, uint32_t id_to);

	void deserialize_sample_names(vector<uint8_t>& v_data);
	void deserialize_contig_names(vector<uint8_t>& v_data, size_t i_sample);
	void deserialize_contig_details(array<vector<uint8_t>, 5>& v_data, size_t i_sample);

	bool prepare_for_compression();
	bool prepare_for_appending_copy();
	bool prepare_for_decompression();

	void zstd_compress(ZSTD_CCtx*& cctx, vector<uint8_t>& v_input, vector<uint8_t>& v_output, int level);
	void zstd_decompress(ZSTD_DCtx*& dctx, vector<uint8_t>& v_input, vector<uint8_t>& v_output, size_t raw_size);

	// Just check
	int get_in_group_id(int pos)
	{
		if ((size_t) pos >= v_in_group_ids.size())
			return -1;
		return v_in_group_ids[pos];
	}
	
	// Check but resize first if necessary
	int read_in_group_id(int pos)
	{
		if ((size_t) pos >= v_in_group_ids.size())
			v_in_group_ids.resize((int)(pos * 1.2), -1);

		return v_in_group_ids[pos];
	}

	void set_in_group_id(int pos, int val)
	{
		if ((size_t) pos >= v_in_group_ids.size())
			v_in_group_ids.resize((int) (pos * 1.2) + 1, -1);

		v_in_group_ids[pos] = val;
	}

	void clear_in_group_ids()
	{
		v_in_group_ids.clear();
	}

	void determine_collection_samples_id()
	{
		if (collection_samples_id >= 0)
			return;

		if(out_archive != nullptr)
			collection_samples_id = out_archive->GetStreamId("collection-samples");
		else
			collection_samples_id = in_archive->GetStreamId("collection-samples");
	}

	void determine_collection_contig_id()
	{
		if (collection_contig_id >= 0)
			return;

		if(out_archive != nullptr)
			collection_contig_id = out_archive->GetStreamId("collection-contigs");
		else
			collection_contig_id = in_archive->GetStreamId("collection-contigs");
	}

	void determnine_collection_details_id()
	{
		if (collection_details_id >= 0)
			return;

		if(out_archive != nullptr)
			collection_details_id = out_archive->GetStreamId("collection-details");
		else
			collection_details_id = in_archive->GetStreamId("collection-details");
	}

	vector<string> split_string(const string& s);
	string encode_split(vector<string>& prev_split, vector<string>& curr_split);
	string decode_split(vector<string>& prev_split, vector<string>& curr_split);

public:
	CCollection_V3() : CCollection() {
		batch_size = 1ull << 20;

		no_threads = 1;

		collection_samples_id = -1;
		collection_contig_id = -1;
		collection_details_id = -1;

		placing_sample_id = ~0u;
		placing_sample_name = "";
	}
	virtual ~CCollection_V3() {
		if (zstd_cctx_samples)	ZSTD_freeCCtx(zstd_cctx_samples);
		if (zstd_cctx_contigs)	ZSTD_freeCCtx(zstd_cctx_contigs);
		for(auto &x : zstd_cctx_details)
			if (x)	ZSTD_freeCCtx(x);

		if (zstd_dctx_samples)	ZSTD_freeDCtx(zstd_dctx_samples);
		if (zstd_dctx_contigs)	ZSTD_freeDCtx(zstd_dctx_contigs);
		for(auto &x : zstd_dctx_details)
			if (x)	ZSTD_freeDCtx(x);
	};

	bool set_archives(shared_ptr<CArchive> _in_archive, shared_ptr<CArchive> _out_archive,
		uint32_t _no_threads, size_t _batch_size, uint32_t _segment_size, uint32_t _kmer_length);

	void complete_serialization();

	bool prepare_for_appending_load_last_batch();

	virtual bool register_sample_contig(const string& sample_name, const string& contig_name);
	
	void reset_prev_sample_name();
	virtual void add_segment_placed(const string &sample_name, const string& contig_name, const uint32_t place, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length);
	void add_segments_placed(vector<segments_to_place_t>& segments_to_place);
	virtual bool get_reference_name(string& reference_name);
	virtual bool get_samples_list(vector<string>& v_samples);
	virtual bool get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names);
	virtual bool get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc_);
	virtual bool get_contig_desc(const string& sample_name, string& contig_name, vector<segment_desc_t>& contig_desc);
	virtual bool is_contig_desc(const string& sample_name, const string& contig_name);
	virtual vector<string> get_samples_for_contig(const string& contig_name);
	virtual size_t get_no_samples();
	virtual int32_t get_no_contigs(const string& sample_name);

	void store_contig_batch(uint32_t id_from, uint32_t id_to);
};

// EOF
#endif
