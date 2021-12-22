#ifndef _LZ_DIFF_H
#define _LZ_DIFF_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include <string>
#include <vector>
#include "../core/utils.h"

using namespace std;

#define USE_SPARSE_HT

class CLZDiff
{
	const uint32_t empty_key32 = ~0u;
	const uint32_t empty_key16 = (uint16_t) ~0u;
	const double max_load_factor = 0.7;
	const uint32_t max_no_tries = 64;
	const uint8_t invalid_symbol = 31;
	const uint8_t N_code = 4;
	const uint8_t N_run_starter_code = 30;
	const uint32_t min_Nrun_len = 4;

#ifdef USE_SPARSE_HT
	const uint32_t hashing_step = 4;
#else
	const uint32_t hashing_step = 1;
#endif

	contig_t reference;
	vector<uint32_t> ht32;
	vector<uint16_t> ht16;
	uint64_t ht_size;
	uint64_t ht_mask;
	uint32_t key_len;
	uint32_t min_match_len;
	bool short_ht_ver;
	bool index_ready;

	void make_index16();
	void make_index32();

	uint64_t get_code(const uint8_t* s) const;
	uint32_t get_Nrun_len(const uint8_t* s, const uint32_t max_len) const;

	void encode_literal(const uint8_t c, contig_t& encoded);
	void encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded);
	void encode_Nrun(const uint32_t len, contig_t& encoded);

	uint32_t coding_cost_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos) const;
	uint32_t coding_cost_Nrun(const uint32_t len) const;
	uint32_t int_len(const uint32_t x) const
	{
		if (x < 10)	return 1;
		if (x < 100)	return 2;
		if (x < 1000)	return 3;
		if (x < 10000)	return 4;
		if (x < 100000)	return 5;
		if (x < 1000000)	return 6;
		if (x < 10000000)	return 7;
		if (x < 100000000)	return 8;
		if (x < 1000000000)	return 9;
		return 10;
	}

	bool is_literal(const contig_t::const_iterator& p) const;
	bool is_Nrun(const contig_t::const_iterator& p) const;
	void decode_literal(contig_t::const_iterator& p, uint8_t &c);
	void decode_match(contig_t::const_iterator& p, uint32_t &ref_pos, uint32_t &len, uint32_t &pred_pos);
	void decode_Nrun(contig_t::const_iterator& p, uint32_t& len);

	bool find_best_match16(uint32_t ht_pos, const uint8_t *s, const uint32_t max_len, const uint32_t no_prev_literals,
		uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd);
	bool find_best_match32(uint32_t ht_pos, const uint8_t *s, const uint32_t max_len, const uint32_t no_prev_literals,
		uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd);
	void append_int(contig_t& text, int64_t x);
	void read_int(contig_t::const_iterator& p, int64_t &x);

	void prepare_gen(const contig_t& _reference);
	void prepare_index();

public:
	CLZDiff(const uint32_t _min_match_len = 18);
	~CLZDiff();

	bool SetMinMatchLen(const uint32_t _min_match_len = 18);
	void Prepare(const contig_t& _reference);
	void Encode(const contig_t& text, contig_t&encoded);
	void GetReference(contig_t& s);
	void Decode(const contig_t &reference, const contig_t& encoded, contig_t& decoded);
	void GetCodingCostVector(const contig_t& text, vector<uint32_t> &v_costs, const bool prefix_costs);
};

// EOF
#endif