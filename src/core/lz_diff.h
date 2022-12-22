#ifndef _LZ_DIFF_H
#define _LZ_DIFF_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include <string>
#include <vector>
#include <array>
#include "../core/utils.h"

using namespace std;

#define USE_SPARSE_HT

// *******************************************************************************************
class CLZDiffBase
{
protected:
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
	uint64_t key_mask;
	uint32_t min_match_len;
	bool short_ht_ver;
	bool index_ready;

	void make_index16();
	void make_index32();

	uint64_t get_code(const uint8_t* s) const;
	uint64_t get_code_skip1(uint64_t code, const uint8_t* s) const;
	uint32_t get_Nrun_len(const uint8_t* s, const uint32_t max_len) const;

	void encode_literal(const uint8_t c, contig_t& encoded);
	void encode_literal_diff(const uint8_t c, const uint8_t r, contig_t& encoded);
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
	void decode_Nrun(contig_t::const_iterator& p, uint32_t& len);

	bool find_best_match16(uint32_t ht_pos, const uint8_t *s, const uint32_t max_len, const uint32_t no_prev_literals,
		uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd);
	bool find_best_match32(uint32_t ht_pos, const uint8_t *s, const uint32_t max_len, const uint32_t no_prev_literals,
		uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd);
	inline void append_int(contig_t& text, int64_t x)
	{
		if (x == 0)
		{
			text.emplace_back('0');
			return;
		}

		if (x < 0)
		{
			text.push_back('-');
			x = -x;
		}

		char tmp[16];
		char* p = tmp + 16;

		for (; x; x /= 10)
		{
			*--p = (uint8_t)('0' + (x % 10));
			x /= 10;
			if (!x)
				break;
			*--p = (uint8_t)('0' + (x % 10));
		}

		int i = text.size();
		text.resize(text.size() + (tmp + 16 - p));

		for (; p != tmp + 16; ++p, ++i)
			text[i] = *p;

		text.insert(text.end(), p, tmp + 16);
	}

	void read_int(contig_t::const_iterator& p, int64_t &x)
	{
		bool is_neg = false;
		x = 0;

		if (*p == '-')
		{
			is_neg = true;
			++p;
		}

		while (*p >= '0' && *p <= '9')
			x = x * 10 + (int64_t)(*p++ - '0');

		if (is_neg)
			x = -x;
	}

	uint32_t compare_fwd(uint8_t* p, uint8_t* q, uint32_t max_len)
	{
		uint32_t len = 0;

		switch (max_len % 4)
		{
		case 3:
			if (*p++ != *q++)
				return len;
			++len;
		case 2:
			if (*p++ != *q++)
				return len;
			++len;
		case 1:
			if (*p++ != *q++)
				return len;
			++len;
		}

		for (; len < max_len; len += 4)
		{
			if (*p++ != *q++)
				return len;
			if (*p++ != *q++)
				return len+1;
			if (*p++ != *q++)
				return len+2;
			if (*p++ != *q++)
				return len+3;
		}

		return len;
	}

	void prepare_gen(const contig_t& _reference);
	void prepare_index();

public:
	CLZDiffBase(const uint32_t _min_match_len = 18);
	virtual ~CLZDiffBase();

	bool SetMinMatchLen(const uint32_t _min_match_len = 18);
	void Prepare(const contig_t& _reference);

	virtual void Encode(const contig_t& text, contig_t&encoded) = 0;
	virtual void Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded) = 0;

	virtual size_t Estimate(const contig_t& text, uint32_t bound = 0) = 0;

	void GetReference(contig_t& s);
	void GetCodingCostVector(const contig_t& text, vector<uint32_t> &v_costs, const bool prefix_costs);
};

// *******************************************************************************************
class CLZDiff_V1 : public CLZDiffBase
{
	void encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded);
	void decode_match(contig_t::const_iterator& p, uint32_t& ref_pos, uint32_t& len, uint32_t& pred_pos);


public:
	CLZDiff_V1(const uint32_t _min_match_len = 18) : CLZDiffBase(_min_match_len)
	{}

	virtual ~CLZDiff_V1() {};

	virtual void Encode(const contig_t& text, contig_t& encoded);
	virtual void Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded);

	virtual size_t Estimate(const contig_t& text, uint32_t bound = 0);
};

// *******************************************************************************************
class CLZDiff_V2 : public CLZDiffBase
{
	void encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded);
	void decode_match(contig_t::const_iterator& p, uint32_t& ref_pos, uint32_t& len, uint32_t& pred_pos);

	uint32_t int_len(int x)
	{
		if (x >= 0)
			return uint_len((uint32_t)x);
		else
			return 1 + uint_len((uint32_t)-x);
	}

	uint32_t uint_len(uint32_t x)
	{
		if (x < 10)
			return 1;
		if (x < 100)
			return 2;
		if (x < 1000)
			return 3;
		if (x < 10000)
			return 4;
		if (x < 100000)
			return 5;
		if (x < 1000000)
			return 6;
		if (x < 10000000)
			return 7;
		return 8;
	}

	uint32_t cost_Nrun(uint32_t x) {
		return 2 + uint_len(x);
	}

	uint32_t cost_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos)
	{
		int dif_pos = (int)ref_pos - (int)pred_pos;

		uint32_t r = int_len(dif_pos);

		if (len != ~0u)
			r += 1 + uint_len(len - min_match_len);

		++r;

		return r;
	}

public:
	CLZDiff_V2(const uint32_t _min_match_len = 18) : CLZDiffBase(_min_match_len)
	{}

	virtual ~CLZDiff_V2() {};

	virtual void Encode(const contig_t& text, contig_t& encoded);
	virtual void Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded);

	virtual size_t Estimate(const contig_t& text, uint32_t bound = ~0u);
};

// EOF
#endif