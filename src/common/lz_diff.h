#ifndef _LZ_DIFF_H
#define _LZ_DIFF_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include <string>
#include <vector>
#include <array>
#include "../common/utils.h"

#include <refresh/string_operations/lib/string_operations.h>

using namespace std;

#define USE_SPARSE_HT

// *******************************************************************************************
class CLZDiffBase
{
protected:
	const uint32_t empty_key32 = ~0u;
	const uint16_t empty_key16 = ~(uint16_t) 0u;
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

	uint64_t get_code(const uint8_t* s) const
	{
		uint64_t x = 0;

		uint32_t i = key_len % 4;

		switch (i)
		{
		case 3:
			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s++;
			[[fallthrough]];
		case 2:
			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s++;
			[[fallthrough]];
		case 1:
			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s++;
		}

		for (; i < key_len; )
		{
			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s;
			++i; ++s;

			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s;
			++i; ++s;

			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s;
			++i; ++s;

			if (*s > 3)
				return ~0ull;
			x = (x << 2) + (uint64_t)*s;
			++i; ++s;
		}

		return x;
	};

	uint64_t get_code_skip1(uint64_t code, const uint8_t* s) const
	{
		s += key_len - 1;

		if (*s > 3)
			return ~0ull;

		code = (code << 2) & key_mask;

		code += *s;

		return code;
	}

	uint32_t get_Nrun_len(const uint8_t* s, const uint32_t max_len) const
	{
		if (*s != N_code || *(s + 1) != N_code || *(s + 2) != N_code)
			return 0;

		uint32_t len;
		for (len = 3; len < max_len && *(s + len) == N_code; ++len)
			;

		return len;
	}

	void encode_literal(const uint8_t c, contig_t& encoded) const
	{
		encoded.push_back('A' + c);
	}

	void encode_literal_diff(const uint8_t c, const uint8_t r, contig_t& encoded) const
	{
		if (r == 0 || (r > 3 || c > 3))
			encoded.push_back(c);
		else
		{
			if (c < r)
				encoded.push_back(3 - c);
			else
				encoded.push_back(c - r);
		}
	}

	void encode_Nrun(const uint32_t len, contig_t& encoded) const
	{
		encoded.emplace_back(N_run_starter_code);		// N-run start marker
		append_int(encoded, len - min_Nrun_len);
		encoded.emplace_back(N_code);					// N-run stop marker
	}

	uint32_t coding_cost_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos) const
	{
		uint32_t r;
		int dif_pos = (int)ref_pos - (int)pred_pos;

		if (dif_pos >= 0)
			r = int_len((uint32_t)dif_pos);
		else
			r = int_len((uint32_t)-dif_pos) + 1;

		r += int_len(len - min_match_len) + 2;

		return r;
	}

	uint32_t coding_cost_Nrun(const uint32_t len) const
	{
		return 2 + int_len(len - min_Nrun_len);
	}

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

	bool is_literal(const contig_t::const_iterator& p) const
	{
		return (*p >= 'A' && *p <= 'A' + 20) || (*p == '!');
	}

	bool is_Nrun(const contig_t::const_iterator& p) const
	{
		return *p == N_run_starter_code;
	}

	void decode_literal(contig_t::const_iterator& p, uint8_t &c) const
	{
		if (*p == '!')
		{
			c = '!';
			++p;
		}
		else
			c = *p++ - 'A';
	}

	void decode_Nrun(contig_t::const_iterator& p, uint32_t& len) const
	{
		int64_t raw_len;

		++p;		// prefix
		read_int(p, raw_len);
		++p;		// suffix

		len = (uint32_t)(raw_len + min_Nrun_len);
	}

	bool find_best_match16(uint32_t ht_pos, const uint8_t *s, const uint32_t max_len, const uint32_t no_prev_literals,
		uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd) const;
	bool find_best_match32(uint32_t ht_pos, const uint8_t *s, const uint32_t max_len, const uint32_t no_prev_literals,
		uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd) const;
	inline void append_int(contig_t& text, int64_t x) const
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

		size_t i = text.size();
		text.resize(text.size() + (tmp + 16 - p));

		for (; p != tmp + 16; ++p, ++i)
			text[i] = *p;

		text.insert(text.end(), p, tmp + 16);
	}

	void read_int(contig_t::const_iterator& p, int64_t &x) const
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

	uint32_t compare_fwd(uint8_t* p, uint8_t* q, uint32_t max_len) const
	{
		return (uint32_t)refresh::matching_length(p, q, max_len);

#if 0
		uint32_t len = 0;

		auto p0 = p;

		switch (max_len % 4)
		{
		case 3:
			if (*p++ != *q++)
				return len;
			++len;
			[[fallthrough]];
		case 2:
			if (*p++ != *q++)
				return len;
			++len;
			[[fallthrough]];
		case 1:
			if (*p++ != *q++)
				return len;
			++len;
		}

		for (; len < max_len; len += 4)
		{
/*			if (*p++ != *q++)
				return len;
			if (*p++ != *q++)
				return len+1;
			if (*p++ != *q++)
				return len+2;
			if (*p++ != *q++)
				return len+3;*/
			int inc = *p == *q;			p += inc;			q += inc;
			inc = *p == *q;				p += inc;			q += inc;
			inc = *p == *q;				p += inc;			q += inc;
			inc = *p == *q;				p += inc;			q += inc;

			if (!inc)
				break;
		}

//		return len;

		return p - p0;
#endif
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

	void AssureIndex();

	void GetReference(contig_t& s);
	void GetCodingCostVector(const contig_t& text, vector<uint32_t> &v_costs, const bool prefix_costs) const;
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

	uint32_t int_len(int x) const
	{
		if (x >= 0)
			return uint_len((uint32_t)x);
		else
			return 1 + uint_len((uint32_t)-x);
	}

	uint32_t uint_len(uint32_t x) const
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

	uint32_t cost_Nrun(uint32_t x) const 
	{
		return 2 + uint_len(x);
	}

	uint32_t cost_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos) const
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