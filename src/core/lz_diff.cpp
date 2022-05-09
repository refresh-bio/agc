// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.1
// Date   : 2022-05-06
// *******************************************************************************************

#include "../core/lz_diff.h"
#include <cmath>
#include <iostream>

// *******************************************************************************************
CLZDiffBase::CLZDiffBase(const uint32_t _min_match_len)
{
	min_match_len = _min_match_len;
	key_len = min_match_len - hashing_step + 1u;
	short_ht_ver = false;
	ht_mask = 0;
	ht_size = 0;
	index_ready = false;
}

// *******************************************************************************************
CLZDiffBase::~CLZDiffBase()
{
}

// *******************************************************************************************
bool CLZDiffBase::SetMinMatchLen(const uint32_t _min_match_len)
{
	if (!reference.empty() || index_ready)
		return false;

	min_match_len = _min_match_len;
	key_len = min_match_len - hashing_step + 1u;

	return true;
}

//#define USE_REVCOMP_REFERENCE
//#define USE_REV_REFERENCE
// *******************************************************************************************
void CLZDiffBase::prepare_gen(const contig_t& _reference)
{
	reference = _reference;

	reference.resize(reference.size() + key_len, invalid_symbol);

#ifdef USE_REVCOMP_REFERENCE
	reference.reserve(2 * reference.size());

	for (auto p = _reference.rbegin(); p != _reference.rend(); ++p)
		if (*p < 4)
			reference.emplace_back(3 - *p);
		else
			reference.emplace_back(*p);

	for (int i = 0; i < key_len; ++i)
		reference.emplace_back(invalid_symbol);
#endif

#ifdef USE_REV_REFERENCE
	reference.reserve(2 * reference.size());

	for (auto p = _reference.rbegin(); p != _reference.rend(); ++p)
		reference.emplace_back(*p);

	for (int i = 0; i < key_len; ++i)
		reference.emplace_back(invalid_symbol);
#endif

	reference.shrink_to_fit();
}

// *******************************************************************************************
void CLZDiffBase::prepare_index()
{
	ht_size = 0;

	uint32_t no_prev_valid = 0;

#ifdef USE_SPARSE_HT
	uint32_t cnt_mod = 0;
	uint32_t key_len_mod = key_len % hashing_step;

	for (auto c : reference)
	{
		if (c < 4)
			++no_prev_valid;
		else
			no_prev_valid = 0;

		if (++cnt_mod == hashing_step)
			cnt_mod = 0;

		if (cnt_mod == key_len_mod && no_prev_valid >= key_len)
			++ht_size;
	}
#else
	for (auto c : reference)
	{
		if (c < 4)
			++no_prev_valid;
		else
			no_prev_valid = 0;

		if (no_prev_valid >= key_len)
			++ht_size;
	}
#endif

	ht_size = (uint64_t)(ht_size / max_load_factor);

	while (ht_size & (ht_size - 1))
		ht_size &= ht_size - 1;

	ht_size <<= 1;
	ht_mask = ht_size - 1;

	if (short_ht_ver)
	{
		ht16.resize(ht_size, empty_key16);
		make_index16();
	}
	else
	{
		ht32.resize(ht_size, empty_key32);
		make_index32();
	}

	index_ready = true;
}

// *******************************************************************************************
void CLZDiffBase::Prepare(const contig_t& _reference)
{
	short_ht_ver = _reference.size() / hashing_step < 65535;

	prepare_gen(_reference);
}

// *******************************************************************************************
void CLZDiffBase::GetCodingCostVector(const contig_t& text, vector<uint32_t>& v_costs, const bool prefix_costs)
{
	if (!index_ready)
		prepare_index();

	v_costs.clear();
	v_costs.reserve(text.size());

	uint32_t text_size = (uint32_t)text.size();
	MurMur64Hash mmh;

	uint32_t i = 0;
	uint32_t pred_pos = 0;

	const uint8_t* text_ptr = text.data();

#ifdef USE_SPARSE_HT
	uint32_t no_prev_literals = 0;
#endif

	for (; i + key_len < text_size; )
	{
		uint64_t x = get_code(text_ptr);
		if (x == ~0ull)
		{
			uint32_t Nrun_len = get_Nrun_len(text_ptr, text_size - i);

			if (Nrun_len >= min_Nrun_len)
			{
				auto tc = coding_cost_Nrun(Nrun_len);
				if (prefix_costs)
				{
					v_costs.emplace_back(tc);
					v_costs.insert(v_costs.end(), Nrun_len - 1, 0);
				}
				else
				{
					v_costs.insert(v_costs.end(), Nrun_len - 1, 0);
					v_costs.emplace_back(tc);
				}

				text_ptr += Nrun_len;
				i += Nrun_len;
#ifdef USE_SPARSE_HT
				no_prev_literals = 0;
#endif
			}
			else
			{
				v_costs.emplace_back(1);
				++i;
				++pred_pos;
				++text_ptr;
#ifdef USE_SPARSE_HT
				++no_prev_literals;
#endif
			}

			continue;
		}

		uint64_t ht_pos = mmh(x) & ht_mask;

		uint32_t len_bck = 0;
		uint32_t len_fwd = 0;
		uint32_t match_pos = 0;
		uint32_t max_len = text_size - i;

		if (short_ht_ver ?
			!find_best_match16((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd) :
			!find_best_match32((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd))
		{
			v_costs.emplace_back(1);
			//			encode_literal(*text_ptr, encoded);
			++i;
			++text_ptr;
			++pred_pos;
#ifdef USE_SPARSE_HT
			++no_prev_literals;
#endif
			continue;
		}
		else
		{
#ifdef USE_SPARSE_HT
			if (len_bck)
			{
				for (uint32_t k = 0; k < len_bck; ++k)
					v_costs.pop_back();
				match_pos -= len_bck;
				pred_pos -= len_bck;
				text_ptr -= len_bck;
				i -= len_bck;
			}
#endif

			auto tc = coding_cost_match(match_pos, len_bck + len_fwd, pred_pos);

			if (prefix_costs)
			{
				v_costs.emplace_back(tc);
				v_costs.insert(v_costs.end(), len_bck + len_fwd - 1, 0);
			}
			else
			{
				v_costs.insert(v_costs.end(), len_bck + len_fwd - 1, 0);
				v_costs.emplace_back(tc);
			}

			pred_pos = match_pos + len_bck + len_fwd;
			i += len_bck + len_fwd;
			text_ptr += len_bck + len_fwd;

#ifdef USE_SPARSE_HT
			no_prev_literals = 0;
#endif
		}
	}

	for (; i < text_size; ++i)
		v_costs.emplace_back(1);
}

// *******************************************************************************************
bool CLZDiffBase::find_best_match16(uint32_t ht_pos, const uint8_t* s, const uint32_t max_len, const uint32_t no_prev_literals,
	uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd)
{
	len_fwd = 0;
	len_bck = 0;

	const uint8_t* ref_ptr = reference.data();

	for (uint32_t i = 0; i < max_no_tries; ++i)
	{
		if (ht16[ht_pos] == empty_key16)
			break;

		uint32_t f_len = 0;
		uint32_t h_pos = ht16[ht_pos] * hashing_step;

		const uint8_t* p = ref_ptr + h_pos;

		for (; f_len < max_len; ++f_len)
			if (s[f_len] != p[f_len])
				break;

		uint32_t b_len = 0;
		for (; b_len < min(no_prev_literals, h_pos); ++b_len)
			if (*(s - b_len - 1) != *(p - b_len - 1))
				break;

		if (b_len + f_len > len_bck + len_fwd)
		{
			len_bck = b_len;
			len_fwd = f_len;
			ref_pos = h_pos;
		}

		ht_pos = (uint32_t)((ht_pos + 1u) & ht_mask);
	}

	return len_bck + len_fwd >= min_match_len;
}

// *******************************************************************************************
bool CLZDiffBase::find_best_match32(uint32_t ht_pos, const uint8_t* s, const uint32_t max_len, const uint32_t no_prev_literals,
	uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd)
{
	len_fwd = 0;
	len_bck = 0;

	const uint8_t* ref_ptr = reference.data();

	for (uint32_t i = 0; i < max_no_tries; ++i)
	{
		if (ht32[ht_pos] == empty_key32)
			break;

		uint32_t f_len = 0;
		uint32_t h_pos = ht32[ht_pos] * hashing_step;
		const uint8_t* p = ref_ptr + h_pos;

		for (; f_len < max_len; ++f_len)
			if (s[f_len] != p[f_len])
				break;

		uint32_t b_len = 0;
		for (; b_len < min(no_prev_literals, h_pos); ++b_len)
			if (*(s - b_len - 1) != *(p - b_len - 1))
				break;

		if (b_len + f_len > len_bck + len_fwd)
		{
			len_bck = b_len;
			len_fwd = f_len;
			ref_pos = h_pos;
		}

		ht_pos = (uint32_t) ((ht_pos + 1u) & ht_mask);
	}

	return len_bck + len_fwd >= min_match_len;
}

// *******************************************************************************************
void CLZDiffBase::append_int(contig_t& text, int64_t x)
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

	contig_t tmp;

	for (; x; x /= 10)
		tmp.emplace_back((uint8_t) ('0' + (x % 10)));

	text.insert(text.end(), tmp.rbegin(), tmp.rend());
}

// *******************************************************************************************
void CLZDiffBase::read_int(contig_t::const_iterator& p, int64_t& x)
{
	bool is_neg = false;
	x = 0;

	if (*p == '-')
	{
		is_neg = true;
		++p;
	}

	while(*p >= '0' && *p <= '9')
		x = x * 10 + (int64_t)(*p++ - '0');

	if (is_neg)
		x = -x;
}

// *******************************************************************************************
void CLZDiffBase::encode_literal(const uint8_t c, contig_t& encoded)
{
	encoded.push_back('A' + c);
}

// *******************************************************************************************
void CLZDiffBase::encode_literal_diff(const uint8_t c, const uint8_t r, contig_t& encoded)
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

// *******************************************************************************************
void CLZDiffBase::encode_Nrun(const uint32_t len, contig_t &encoded)
{
	encoded.emplace_back(N_run_starter_code);		// N-run start marker
	append_int(encoded, len - min_Nrun_len);
	encoded.emplace_back(N_code);					// N-run stop marker
}

// *******************************************************************************************
uint32_t CLZDiffBase::coding_cost_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos) const
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

// *******************************************************************************************
uint32_t CLZDiffBase::coding_cost_Nrun(const uint32_t len) const
{
	return 2 + int_len(len - min_Nrun_len);
}

// *******************************************************************************************
uint32_t CLZDiffBase::get_Nrun_len(const uint8_t* s, const uint32_t max_len) const
{
	if (*s != N_code || *(s + 1) != N_code || *(s + 2) != N_code)
		return 0;

	uint32_t len;
	for (len = 3; len < max_len && *(s + len) == N_code; ++len)
		;

	return len;
}

// *******************************************************************************************
bool CLZDiffBase::is_literal(const contig_t::const_iterator& p) const
{
	return (*p >= 'A' && *p <= 'A' + 20) || (*p == '!');
}

// *******************************************************************************************
bool CLZDiffBase::is_Nrun(const contig_t::const_iterator& p) const
{
	return *p == N_run_starter_code;
}

// *******************************************************************************************
void CLZDiffBase::decode_literal(contig_t::const_iterator& p, uint8_t &c)
{
	if (*p == '!')
	{
		c = '!';
		++p;
	}
	else
		c = *p++ - 'A';
}

// *******************************************************************************************
void CLZDiffBase::decode_Nrun(contig_t::const_iterator& p, uint32_t& len)
{
	int64_t raw_len;

	++p;		// prefix
	read_int(p, raw_len);
	++p;		// suffix

	len = (uint32_t) (raw_len + min_Nrun_len);
}

// *******************************************************************************************
uint64_t CLZDiffBase::get_code(const uint8_t* s) const
{
	uint64_t x = 0;

	for (uint32_t i = 0; i < key_len; ++i, ++s)
	{
		if (*s > 3)
			return ~0ull;
		x = (x << 2) + (uint64_t) *s;
	}

	return x;
}

// *******************************************************************************************
void CLZDiffBase::make_index16()
{
	uint32_t ref_size = (uint32_t)reference.size();
	MurMur64Hash mmh;

	const uint8_t* ptr = reference.data();

#ifdef USE_SPARSE_HT
	for (uint32_t i = 0; i + key_len < ref_size; i += hashing_step, ptr += hashing_step)
#else
	for (uint32_t i = 0; i + key_len < ref_size; ++i, ++ptr)
#endif
	{
		uint64_t x = get_code(ptr);
		if (x == ~0ull)
			continue;
		uint64_t pos = mmh(x) & ht_mask;

		for (uint32_t j = 0; j < max_no_tries; ++j)
			if (ht16[(pos + j) & ht_mask] == empty_key16)
			{
				ht16[(pos + j) & ht_mask] = i / hashing_step;
				break;
			}
	}
}

// *******************************************************************************************
void CLZDiffBase::make_index32()
{
	uint32_t ref_size = (uint32_t)reference.size();
	MurMur64Hash mmh;

	const uint8_t* ptr = reference.data();

#ifdef USE_SPARSE_HT
	for (uint32_t i = 0; i + key_len < ref_size; i += hashing_step, ptr += hashing_step)
#else
		for (uint32_t i = 0; i + key_len < ref_size; ++i, ++ptr)
#endif
	{
		uint64_t x = get_code(ptr);
		if (x == ~0ull)
			continue;
		uint64_t pos = mmh(x) & ht_mask;

		for (uint32_t j = 0; j < max_no_tries; ++j)
			if (ht32[(pos + j) & ht_mask] == empty_key32)
			{
				ht32[(pos + j) & ht_mask] = i / hashing_step;
				break;
			}
	}
}

// *******************************************************************************************
void CLZDiffBase::GetReference(contig_t& s)
{
	if (reference.empty())
		s.clear();
	else
		s.assign(reference.begin(), reference.begin() + (reference.size() - key_len));
}


// *******************************************************************************************
//
// *******************************************************************************************
void CLZDiff_V1::encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded)
{
	int dif_pos = (int)ref_pos - (int)pred_pos;

	append_int(encoded, dif_pos);
	encoded.emplace_back(',');
	append_int(encoded, len - min_match_len);

	encoded.emplace_back('.');
}

// *******************************************************************************************
void CLZDiff_V1::decode_match(contig_t::const_iterator& p, uint32_t& ref_pos, uint32_t& len, uint32_t& pred_pos)
{
	int64_t raw_pos;
	int64_t raw_len;

	read_int(p, raw_pos);
	++p;		// ','

	ref_pos = (uint32_t)(raw_pos + (int64_t)pred_pos);

	if (*p == '.')
		len = ~0u;
	else
	{
		read_int(p, raw_len);
		len = (uint32_t)(raw_len + min_match_len);
	}

	++p;		// '.'
}

// *******************************************************************************************
void CLZDiff_V1::Encode(const contig_t& text, contig_t& encoded)
{
	if (!index_ready)
		prepare_index();

	uint32_t text_size = (uint32_t)text.size();

	encoded.clear();

#ifdef IMPROVED_LZ_ENCODING
	if (text_size == reference.size() - key_len)
		if (equal(text.begin(), text.end(), reference.begin()))
			return;														// equal sequences
#endif

	MurMur64Hash mmh;

	uint32_t i = 0;
	uint32_t pred_pos = 0;


	const uint8_t* text_ptr = text.data();

#ifdef USE_SPARSE_HT
//	const uint32_t thr_literal_len = 1000;

	uint32_t no_prev_literals = 0;
#endif

	for (; i + key_len < text_size; )
	{
		uint64_t x = get_code(text_ptr);

		if (x == ~0ull)
		{
			uint32_t Nrun_len = get_Nrun_len(text_ptr, text_size - i);

			if (Nrun_len >= min_Nrun_len)
			{
				encode_Nrun(Nrun_len, encoded);
				text_ptr += Nrun_len;
				i += Nrun_len;
#ifdef USE_SPARSE_HT
				/*				if (no_prev_literals > thr_literal_len)
									cerr << "literal len: " + to_string(no_prev_literals) + "    \n";*/
				no_prev_literals = 0;
#endif
			}
			else
			{
				encode_literal(*text_ptr, encoded);

				++i;
				++pred_pos;
				++text_ptr;
#ifdef USE_SPARSE_HT
				++no_prev_literals;
#endif
			}

			continue;
		}

		uint64_t ht_pos = mmh(x) & ht_mask;

		uint32_t len_bck = 0;
		uint32_t len_fwd = 0;
		uint32_t match_pos = 0;
		uint32_t max_len = text_size - i;

		if (short_ht_ver ?
			!find_best_match16((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd) :
			!find_best_match32((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd))
		{
			encode_literal(*text_ptr, encoded);

			++i;
			++text_ptr;
			++pred_pos;
#ifdef USE_SPARSE_HT
			++no_prev_literals;
#endif
			continue;
		}
		else
		{
#ifdef USE_SPARSE_HT
			if (len_bck)
			{
				for (uint32_t k = 0; k < len_bck; ++k)
					encoded.pop_back();
				match_pos -= len_bck;
				pred_pos -= len_bck;
				text_ptr -= len_bck;
				i -= len_bck;
			}
#endif

			encode_match(match_pos, len_bck + len_fwd, pred_pos, encoded);

			pred_pos = match_pos + len_bck + len_fwd;
			i += len_bck + len_fwd;
			text_ptr += len_bck + len_fwd;

#ifdef USE_SPARSE_HT
			/*			if (no_prev_literals > thr_literal_len)
							cerr << "litera len: " + to_string(no_prev_literals) + "     \n";*/

			no_prev_literals = 0;
#endif
		}
	}

	for (; i < text_size; ++i)
		encode_literal(text[i], encoded);
}

// *******************************************************************************************
void CLZDiff_V1::Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded)
{
	uint8_t c;
	uint32_t ref_pos, len;
	uint32_t pred_pos = 0;

	decoded.clear();

	for (auto p = encoded.begin(); p != encoded.end(); )
	{
		if (is_literal(p))
		{
			decode_literal(p, c);
			decoded.emplace_back(c);
			++pred_pos;
		}
		else if (is_Nrun(p))
		{
			decode_Nrun(p, len);
			decoded.insert(decoded.end(), len, N_code);
		}
		else
		{
			decode_match(p, ref_pos, len, pred_pos);
			decoded.insert(decoded.end(), reference.begin() + ref_pos, reference.begin() + ref_pos + len);
			pred_pos = ref_pos + len;
		}
	}
}


// *******************************************************************************************
//
// *******************************************************************************************
void CLZDiff_V2::encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded)
{
	int dif_pos = (int)ref_pos - (int)pred_pos;

	append_int(encoded, dif_pos);
	if (len != ~0u)
	{
		encoded.emplace_back(',');
		append_int(encoded, len - min_match_len);
	}

	encoded.emplace_back('.');
}

// *******************************************************************************************
void CLZDiff_V2::decode_match(contig_t::const_iterator& p, uint32_t& ref_pos, uint32_t& len, uint32_t& pred_pos)
{
	int64_t raw_pos;
	int64_t raw_len;

	read_int(p, raw_pos);
	ref_pos = (uint32_t)(raw_pos + (int64_t)pred_pos);

	if(*p == ',')
	{
		++p;
		read_int(p, raw_len);
		len = (uint32_t)(raw_len + min_match_len);
		++p;		// '.'
	}
	else
	{
		len = ~0u;
		++p;		// '.'
	}
}

// *******************************************************************************************
void CLZDiff_V2::Encode(const contig_t& text, contig_t& encoded)
{
	if (!index_ready)
		prepare_index();

	uint32_t text_size = (uint32_t)text.size();

	encoded.clear();

	if (text_size == reference.size() - key_len)
		if (equal(text.begin(), text.end(), reference.begin()))
			return;														// equal sequences

	MurMur64Hash mmh;

	uint32_t i = 0;
	uint32_t pred_pos = 0;

	const uint8_t* text_ptr = text.data();

#ifdef USE_SPARSE_HT
//	const uint32_t thr_literal_len = 1000;

	uint32_t no_prev_literals = 0;
#endif

	for (; i + key_len < text_size; )
	{
		uint64_t x = get_code(text_ptr);

		if (x == ~0ull)
		{
			uint32_t Nrun_len = get_Nrun_len(text_ptr, text_size - i);

			if (Nrun_len >= min_Nrun_len)
			{
				encode_Nrun(Nrun_len, encoded);
				text_ptr += Nrun_len;
				i += Nrun_len;
#ifdef USE_SPARSE_HT
				/*				if (no_prev_literals > thr_literal_len)
									cerr << "literal len: " + to_string(no_prev_literals) + "    \n";*/
				no_prev_literals = 0;
#endif
			}
			else
			{
				encode_literal(*text_ptr, encoded);

				++i;
				++pred_pos;
				++text_ptr;
#ifdef USE_SPARSE_HT
				++no_prev_literals;
#endif
			}

			continue;
		}

		uint64_t ht_pos = mmh(x) & ht_mask;

		uint32_t len_bck = 0;
		uint32_t len_fwd = 0;
		uint32_t match_pos = 0;
		uint32_t max_len = text_size - i;

		if (short_ht_ver ?
			!find_best_match16((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd) :
			!find_best_match32((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd))
		{
			encode_literal(*text_ptr, encoded);

			++i;
			++text_ptr;
			++pred_pos;
#ifdef USE_SPARSE_HT
			++no_prev_literals;
#endif
			continue;
		}
		else
		{
#ifdef USE_SPARSE_HT
			if (len_bck)
			{
				for (uint32_t k = 0; k < len_bck; ++k)
					encoded.pop_back();
				match_pos -= len_bck;
				pred_pos -= len_bck;
				text_ptr -= len_bck;
				i -= len_bck;
			}
#endif

			if (match_pos == pred_pos)
			{
				uint32_t e_size = encoded.size();
				for (uint32_t i = 1; i < e_size && i < match_pos; ++i)
				{
					if (encoded[e_size - i] < 'A' || encoded[e_size - i] > 'Z')
						break;
					if (encoded[e_size - i] - 'A' == reference[match_pos - i])
						encoded[e_size - i] = '!';
				}
			}

			if (i + len_bck + len_fwd == text_size && match_pos + len_bck + len_fwd == reference.size() - key_len)			// is match to end of sequence?
				encode_match(match_pos, ~0u, pred_pos, encoded);
			else
				encode_match(match_pos, len_bck + len_fwd, pred_pos, encoded);

			pred_pos = match_pos + len_bck + len_fwd;
			i += len_bck + len_fwd;
			text_ptr += len_bck + len_fwd;

#ifdef USE_SPARSE_HT
			/*			if (no_prev_literals > thr_literal_len)
							cerr << "litera len: " + to_string(no_prev_literals) + "     \n";*/

			no_prev_literals = 0;
#endif
		}
	}

	for (; i < text_size; ++i)
		encode_literal(text[i], encoded);
}

// *******************************************************************************************
void CLZDiff_V2::Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded)
{
	uint8_t c;
	uint32_t ref_pos, len;
	uint32_t pred_pos = 0;

	decoded.clear();

	for (auto p = encoded.begin(); p != encoded.end(); )
	{
		if (is_literal(p))
		{
			decode_literal(p, c);

			if (c == '!')
				c = reference[pred_pos];
			decoded.emplace_back(c);
			++pred_pos;
		}
		else if (is_Nrun(p))
		{
			decode_Nrun(p, len);
			decoded.insert(decoded.end(), len, N_code);
		}
		else
		{
			decode_match(p, ref_pos, len, pred_pos);
			
			if (len == ~0u)
				len = reference.size() - ref_pos;

			decoded.insert(decoded.end(), reference.begin() + ref_pos, reference.begin() + ref_pos + len);
			pred_pos = ref_pos + len;
		}
	}
}

// EOL