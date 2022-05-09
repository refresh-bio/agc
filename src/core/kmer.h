#ifndef _KMER_H
#define _KMER_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.1
// Date   : 2022-05-06
// *******************************************************************************************

#include "../core/defs.h"
#include <algorithm>
#include <bitset>
#include "../core/utils.h"

enum class kmer_mode_t {direct, rev_comp, canonical};

class CKmer
{
	uint64_t kmer_dir;
	uint64_t kmer_rc;
	uint32_t cur_size;
	uint32_t max_size;
	kmer_mode_t variant;
	uint64_t mask;
	uint32_t shift;
	uint64_t kernel_mask;
	uint32_t kernel_shift;

	// *******************************************************************************************
	inline void insert_direct(uint64_t symbol) {
		if (cur_size == max_size)
		{
			kmer_dir <<= 2;
			kmer_dir += symbol << shift;
		}
		else
		{
			++cur_size;
			kmer_dir += symbol << (64 - 2 * cur_size);
		}
	}

	// *******************************************************************************************
	inline void insert_direct_zero() {
		if (cur_size == max_size)
			kmer_dir <<= 2;
		else
			++cur_size;
	}

	// *******************************************************************************************
	inline void insert_rev_comp(uint64_t symbol) {
		kmer_rc >>= 2;
		kmer_rc += reverse_complement(symbol) << 62;
		kmer_rc &= mask;

		if (cur_size < max_size)
			++cur_size;
	}

	// *******************************************************************************************
	inline void insert_rev_comp_zero() {
		kmer_rc >>= 2;
		kmer_rc += reverse_complement(0ull) << 62;
		kmer_rc &= mask;

		if (cur_size < max_size)
			++cur_size;
	}
	
	// *******************************************************************************************
	inline void insert_canonical(uint64_t symbol) {
		// rev. comp. code
		kmer_rc >>= 2;
		kmer_rc += reverse_complement(symbol) << 62;
		kmer_rc &= mask;

		// direct code
		if (cur_size == max_size)
		{
			kmer_dir <<= 2;
			kmer_dir += symbol << shift;
		}
		else
		{
			++cur_size;
			kmer_dir += symbol << (64 - 2 * cur_size);
		}
	}

	// *******************************************************************************************
	inline void insert_canonical_zero() {
		// rev. comp. code
		kmer_rc >>= 2;
		kmer_rc += reverse_complement(0ull) << 62;
		kmer_rc &= mask;

		// direct code
		if (cur_size == max_size)
		{
			kmer_dir <<= 2;
		}
		else
		{
			++cur_size;
		}
	}

	// *******************************************************************************************
	inline void insert_front_direct(uint64_t symbol) {
		if (cur_size < max_size)
		{
			kmer_dir >>= 2;
			kmer_dir += symbol << 62;
			++cur_size;
		}
	}

	// *******************************************************************************************
	inline void insert_front_rev_comp(uint64_t symbol) {
		if (cur_size < max_size)
		{
			kmer_rc += reverse_complement(symbol) << (62 - 2 * cur_size);

			++cur_size;
		}
	}

	// *******************************************************************************************
	inline void insert_front_canonical(uint64_t symbol) {
		if (cur_size < max_size)
		{
			kmer_dir >>= 2;
			kmer_dir += symbol << 62;

			kmer_rc += reverse_complement(symbol) << (62 - 2 * cur_size);

			++cur_size;
		}
	}

	// *******************************************************************************************
	inline void replace_direct(uint64_t symbol, uint32_t pos) {
		uint32_t sym_shift = 62 - 2 * pos;
		uint64_t sym_mask = ~(3ull << sym_shift);

		kmer_dir &= sym_mask;
		kmer_dir += symbol << sym_shift;
	}

	// *******************************************************************************************
	inline void replace_direct_last(uint64_t symbol) {
		uint32_t sym_shift = 64 - 2 * cur_size;
		uint64_t m = ~(3ull << sym_shift);
		kmer_dir &= m;
		kmer_dir += symbol << sym_shift;
	}

	// *******************************************************************************************
	inline void replace_rev_comp(uint64_t symbol, uint32_t pos) {
		uint32_t sym_shift = 64 - 2 * cur_size + 2 * pos;
		uint64_t sym_mask = ~(3ull << sym_shift);

		kmer_rc &= sym_mask;
		kmer_rc += reverse_complement(symbol) << sym_shift;
	}

	// *******************************************************************************************
	inline void replace_rev_comp_last(uint64_t symbol) {
		kmer_rc <<= 2;
		kmer_rc >>= 2;
		kmer_rc += reverse_complement(symbol) << 62;
	}

public:
	// *******************************************************************************************
	CKmer() : kernel_shift(0)
	{
		max_size = 0;
	}

	// *******************************************************************************************
	CKmer(uint32_t _max_size, kmer_mode_t _variant) : kernel_shift(0) {
		Reset(_max_size, _variant);
	}

	// *******************************************************************************************
	CKmer(uint64_t _kmer_dir, uint64_t _kmer_rc, uint32_t _max_size, kmer_mode_t _variant) {
		kmer_dir = _kmer_dir;
		kmer_rc = _kmer_rc;
		max_size = _max_size;
		variant = _variant;
		cur_size = _max_size;

		shift = 64 - 2 * max_size;
		mask = (~0ull) << shift;

#ifdef KMER_MARGIN_2_SYMBOLS
		kernel_mask = (1ull << (2 * max_size - 8)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 4;
		kernel_mask <<= kernel_shift;
#else
		kernel_mask = (1ull << (2 * max_size - 4)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 2;
		kernel_mask <<= kernel_shift;
#endif
	}

	// *******************************************************************************************
	CKmer(uint64_t _kmer_dir, uint64_t _kmer_rc, uint32_t _max_size, uint32_t _cur_size, kmer_mode_t _variant) {
		kmer_dir = _kmer_dir;
		kmer_rc = _kmer_rc;
		max_size = _max_size;
		variant = _variant;
		cur_size = _cur_size;

		shift = 64 - 2 * max_size;
		mask = (~0ull) << shift;

#ifdef KMER_MARGIN_2_SYMBOLS
		kernel_mask = (1ull << (2 * max_size - 8)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 4;
		kernel_mask <<= kernel_shift;
#else
		kernel_mask = (1ull << (2 * max_size - 4)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 2;
		kernel_mask <<= kernel_shift;
#endif
	}

	// *******************************************************************************************
	void Reset() {
		kmer_dir = 0ull;
		kmer_rc = 0ull;
		cur_size = 0;
	}

	// *******************************************************************************************
	void ResetFromCan(CKmer &_kmer_can, kmer_mode_t _variant)
	{
		variant = _variant;
		max_size = _kmer_can.max_size;
		cur_size = _kmer_can.cur_size;
		
		if (variant == kmer_mode_t::direct)
		{
			kmer_dir = _kmer_can.kmer_dir;
			kmer_rc = 0ull;
		}
		else
		{
			kmer_dir = 0ull;
			kmer_rc = _kmer_can.kmer_rc;
		}

		shift = 64 - 2 * max_size;
		mask = (~0ull) << shift;

#ifdef KMER_MARGIN_2_SYMBOLS
		kernel_mask = (1ull << (2 * max_size - 8)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 4;
		kernel_mask <<= kernel_shift;
#else
		kernel_mask = (1ull << (2 * max_size - 4)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 2;
		kernel_mask <<= kernel_shift;
#endif
	}
	
	// *******************************************************************************************
	void Reset(uint32_t _max_size, kmer_mode_t _variant) {
		max_size = _max_size;
		variant = _variant;
		kmer_dir = 0ull;
		kmer_rc = 0ull;
		cur_size = 0;

		shift = 64 - 2 * max_size;
		mask = (~0ull) << shift;

#ifdef KMER_MARGIN_2_SYMBOLS
		kernel_mask = (1ull << (2 * max_size - 8)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 4;
		kernel_mask <<= kernel_shift;
#else
		kernel_mask = (1ull << (2 * max_size - 4)) - 1ull;
		kernel_shift = 64 - 2 * max_size + 2;
		kernel_mask <<= kernel_shift;
#endif
	}

	// *******************************************************************************************
	void insert(uint64_t symbol) {
		if (variant == kmer_mode_t::direct)
			insert_direct(symbol);
		else if (variant == kmer_mode_t::rev_comp)
			insert_rev_comp(symbol);
		else
			insert_canonical(symbol);
	}

	// *******************************************************************************************
	void insert_zero() {
		if (variant == kmer_mode_t::direct)
			insert_direct_zero();
		else if (variant == kmer_mode_t::rev_comp)
			insert_rev_comp_zero();
		else
			insert_canonical_zero();
	}
	
	// *******************************************************************************************
	void insert_front(uint64_t symbol) {
		if (variant == kmer_mode_t::direct)
			insert_front_direct(symbol);
		else if (variant == kmer_mode_t::rev_comp)
			insert_front_rev_comp(symbol);
		else
			insert_front_canonical(symbol);
	}
	
	// *******************************************************************************************
	void replace(uint64_t symbol, uint32_t pos) {
		if (variant != kmer_mode_t::rev_comp)
			replace_direct(symbol, pos);
		if (variant != kmer_mode_t::direct)
			replace_rev_comp(symbol, pos);
	}

	// *******************************************************************************************
	void replace_last(uint64_t symbol) {
		if (variant != kmer_mode_t::rev_comp)
			replace_direct_last(symbol);
		if (variant != kmer_mode_t::direct)
			replace_rev_comp_last(symbol);
	}

	// *******************************************************************************************
	uint64_t data() const {
		if (variant == kmer_mode_t::direct)
			return kmer_dir;
		else if (variant == kmer_mode_t::rev_comp)
			return kmer_rc;
		else
			return min(kmer_dir, kmer_rc);
	}

	// *******************************************************************************************
	uint64_t data_dir() const {
		if (variant == kmer_mode_t::direct || variant == kmer_mode_t::canonical)
			return kmer_dir;

		return 0;
	}

	// *******************************************************************************************
	uint64_t data_rc() const {
		if (variant == kmer_mode_t::rev_comp || variant == kmer_mode_t::canonical)
			return kmer_rc;

		return 0;
	}

	// *******************************************************************************************
	uint64_t data_normalized() const {
		if (variant != kmer_mode_t::canonical)
			return 0;

		uint64_t kernel_dir = kmer_dir & kernel_mask;
		uint64_t kernel_rc = kmer_rc & kernel_mask;

		if (kernel_dir < kernel_rc)
			return kmer_dir;
		else
			return kmer_rc;
	}

	// *******************************************************************************************
	bool is_normalized_dir() const {
		uint64_t kernel_dir = kmer_dir & kernel_mask;
		uint64_t kernel_rc = kmer_rc & kernel_mask;

		return kernel_dir < kernel_rc;
	}

	// *******************************************************************************************
	uint64_t data_aligned() const {
		if (variant == kmer_mode_t::direct)
			return kmer_dir >> (64 - 2 * cur_size);
		else if (variant == kmer_mode_t::rev_comp)
			return kmer_rc >> (64 - 2 * cur_size);
		else
			return min((kmer_dir >> (64 - 2 * cur_size)), kmer_rc >> (64 - 2 * cur_size));
	}

	// *******************************************************************************************
	uint64_t data_aligned_dir() const {
		return kmer_dir >> (64 - 2 * cur_size);
	}

	// *******************************************************************************************
	uint64_t data_aligned_rc() const {
		return kmer_rc >> (64 - 2 * cur_size);
	}
	
	// *******************************************************************************************
	uint64_t kernel_canonical() const
	{
		uint64_t kernel_dir = (kmer_dir & kernel_mask) >> kernel_shift;
		uint64_t kernel_rc = (kmer_rc & kernel_mask) >> kernel_shift;

		return min(kernel_dir, kernel_rc);
	}

	// *******************************************************************************************
	uint64_t kernel_canonical_plus1()
	{
		uint64_t kernel_dir = ((kmer_dir << 2) & kernel_mask) >> kernel_shift;
		uint64_t kernel_rc = ((kmer_rc >> 2) & kernel_mask) >> kernel_shift;

		return min(kernel_dir, kernel_rc);
	}

	// *******************************************************************************************
	uint64_t kernel_canonical_plus2()
	{
		uint64_t kernel_dir = ((kmer_dir << 4) & kernel_mask) >> kernel_shift;
		uint64_t kernel_rc = ((kmer_rc >> 4) & kernel_mask) >> kernel_shift;

		return min(kernel_dir, kernel_rc);
	}

	// *******************************************************************************************
	bool operator==(const CKmer&x) const {
		if (variant != kmer_mode_t::rev_comp)
			return kmer_dir == x.kmer_dir;
		else
			return kmer_rc == x.kmer_rc;
	}

	// *******************************************************************************************
	bool operator!=(const CKmer&x) const {
		if (variant != kmer_mode_t::rev_comp)
			return kmer_dir != x.kmer_dir;
		else
			return kmer_rc != x.kmer_rc;
	}

	// *******************************************************************************************
	bool cmp_symbol(const CKmer &x, uint32_t pos) {
		return get_symbol(pos) == x.get_symbol(pos);
	}

	// *******************************************************************************************
	bool cmp_dir_rc(const CKmer &x) {
		if (variant == x.variant)
			return false;

		if (cur_size != x.cur_size)
			return false;

		for (uint32_t i = 0; i < cur_size; ++i)
			if (get_symbol(i) != reverse_complement(x.get_symbol(i)))
				return false;

		return true;
	}

	// *******************************************************************************************
	uint64_t get_symbol(uint32_t pos) const {
		uint32_t sym_shift;

		if (variant != kmer_mode_t::rev_comp)
		{
			sym_shift = 62 - 2 * pos;
			return (kmer_dir >> sym_shift) & 3;
		}
		else
		{
			sym_shift = 64 - 2 * cur_size + 2 * pos;
			return (kmer_rc >> sym_shift) & 3;
		}
	}

	// *******************************************************************************************
	uint64_t get_prefix(uint32_t len) const
	{
		if (variant != kmer_mode_t::rev_comp)
		{
			uint32_t prefix_shift = (64 - 2 * len);
			return kmer_dir >> prefix_shift;
		}
		else
			;

		return 0;
	}

	// *******************************************************************************************
	void shorten(uint32_t len)
	{
		if (variant != kmer_mode_t::rev_comp)
		{
			kmer_dir <<= 2ull * (cur_size - len);
			cur_size = len;
		}

		if (variant != kmer_mode_t::direct)
		{
			cur_size = len;
			uint64_t loc_mask = (~0ull) << (64 - 2 * len);
			kmer_rc &= loc_mask;
		}
	}

	// *******************************************************************************************
	bool is_full() const
	{
		return cur_size == max_size;
	}

	// *******************************************************************************************
	bool is_almost_full(uint32_t margin) const
	{
		return cur_size + margin >= max_size;
	}

	// *******************************************************************************************
	bool is_dir_oriented() const
	{
		if (variant == kmer_mode_t::canonical)
			return kmer_dir <= kmer_rc;

		return false;
	}

	// *******************************************************************************************
	void swap_dir_rc()
	{
		if (variant == kmer_mode_t::canonical)
		{
			auto t = kmer_dir;
			kmer_dir = kmer_rc;
			kmer_rc = t;
		}
	}

	// *******************************************************************************************
	uint32_t get_cur_size() const
	{
		return cur_size;
	}

	// ************************************************************************************
	uint32_t get_max_size() const
	{
		return max_size;
	}

	// *******************************************************************************************
	void do_rev_comp()
	{
		swap(kmer_dir, kmer_rc);
	}
};

// EOF
#endif