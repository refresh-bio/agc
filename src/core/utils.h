#ifndef _UTILS_H
#define _UTILS_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include "../core/defs.h"
#include <mutex>
#include <nmmintrin.h>
#include <string>
#include <random>
#include <vector>
#include <algorithm>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <cstddef>

using namespace std;

#include <chrono>
#include <cstdint>

// **********************************************************************************
template<typename T>
constexpr T dna_code(const T x) 
{
	switch (x)
	{
	case (T) 'A':	return (T) 0; break;
	case (T) 'C':	return (T) 1; break;
	case (T) 'G':	return (T) 2; break;
	case (T) 'T':	return (T) 3; break;
	}

	return (T)4;
}

// **********************************************************************************
template<typename T>
constexpr T reverse_complement(const T x)
{
	switch (x)
	{
	case dna_code('A'):	return (T)dna_code('T'); break;
	case dna_code('C'):	return (T)dna_code('G'); break;
	case dna_code('G'):	return (T)dna_code('C'); break;
	case dna_code('T'):	return (T)dna_code('A'); break;
	}

	return (T) 4;
}

// **********************************************************************************
template<typename T>
constexpr T complement(const T x)
{
	switch (x)
	{
	case dna_code('A'):	return (T)dna_code('T'); break;
	case dna_code('C'):	return (T)dna_code('G'); break;
	case dna_code('G'):	return (T)dna_code('C'); break;
	case dna_code('T'):	return (T)dna_code('A'); break;
	}

	return (T)4;
}

// **********************************************************************************
constexpr uint8_t reverse_complement_alhpa(const uint8_t x)
{
	switch (x)
	{
	case 'A': return 'T'; break;
	case 'C': return 'G'; break;
	case 'G': return 'C'; break;
	case 'T': return 'A'; break;
	}

	return 'N';
}

// **********************************************************************************
template<typename T>
constexpr uint32_t pop_count(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x &= x - 1;

	return r;
}

// **********************************************************************************
template<typename T>
constexpr bool is_power_2(const T x)
{
	return (x & (x - (T)1)) == 0;
}

// **********************************************************************************
constexpr uint64_t ilog2(uint64_t x)
{
	uint64_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

// **********************************************************************************
constexpr uint64_t no_bytes(uint64_t x)
{
	uint64_t r = 1;

	x >>= 8;

	for (; x; ++r)
		x >>= 8;

	return r;
}

// **********************************************************************************
constexpr uint64_t zigzag_encode(int64_t x)
{
	if (x >= 0)
		return (uint64_t)(2 * x);
	else
		return (uint64_t)(2 * (-x) - 1);
}

// **********************************************************************************
constexpr int64_t zigzag_decode(uint64_t x)
{
	if (x & 1)
		return -(int64_t) (x + 1) / 2;
	else
		return (int64_t)(x / 2);
}

// **********************************************************************************
constexpr uint64_t zigzag_encode(uint64_t x_curr, uint64_t x_prev)
{
	if (x_curr < x_prev)
		return 2 * (x_prev - x_curr) - 1u;

	if (x_curr < 2 * x_prev)
		return 2 * (x_curr - x_prev);

	return x_curr;
}

// **********************************************************************************
constexpr uint64_t zigzag_decode(uint64_t x_val, uint64_t x_prev)
{
	if (x_val >= 2 * x_prev)
		return x_val;

	if (x_val & 1)
//		return (2 * x_prev - x_val - 1u) / 2;
		return (2 * x_prev - x_val) / 2;			// optimization (-1 is unnecessary due to /2)

	return (x_val + 2 * x_prev) / 2;
}

// *****************************************************************************************
string ss_prefix(uint32_t archive_version);
string ss_base(uint32_t archive_version, uint32_t n);
string ss_ref_name(uint32_t archive_version, uint32_t n);
string ss_delta_name(uint32_t archive_version, uint32_t n);
string ss_ref_ext(uint32_t archive_version);
string ss_delta_ext(uint32_t archive_version);
string int_to_hex(uint32_t n);
string int_to_base64(uint32_t n);


// *****************************************************************************************
//
class CBarrier
{
public:
	CBarrier(const CBarrier&) = delete;
	CBarrier& operator=(const CBarrier&) = delete;
	explicit CBarrier(unsigned int count) :
		m_count(count), m_generation(0),
		m_count_reset_value(count)
	{
	}
	void arrive_and_wait()
	{
		std::unique_lock< std::mutex > lock(m_mutex);
		unsigned int gen = m_generation;
		if (--m_count == 0)
		{
			m_generation++;
			m_count = m_count_reset_value;
			m_cond.notify_all();
			return;
		}
		while (gen == m_generation)
			m_cond.wait(lock);
	}
private:
	std::mutex m_mutex;
	std::condition_variable m_cond;
	unsigned int m_count;
	unsigned int m_generation;
	unsigned int m_count_reset_value;
};

// *****************************************************************************************
//
class CAtomicBarrier
{
public:
	CAtomicBarrier(const CAtomicBarrier&) = delete;
	CAtomicBarrier& operator=(const CAtomicBarrier&) = delete;
	explicit CAtomicBarrier(int32_t count) :
		a_count(count - 1), a_generation(0),
		count_reset_value(count - 1)
	{
	}

	void arrive_and_wait()
	{
		int32_t old_generation = a_generation;

		if (!a_count.fetch_sub(1, memory_order_relaxed))
		{
			a_count = count_reset_value;
			++a_generation;
		}

		while (a_generation == old_generation)
			;
	}
private:
	atomic<int32_t> a_count;
	atomic<int32_t> a_generation;
	int32_t count_reset_value;
};

// **********************************************************************************
struct MurMur32Hash
{
	std::size_t operator()(uint32_t h) const noexcept
	{
		h ^= h >> 16;
		h *= 0x85ebca6b;
		h ^= h >> 13;
		h *= 0xc2b2ae35;
		h ^= h >> 16;

		return h;
	}
};

// **********************************************************************************
// MurMurHash3
struct MurMur64Hash
{
	std::size_t operator()(size_t h) const noexcept
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53L;
		h ^= h >> 33;

		return h;
	}
};

// **********************************************************************************
struct hash_pair {
	template <typename T, typename S>
	size_t operator()(const pair<T, S>& x) const
	{
		return hash<T>{}(x.first) ^ hash<S>{}(x.second);
	}
};

template <>
struct std::hash<pair<uint64_t, uint64_t>>
{
	std::size_t operator()(const pair<uint64_t, uint64_t>& k) const
	{
		using std::size_t;
		using std::hash;

		return (hash<uint64_t>()(k.first)) ^ (hash<uint64_t>()(k.second));
	}
};

// **********************************************************************************
class bloom_set_t {
//	const uint32_t no_hashes = 2;
	const uint32_t no_hashes = 3;

	MurMur64Hash mmh;

	vector<uint64_t> arr;

	size_t no_elements;
	size_t allocated;
	size_t mask;
	uint32_t mask_shift;

	uint64_t normalize_size(uint64_t size)
	{
		size *= no_hashes;
		size *= 2;

		while (size & (size - 1))
			size &= size - 1;

		return max((uint64_t) 256, size * 2);
	}

	void allocate(size_t size)
	{
		arr.clear();
		no_elements = 0;

		allocated = normalize_size(size);
		
		arr.resize(allocated / 64, 0);

		mask_shift = 6 * no_hashes;
		mask = (allocated / 64 - 1) << mask_shift;
	}

	void insert_impl(uint64_t x)
	{
		uint64_t h = mmh(x);

		uint64_t pos = (h & mask) >> mask_shift;

		arr[pos] |=	(1ull << (h & 63)) | (1ull << ((h >> 6) & 63)) | (1ull << ((h >> 12) & 63));
//		arr[pos] |=	(1ull << (h & 63)) | (1ull << ((h >> 6) & 63));

		++no_elements;
	}

public:
	bloom_set_t(size_t size = 64)
	{
		allocate(size);
	}

	void resize(size_t size)
	{
		allocate(size);
	}

	template<typename Iter>
	void insert(Iter begin, Iter end)
	{
		for (auto p = begin; p != end; ++p)
			insert_impl(*p);
	}

	void insert(uint64_t x)
	{
		insert_impl(x);
	}

	bool check(uint64_t x)
	{
		uint64_t h = mmh(x);

		uint64_t pos = (h & mask) >> mask_shift;

		return (arr[pos] & (1ull << (h & 63))) && (arr[pos] & (1ull << ((h >> 6) & 63))) && (arr[pos] & (1ull << ((h >> 12) & 63)));
//		return (arr[pos] & (1ull << (h & 63))) && (arr[pos] & (1ull << ((h >> 6) & 63)));
	}

	double filling_factor()
	{
		return (double)no_hashes * no_elements / allocated;
	}
};

// EOF
#endif