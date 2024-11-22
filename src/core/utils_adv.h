#ifndef _UTILS_ADV_H
#define _UTILS_ADV_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "../common/defs.h"
#include <mutex>
#if defined(ARCH_X64)
#include <nmmintrin.h>
#elif defined(ARCH_ARM)
#include <arm_neon.h>
#endif
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

		//		m_cond.wait(lock, [&] {return gen != m_generation; });
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

// *****************************************************************************************
//
class CAtomicBarrierWithIncrementing
{
	std::atomic<int32_t> a_count;
	std::atomic<int32_t> a_generation;
	int32_t count_reset_value;

public:
	CAtomicBarrierWithIncrementing(const CAtomicBarrierWithIncrementing&) = delete;
	CAtomicBarrierWithIncrementing& operator=(const CAtomicBarrierWithIncrementing&) = delete;
	explicit CAtomicBarrierWithIncrementing(int32_t count) :
		a_count(count - 1), a_generation(0),
		count_reset_value(count - 1)
	{
	}

	void arrive_and_wait()
	{
		int32_t old_generation = a_generation.load();

		if (!a_count.fetch_sub(1, memory_order_relaxed))
		{
			a_count = count_reset_value;
			++a_generation;
			a_generation.notify_all();
			return;
		}

		a_generation.wait(old_generation);
	}

	bool try_increment(int32_t inc = 1)
	{
		//		return false;

		auto new_count = a_count.fetch_add(inc, memory_order_relaxed) + inc;

		if (new_count <= count_reset_value)
			return true;

		a_count.fetch_sub(inc, memory_order_relaxed);
		return false;
	}

	int32_t try_increment_max(int32_t inc_req)
	{
		//		return 0;

		inc_req = min(inc_req, count_reset_value);

		auto new_count = a_count.fetch_add(inc_req, memory_order_relaxed) + inc_req;

		if (new_count <= count_reset_value)
			return inc_req;

		int32_t to_many = new_count - count_reset_value;

		if (to_many >= inc_req)		// Can happen if other thread is also in this function at the same time
			to_many = inc_req;

		int32_t no_ext = inc_req - to_many;

		a_count.fetch_sub(to_many, memory_order_relaxed);

		return no_ext;
	}

	void decrement(int32_t dec = 1)
	{
		//		return;

		a_count.fetch_sub(dec);
	}
};

// **********************************************************************************
class bloom_set_t {
	//	const uint32_t no_hashes = 2;
	const uint32_t no_hashes = 3;

	MurMur64Hash mmh;

	//	vector<uint64_t> arr;
	uint64_t* arr = nullptr;
	uint64_t* raw_arr = nullptr;

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

		return max((uint64_t)256, size * 2);
	}

	void allocate(size_t size)
	{
		if (raw_arr)
			delete[] raw_arr;
		//		arr.clear();
		no_elements = 0;

		allocated = normalize_size(size);

		//		arr.resize(allocated / 64, 0);
		raw_arr = new uint64_t[allocated / 64 + 7];
		arr = raw_arr;
		while (((uint64_t)arr) % 64 != 0)
			++arr;

		fill_n(arr, allocated / 64, 0ull);

		mask_shift = 6 * no_hashes;
		mask = (allocated / 64 - 1) << mask_shift;
	}

	void insert_impl(uint64_t x)
	{
		uint64_t h = mmh(x);

		uint64_t pos = (h & mask) >> mask_shift;

		arr[pos] |= (1ull << (h & 63)) | (1ull << ((h >> 6) & 63)) | (1ull << ((h >> 12) & 63));
		//		arr[pos] |=	(1ull << (h & 63)) | (1ull << ((h >> 6) & 63));

		++no_elements;
	}

public:
	bloom_set_t(size_t size = 64)
	{
		allocate(size);
	}

	~bloom_set_t()
	{
		if (raw_arr)
			delete[] raw_arr;
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

// EOF
#endif