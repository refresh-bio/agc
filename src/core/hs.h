#ifndef _HS_H
#define _HS_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-04-05
// *******************************************************************************************

#include <mmintrin.h>
#include <cstdint>
#include <xmmintrin.h>
#include <cstddef>

#include <algorithm>
#include <vector>
#include <utility>
#include <iterator>

#if defined(_MSC_VER)  /* Visual Studio */
#define FORCE_INLINE __forceinline
#define NO_INLINE __declspec(noinline)
#elif defined(__GNUC__)
#define FORCE_INLINE __inline__ __attribute__((always_inline, unused))
#define NO_INLINE __attribute__((noinline))
#else
#define FORCE_INLINE
#define NO_INLINE
#endif

// *******************************************************************************************
// *** Global const iterator
template<typename HashSetLP>
class const_hash_set_lp_iterator
{
	friend HashSetLP;

public:
	using key_type = typename HashSetLP::key_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = key_type*;
	using reference = key_type&;
	using vector_iterator_type = typename HashSetLP::VectorType::const_iterator;
	using key_equal = typename HashSetLP::key_equal;

	// *******************************************************************************************
	const_hash_set_lp_iterator() = default;

	// *******************************************************************************************
	const_hash_set_lp_iterator(HashSetLP* _p_hs)
	{
		empty_key = _p_hs->empty_key;
		iter = _p_hs->data.end();
		iter_end = _p_hs->data.end();
		compare = _p_hs->compare;
	}

	// *******************************************************************************************
	const_hash_set_lp_iterator(HashSetLP* _p_hs, vector_iterator_type _iter)
	{
		iter = _iter;
		iter_end = _p_hs->data.end();
		empty_key = _p_hs->empty_key;
		compare = _p_hs->compare;
	}

	// *******************************************************************************************
	const key_type& operator*() const
	{
		return *iter;
	}

	// *******************************************************************************************
	const key_type* operator->() const
	{
		return &(*iter);
	}

	// *******************************************************************************************
	const_hash_set_lp_iterator<HashSetLP>& operator++()
	{
		increment();
		return *this;
	}

	// *******************************************************************************************
	const_hash_set_lp_iterator<HashSetLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	// *******************************************************************************************
	bool operator==(const const_hash_set_lp_iterator<HashSetLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	// *******************************************************************************************
	bool operator!=(const const_hash_set_lp_iterator<HashSetLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	vector_iterator_type iter;
	vector_iterator_type iter_end;
	key_type empty_key;
	key_equal compare;

	// *******************************************************************************************
	void increment()
	{
		for (++iter; iter != iter_end && compare(*iter, empty_key); ++iter)
			;
	}
};

// *******************************************************************************************
// *** Global iterator
template <typename HashSetLP>
class hash_set_lp_iterator : public const_hash_set_lp_iterator<HashSetLP>
{
	friend HashSetLP;

public:
	using key_type = typename const_hash_set_lp_iterator<HashSetLP>::key_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = key_type*;
	using reference = key_type&;
	using vector_iterator_type = typename HashSetLP::VectorType::iterator;

	// *******************************************************************************************
	hash_set_lp_iterator() = default;

	// *******************************************************************************************
	hash_set_lp_iterator(HashSetLP* _p_hs) : const_hash_set_lp_iterator<HashSetLP>(_p_hs)
	{
	}

	// *******************************************************************************************
	hash_set_lp_iterator(HashSetLP* _p_hs, vector_iterator_type _iter) : const_hash_set_lp_iterator<HashSetLP>(_p_hs, _iter)
	{
	}

	// *******************************************************************************************
	key_type& operator*()
	{
		return const_cast<key_type&>(*this->iter);
	}

	// *******************************************************************************************
	key_type* operator->()
	{
		return const_cast<key_type*>(&(*this->iter));
	}
};

// *******************************************************************************************
// *** Hash set with linear probing 
template<typename Key_t,
	typename Compare_t = std::equal_to<Key_t>,
	typename Hash_t = std::hash<Key_t>>
	class hash_set_lp {
	public:
		using key_type = Key_t;
		using hasher = Hash_t;
		using key_equal = Compare_t;
		using reference = Key_t&;
		using const_reference = const Key_t&;
		using size_type = size_t;
		using difference_type = ptrdiff_t;
		using hash_set_lp_type = hash_set_lp<Key_t, Compare_t, Hash_t>;
		using iterator = hash_set_lp_iterator<hash_set_lp_type>;
		using const_iterator = const_hash_set_lp_iterator<hash_set_lp_type>;

	private:
		using VectorType = std::vector<key_type>;

	private:
		Key_t empty_key;

		Compare_t compare;
		Hash_t hash;

		double max_fill_factor;
		size_t no_elements;
		std::vector<key_type> data;
		size_t allocated;
		size_t size_when_restruct;
		size_t allocated_mask;

	public:
		friend class hash_set_lp_iterator<hash_set_lp_type>;
		friend class const_hash_set_lp_iterator<hash_set_lp_type>;

		// *******************************************************************************************
		virtual ~hash_set_lp() = default;

		// *******************************************************************************************
		explicit hash_set_lp(const Key_t _empty_key = Key_t(), size_t _init_reserved = 16, double _max_fill_factor = 0.7,
//			const Compare_t& _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			const Compare_t _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			empty_key(_empty_key),
			compare(_compare),
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);
		}

		// *******************************************************************************************
		hash_set_lp(const hash_set_lp<Key_t, Compare_t, Hash_t>& src) = default;

		// *******************************************************************************************
		hash_set_lp(hash_set_lp<Key_t, Compare_t, Hash_t>&& src) = default;

		// *******************************************************************************************
		template <typename InputIterator>
		hash_set_lp(InputIterator first, InputIterator last,
			const Key_t _empty_key = Key_t(), size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Compare_t& _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			empty_key(_empty_key),
			compare(_compare),
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);

			for (auto p = first; p != last; ++p)
				insert_fast(*p);
		}

		// *******************************************************************************************
		explicit hash_set_lp(std::initializer_list<key_type> il,
			const Key_t _empty_key = Key_t(), size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Compare_t& _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			empty_key(_empty_key),
			compare(_compare),
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);

			for (auto& x : il)
				insert_fast(x);
		}


		// *******************************************************************************************
		hash_set_lp_type& operator=(
			const hash_set_lp<Key_t, Compare_t, Hash_t>& rhs)
		{
			if (this != &rhs)
			{
				data.clear();
				data.shrink_to_fit();

				compare = rhs.compare;
				hash = rhs.hash;
				data = rhs.data;
				no_elements = rhs.no_elements;
				allocated = rhs.allocated;
				size_when_restruct = rhs.size_when_restruct;
				allocated_mask = rhs.allocated_mask;
				empty_key = rhs.empty_key;
			}

			return *this;
		}

		// *******************************************************************************************
		hash_set_lp_type& operator=(
			hash_set_lp<Key_t, Compare_t, Hash_t>&& rhs)
		{
			if (this != &rhs)
			{
				compare = rhs.compare;
				hash = rhs.hash;
				data = move(rhs.data);
				no_elements = rhs.no_elements;
				allocated = rhs.allocated;
				size_when_restruct = rhs.size_when_restruct;
				allocated_mask = rhs.allocated_mask;
				empty_key = rhs.empty_key;
			}

			return *this;
		}

		//	hash_set_type& operator=(std::initializer_list<key_type> il);

// *******************************************************************************************
		void swap(hash_set_lp& x)
		{
			swap(compare, x.compare);
			swap(hash, x.hash);

			swap(max_fill_factor, x.max_fill_factor);
			swap(no_elements, x.no_elements);
			swap(data, x.data);
			swap(allocated, x.allocated);
			swap(size_when_restruct, x.size_when_restruct);
			swap(allocated_mask, x.allocated_mask);

			swap(empty_key, x.empty_key);
		}

		// *******************************************************************************************
		void clear()
		{
			fill(data.begin(), data.end(), empty_key);
			no_elements = 0;
		}

		// *******************************************************************************************
		iterator begin()
		{
			if (!no_elements)
				return end();

			auto p = hash_set_lp_iterator<hash_set_lp_type>(this, data.begin());
			if (*p == empty_key)
				++p;

			return p;
		}

		// *******************************************************************************************
		iterator end()
		{
			return hash_set_lp_iterator<hash_set_lp_type>(this);
		}

		// *******************************************************************************************
		const_iterator begin() const
		{
			return const_cast<hash_set_lp_type*>(this)->begin();
		}

		// *******************************************************************************************
		const_iterator end() const
		{
			return const_cast<hash_set_lp_type*>(this)->end();
		}

		// *******************************************************************************************
		const_iterator cbegin() const
		{
			return begin();
		}

		// *******************************************************************************************
		const_iterator cend() const
		{
			return end();
		}

		// *******************************************************************************************
		bool empty()
		{
			return no_elements == 0;
		}

		// *******************************************************************************************
		size_type size()
		{
			return no_elements;
		}

		// *******************************************************************************************
		size_type allocated_size()
		{
			return allocated;
		}

		// *******************************************************************************************
		size_type max_size()
		{
			return data.max_size();
		}

		// *******************************************************************************************
		size_t reserve(size_t _requested_reserve)
		{
			if (_requested_reserve < allocated)
				return allocated;

			restruct(_requested_reserve);
		}

		//	T& operator[](const key_type& k);

		// *******************************************************************************************
		template <typename InputIterator>
		void insert(InputIterator first, InputIterator last)
		{
			for (auto p = first; p != last; ++p)
				insert_fast(*p);
		}

		// *******************************************************************************************
		void insert(std::initializer_list<key_type> il)
		{
			for (auto& x : il)
				insert_fast(x);
		}

		// *******************************************************************************************
		std::pair<iterator, bool> insert(const key_type& x)
		{
			if (no_elements >= size_when_restruct)
				restruct(allocated * 2);

			size_t h = hash(x) & allocated_mask;

			if (!compare(data[h], empty_key))
			{
				if (compare(data[h], x))
					return std::make_pair(iterator(this, data.begin() + h), false);

				do
				{
					h = (h + 1) & allocated_mask;

					if (compare(data[h], x))
						return std::make_pair(iterator(this, data.begin() + h), false);
				} while (data[h] != empty_key);
			}

			++no_elements;

			data[h] = x;

			return std::make_pair(iterator(this, data.begin() + h), true);
		}
			
		// *******************************************************************************************
		bool insert_fast(const key_type& x)
		{
			if (no_elements >= size_when_restruct)
				restruct(allocated * 2);

			size_t h = hash(x) & allocated_mask;

			if (!compare(data[h], empty_key))
			{
				if (compare(data[h], x))
					return false;

				do
				{
					h = (h + 1) & allocated_mask;

					if (compare(data[h], x))
						return false;
				} while (data[h] != empty_key);
			}

			++no_elements;

			data[h] = x;

			return true;
		}

		// *******************************************************************************************
		iterator find(const key_type& key)
		{
			// !!! For unknown reasons turning off inline here speeds things up significantly, but only in "find", it does not affect "check"
			size_t pos = _find_noinline(key, hash(key) & allocated_mask);

			if (compare(data[pos], key))
				return iterator(this, data.begin() + pos);
			else
				return iterator(this);
		}

		// *******************************************************************************************
		bool check(const key_type& key)
		{
			return _check_noinline(key, hash(key) & allocated_mask);

/*			size_t pos = _find(key, hash(key) & allocated_mask);

			return compare(data[pos], key);*/
		}

		// *******************************************************************************************
		size_type count(const key_type& key) const
		{
			size_t pos = _find(key, hash(key) & allocated_mask);

			return (size_t)compare(data[pos].key, key);
		}

		// *******************************************************************************************
		key_equal key_eq() const
		{
			return compare;
		}

		// *******************************************************************************************
		hasher hash_function() const
		{
			return hash;
		}

		// *******************************************************************************************
		void prefetch(const Key_t& key)
		{
			size_t h = hash(key) & allocated_mask;

#ifdef _WIN32
			_mm_prefetch((const char*)(data.data() + h), _MM_HINT_T0);
#else
			__builtin_prefetch(data.data() + h);
#endif
		}

	private:
		// *******************************************************************************************
		void construct(size_t _init_reserved, double _max_fill_factor)
		{
			max_fill_factor = _max_fill_factor;

			if (max_fill_factor > 0.99)
				max_fill_factor = 0.99;
			else if (max_fill_factor < 0.01)
				max_fill_factor = 0.1;

			_reserve(_init_reserved);

			size_when_restruct = (size_t)(allocated * max_fill_factor);
			if (size_when_restruct == allocated)
				--size_when_restruct;
		}

		// *******************************************************************************************
		void _reserve(size_t requested_allocated)
		{
			allocated = requested_allocated;

			if (allocated < 8)
				allocated = 8;

			// Round up to the power of 2
			if ((allocated & (allocated - 1)))
			{
				while ((allocated & (allocated - 1)))
					allocated &= allocated - 1;
				allocated *= 2;
			}

			allocated_mask = allocated - 1ull;
			size_when_restruct = (size_t)(allocated * max_fill_factor);
			if (size_when_restruct == allocated)
				--size_when_restruct;

			data.resize(allocated, empty_key);

			no_elements = 0;
		}

		// *******************************************************************************************
		void restruct(size_t new_allocated)
		{
			std::vector<key_type> old_data;
			old_data = move(data);
			size_t old_allocated = allocated;

			_reserve(new_allocated);

			for (size_t i = 0; i < old_allocated; ++i)
				if (!compare(old_data[i], empty_key))
					insert_fast(old_data[i]);
		}

		// *******************************************************************************************
		size_t _find(const Key_t key, size_t pos)
		{
			if (compare(data[pos], empty_key) || compare(data[pos], key))
				return pos;

			for (++pos; pos < allocated; ++pos)
				if (compare(data[pos], empty_key) || compare(data[pos], key))
					return pos;

			for (pos = 0; pos < allocated; ++pos)
				if (compare(data[pos], empty_key) || compare(data[pos], key))
					return pos;

			return pos;
		}

		// *******************************************************************************************
		NO_INLINE size_t _find_noinline(const Key_t key, size_t pos)
		{
			if (compare(data[pos], empty_key) || compare(data[pos], key))
				return pos;

			for (++pos; pos < allocated; ++pos)
				if (compare(data[pos], empty_key) || compare(data[pos], key))
					return pos;

			for (pos = 0; pos < allocated; ++pos)
				if (compare(data[pos], empty_key) || compare(data[pos], key))
					return pos;

			return pos;
		}

		// *******************************************************************************************
		bool _check(const Key_t key, size_t pos)
		{
			if (compare(data[pos], empty_key))
				return false;
			if (compare(data[pos], key))
				return true;

			for (++pos; pos < allocated; ++pos)
			{
				if (compare(data[pos], empty_key))
					return false;
				if (compare(data[pos], key))
					return true;
			}

			for (pos = 0; pos < allocated; ++pos)
			{
				if (compare(data[pos], empty_key))
					return false;
				if (compare(data[pos], key))
					return true;
			}

			return false;
		}

		// *******************************************************************************************
		NO_INLINE bool _check_noinline(const Key_t key, size_t pos)
		{
			if (compare(data[pos], empty_key))
				return false;
			if (compare(data[pos], key))
				return true;

			for (++pos; pos < allocated; ++pos)
			{
				if (compare(data[pos], empty_key))
					return false;
				if (compare(data[pos], key))
					return true;
			}

			for (pos = 0; pos < allocated; ++pos)
			{
				if (compare(data[pos], empty_key))
					return false;
				if (compare(data[pos], key))
					return true;
			}

			return false;
		}
};

// EOF
#endif