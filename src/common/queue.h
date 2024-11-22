#ifndef _QUEUE_H
#define _QUEUE_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include <map>
#include <mutex>
#include <condition_variable>

using namespace std;

// *******************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CBoundedQueue
{
	typedef list<pair<T, size_t>> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;
	size_t current_cost;
	size_t max_cost;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;
	condition_variable cv_queue_full;

public:
	typename queue_t::iterator q_it;

	// *******************************************************************************************
	CBoundedQueue(const int _n_producers, const size_t _max_cost)
	{
		max_cost = _max_cost;

		Restart(_n_producers);
	};

	// *******************************************************************************************
	~CBoundedQueue()
	{};

	// *******************************************************************************************
	void Restart(const int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
		current_cost = 0;
	}

	// *******************************************************************************************
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *******************************************************************************************
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *******************************************************************************************
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Push(const T data, const size_t cost)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_full.wait(lck, [this] {return current_cost < max_cost; });

		bool was_empty = n_elements == 0;
		q.emplace_back(data, cost);
		++n_elements;
		current_cost += cost;

		if (was_empty)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Emplace(T&& data, const size_t cost)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_full.wait(lck, [this] {return current_cost < max_cost; });

		bool was_empty = n_elements == 0;
		q.emplace_back(move(data), cost);
		++n_elements;
		current_cost += cost;

		if (was_empty)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	bool Pop(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_producers; });

		if (n_elements == 0)
			return false;

		data = move(q.front().first);
		size_t cost = q.front().second;

		q.pop_front();
		--n_elements;
		current_cost -= cost;

		if (n_elements == 0)
			cv_queue_empty.notify_all();

		cv_queue_full.notify_all();

		return true;
	}

	// *******************************************************************************************
	pair<uint32_t, size_t> GetSize()
	{
		unique_lock<mutex> lck(mtx);

		return make_pair(n_elements, current_cost);
	}
};

// *******************************************************************************************
template<typename T> class CBoundedPQueue
{
	typedef multimap<pair<size_t, size_t>, T> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;
	size_t current_cost;
	size_t max_cost;

	mutable mutex mtx;							
	condition_variable cv_queue_empty;
	condition_variable cv_queue_full;

public:
	typename queue_t::iterator q_it;

	enum class result_t { empty, completed, normal };

	// *******************************************************************************************
	CBoundedPQueue(const int _n_producers, const size_t _max_cost)
	{
		current_cost = 0;
		max_cost = _max_cost;

		Restart(_n_producers);
	};

	// *******************************************************************************************
	~CBoundedPQueue()
	{};

	// *******************************************************************************************
	void Restart(const int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
	}

	// *******************************************************************************************
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *******************************************************************************************
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *******************************************************************************************
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Push(const T data, const size_t priority, const size_t cost)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_full.wait(lck, [this] {return current_cost < max_cost; });

		bool was_empty = n_elements == 0;
		q.emplace(make_pair<priority, cost>, data);
		++n_elements;
		current_cost += cost;

		if (was_empty)
			cv_queue_empty.notify_all();
//		cv_queue_empty.notify_one();
	}

	// *******************************************************************************************
	void Emplace(T&& data, const size_t priority, const size_t cost)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_full.wait(lck, [this] {return current_cost < max_cost; });

		bool was_empty = n_elements == 0;
		q.emplace(make_pair(priority, cost), move(data));
		++n_elements;
		current_cost += cost;

		if (was_empty)
			cv_queue_empty.notify_all();
//		cv_queue_empty.notify_one();
	}

	// *******************************************************************************************
	void EmplaceNoLock(T&& data, const size_t priority, const size_t cost)
	{
//		unique_lock<mutex> lck(mtx);
//		cv_queue_full.wait(lck, [this] {return current_cost < max_cost; });

		bool was_empty = n_elements == 0;
		q.emplace(make_pair(priority, cost), move(data));
		++n_elements;
		current_cost += cost;

		if (was_empty)
			cv_queue_empty.notify_all();
//		cv_queue_empty.notify_one();
	}

	// *******************************************************************************************
	void EmplaceManyNoCost(T&& data, const size_t priority, size_t n_items)
	{
		unique_lock<mutex> lck(mtx);

//		bool was_empty = n_elements == 0;
		for(size_t i = 0; i < n_items; ++i)
			q.emplace(make_pair(priority, 0), move(data));
		n_elements += n_items;

		cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
//	bool PopLarge(T& data)
	CBoundedPQueue::result_t PopLarge(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_producers; });

		if (n_elements == 0)
		{
			return n_producers ? result_t::empty : result_t::completed;
		}
//			return false;

//		data = move(q.rbegin()->second);

		data.swap(q.rbegin()->second);

		size_t cost = q.rbegin()->first.second;

		//		q.pop();
		q.erase(--q.end());
		--n_elements;
		current_cost -= cost;

		if (n_elements == 0)
			cv_queue_empty.notify_all();

		cv_queue_full.notify_all();

//		return true;
		return result_t::normal;
	}

	// *******************************************************************************************
	bool PopSmall(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_producers; });

		if (n_elements == 0)
			return false;

		data = move(q.begin()->second);
		size_t cost = q.begin()->first.second;

		q.erase(q.begin());
		--n_elements;
		current_cost -= cost;

		if (n_elements == 0)
			cv_queue_empty.notify_all();

		cv_queue_full.notify_all();

		return true;
	}

	// *******************************************************************************************
	pair<uint32_t, size_t> GetSize()
	{
		unique_lock<mutex> lck(mtx);

		return make_pair(n_elements, current_cost);
	}
};

// *******************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CPriorityQueue
{
	typedef map<size_t, T> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;
	size_t current_priority;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	typename queue_t::iterator q_it;

	// *******************************************************************************************
	CPriorityQueue(const int _n_producers)
	{
		current_priority = 0;

		Restart(_n_producers);
	};

	// *******************************************************************************************
	~CPriorityQueue()
	{};

	// *******************************************************************************************
	void Restart(const int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
		current_priority = 0;
	}

	// *******************************************************************************************
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *******************************************************************************************
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *******************************************************************************************
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Push(const T data, const size_t priority)
	{
		unique_lock<mutex> lck(mtx);

		bool was_empty = n_elements == 0;
		q.emplace(priority, data);
		++n_elements;

//		if (was_empty)
		cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	void Emplace(const size_t priority, T&& data)
	{
		unique_lock<mutex> lck(mtx);

//		bool was_empty = n_elements == 0;
		q.emplace(priority, move(data));
		++n_elements;

		cv_queue_empty.notify_all();
	}

	// *******************************************************************************************
	bool Pop(T& data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this] {return (!this->q.empty() && current_priority == q.begin()->first) || !this->n_producers; });

		if (n_elements == 0)
			return false;

		data = move(q.begin()->second);

		q.erase(q.begin());
		--n_elements;
		++current_priority;

		cv_queue_empty.notify_all();

		return true;
	}

	// *******************************************************************************************
	size_t GetSize()
	{
		unique_lock<mutex> lck(mtx);

		return n_elements;
	}
};

// EOF
#endif