#ifndef _AGC_COMPRESSOR_H
#define _AGC_COMPRESSOR_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "../common/agc_basic.h"
#include "../core/genome_io.h"
#include "../core/hs.h"
#include "../core/kmer.h"
#include "../common/utils.h"
#include "../core/utils_adv.h"

#include <list>
#include <set>
#include <map>
#include <future>

// *******************************************************************************************
class CBufferedSegPart
{
	struct seg_part_t {
		uint64_t kmer1;
		uint64_t kmer2;
		string sample_name;
		string contig_name;
		contig_t seg_data;
		bool is_rev_comp;
		uint32_t seg_part_no;

		seg_part_t(const uint64_t _kmer1, const uint64_t _kmer2,
			const string& _sample_name, const string& _contig_name, contig_t& _seg_data, bool _is_rev_comp, uint32_t _seg_part_no) :
			kmer1(_kmer1),
			kmer2(_kmer2),
			sample_name(_sample_name),
			contig_name(_contig_name),
			seg_data(move(_seg_data)),
			is_rev_comp(_is_rev_comp),
			seg_part_no(_seg_part_no)
		{};

		seg_part_t() :
			kmer1(~0ull),
			kmer2(~0ull),
			sample_name{},
			contig_name{},
			seg_data{},
			is_rev_comp(false),
			seg_part_no{ 0 }
		{};

		seg_part_t(const seg_part_t& rhs)
		{
			kmer1 = rhs.kmer1;
			kmer2 = rhs.kmer2;
			sample_name = rhs.sample_name;
			contig_name = rhs.contig_name;
			seg_data = rhs.seg_data;
			is_rev_comp = rhs.is_rev_comp;
			seg_part_no = rhs.seg_part_no;
		}

		seg_part_t(seg_part_t&& rhs) noexcept
		{
			kmer1 = rhs.kmer1;
			kmer2 = rhs.kmer2;
			sample_name = move(rhs.sample_name);
			contig_name = move(rhs.contig_name);
			seg_data = move(rhs.seg_data);
			is_rev_comp = rhs.is_rev_comp;
			seg_part_no = rhs.seg_part_no;
		}

		seg_part_t& operator=(const seg_part_t& rhs) {
			if (this == &rhs)
				return *this;

			kmer1 = rhs.kmer1;
			kmer2 = rhs.kmer2;
			sample_name = rhs.sample_name;
			contig_name = rhs.contig_name;
			seg_data = rhs.seg_data;
			is_rev_comp = rhs.is_rev_comp;
			seg_part_no = rhs.seg_part_no;

			return *this;
		}

		seg_part_t& operator=(seg_part_t&& rhs) noexcept
		{
			if (this == &rhs)
				return *this;

			kmer1 = rhs.kmer1;
			kmer2 = rhs.kmer2;
			sample_name = move(rhs.sample_name);
			contig_name = move(rhs.contig_name);
			seg_data = move(rhs.seg_data);
			is_rev_comp = rhs.is_rev_comp;
			seg_part_no = rhs.seg_part_no;

			return *this;
		}

		bool operator<(const struct seg_part_t& x) const
		{
			if (sample_name != x.sample_name)
				return sample_name < x.sample_name;
			if (contig_name != x.contig_name)
				return contig_name < x.contig_name;
			return seg_part_no < x.seg_part_no;
		}
	};

public:
	// *******************************************************************************************
	struct kk_seg_part_t {
		uint64_t kmer1;
		uint64_t kmer2;
		string sample_name;
		string contig_name;
		contig_t seg_data;
		bool is_rev_comp;
		uint32_t seg_part_no;

		kk_seg_part_t(const uint64_t _kmer1, const uint64_t _kmer2, const string& _sample_name, const string& _contig_name, contig_t& _seg_data, bool _is_rev_comp, uint32_t _seg_part_no) :
			kmer1(_kmer1),
			kmer2(_kmer2),
			sample_name(_sample_name),
			contig_name(_contig_name),
			seg_data(move(_seg_data)),
			is_rev_comp(_is_rev_comp),
			seg_part_no(_seg_part_no)
		{};

		kk_seg_part_t() :
			kmer1{},
			kmer2{},
			sample_name{},
			contig_name{},
			seg_data{},
			is_rev_comp(false),
			seg_part_no{ 0 }
		{};

		kk_seg_part_t(kk_seg_part_t&&) = default;
		kk_seg_part_t& operator=(const kk_seg_part_t&) = default;
		kk_seg_part_t(const kk_seg_part_t&) = default;

		bool operator<(const struct kk_seg_part_t& x) const
		{
			if (sample_name != x.sample_name)
				return sample_name < x.sample_name;
			if (contig_name != x.contig_name)
				return contig_name < x.contig_name;
			return seg_part_no < x.seg_part_no;
		}
	};

private:
	// *******************************************************************************************
	struct list_seg_part_t {
		mutex mtx;

		vector<seg_part_t> l_seg_part;
		size_t virt_begin = 0;

		list_seg_part_t() = default;
		list_seg_part_t(const list_seg_part_t& x)
		{
			l_seg_part = x.l_seg_part;
			virt_begin = x.virt_begin;
		};

		list_seg_part_t(list_seg_part_t&& x) noexcept
		{
			l_seg_part = move(x.l_seg_part);
			virt_begin = x.virt_begin;
		};

		~list_seg_part_t() = default;

		list_seg_part_t(const seg_part_t& seg_part)
		{
			l_seg_part.push_back(seg_part);
			virt_begin = 0;
		}

		list_seg_part_t& operator=(const list_seg_part_t& x)
		{
			if (&x != this)
			{
				l_seg_part = x.l_seg_part;
				virt_begin = x.virt_begin;
			}

			return *this;
		}

		list_seg_part_t& operator=(list_seg_part_t&& x) noexcept
		{
			if (&x != this)
			{
				l_seg_part = move(x.l_seg_part);
				virt_begin = x.virt_begin;
			}

			return *this;
		}

		void append(seg_part_t&& seg_part)
		{
			lock_guard<mutex> lck(mtx);
			l_seg_part.emplace_back(move(seg_part));
		}

		void emplace(uint64_t kmer1, uint64_t kmer2, const string& sample_name, const string& contig_name, contig_t& seg_data, bool is_rev_comp, uint32_t seg_part_no)
		{
			lock_guard<mutex> lck(mtx);
			l_seg_part.emplace_back(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no);
		}

		void append_no_lock(seg_part_t& seg_part)
		{
			l_seg_part.emplace_back(move(seg_part));
		}

		void sort()
		{
			std::sort(l_seg_part.begin(), l_seg_part.end());
		}

		void clear()
		{
			l_seg_part.clear();
			virt_begin = 0;
		}

		bool empty()
		{
			return virt_begin >= l_seg_part.size();
		}

		bool pop(seg_part_t& seg_part)
		{
			if (virt_begin >= l_seg_part.size())
			{
				virt_begin = 0;
				l_seg_part.clear();

				return false;
			}

			seg_part = move(l_seg_part[virt_begin]);
			++virt_begin;

			return true;
		}

		bool pop(uint64_t& kmer1, uint64_t& kmer2, string& sample_name, string& contig_name, contig_t& seg_data, bool& is_rev_comp, uint32_t& seg_part_no)
		{
			if (virt_begin >= l_seg_part.size())
			{
				virt_begin = 0;
				l_seg_part.clear();

				return false;
			}

			auto& x = l_seg_part[virt_begin];

			kmer1 = x.kmer1;
			kmer2 = x.kmer2;
			sample_name = move(x.sample_name);
			contig_name = move(x.contig_name);
			seg_data = move(x.seg_data);
			is_rev_comp = x.is_rev_comp;
			seg_part_no = x.seg_part_no;

			++virt_begin;

			return true;
		}

		uint32_t size()
		{
			return (uint32_t)l_seg_part.size();
		}
	};

	vector<list_seg_part_t> vl_seg_part;

	set<kk_seg_part_t> s_seg_part;
	mutex mtx;

	atomic<int32_t> a_v_part_id;

public:
	static const int32_t part_id_step = 1;

	CBufferedSegPart(uint32_t no_raw_groups)
	{
		vl_seg_part.resize(no_raw_groups);
	}

	~CBufferedSegPart() = default;

	void resize(uint32_t no_groups)
	{
		vl_seg_part.resize(no_groups);
	}

	void add_known(uint32_t group_id, uint64_t kmer1, uint64_t kmer2, const string& sample_name, const string& contig_name, contig_t&& seg_data, bool is_rev_comp, uint32_t seg_part_no)
	{
		// !!! TODO: use move() here?
		vl_seg_part[group_id].emplace(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no);		// internal mutex
	}

	void add_new(uint64_t kmer1, uint64_t kmer2, const string& sample_name, const string& contig_name, contig_t& seg_data, bool is_rev_comp, uint32_t seg_part_no)
	{
		lock_guard<mutex> lck(mtx);
		// !!! TODO: use move() here?
		s_seg_part.emplace(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no);
	}

	void sort_known(uint32_t nt)
	{
		lock_guard<mutex> lck(mtx);
		const uint32_t min_parts_in_job = 16;

		uint64_t n_seg_part = vl_seg_part.size();

		nt = clamp<uint32_t>(n_seg_part / min_parts_in_job, 1u, nt);

		auto vl_seg_part_begin = vl_seg_part.begin();

		atomic<uint64_t> seg_part_id = 0;
		uint64_t job_step = n_seg_part / (16 * nt);
		if (job_step < 16)
			job_step = 16;

		auto job = [&seg_part_id, job_step, n_seg_part, vl_seg_part_begin] {
			while (true)
			{
				uint64_t j_from = seg_part_id.fetch_add(job_step);

				if (j_from >= n_seg_part)
					break;

				uint64_t j_to = min(j_from + job_step, n_seg_part);

				auto p = vl_seg_part_begin + j_from;

				for (uint64_t j = j_from; j < j_to; ++j, ++p)
					p->sort();
			}
			};

		vector<future<void>> v_fut;

		v_fut.reserve(nt);

		for (uint64_t i = 0; i < nt - 1; ++i)
			v_fut.emplace_back(async(job));

		job();

		for (auto& f : v_fut)
			f.wait();
	}

	const set<kk_seg_part_t>& get_seg_parts() const
	{
		return s_seg_part;
	}

	uint32_t process_new()
	{
		lock_guard<mutex> lck(mtx);

		map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
		uint32_t group_id = (uint32_t)vl_seg_part.size();

		// Assign group ids to new segments
		for (const auto& x : s_seg_part)
		{
			auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

			if (p == m_kmers.end())
				m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
		}

		uint32_t no_new = group_id - (uint32_t)vl_seg_part.size();

		if (vl_seg_part.capacity() < group_id)
			vl_seg_part.reserve((uint64_t)(group_id * 1.2));
		vl_seg_part.resize(group_id);

		for (auto& x : s_seg_part)
		{
			add_known(m_kmers[make_pair(x.kmer1, x.kmer2)], x.kmer1, x.kmer2,
				x.sample_name, x.contig_name, const_cast<contig_t&&>(x.seg_data), x.is_rev_comp, x.seg_part_no);
		}

		s_seg_part.clear();

		return no_new;
	}

	void distribute_segments(uint32_t src_id, uint32_t dest_id_from, uint32_t dest_id_to)
	{
		uint32_t no_in_src = vl_seg_part[src_id].size();
		uint32_t dest_id_curr = dest_id_from;

		seg_part_t seg_part;

		for (uint32_t i = 0; i < no_in_src; ++i)
		{
			if (dest_id_curr != src_id)
			{
				vl_seg_part[src_id].pop(seg_part);
				vl_seg_part[dest_id_curr].append_no_lock(seg_part);
			}

			if (++dest_id_curr == dest_id_to)
				dest_id_curr = dest_id_from;
		}
	}

	atomic<uint64_t> clear_idx{};
	const uint64_t job_clear_step = 64;

	void clear_job()
	{
		uint64_t n_seg_part = vl_seg_part.size();
		auto vl_seg_part_begin = vl_seg_part.begin();

		while (true)
		{
			uint64_t loc_idx = clear_idx.fetch_add(job_clear_step);
			uint64_t upp_idx = min<uint64_t>(loc_idx + job_clear_step, n_seg_part);

			auto p = vl_seg_part_begin + loc_idx;

			for (; loc_idx < upp_idx; ++loc_idx, ++p)
				//					vl_seg_part[loc_idx].clear();
				p->clear();

			if (loc_idx >= n_seg_part)
				break;
		}
	}

	void clear(uint32_t nt)
	{
		lock_guard<mutex> lck(mtx);
		const uint32_t min_parts_in_job = 512;

		uint64_t n_seg_part = vl_seg_part.size();

		nt = clamp<uint32_t>(n_seg_part / min_parts_in_job, 1u, nt);

		uint64_t job_step = n_seg_part / (64 * nt);
		if (job_step < 64)
			job_step = 64;

		atomic<uint64_t> idx = 0;

		auto vl_seg_part_begin = vl_seg_part.begin();

		auto job = [&idx, job_step, n_seg_part, vl_seg_part_begin] {
			while (true)
			{
				uint64_t loc_idx = idx.fetch_add(job_step);

				if (loc_idx >= n_seg_part)
					break;

				uint64_t upp_idx = min<uint64_t>(loc_idx + job_step, n_seg_part);

				auto p = vl_seg_part_begin + loc_idx;

				for (; loc_idx < upp_idx; ++loc_idx, ++p)
					p->clear();
			}
			};

		vector<future<void>> v_fut;

		v_fut.reserve(nt);

		for (uint64_t i = 0; i < nt - 1; ++i)
			v_fut.emplace_back(async(job));

		s_seg_part.clear();
		job();

		for (auto& f : v_fut)
			f.wait();
	}

	void restart_read_vec()
	{
		lock_guard<mutex> lck(mtx);

		a_v_part_id = (int32_t)vl_seg_part.size() - 1;
	}

	int get_vec_id()
	{
//		return a_v_part_id.fetch_sub(1);
		return a_v_part_id.fetch_sub(part_id_step);
	}

	int get_no_parts()
	{
		return (int)vl_seg_part.size();
	}

	bool is_empty_part(int group_id)
	{
		return group_id < 0 || vl_seg_part[group_id].empty();
	}

	bool get_part(int group_id, uint64_t &kmer1, uint64_t &kmer2, string& sample_name, string& contig_name, contig_t& seg_data, bool& is_rev_comp, uint32_t& seg_part_no)
	{
		return vl_seg_part[group_id].pop(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no);
	}
};

// *******************************************************************************************
// Class supporting decompression and compresion of AGC files
class CAGCCompressor : public CAGCBasic
{
	enum class splitter_orientation_t { unknown, direct, rev_comp };

	using hash_set_t = hash_set_lp<uint64_t, equal_to<uint64_t>, MurMur64Hash>;

	// *******************************************************************************************
	struct hash_pair {
		template <class T1, class T2>
		size_t operator()(const pair<T1, T2>& p) const
		{
			auto hash1 = hash<T1>{}(p.first);
			auto hash2 = hash<T2>{}(p.second);
			return hash1 ^ hash2;
		}
	};

	// *******************************************************************************************
	struct splitter_desc_t
	{
		splitter_orientation_t orientation;
		int32_t as_front_id;
		int32_t as_back_id;

		splitter_desc_t(const splitter_orientation_t _orientation = splitter_orientation_t::unknown, const int32_t _as_front_id = -1, const int32_t _as_back_id = -1) :
			orientation(_orientation), as_front_id(_as_front_id), as_back_id(_as_back_id)
		{}
	};

	// *******************************************************************************************
	class kmer_filter_t
	{
		uint64_t thr = 0ull;
		MurMur64Hash mmh;
		uint64_t rnd = 0xD73F8BF11046C40Eull;

	public:
		kmer_filter_t(double fraction = 0.0)
		{
			reset(fraction);
		}

		void reset(double fraction)
		{
			if (fraction == 0.0)
				thr = 0ull;
			else
				thr = (uint64_t)(((double)~0ull) * fraction);
		}

		explicit operator bool() const
		{
			return thr != 0ull;
		}

		bool operator()(uint64_t x) const 
		{
			return (mmh(x) ^ rnd) < thr;
		}
	};

//	using my_barrier = CBarrier;
//	using my_barrier = CAtomicBarrier;
	using my_barrier = CAtomicBarrierWithIncrementing;
//	using my_barrier = barrier<>;

	const pair<uint64_t, uint64_t> pk_empty = make_pair(~0ull, ~0ull);

	shared_mutex seg_map_mtx;
	shared_mutex seg_vec_mtx;

	string out_archive_name;
	size_t no_samples_in_archive;

	vector<string> v_file_names;

	bool concatenated_genomes;
	bool adaptive_compression;

	shared_ptr<CArchive> out_archive;															// internal mutexes

	vector<uint64_t> v_candidate_kmers;
	vector<uint64_t> v_duplicated_kmers;
	uint64_t v_candidate_kmers_offset = 0;

	hash_set_t hs_splitters{ ~0ull, 16ull, 0.4, equal_to<uint64_t>{}, MurMur64Hash{} };			// only reads after init - no need to lock
	bloom_set_t bloom_splitters;

	unordered_map<pair<uint64_t, uint64_t>, int32_t, MurMurPair64Hash> map_segments;			// shared_mutex (seg_map_mtx)
	unordered_map<uint64_t, vector<uint64_t>, MurMur64Hash> map_segments_terminators;			// shared_mutex (seg_map_mtx)
	vector<shared_ptr<CSegment>> v_segments;													// shared_mutex to vector (seg_vec_mtx) + internal mutexes in stored objects

	mutex mtx_fallback_map;
	kmer_filter_t fallback_filter;
	double fallback_frac = 0.0;
	vector<vector<array<uint64_t, 4>>> vv_fallback_minimizers;
	unordered_map<uint64_t, vector<pair<uint64_t, uint64_t>>> map_fallback_minimizers;			// mtx_fallback_map

	uint32_t no_segments;
	atomic<uint32_t> id_segment = 0;

	const size_t contig_part_size = 512 << 10;

	CBufferedSegPart buffered_seg_part{ no_raw_groups };

	atomic<size_t> processed_bases;
	atomic<uint64_t> a_part_id;
	uint32_t processed_samples{ 0 };

	vector<vector<uint64_t>> vv_splitters;

	vector<tuple<string, string, contig_t>> v_raw_contigs;
	mutex mtx_raw_contigs;

	enum class contig_processing_stage_t {unknown, all_contigs, new_splitters, hard_contigs, registration};

	using task_t = tuple<contig_processing_stage_t, string, string, contig_t>;
	
	shared_ptr<CBoundedPQueue<task_t>> pq_contigs_desc;											// internal mutexes
	shared_ptr<CBoundedPQueue<task_t>> pq_contigs_desc_aux;										// internal mutexes
	shared_ptr<CBoundedPQueue<task_t>> pq_contigs_desc_working;									// internal mutexes

	unique_ptr<CBoundedQueue<tuple<string, string, contig_t>>> q_contigs_desc;					// internal mutexes
	unique_ptr<CBoundedPQueue<contig_t>> pq_contigs_raw;										// internal mutexes
	unique_ptr<CBoundedQueue<contig_t>> q_contigs_data;											// internal mutexes

	vector<ZSTD_CCtx*> v_cctx;
	vector<ZSTD_DCtx*> v_dctx;
	ZSTD_DCtx* zstd_dctx_for_fallback = nullptr;

	bool compress_contig(contig_processing_stage_t contig_processing_stage, string sample_name, string id, contig_t& contig, 
		ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx, uint32_t thread_id, my_barrier &bar);
	pair_segment_desc_t add_segment(const string &sample_name, const string &contig_name, uint32_t seg_part_no,
		contig_t &&segment, CKmer kmer_front, CKmer kmer_back, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx, uint32_t thread_id, my_barrier& bar);
	void register_segments(uint32_t n_t);
	void store_segments(ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx);

	pair<pair<uint64_t, uint64_t>, bool> find_cand_segment_with_one_splitter(CKmer kmer, contig_t& segment_dir, contig_t& segment_rc, ZSTD_DCtx* zstd_dctx, my_barrier& bar);
	pair<uint64_t, uint32_t> find_cand_segment_with_missing_middle_splitter(CKmer kmer_front, CKmer kmer_back, contig_t& segment_dir, contig_t& segment_rc, ZSTD_DCtx* zstd_dctx, my_barrier& bar);
	pair<pair<uint64_t, uint64_t>, bool> find_cand_segment_using_fallback_minimizers(contig_t& segment, uint64_t max_val);

	contig_t get_part(const contig_t& contig, uint64_t pos, uint64_t len);
	void preprocess_raw_contig(contig_t& ctg);
	void find_new_splitters(contig_t& ctg, uint32_t thread_id);

	void add_fallback_kmers(vector<uint64_t>::iterator first, vector<uint64_t>::iterator last);
	void add_fallback_mapping(uint64_t splitter1, uint64_t splitter2, vector<pair<uint64_t, bool>>& cand_fallback_kmers);
	void add_fallback_mapping(uint64_t splitter1, uint64_t splitter2, uint64_t kmer, bool is_dir_oriented);
	void add_fallback_mapping(uint64_t splitter1, uint64_t splitter2, const contig_t &segment);

	// *******************************************************************************************
	void append(vector<uint8_t>& data, uint32_t num)
	{
		for (int i = 0; i < 4; ++i)
		{
			data.emplace_back(num & 0xff);
			num >>= 8;
		}
	}

	// *******************************************************************************************
	void append64(vector<uint8_t>& data, uint64_t num)
	{
		for (int i = 0; i < 8; ++i)
		{
			data.emplace_back(num & 0xff);
			num >>= 8;
		}
	}

	// *******************************************************************************************
	void append(vector<uint8_t>& data, const string& str)
	{
		data.insert(data.end(), str.begin(), str.end());
		data.emplace_back(0);
	}

	// *******************************************************************************************
	bool close_compression(const uint32_t no_threads);

	void prepare_compressing_stuctures(const uint32_t n_t);
	void release_compressing_stuctures();

	void compressing_stage1_job(uint32_t thread_id, uint32_t n_t);
	void compressing_stage2_job(uint32_t thread_id, uint32_t n_t);

	void start_compressing_threads(vector<thread> &v_threads, my_barrier &bar, const uint32_t n_t);
	void start_finalizing_threads(vector<thread>& v_threads, const uint32_t n_t);
	void start_splitter_finding_threads(vector<thread>& v_threads, const uint32_t n_t, const vector<uint64_t>::iterator v_begin, const vector<uint64_t>::iterator v_end, vector<vector<uint64_t>>& v_splitters);
	void start_kmer_collecting_threads(vector<thread>& v_threads, const uint32_t n_t, vector<uint64_t>& v_kmers, const size_t extra_items);

	void store_metadata_impl_v1(uint32_t no_threads);
	void store_metadata_impl_v2(uint32_t no_threads);
	void store_metadata_impl_v3(uint32_t no_threads);

	void store_metadata(uint32_t no_threads);
	void appending_init();
	bool determine_splitters(const string& reference_file_name, const size_t segment_size, const uint32_t no_threads);
	bool count_kmers(vector<pair<string, vector<uint8_t>>>& v_contig_data, const uint32_t no_threads);

	void remove_non_singletons(vector<uint64_t>& vec, size_t virtual_begin);
	void remove_non_singletons(vector<uint64_t>& vec, vector<uint64_t>& v_duplicated, size_t virtual_begin);

	void enumerate_kmers(contig_t& ctg, vector<uint64_t> &vec);
	void find_splitters_in_contig(contig_t& ctg, const vector<uint64_t>::iterator v_begin, const vector<uint64_t>::iterator v_end, vector<uint64_t>& v_splitters, vector<array<uint64_t, 4>>&v_fallbacks);

	void store_file_type_info();

	void build_candidate_kmers_from_archive(const uint32_t n_t);

public:
	CAGCCompressor();
	~CAGCCompressor();

	bool Create(const string& _file_name, const uint32_t _pack_cardinality, const uint32_t _kmer_length, const string& reference_file_name, const uint32_t _segment_size,
		const uint32_t _min_match_len, const bool _concatenated_genomes, const bool _adaptive_compression, const uint32_t _verbosity, const uint32_t _no_threads, double _fallback_frac);
	bool Append(const string& _in_archive_fn, const string& _out_archive_fn, const uint32_t _verbosity, const bool _prefetch_archive, const bool _concatenated_genomes, const bool _adaptive_compression,
		const uint32_t no_threads, double _fallback_frac);

	void AddCmdLine(const string& cmd_line);

	bool Close(const uint32_t no_threads = 1);

	bool AddSampleFiles(vector<pair<string, string>> _v_sample_file_name, const uint32_t _no_threads);
};

// EOF
#endif