#ifndef _AGC_COMPRESSOR_H
#define _AGC_COMPRESSOR_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include "../core/agc_basic.h"
#include "../core/genome_io.h"
#include "../core/hs.h"
#include "../core/kmer.h"
#include "../core/utils.h"
#include <list>
#include <set>
#include <map>

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

		seg_part_t(seg_part_t&&) = default;
		seg_part_t& operator=(const seg_part_t&) = default;
		seg_part_t(const seg_part_t&) = default;

		bool operator<(const struct seg_part_t& x) const
		{
			if (sample_name != x.sample_name)
				return sample_name < x.sample_name;
			if (contig_name != x.contig_name)
				return contig_name < x.contig_name;
			return seg_part_no < x.seg_part_no;
		}
	};

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

	// *******************************************************************************************
	struct list_seg_part_t {
		mutex mtx;

		list<seg_part_t> l_seg_part;

		list_seg_part_t() = default;
		list_seg_part_t(const list_seg_part_t& x)
		{
			l_seg_part = x.l_seg_part;
		};
		~list_seg_part_t() = default;

		list_seg_part_t(const seg_part_t& seg_part)
		{
			l_seg_part.push_back(seg_part);
		}

		void append(seg_part_t&& seg_part)
		{
			lock_guard<mutex> lck(mtx);
			l_seg_part.emplace_back(move(seg_part));
		}

		void append_no_lock(seg_part_t& seg_part)
		{
			l_seg_part.emplace_back(move(seg_part));
		}

		void sort()
		{
			l_seg_part.sort();
		}

		void clear()
		{
			l_seg_part.clear();
		}

		bool empty()
		{
			return l_seg_part.empty();
		}

		bool pop(seg_part_t& seg_part)
		{
			if (l_seg_part.empty())
				return false;

			seg_part = move(l_seg_part.front());

			l_seg_part.pop_front();
			return true;
		}

		uint32_t size()
		{
			return (uint32_t) l_seg_part.size();
		}
	};

	vector<list_seg_part_t> vl_seg_part;

	set<kk_seg_part_t> s_seg_part;
	mutex mtx;

	atomic<uint32_t> a_v_part_id;
	atomic<uint32_t> a_l_part_id;

public:
	CBufferedSegPart(uint32_t no_raw_groups)
	{
		vl_seg_part.resize(no_raw_groups);
	}

	~CBufferedSegPart() = default;

	void resize(uint32_t no_groups)
	{
		vl_seg_part.resize(no_groups);
	}

	void add_known(uint32_t group_id, const uint64_t kmer1, const uint64_t kmer2, const string& sample_name, const string& contig_name, contig_t& seg_data, bool is_rev_comp, uint32_t seg_part_no)
	{
		vl_seg_part[group_id].append(seg_part_t(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no));		// internal mutex
	}

	void add_new(uint64_t kmer1, uint64_t kmer2, const string& sample_name, const string& contig_name, contig_t& seg_data, bool is_rev_comp, uint32_t seg_part_no)
	{
		lock_guard<mutex> lck(mtx);
		s_seg_part.emplace(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no);
	}

	void sort_known()
	{
		lock_guard<mutex> lck(mtx);

		for (auto& v : vl_seg_part)
			v.sort();
	}

	uint32_t process_new()
	{
		lock_guard<mutex> lck(mtx);

		map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
		uint32_t group_id = (uint32_t) vl_seg_part.size();

		// Assign group ids to new segments
		for (const auto& x : s_seg_part)
		{
			auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

			if (p == m_kmers.end())
				m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
		}

		uint32_t no_new = group_id - vl_seg_part.size();

		vl_seg_part.resize(group_id);
		
		for (auto& x : s_seg_part)
		{
			add_known(m_kmers[make_pair(x.kmer1, x.kmer2)], x.kmer1, x.kmer2,
				x.sample_name, x.contig_name, const_cast<contig_t&>(x.seg_data), x.is_rev_comp, x.seg_part_no);
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

	void clear()
	{
		lock_guard<mutex> lck(mtx);

		s_seg_part.clear();

		for (auto& x : vl_seg_part)
			x.clear();
	}

	void restart_read_vec()
	{
		a_v_part_id = 0;
	}

	int get_vec_id()
	{
		int part_id = (int)a_v_part_id.fetch_add(1);

		if (part_id >= (int) vl_seg_part.size())
			return -1;
		else
			return part_id;
	}

	bool get_part(uint32_t group_id, uint64_t &kmer1, uint64_t &kmer2, string& sample_name, string& contig_name, contig_t& seg_data, bool& is_rev_comp, uint32_t& seg_part_no)
	{
		if (group_id >= vl_seg_part.size() || vl_seg_part[group_id].empty())
			return false;

		seg_part_t seg_part;
		vl_seg_part[group_id].pop(seg_part);

		kmer1 = seg_part.kmer1;
		kmer2 = seg_part.kmer2;
		sample_name = move(seg_part.sample_name);
		contig_name = move(seg_part.contig_name);
		seg_data = move(seg_part.seg_data);
		is_rev_comp = seg_part.is_rev_comp;
		seg_part_no = seg_part.seg_part_no;

		return true;
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

//	using my_barrier = CBarrier;
	using my_barrier = CAtomicBarrier;

	shared_mutex seg_map_mtx;
	shared_mutex seg_vec_mtx;

	string out_archive_name;
	size_t no_samples_in_archive;

	vector<string> v_file_names;

	bool concatenated_genomes;
	uint32_t segment_size;

	shared_ptr<CArchive> out_archive;															// internal mutexes
	hash_set_t hs_splitters{ ~0ull, 16ull, 0.5, equal_to<uint64_t>{}, MurMur64Hash{} };			// only reads after init - no need to lock

	map<pair<uint64_t, uint64_t>, int32_t> map_segments;										// shared_mutex (seg_map_mtx)
	map<uint64_t, vector<uint64_t>> map_segments_terminators;									// shared_mutex (seg_map_mtx)
	vector<shared_ptr<CSegment>> v_segments;													// shared_mutex to vector (seg_vec_mtx) + internal mutexes in stored objects

	uint32_t no_segments;
	atomic<uint32_t> id_segment = 0;

	const size_t contig_part_size = 512 << 10;

	CBufferedSegPart buffered_seg_part{ no_raw_groups };
	bool reproducibility_mode = false;

	atomic<size_t> processed_bases;
	atomic<uint64_t> a_part_id;

	unique_ptr<CBoundedPQueue<tuple<string, string, contig_t>>> pq_contigs_desc;				// internal mutexes
	unique_ptr<CBoundedQueue<tuple<string, string, contig_t>>> q_contigs_desc;					// internal mutexes
	unique_ptr<CBoundedPQueue<contig_t>> pq_contigs_raw;										// internal mutexes
	unique_ptr<CBoundedQueue<contig_t>> q_contigs_data;											// internal mutexes

	bool compress_contig(string sample_name, string id, contig_t& contig, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx);
	bool compress_contig_rep(string sample_name, string id, contig_t& contig, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx);
	pair_segment_desc_t add_segment(const string &sample_name, const string &contig_name, uint32_t seg_part_no,
		contig_t segment, CKmer kmer_front, CKmer kmer_back, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx);

	pair<pair<uint64_t, uint64_t>, bool> find_cand_segment_with_one_splitter(CKmer kmer, contig_t& segment_dir, contig_t& segment_rc, ZSTD_DCtx* zstd_dctx);
	pair<uint64_t, uint32_t> find_cand_segment_with_missing_middle_splitter(CKmer kmer_front, CKmer kmer_back, contig_t& segment_dir, contig_t& segment_rc, ZSTD_DCtx* zstd_dctx);

	contig_t get_part(const contig_t& contig, uint64_t pos, uint64_t len);
	void preprocess_raw_contig(contig_t& ctg);

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
	bool close_compression(const uint32_t no_threads, const bool buffered);

	void start_compressing_threads(vector<thread> &v_threads, const uint32_t n_t);
	void start_compressing_rep_threads(vector<thread> &v_threads, my_barrier &bar, const uint32_t n_t);
	void start_finalizing_threads(vector<thread>& v_threads, const uint32_t n_t, bool buffered);
	void start_splitter_finding_threads(vector<thread>& v_threads, const uint32_t n_t, const vector<uint64_t>::iterator v_begin, const vector<uint64_t>::iterator v_end, vector<vector<uint64_t>>& v_splitters);
	void start_kmer_collecting_threads(vector<thread>& v_threads, const uint32_t n_t, vector<uint64_t>& v_kmers, const size_t extra_items);

	void store_metadata();
	void appending_init();
	bool determine_splitters(const string& reference_file_name, const size_t segment_size, const uint32_t no_threads);

	void store_file_type_info();

public:
	CAGCCompressor();
	~CAGCCompressor();

	bool Create(const string& _file_name, const uint32_t _pack_cardinality, const uint32_t _kmer_length, const string& reference_file_name, const uint32_t _segment_size,
		const uint32_t _min_match_len, const bool _concatenated_genomes, const uint32_t _verbosity, const uint32_t _no_threads);
	bool Append(const string& _in_archive_fn, const string& _out_archive_fn, const uint32_t _verbosity, const bool _prefetch_archive, const bool _concatenated_genomes);

	void AddCmdLine(const string& cmd_line);

	bool Close(const uint32_t no_threads = 1);

	bool AddSampleFile(const string& _sample_name, const string& _file_name, const uint32_t _no_threads);
	bool AddSampleFiles(vector<pair<string, string>> _v_sample_file_name, const uint32_t _no_threads);
	bool AddSampleFilesRep(vector<pair<string, string>> _v_sample_file_name, const uint32_t _no_threads);
};

// EOF
#endif