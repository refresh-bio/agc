// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include <numeric>
#include <memory>
#include <filesystem>
#include "../core/agc_compressor.h"
#include "../core/agc_decompressor.h"

#ifndef _DEBUG
#include "../../libs/raduls.h"
#endif

using namespace std;

// *******************************************************************************************
CAGCCompressor::CAGCCompressor() : CAGCBasic()
{
    no_segments = 0;

    no_samples_in_archive = 0;

    concatenated_genomes = false;
    adaptive_compression = false;
    segment_size = 0;

    m_file_type_info["producer"] = "agc";
    m_file_type_info["producer_version_major"] = to_string(AGC_VER_MAJOR);
    m_file_type_info["producer_version_minor"] = to_string(AGC_VER_MINOR);
    m_file_type_info["producer_version_build"] = AGC_VER_BUILD;
    m_file_type_info["file_version_major"] = to_string(AGC_FILE_MAJOR);
    m_file_type_info["file_version_minor"] = to_string(AGC_FILE_MINOR);
    m_file_type_info["comment"] = AGC_VERSION;

    map_segments.max_load_factor(1);
    map_segments_terminators.max_load_factor(1);
}

// *******************************************************************************************
CAGCCompressor::~CAGCCompressor()
{
    if (working_mode == working_mode_t::compression)
        close_compression(1);
    else if (working_mode == working_mode_t::appending)
        close_compression(1);
}

// *******************************************************************************************
void CAGCCompressor::store_metadata_impl_v1(uint32_t no_threads)
{
    vector<uint8_t> coll_desc_ser;

    auto collection_desc_v1 = dynamic_pointer_cast<CCollection_V1>(collection_desc);
    collection_desc_v1->serialize(coll_desc_ser, false);

    vector<uint8_t> zstd_coll_desc_ser;
    uint64_t coll_bound_size = ZSTD_compressBound(coll_desc_ser.size());

    zstd_coll_desc_ser.resize(coll_bound_size);
    //    auto c_size = ZSTD_compress(zstd_coll_desc_ser.data(), zstd_coll_desc_ser.size(), coll_desc_ser.data(), coll_desc_ser.size(), 13);
    auto c_size = ZSTD_compress(zstd_coll_desc_ser.data(), zstd_coll_desc_ser.size(), coll_desc_ser.data(), coll_desc_ser.size(), 19);

    zstd_coll_desc_ser.resize(c_size);

    auto collection_desc_id = out_archive->RegisterStream("collection-desc");
    out_archive->AddPart(collection_desc_id, zstd_coll_desc_ser, coll_desc_ser.size());
}

// *******************************************************************************************
void CAGCCompressor::store_metadata_impl_v2(uint32_t no_threads)
{
    vector<uint8_t> coll_data_main;
    vector<vector<uint8_t>> coll_data_details;
    vector<vector<uint8_t>> ztsd_coll_data_details;

    auto collection_desc_v2 = dynamic_pointer_cast<CCollection_V2>(collection_desc);

    collection_desc_v2->serialize(coll_data_main, coll_data_details, false, pack_cardinality * 5);    // Consider larger multiplier

    ztsd_coll_data_details.resize(coll_data_details.size());

    auto collection_desc_main_id = out_archive->RegisterStream("collection-main");
    auto collection_desc_details_id = out_archive->RegisterStream("collection-details");

    unique_ptr<CBoundedQueue<int32_t>> pq_coll_part = make_unique<CBoundedQueue<int32_t>>(1, 1);

    // Main part
    for (int i = 0; i < (int)coll_data_details.size(); ++i)
        pq_coll_part->Push(i, 0);

    pq_coll_part->MarkCompleted();

    // Details compression
    vector<thread> v_thread;
    v_thread.reserve(no_threads);

    for (uint32_t i = 0; i < max(1u, no_threads - 1u); ++i)
        v_thread.emplace_back([&coll_data_details, &ztsd_coll_data_details, &pq_coll_part]
            {
                ZSTD_CCtx* zstd_cctx = ZSTD_createCCtx();
                int part_id;

                while (!pq_coll_part->IsCompleted())
                {
                    if (!pq_coll_part->Pop(part_id))
                        break;

                    uint64_t coll_bound_size = ZSTD_compressBound(coll_data_details[part_id].size());

                    ztsd_coll_data_details[part_id].resize(coll_bound_size);

                    auto c_size = ZSTD_compressCCtx(zstd_cctx, ztsd_coll_data_details[part_id].data(), ztsd_coll_data_details[part_id].size(), coll_data_details[part_id].data(), coll_data_details[part_id].size(), 19);
                    ztsd_coll_data_details[part_id].resize(c_size);
                }

                ZSTD_freeCCtx(zstd_cctx);
            });

    // Main part compression
    vector<uint8_t> zstd_coll_data;
    uint64_t coll_bound_size = ZSTD_compressBound(coll_data_main.size());
    zstd_coll_data.resize(coll_bound_size);
    //        auto c_size = ZSTD_compress(zstd_coll_data.data(), zstd_coll_data.size(), coll_data_main.data(), coll_data_main.size(), 19);
    auto c_size = ZSTD_compress(zstd_coll_data.data(), zstd_coll_data.size(), coll_data_main.data(), coll_data_main.size(), 15);
    zstd_coll_data.resize(c_size);

    join_threads(v_thread);

    out_archive->AddPart(collection_desc_main_id, zstd_coll_data, coll_data_main.size());

    // Details storage
    for (uint32_t i = 0; i < coll_data_details.size(); ++i)
        out_archive->AddPart(collection_desc_details_id, ztsd_coll_data_details[i], coll_data_details[i].size());
}

// *******************************************************************************************
void CAGCCompressor::store_metadata_impl_v3(uint32_t no_threads)
{
    // Nothing to do here
}

// *******************************************************************************************
void CAGCCompressor::store_metadata(uint32_t no_threads)
{
    if (archive_version < 2000)
        store_metadata_impl_v1(no_threads);
    else if (archive_version < 3000)
        store_metadata_impl_v2(no_threads);
    else
        store_metadata_impl_v3(no_threads);

    uint64_t total_size_ref = 0;
    uint64_t total_size_delta = 0;
    uint64_t no_only_ref_segments = 0;
    uint64_t total_size_only_ref = 0;

    for (uint32_t i = 0; i < no_segments; ++i)
    {
        total_size_ref += out_archive->GetStreamPackedSize(out_archive->GetStreamId(ss_ref_name(archive_version, i)));
        total_size_delta += out_archive->GetStreamPackedSize(out_archive->GetStreamId(ss_delta_name(archive_version, i)));

        if (out_archive->GetNoParts(out_archive->GetStreamId(ss_delta_name(archive_version, i))) == 0)
        {
            ++no_only_ref_segments;
            total_size_only_ref += out_archive->GetStreamPackedSize(out_archive->GetStreamId(ss_ref_name(archive_version, i)));
        }
    }

    uint64_t total_size_raw = 0;
    
    for(uint32_t i = 0; i < no_raw_groups; ++i)
        total_size_raw += out_archive->GetStreamPackedSize(out_archive->GetStreamId(ss_delta_name(archive_version, i)));

    compression_params_t compression_params{ kmer_length, min_match_len, pack_cardinality, segment_size };

    vector<uint8_t> v_params;
    append(v_params, compression_params.kmer_length);
    append(v_params, compression_params.min_match_len);
    append(v_params, compression_params.pack_cardinality);

    if(archive_version >= 2000)
        append(v_params, compression_params.segment_size);

    auto compression_params_id = out_archive->RegisterStream("params");
    out_archive->AddPart(compression_params_id, v_params);

    vector<uint8_t> v_tmp;
    vector<uint64_t> v_hs_splitters;
    for (auto x : hs_splitters)
        v_hs_splitters.emplace_back(x);
    sort(v_hs_splitters.begin(), v_hs_splitters.end());

    for (auto x : v_hs_splitters)
        append64(v_tmp, x);
    auto splitters_id = out_archive->RegisterStream("splitters");
    out_archive->AddPart(splitters_id, v_tmp, hs_splitters.size());

    v_tmp.clear();
    for (auto& x : map_segments)
    {
        append64(v_tmp, x.first.first);
        append64(v_tmp, x.first.second);
        append(v_tmp, x.second);
    }
    auto map_segments_id = out_archive->RegisterStream("segment-splitters");
    out_archive->AddPart(map_segments_id, v_tmp, map_segments.size());

    uint32_t no_segments_one_side = 0;
    for (auto& x : map_segments)
        if (x.first.first == ~0ull || x.first.second == ~0ull)
            ++no_segments_one_side;

    if (verbosity > 0 && is_app_mode)
    {
        cerr << endl;
        cerr << "*** Component sizes ***" << endl;
        cerr << "Reference sequences    : " << total_size_ref << endl;
        cerr << "   (only ref)          : " << total_size_only_ref << endl;
        cerr << "Raw sequences          : " << total_size_raw << endl;
        cerr << "Delta sequences        : " << total_size_delta - total_size_raw << endl;
        cerr << "Params                 : " << out_archive->GetStreamPackedSize(out_archive->GetStreamId("params")) << endl;
        cerr << "Splitters              : " << out_archive->GetStreamPackedSize(out_archive->GetStreamId("splitters")) << endl;
        cerr << "Segment splitters      : " << out_archive->GetStreamPackedSize(out_archive->GetStreamId("segment-splitters")) << endl;

        if(archive_version < 2000)
            cerr << "Collection desc.       : " << 
                out_archive->GetStreamPackedSize(out_archive->GetStreamId("collection-desc")) << endl;
        else if(archive_version < 3000)
            cerr << "Collection desc.       : " << 
                out_archive->GetStreamPackedSize(out_archive->GetStreamId("collection-main")) + 
                out_archive->GetStreamPackedSize(out_archive->GetStreamId("collection-details")) << endl;
        else
            cerr << "Collection desc.       : " << 
                out_archive->GetStreamPackedSize(out_archive->GetStreamId("collection-samples")) +
                out_archive->GetStreamPackedSize(out_archive->GetStreamId("collection-contigs")) +
                out_archive->GetStreamPackedSize(out_archive->GetStreamId("collection-details")) << endl;

        cerr << "*** Stats ***" << endl;
        cerr << "No. segments           : " << no_segments << endl;
        cerr << "No. one-side segments  : " << no_segments_one_side << endl;
        cerr << "No. only ref. segments : " << no_only_ref_segments << endl;
    }
}

// *******************************************************************************************
void CAGCCompressor::store_file_type_info()
{
    vector<uint8_t> v_data;

    for (auto& x : m_file_type_info)
    {
        append(v_data, x.first);
        append(v_data, x.second);
    }

    auto s_id = out_archive->RegisterStream("file_type_info");

    out_archive->AddPart(s_id, v_data, m_file_type_info.size());
}

// *******************************************************************************************
void CAGCCompressor::appending_init()
{
    no_segments = 0;

    if(archive_version >= 3000)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->prepare_for_appending_load_last_batch();

    while (true)
    {
        auto ref_stream_id = in_archive->GetStreamId(ss_ref_name(archive_version, no_segments));
        auto delta_stream_id = in_archive->GetStreamId(ss_delta_name(archive_version, no_segments));

        if (ref_stream_id < 0 && delta_stream_id < 0)
            break;

        v_segments.emplace_back(make_shared<CSegment>(ss_base(archive_version, no_segments), in_archive, out_archive, pack_cardinality, min_match_len, concatenated_genomes, archive_version));
        v_segments.back()->appending_init();

        ++no_segments;
    }

    auto vss = v_segments.size();
    if (!is_power_2(vss))
    {
        while (!is_power_2(vss))
            vss &= vss - 1;
        vss *= 2;
        v_segments.resize(vss);
    }

    auto splitters_stream_id = in_archive->GetStreamId("splitters");
    vector<uint8_t> v_tmp;
    uint64_t no_splitters;
    uint64_t x1, x2;
    uint32_t x3;

    in_archive->GetPart(splitters_stream_id, v_tmp, no_splitters);

    auto p = v_tmp.begin();

    hs_splitters.clear();
    for (uint32_t i = 0; i < no_splitters; ++i)
    {
        read64(p, x1);
        hs_splitters.insert(x1);
    }

    bloom_splitters.resize((uint64_t) (no_splitters / 0.25));
    bloom_splitters.insert(hs_splitters.begin(), hs_splitters.end());

    auto map_segments_id = in_archive->GetStreamId("segment-splitters");
    uint64_t no_stored_segment_maps;
    in_archive->GetPart(map_segments_id, v_tmp, no_stored_segment_maps);

    p = v_tmp.begin();
    map_segments[make_pair(~0ull, ~0ull)] = 0;

    for (uint32_t i = 0; i < no_stored_segment_maps; ++i)
    {
        read64(p, x1);
        read64(p, x2);
        read(p, x3);

        map_segments[make_pair(x1, x2)] = x3;

        if (x1 != ~0ull && x2 != ~0ull)
        {
            map_segments_terminators[x1].emplace_back(x2);
            if (x1 != x2)
                map_segments_terminators[x2].emplace_back(x1);
        }
    }

    buffered_seg_part.resize(no_segments);

    for (auto& term : map_segments_terminators)
        sort(term.second.begin(), term.second.end());
}

// *******************************************************************************************
bool CAGCCompressor::determine_splitters(const string& reference_file_name, const size_t segment_size, const uint32_t no_threads)
{
    CGenomeIO gio;
    if (!gio.Open(reference_file_name, false))
    {
        if (is_app_mode)
            cerr << "Cannot open file: " << reference_file_name << endl;
        return false;
    }

    string id;
    contig_t contig;

#ifndef _DEBUG
    v_candidate_kmers.resize(gio.FileSize() + raduls::ALIGNMENT / 8 + contig_part_size * (no_threads + 1) + 1, ~0ull);

    auto alignment_shift = ((uint64_t)v_candidate_kmers.data()) % raduls::ALIGNMENT;
//    auto extra_items = raduls::ALIGNMENT / 8 - alignment_shift / 8;
    v_candidate_kmers_offset = raduls::ALIGNMENT / 8 - alignment_shift / 8;
#else
    v_candidate_kmers.resize(gio.FileSize() + 256 / 8 + contig_part_size * (no_threads + 1) + 1, ~0ull);

    auto alignment_shift = ((uint64_t)v_candidate_kmers.data()) % 256;
    //    auto extra_items = raduls::ALIGNMENT / 8 - alignment_shift / 8;
    v_candidate_kmers_offset = 256 / 8 - alignment_shift / 8;
#endif

    if (verbosity > 0 && is_app_mode)
        cerr << "Gathering reference k-mers\n";

    q_contigs_data = make_unique<CBoundedQueue<contig_t>>(1, contig_part_size * no_threads * 3);
    vector<thread> v_threads;

    start_kmer_collecting_threads(v_threads, no_threads, v_candidate_kmers, v_candidate_kmers_offset);

    while (gio.ReadContigRaw(id, contig))
    {
        preprocess_raw_contig(contig);

        size_t c_size = contig.size();
        size_t start_pos;

        for (start_pos = 0; start_pos + (kmer_length - 1) < c_size; )
        {
            contig_t part = get_part(contig, start_pos, contig_part_size);

            start_pos += part.size() - (kmer_length - 1);
            q_contigs_data->Emplace(move(part), part.size());
        }
    }

    q_contigs_data->MarkCompleted();

    join_threads(v_threads);

    q_contigs_data.release();

    // Sort k-mers
    if (verbosity > 0 && is_app_mode)
        cerr << "Determination of splitters\n";

#ifdef _DEBUG
    sort(v_candidate_kmers.begin(), v_candidate_kmers.end());
#else
       raduls::RadixSortMSD((uint8_t*)(v_candidate_kmers.data() + v_candidate_kmers_offset), nullptr, v_candidate_kmers.size() - v_candidate_kmers_offset, 8, 8, no_threads);
#endif

    if(adaptive_compression)
        remove_non_singletons(v_candidate_kmers, v_duplicated_kmers, v_candidate_kmers_offset);
    else
        remove_non_singletons(v_candidate_kmers, v_candidate_kmers_offset);

    if (verbosity > 1 && is_app_mode)
        cerr << "No. of singletons: " << v_candidate_kmers.size() - v_candidate_kmers_offset << endl;

    auto v_begin = v_candidate_kmers.begin() + v_candidate_kmers_offset;
    auto v_end = v_candidate_kmers.end();

    gio.Close();
    
    if (!gio.Open(reference_file_name, false))
    {
        if (is_app_mode)
            cerr << "Cannot open file: " << reference_file_name << endl;
        return false;
    }

    // Determine splitters
    pq_contigs_raw = make_unique<CBoundedPQueue<contig_t>>(1, 4ull << 30);

    vv_splitters.resize(no_threads);
    start_splitter_finding_threads(v_threads, no_threads, v_begin, v_end, vv_splitters);

    while (gio.ReadContigRaw(id, contig))
        pq_contigs_raw->Emplace(move(contig), 0, contig.size());

    pq_contigs_raw->MarkCompleted();

    join_threads(v_threads);

    pq_contigs_raw.release();

    if (!adaptive_compression)
    {
        v_candidate_kmers.clear();
        v_candidate_kmers.shrink_to_fit();
        v_candidate_kmers_offset = 0;
    }

    hs_splitters.clear();

    for (auto& vec : vv_splitters)
    {
        for (auto d : vec)
            hs_splitters.insert(d);

        vec.clear();
        vec.shrink_to_fit();
    }

    bloom_splitters.resize((uint64_t)(hs_splitters.size() / 0.25));
    bloom_splitters.insert(hs_splitters.begin(), hs_splitters.end());

    gio.Close();

    if (verbosity > 1 && is_app_mode)
        cerr << "No. of splitters: " << hs_splitters.size() << endl;

    return true;
}

// *******************************************************************************************
bool CAGCCompressor::count_kmers(vector<pair<string, vector<uint8_t>>>& v_contig_data, const uint32_t no_threads)
{
    size_t tot_contig_len = 0;

    for (auto& cd : v_contig_data)
        tot_contig_len += cd.second.size();

#ifndef _DEBUG
    v_candidate_kmers.resize(tot_contig_len + raduls::ALIGNMENT / 8 + contig_part_size * (no_threads + 1) + 1, ~0ull);

    auto alignment_shift = ((uint64_t)v_candidate_kmers.data()) % raduls::ALIGNMENT;
    //    auto extra_items = raduls::ALIGNMENT / 8 - alignment_shift / 8;
    v_candidate_kmers_offset = raduls::ALIGNMENT / 8 - alignment_shift / 8;
#else
    v_candidate_kmers.resize(tot_contig_len + 256 / 8 + contig_part_size * (no_threads + 1) + 1, ~0ull);

    auto alignment_shift = ((uint64_t)v_candidate_kmers.data()) % 256;
    //    auto extra_items = raduls::ALIGNMENT / 8 - alignment_shift / 8;
    v_candidate_kmers_offset = 256 / 8 - alignment_shift / 8;
#endif

    if (verbosity > 0 && is_app_mode)
        cerr << "Gathering reference k-mers\n";

    q_contigs_data = make_unique<CBoundedQueue<contig_t>>(1, contig_part_size * no_threads * 3);
    vector<thread> v_threads;

    start_kmer_collecting_threads(v_threads, no_threads, v_candidate_kmers, v_candidate_kmers_offset);

    for (auto& cd : v_contig_data)
    {
        preprocess_raw_contig(cd.second);

        size_t c_size = cd.second.size();
        size_t start_pos;

        for (start_pos = 0; start_pos + (kmer_length - 1) < c_size; )
        {
            contig_t part = get_part(cd.second, start_pos, contig_part_size);

            start_pos += part.size() - (kmer_length - 1);
            q_contigs_data->Emplace(move(part), part.size());
        }
    }

    q_contigs_data->MarkCompleted();

    join_threads(v_threads);

    q_contigs_data.release();

    // Sort k-mers
    if (verbosity > 0 && is_app_mode)
        cerr << "Determination of splitters\n";

#ifdef _DEBUG
    sort(v_candidate_kmers.begin(), v_candidate_kmers.end());
#else
    raduls::RadixSortMSD((uint8_t*)(v_candidate_kmers.data() + v_candidate_kmers_offset), nullptr, v_candidate_kmers.size() - v_candidate_kmers_offset, 8, 8, no_threads);
#endif

    remove_non_singletons(v_candidate_kmers, v_duplicated_kmers, v_candidate_kmers_offset);

    if (verbosity > 1 && is_app_mode)
        cerr << "No. of singletons: " << v_candidate_kmers.size() - v_candidate_kmers_offset << endl;

    return true;
}

// *******************************************************************************************
void CAGCCompressor::enumerate_kmers(contig_t& ctg, vector<uint64_t>& vec)
{
    CKmer kmer(kmer_length, kmer_mode_t::canonical);

    vec.clear();

    if (ctg.size() < kmer_length)
        return;

    vec.reserve(ctg.size() + 1 - kmer_length);

    kmer.Reset();

    for (auto x : ctg)
    {
        if (x > 3)
            kmer.Reset();
        else
        {
            kmer.insert(x);

            if (kmer.is_full())
                vec.emplace_back(kmer.data());
        }
    }
}

// *******************************************************************************************
void CAGCCompressor::remove_non_singletons(vector<uint64_t>& vec, size_t virtual_begin)
{
    // Remove non-singletons
    size_t curr_end = virtual_begin;
    for (size_t i = virtual_begin; i < vec.size();)
    {
        size_t j;
        for (j = i + 1; j < vec.size() && vec[i] == vec[j]; ++j)
            ;

        if (i + 1 == j)
            vec[curr_end++] = vec[i];

        i = j;
    }

    vec.resize(curr_end);
}

// *******************************************************************************************
void CAGCCompressor::remove_non_singletons(vector<uint64_t>& vec, vector<uint64_t>& v_duplicated, size_t virtual_begin)
{
    v_duplicated.clear();

    size_t curr_end = virtual_begin;
    for (size_t i = virtual_begin; i < vec.size();)
    {
        size_t j;
        for (j = i + 1; j < vec.size() && vec[i] == vec[j]; ++j)
            ;

        if (i + 1 == j)
            vec[curr_end++] = vec[i];
        else
            v_duplicated.emplace_back(vec[i]);

        i = j;
    }

    vec.resize(curr_end);
}

// *******************************************************************************************
void CAGCCompressor::start_kmer_collecting_threads(vector<thread> &v_threads, const uint32_t n_t, vector<uint64_t>& v_kmers, const size_t extra_items)
{
    v_threads.clear();
    v_threads.reserve(n_t);

    a_part_id = 0;

    for (uint32_t i = 0; i < n_t; ++i)
        v_threads.emplace_back([&, extra_items] {

        CKmer kmer(kmer_length, kmer_mode_t::canonical);

        uint64_t curr_part_id = 0;
        uint64_t j = contig_part_size;
        uint64_t part_start_pos = extra_items;

        while (!q_contigs_data->IsCompleted())
        {
            contig_t task;

            if (!q_contigs_data->Pop(task))
                continue;

            kmer.Reset();

            for (auto x : task)
            {
                if (x > 3)
                    kmer.Reset();
                else
                {
                    kmer.insert(x);

                    if (kmer.is_full())
                    {
                        if (j == contig_part_size)
                        {
                            curr_part_id = atomic_fetch_add(&a_part_id, 1);
                            part_start_pos = extra_items + curr_part_id * contig_part_size;
                            j = 0;

                            if (part_start_pos + contig_part_size >= v_kmers.size() && is_app_mode)
                                cerr << "Problem with v_kmers\n";
                        }

                        v_kmers[part_start_pos + j++] = kmer.data();
                    }
                }
            }
        }

        });
}

// *******************************************************************************************
void CAGCCompressor::find_splitters_in_contig(contig_t& ctg, const vector<uint64_t>::iterator v_begin, const vector<uint64_t>::iterator v_end, vector<uint64_t>& v_splitters)
{
    // Initialization to large value to add 1st candidate k-mer
    uint64_t current_len = segment_size;
    vector<uint64_t> v_recent_kmers;
    CKmer kmer(kmer_length, kmer_mode_t::canonical);

    kmer.Reset();

    for (auto x : ctg)
    {
        if (x > 3)
            kmer.Reset();
        else
        {
            kmer.insert(x);

            if (kmer.is_full())
            {
                v_recent_kmers.emplace_back(kmer.data());

                if (current_len >= segment_size)
                {
                    uint64_t d = kmer.data();

                    if (binary_search(v_begin, v_end, d))
                    {
                        v_splitters.emplace_back(d);
                        current_len = 0;
                        kmer.Reset();
                        v_recent_kmers.clear();
                    }
                }
            }
        }

        ++current_len;
    }

    // Try add the rightmost candidate k-mer
    for (auto p = v_recent_kmers.rbegin(); p != v_recent_kmers.rend(); ++p)
        if (binary_search(v_begin, v_end, *p))
        {
            v_splitters.emplace_back(*p);
            break;
        }
}

// *******************************************************************************************
void CAGCCompressor::build_candidate_kmers_from_archive(const uint32_t n_t)
{
    CAGCDecompressor agc_decompressor(false);

    agc_decompressor.AssignArchive(*this);

    string ref_sample_name;

    if (!collection_desc->get_reference_name(ref_sample_name))
        return;

    vector<pair<string, vector<uint8_t>>> v_contig_data;

    agc_decompressor.GetSampleSequences(ref_sample_name, v_contig_data, n_t);

    count_kmers(v_contig_data, n_t);

    vv_splitters.resize(n_t);
}

// *******************************************************************************************
void CAGCCompressor::start_splitter_finding_threads(vector<thread>& v_threads, const uint32_t n_t, 
    const vector<uint64_t>::iterator v_begin, const vector<uint64_t>::iterator v_end, vector<vector<uint64_t>>& v_splitters)
{
    v_threads.clear();
    v_threads.reserve(n_t);

    for (uint32_t i = 0; i < n_t; ++i)
        v_threads.emplace_back([&, i, v_begin, v_end] {

        uint32_t thread_id = i;

        while (!pq_contigs_raw->IsCompleted())
        {
            contig_t task;

            if (!pq_contigs_raw->PopLarge(task))
                continue;

            preprocess_raw_contig(task);

            find_splitters_in_contig(task, v_begin, v_end, v_splitters[thread_id]);
        }

        });
}

// *******************************************************************************************
void CAGCCompressor::start_finalizing_threads(vector<thread>& v_threads, const uint32_t n_t)
{
    v_threads.clear();
    v_threads.reserve(n_t);

    id_segment = 0;

    for (uint32_t i = 0; i < n_t; ++i)
        v_threads.emplace_back([&] {
        auto zstd_ctx = ZSTD_createCCtx();

        while (true)
        {
            uint32_t j = atomic_fetch_add(&id_segment, 1);

            if (j >= no_segments)
                break;

            v_segments[j]->finish(zstd_ctx);
            v_segments[j].reset();
        }

        ZSTD_freeCCtx(zstd_ctx);
        });
}

// *******************************************************************************************
void CAGCCompressor::preprocess_raw_contig(contig_t& ctg)
{
    size_t len = ctg.size();
    size_t in_pos = 0;
    size_t out_pos = 0;

    uint8_t c;

    switch (len % 4)
    {
    case 3:
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
    case 2:
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
    case 1:
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
    }

    for (; in_pos < len;)
    {
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
        c = ctg[in_pos++];
        if (c >> 6)                          // (c >= 64)
            ctg[out_pos++] = cnv_num[c];
    }

    ctg.resize(out_pos);
//    ctg.shrink_to_fit();
}

// *******************************************************************************************
void CAGCCompressor::register_segments(uint32_t n_t)
{
    buffered_seg_part.sort_known(n_t);           // TODO: sort in parallel by many threads
    uint32_t no_new = buffered_seg_part.process_new();

    for (uint32_t i = 0; i < no_new; ++i)
    {
        out_archive->RegisterStream(ss_ref_name(archive_version, no_segments + i));
        out_archive->RegisterStream(ss_delta_name(archive_version, no_segments + i));
    }

    no_segments += no_new;

    if (no_segments > v_segments.size())
        v_segments.resize(no_segments);

    buffered_seg_part.distribute_segments(0, 0, no_raw_groups);

    buffered_seg_part.restart_read_vec();
}

// *******************************************************************************************
void CAGCCompressor::store_segments(ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    string sample_name;
    string contig_name;
    contig_t seg_data;
    bool is_rev_comp;
    uint32_t seg_part_no;
    uint64_t kmer1;
    uint64_t kmer2;

    while (true)
    {
        int block_group_id = buffered_seg_part.get_vec_id();
        int in_group_id;

        if (block_group_id < 0)
            break;

        for(int group_id = block_group_id; group_id > block_group_id - 10; --group_id)
            if(!buffered_seg_part.is_empty_part(group_id))
                while (buffered_seg_part.get_part(group_id, kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no))
                {
/*                    if (contig_name == "cluster3_contig_100")
                        cout << "!";*/

                    if (v_segments[group_id] == nullptr)
                    {
                        v_segments[group_id] = make_shared<CSegment>(ss_base(archive_version, group_id), nullptr, out_archive, pack_cardinality, min_match_len, concatenated_genomes, archive_version);

                        seg_map_mtx.lock();

                        auto p = map_segments.find(make_pair(kmer1, kmer2));
                        if (p == map_segments.end())
                            map_segments[make_pair(kmer1, kmer2)] = group_id;
                        else if (p->second > group_id)
                            p->second = group_id;

                        if (kmer1 != ~0ull && kmer2 != ~0ull)
                        {
                            map_segments_terminators[kmer1].push_back(kmer2);
                            sort(map_segments_terminators[kmer1].begin(), map_segments_terminators[kmer1].end());

                            if (kmer1 != kmer2)
                            {
                                map_segments_terminators[kmer2].push_back(kmer1);
                                sort(map_segments_terminators[kmer2].begin(), map_segments_terminators[kmer2].end());
                            }
                        }

                        seg_map_mtx.unlock();
                    }

                    if (group_id < (int)no_raw_groups)
                    {
                        in_group_id = v_segments[group_id]->add_raw(seg_data, zstd_cctx, zstd_dctx);
                    }
                    else
                        in_group_id = v_segments[group_id]->add(seg_data, zstd_cctx, zstd_dctx);

                    collection_desc->add_segment_placed(sample_name, contig_name, seg_part_no, group_id, in_group_id, is_rev_comp, (uint32_t)seg_data.size());
                }
    }
}

// *******************************************************************************************
// Start compressing threads
void CAGCCompressor::start_compressing_threads(vector<thread>& v_threads, my_barrier &bar, const uint32_t n_t)
{
    v_threads.clear();

    for (uint32_t i = 0; i < n_t; ++i)
    {
        v_threads.emplace_back([&, i, n_t]() {
            auto zstd_cctx = ZSTD_createCCtx();
            auto zstd_dctx = ZSTD_createDCtx();

            uint32_t thread_id = i;

            while (!pq_contigs_desc_working->IsCompleted())
            {
                task_t task;

                if (!pq_contigs_desc_working->PopLarge(task))
                    continue;

                if (get<0>(task) == contig_processing_stage_t::registration)
                {
                    // synchronization token received
                    bar.arrive_and_wait();
                    if (thread_id == 0)
                        register_segments(n_t);
                    bar.arrive_and_wait();

                    store_segments(zstd_cctx, zstd_dctx);

                    bar.arrive_and_wait();

                    if (thread_id == 0)
                    {
                        buffered_seg_part.clear(max(1u, n_t-1));

                        if (n_t == 1)
                        {
                            if (!concatenated_genomes)
                                ++processed_samples;
                            else
                            {
                                processed_samples = processed_samples / pack_cardinality * pack_cardinality + pack_cardinality;

                                auto max_ps = dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();
                                if (max_ps < processed_samples)
                                    processed_samples = max_ps;
                            }

                            if (archive_version >= 3000 && processed_samples % pack_cardinality == 0)
                                dynamic_pointer_cast<CCollection_V3>(collection_desc)->store_contig_batch(processed_samples - pack_cardinality, processed_samples);

                            out_archive->FlushOutBuffers();
                        }

                        // !!! ???
                        if(adaptive_compression)
                            pq_contigs_desc_working = pq_contigs_desc;
                    }
                    else if (thread_id == 1)
                    {
                        if (!concatenated_genomes)
                            ++processed_samples;
                        else
                        {
                            processed_samples = processed_samples / pack_cardinality * pack_cardinality + pack_cardinality;

                            auto max_ps = dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();
                            if (max_ps < processed_samples)
                                processed_samples = max_ps;
                        }

                        if (archive_version >= 3000 && processed_samples % pack_cardinality == 0)
                            dynamic_pointer_cast<CCollection_V3>(collection_desc)->store_contig_batch(processed_samples - pack_cardinality, processed_samples);

                        out_archive->FlushOutBuffers();
                    }

                    bar.arrive_and_wait();

                    continue;
                }

                if (get<0>(task) == contig_processing_stage_t::new_splitters)
                {
                    bar.arrive_and_wait();
                    if (thread_id == 0)
                    {
                        // Add new splitters
                        for(auto &v : vv_splitters)
                        { 
                            for (auto& x : v)
                            {
                                hs_splitters.insert_fast(x);
                                bloom_splitters.insert(x);
                            }

                            v.clear();
                        }

                        if (bloom_splitters.filling_factor() > 0.3)
                        {
                            bloom_splitters.resize((uint64_t)(hs_splitters.size() / 0.25));
                            bloom_splitters.insert(hs_splitters.begin(), hs_splitters.end());
                        }

                        for (auto& x : v_raw_contigs)
                            pq_contigs_desc_aux->Emplace(make_tuple(contig_processing_stage_t::hard_contigs, get<0>(x), get<1>(x), move(get<2>(x))), 1, get<2>(x).size());

                        v_raw_contigs.clear();

                        for(uint32_t i = 0; i < n_t; ++i)
                            pq_contigs_desc_aux->Emplace(make_tuple(contig_processing_stage_t::registration, "", "", contig_t()), 0, 0);

                        pq_contigs_desc_working = pq_contigs_desc_aux;
                    }

                    bar.arrive_and_wait();

                    continue;
                }

                if (get<0>(task) == contig_processing_stage_t::all_contigs)
                {
                    preprocess_raw_contig(get<3>(task));
                }

                size_t ctg_size = get<3>(task).size();

                if (compress_contig(get<0>(task), get<1>(task), get<2>(task), get<3>(task), zstd_cctx, zstd_dctx, thread_id))
                {
                    processed_bases += ctg_size;

                    if (verbosity > 0 && is_app_mode)
                        cerr << "Compressed: " + to_string(processed_bases / 1'000'000) + " Mb\r";
                    fflush(stdout);
                }
                else
                {
                    lock_guard<mutex> lck(mtx_raw_contigs);
                    v_raw_contigs.emplace_back(get<1>(task), get<2>(task), move(get<3>(task)));
                }

                get<3>(task).clear();
                get<3>(task).shrink_to_fit();
            }

            ZSTD_freeCCtx(zstd_cctx);
            ZSTD_freeDCtx(zstd_dctx);
            });
    }
}

// *******************************************************************************************
pair_segment_desc_t CAGCCompressor::add_segment(const string& sample_name, const string& contig_name, uint32_t seg_part_no,
    contig_t segment, CKmer kmer_front, CKmer kmer_back, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    pair<uint64_t, uint64_t> pk, pk2(~0ull, ~0ull);
    contig_t segment_rc;
    contig_t segment2, segment2_rc;
    bool store_rc = false;
    bool store2_rc = false;

    int segment_id{ -1 }, segment_id2{ -1 };

    if (!kmer_front.is_full() && !kmer_back.is_full())
    {
        // No terminal splitters present

        pk = make_pair(~0ull, ~0ull);
    }
    else if (kmer_front.is_full() && kmer_back.is_full())
    {
        // Both terminal splitters present

        if (kmer_front.data() < kmer_back.data())
            pk = make_pair(kmer_front.data(), kmer_back.data());
        else
        {
            pk = make_pair(kmer_back.data(), kmer_front.data());
            reverse_complement_copy(segment, segment_rc);
            store_rc = true;
        }
    }
    else if (kmer_front.is_full())
    {
        CKmer kmer = kmer_front;
        reverse_complement_copy(segment, segment_rc);

        tie(pk, store_rc) = find_cand_segment_with_one_splitter(kmer, segment, segment_rc, zstd_dctx);
    }
    else if (kmer_back.is_full())
    {
        CKmer kmer = kmer_back;
        kmer.swap_dir_rc();
        reverse_complement_copy(segment, segment_rc);
        bool store_dir;

        tie(pk, store_dir) = find_cand_segment_with_one_splitter(kmer, segment_rc, segment, zstd_dctx);
        store_rc = !store_dir;
    }

    auto p = map_segments.find(pk);

    // There is no such a segment terminated by the splitters, so let's try to check if the segment can be splitted into two segments
    if (!concatenated_genomes &&
        p == map_segments.end() &&
        pk.first != ~0ull && pk.second != ~0ull &&
        map_segments_terminators.count(pk.first) && map_segments_terminators.count(pk.second))
    {
        if (segment_rc.empty())
            reverse_complement_copy(segment, segment_rc);

        if (kmer_front.data() == kmer_back.data())
        {
            if (!kmer_front.is_dir_oriented())
                store_rc = true;
        }
        else
        {
            CKmer kmer1 = kmer_front;
            CKmer kmer2 = kmer_back;
            bool use_rc = false;

            if (kmer1.data() > kmer2.data())
            {
                swap(kmer1, kmer2);
                use_rc = true;
                kmer1.swap_dir_rc();
                kmer2.swap_dir_rc();
            }

            auto split_match = find_cand_segment_with_missing_middle_splitter(kmer1, kmer2, use_rc ? segment_rc : segment, use_rc ? segment : segment_rc, zstd_dctx);

            if (split_match.first != ~0ull)
            {
                uint32_t left_size = split_match.second;
                uint32_t right_size = (uint32_t)(segment.size() - split_match.second);

                if (left_size == 0)
                {
                    if (split_match.first < kmer2.data())
                        store_rc = use_rc;
                    else
                        store_rc = !use_rc;
                    pk = minmax(split_match.first, kmer2.data());
                }
                else if (right_size == 0)
                {
                    if (kmer1.data() < split_match.first)
                        store_rc = use_rc;
                    else
                        store_rc = !use_rc;
                    pk = minmax(kmer1.data(), split_match.first);
                }
                else
                {
                    if (use_rc)
                        swap(left_size, right_size);

                    // Split segment into 2 parts (with overlap of size kmer_length)
                    uint32_t seg2_start_pos = left_size - kmer_length / 2;
                    segment2.assign(segment.begin() + seg2_start_pos, segment.end());

                    segment.resize((size_t)seg2_start_pos + kmer_length);

                    if (kmer_front.data() < split_match.first)
                    {
                        store_rc = false;
                        pk = make_pair(kmer_front.data(), split_match.first);
                    }
                    else
                    {
                        store_rc = true;
                        reverse_complement_copy(segment, segment_rc);
                        pk = make_pair(split_match.first, kmer_front.data());
                    }

                    segment_id = map_segments.at(pk);          // must exists

                    if (split_match.first < kmer_back.data())
                    {
                        store2_rc = false;
                        pk2 = make_pair(split_match.first, kmer_back.data());
                    }
                    else
                    {
                        store2_rc = true;
                        reverse_complement_copy(segment2, segment2_rc);
                        pk2 = make_pair(kmer_back.data(), split_match.first);
                    }

                    segment_id2 = map_segments.at(pk2);         // must exists
                }
            }
        }

        p = map_segments.find(pk);
    }
    
    uint32_t segment_size = (uint32_t) segment.size();
    uint32_t segment2_size = (uint32_t) segment2.size();

    if (p == map_segments.end())
    {
        buffered_seg_part.add_new(pk.first, pk.second, sample_name, contig_name, store_rc ? segment_rc : segment, store_rc, seg_part_no);
    }
    else
    {
        if (segment_id2 == -1)
            segment_id = p->second;

        buffered_seg_part.add_known(segment_id, ~0ull, ~0ull, sample_name, contig_name, store_rc ? move(segment_rc) : move(segment), store_rc, seg_part_no);

        if (segment_id2 >= 0)
            buffered_seg_part.add_known(segment_id2, ~0ull, ~0ull, sample_name, contig_name, store2_rc ? move(segment2_rc) : move(segment2), store2_rc, seg_part_no + 1);
    }

    return pair_segment_desc_t(segment_desc_t(segment_id, 0, store_rc, segment_size), segment_desc_t(segment_id2, 0, store2_rc, segment2_size), segment_id2 >= 0);
}

// *******************************************************************************************
pair<uint64_t, uint32_t> CAGCCompressor::find_cand_segment_with_missing_middle_splitter(CKmer kmer_front, CKmer kmer_back, contig_t& segment_dir, contig_t& segment_rc, ZSTD_DCtx* zstd_dctx)
{
    vector<uint64_t> shared_splitters;

    shared_splitters.resize(map_segments_terminators[kmer_front.data()].size());

    auto p_front = map_segments_terminators.find(kmer_front.data());
    auto p_back = map_segments_terminators.find(kmer_back.data());

    auto p_shared = set_intersection(
        p_front->second.begin(), p_front->second.end(),
        p_back->second.begin(), p_back->second.end(),
        shared_splitters.begin());

    shared_splitters.erase(p_shared, shared_splitters.end());
    shared_splitters.erase(remove(shared_splitters.begin(), shared_splitters.end(), ~0ull), shared_splitters.end());

    if (shared_splitters.empty())
        return make_pair(~0ull, 0);     // unable to find middle splitter

    uint64_t middle = shared_splitters.front();     // Take 1st shared k-mer and ignore the rest (at least in this version; maybe in the future it is worth to investigate the remaining ones)

    vector<uint32_t> v_costs1, v_costs2;

    auto segment_id1 = map_segments[minmax(kmer_front.data(), middle)];
    auto segment_id2 = map_segments[minmax(middle, kmer_back.data())];


    auto seg1 = v_segments[segment_id1];
    auto seg2 = v_segments[segment_id2];
    
    if (kmer_front.data() < middle)
        seg1->get_coding_cost(segment_dir, v_costs1, true, zstd_dctx);
    else
    {
        seg1->get_coding_cost(segment_rc, v_costs1, false, zstd_dctx);
        reverse(v_costs1.begin(), v_costs1.end());
    }

    if (middle < kmer_back.data())
    {
        seg2->get_coding_cost(segment_dir, v_costs2, false, zstd_dctx);
        reverse(v_costs2.begin(), v_costs2.end());
    }
    else
        seg2->get_coding_cost(segment_rc, v_costs2, true, zstd_dctx);

    partial_sum(v_costs1.begin(), v_costs1.end(), v_costs1.begin());
    partial_sum(v_costs2.begin(), v_costs2.end(), v_costs2.begin());

    reverse(v_costs2.begin(), v_costs2.end());

    uint32_t best_sum = ~0u;
    uint32_t best_pos = 0;

    for (uint32_t i = 0; i < v_costs1.size(); ++i)
    {
        uint32_t cs = v_costs1[i] + v_costs2[i];
        if (cs < best_sum)
        {
            best_sum = cs;
            best_pos = i;
        }
    }

    if (best_pos < kmer_length + 1u)
        best_pos = 0;
    if ((size_t)best_pos + kmer_length + 1u > v_costs1.size())
        best_pos = (uint32_t)v_costs1.size();

    return make_pair(middle, best_pos);
}

// *******************************************************************************************
pair<pair<uint64_t, uint64_t>, bool> CAGCCompressor::find_cand_segment_with_one_splitter(CKmer kmer, contig_t& segment_dir, contig_t& segment_rc, ZSTD_DCtx* zstd_dctx)
{
    const pair<uint64_t, uint64_t> empty_pk(~0ull, ~0ull);
    pair<uint64_t, uint64_t> best_pk(~0ull, ~0ull);
    uint64_t best_estim_size = segment_dir.size() < 16 ? segment_dir.size() : segment_dir.size() - 16u;
    bool is_best_rc = false;

    vector<tuple<uint64_t, uint64_t, bool, shared_ptr<CSegment>>> v_candidates;

    auto p = map_segments_terminators.find(kmer.data());
    if (p == map_segments_terminators.end())
    {
        if (kmer.is_dir_oriented())
            best_pk = make_pair(kmer.data(), ~0ull);
        else
        {
            best_pk = make_pair(~0ull, kmer.data());
            is_best_rc = true;
        }

        return make_pair(best_pk, is_best_rc);
    }

    v_candidates.reserve(p->second.size());

    for (auto cand_kmer : p->second)
    {
        pair<uint64_t, uint64_t> cand_pk;

        v_candidates.emplace_back();
        auto& ck = v_candidates.back();

        if (cand_kmer < kmer.data())
        {
            cand_pk = make_pair(cand_kmer, kmer.data());
            get<0>(ck) = cand_kmer;
            get<1>(ck) = kmer.data();
            get<2>(ck) = true;
        }
        else
        {
            cand_pk = make_pair(kmer.data(), cand_kmer);
            get<0>(ck) = kmer.data();
            get<1>(ck) = cand_kmer;
            get<2>(ck) = false;
        }

        get<3>(ck) = v_segments[map_segments[cand_pk]];
    }

    int64_t segment_size = (int64_t) segment_dir.size();
    stable_sort(v_candidates.begin(), v_candidates.end(), [segment_size](const auto& x, const auto& y) {
        int64_t x_size = get<3>(x)->get_ref_size();
        int64_t y_size = get<3>(y)->get_ref_size();

        if (abs(segment_size - x_size) < abs(segment_size - y_size))
            return true;
        if (abs(segment_size - x_size) > abs(segment_size - y_size))
            return false;
        return x_size < y_size;
        });

    for (auto& candidate : v_candidates)
    {
        auto estim_size = get<3>(candidate)->estimate(get<2>(candidate) ? segment_rc : segment_dir, best_estim_size, zstd_dctx);
        auto cand_pk = make_pair(get<0>(candidate), get<1>(candidate));

        if (estim_size < best_estim_size || 
            (estim_size == best_estim_size && cand_pk < best_pk) ||
            (estim_size == best_estim_size && cand_pk == best_pk && !get<2>(candidate)))
        {
            best_estim_size = estim_size;
            best_pk = cand_pk;
            is_best_rc = get<2>(candidate);
        }
    }

    if (best_pk == empty_pk)
    {
        if (kmer.is_dir_oriented())
            best_pk = make_pair(kmer.data(), ~0ull);
        else
        {
            best_pk = make_pair(~0ull, kmer.data());
            is_best_rc = true;
        }
    }

    return make_pair(best_pk, is_best_rc);
}

// *******************************************************************************************
bool CAGCCompressor::compress_contig(contig_processing_stage_t contig_processing_stage, string sample_name, string id, contig_t& contig, 
    ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx, uint32_t thread_id)
{
    CKmer kmer(kmer_length, kmer_mode_t::canonical);

    uint64_t pos = 0;
    uint64_t split_pos = 0;
    CKmer split_kmer(kmer_length, kmer_mode_t::canonical);
    uint32_t seg_part_no = 0;

    for (auto x : contig)
    {
        if (x >> 2)         // x > 3
            kmer.Reset();
        else
        {
            kmer.insert_canonical(x);           // a bit faster than insert()

            if (kmer.is_full())
            {
                uint64_t d = kmer.data_canonical();     // a bit faster than data()
                if (bloom_splitters.check(d) && hs_splitters.check(d))
                {
                    auto seg_id = add_segment(sample_name, id, seg_part_no,
                        get_part(contig, split_pos, pos + 1 - split_pos), split_kmer, kmer, zstd_cctx, zstd_dctx);

                    ++seg_part_no;

                    if (seg_id.contains_second)
                        ++seg_part_no;

                    split_pos = pos + 1 - kmer_length;
                    split_kmer = kmer;
                    kmer.Reset();
                }
            }
        }

        ++pos;
    }

    if (adaptive_compression && contig_processing_stage == contig_processing_stage_t::all_contigs && split_kmer == CKmer(kmer_length, kmer_mode_t::canonical))
    {
        if(contig.size() >= segment_size)
            find_new_splitters(contig, thread_id);

        return false;
    }

    if (split_pos < contig.size())
        add_segment(sample_name, id, seg_part_no,
            get_part(contig, split_pos, contig.size() - split_pos), split_kmer, CKmer(kmer_length, kmer_mode_t::canonical), zstd_cctx, zstd_dctx);

    return true;
}

// *******************************************************************************************
void CAGCCompressor::find_new_splitters(contig_t& ctg, uint32_t thread_id)
{
    vector<uint64_t> v_contig_kmers;
    vector<uint64_t> v_tmp;

    enumerate_kmers(ctg, v_contig_kmers);
    sort(v_contig_kmers.begin(), v_contig_kmers.end());
    remove_non_singletons(v_contig_kmers, 0);

    v_tmp.resize(v_contig_kmers.size());

    // Exclude k-mers in reference genome - singletons
    auto p_end = set_difference(v_contig_kmers.begin(), v_contig_kmers.end(),
        v_candidate_kmers.begin() + v_candidate_kmers_offset, v_candidate_kmers.end(),
        v_tmp.begin());

    v_tmp.erase(p_end, v_tmp.end());

    // Exclude k-mers in reference genome - duplicated
    p_end = set_difference(v_tmp.begin(), v_tmp.end(),
        v_duplicated_kmers.begin(), v_duplicated_kmers.end(),
        v_contig_kmers.begin());

    v_contig_kmers.erase(p_end, v_contig_kmers.end());

    find_splitters_in_contig(ctg, v_contig_kmers.begin(), v_contig_kmers.end(), vv_splitters[thread_id]);
}

// *******************************************************************************************
contig_t CAGCCompressor::get_part(const contig_t& contig, uint64_t pos, uint64_t len)
{
    if (pos + len < contig.size())
        return contig_t(contig.begin() + pos, contig.begin() + pos + len);
    else
        return contig_t(contig.begin() + pos, contig.end());
}

// *******************************************************************************************
bool CAGCCompressor::close_compression(const uint32_t no_threads)
{
    if (working_mode != working_mode_t::compression && working_mode != working_mode_t::appending)
        return false;

    vector<thread> v_threads;

    start_finalizing_threads(v_threads, no_threads);
    join_threads(v_threads);

    out_archive->FlushOutBuffers();

    store_metadata(no_threads);

    if(archive_version >= 3000)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->complete_serialization();

    store_file_type_info();

    return true;
}

// *******************************************************************************************
// Add sample files
bool CAGCCompressor::AddSampleFiles(vector<pair<string, string>> _v_sample_file_name, const uint32_t no_threads)
{
    if (_v_sample_file_name.empty())
        return true;

    processed_bases = 0;

    size_t queue_capacity = max(2ull << 30, no_threads * (192ull << 20));

    pq_contigs_desc = make_shared<CBoundedPQueue<task_t>>(1, queue_capacity);
    pq_contigs_desc_aux = make_shared<CBoundedPQueue<task_t>>(1, ~0ull);
    pq_contigs_desc_working = pq_contigs_desc;

    uint32_t no_workers = (no_threads < 8) ? no_threads : no_threads - 1;

    vector<thread> v_threads;
    v_threads.reserve((size_t)no_workers);

    my_barrier bar(no_workers);

    start_compressing_threads(v_threads, bar, no_workers);

    // Reading Input
    pair<string, string> sf;
    CGenomeIO gio;
    string id;
    contig_t contig;
    size_t sample_priority = ~0ull;
    size_t cnt_contigs_in_sample = 0;
    const size_t max_no_contigs_before_synchronization = pack_cardinality;
    const size_t min_size_before_synchronization = 1ull << 30;

    if (archive_version >= 3000 && in_archive != nullptr)
        processed_samples = dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();
    else
        processed_samples = 0;

    if (concatenated_genomes)
        cnt_contigs_in_sample = processed_samples % pack_cardinality;

    for(auto sf : _v_sample_file_name)
    {
        if (archive_version >= 3000)
            dynamic_pointer_cast<CCollection_V3>(collection_desc)->reset_prev_sample_name();

        if (!gio.Open(sf.second, false))
        {
            cerr << "Cannot open file: " << sf.second << endl;
            continue;
        }
        else
        {
            ;
        }

        while (gio.ReadContigRaw(id, contig))
        {
            if (concatenated_genomes)
            {
                if (!collection_desc->register_sample_contig("", id))
                    cerr << "Error: Pair sample_name:contig_name " << id << ":" << id << " is already in the archive!\n";
                else
                {
                    pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs, "", id, move(contig)), sample_priority, contig.size());

                    if (++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
                    {
                        // Send synchronization tokens
                        for (uint32_t i = 0; i < no_workers; ++i)
                            pq_contigs_desc->Emplace(make_tuple(
                                adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration, "", "", contig_t()), sample_priority, 0);

                        cnt_contigs_in_sample = 0;
                        --sample_priority;
                    }
                }
            }
            else
            {
                if (collection_desc->register_sample_contig(sf.first, id))
                    pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs, sf.first, id, move(contig)), sample_priority, contig.size());
                else
                    cerr << "Error: Pair sample_name:contig_name " << sf.first << ":" << id << " is already in the archive!\n";
            }
        }

        if (!concatenated_genomes)
        {
            // Send synchronization tokens
            for (uint32_t i = 0; i < no_workers; ++i)
                pq_contigs_desc->Emplace(make_tuple(
                    adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
                    "", "", contig_t()), sample_priority, 0);

            --sample_priority;
        }

        gio.Close();
    }

    if (concatenated_genomes)// && ++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
    {
        // Send synchronization tokens
        for (uint32_t i = 0; i < no_workers; ++i)
            pq_contigs_desc->Emplace(make_tuple(
                adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration, "", "", contig_t()), sample_priority, 0);

        cnt_contigs_in_sample = 0;
        --sample_priority;
    }


    pq_contigs_desc->MarkCompleted();

    join_threads(v_threads);

    if(concatenated_genomes)
        processed_samples = dynamic_pointer_cast<CCollection_V3>(collection_desc)->get_no_samples();

    if (archive_version >= 3000 && processed_samples % pack_cardinality != 0)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->store_contig_batch((processed_samples / pack_cardinality) * pack_cardinality, processed_samples);

    out_archive->FlushOutBuffers();

    pq_contigs_desc.reset();
    pq_contigs_desc_aux.reset();
    pq_contigs_desc_working.reset();

    no_samples_in_archive += _v_sample_file_name.size();

    return true;
}

// *******************************************************************************************
bool CAGCCompressor::Create(const string& _file_name, const uint32_t _pack_cardinality, const uint32_t _kmer_length, const string& reference_file_name, const uint32_t _segment_size,
    const uint32_t _min_match_len, const bool _concatenated_genomes, const bool _adaptive_compression, const uint32_t _verbosity, const uint32_t no_threads)
{
    if (working_mode != working_mode_t::none)
        return false;

    pack_cardinality = _pack_cardinality;
    kmer_length = _kmer_length;
    min_match_len = _min_match_len;
    concatenated_genomes = _concatenated_genomes;
    adaptive_compression = _adaptive_compression;
    segment_size = _segment_size;
    verbosity = _verbosity;

    if (!determine_splitters(reference_file_name, _segment_size, no_threads))
    {
        working_mode = working_mode_t::none;
        return false;
    }

    out_archive = make_shared<CArchive>(false, 32 << 20);
    if (!out_archive->Open(_file_name))
    {
        working_mode = working_mode_t::none;
        return false;
    }
    
    working_mode = working_mode_t::compression;

    if (archive_version >= 3000 && archive_version < 4000)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->set_archives(nullptr, out_archive, no_threads, pack_cardinality, segment_size, kmer_length);

    no_samples_in_archive = 0;

    map_segments[std::make_pair(~0ull, ~0ull)] = 0;

    v_segments.resize(no_raw_groups);

    for (no_segments = 0; no_segments < no_raw_groups; ++no_segments)
    {
        contig_t empty_ctg{ 0x7f };

        out_archive->RegisterStream(ss_delta_name(archive_version, no_segments));

        v_segments[no_segments] = make_shared<CSegment>(ss_base(archive_version, no_segments), nullptr, out_archive, pack_cardinality, min_match_len, concatenated_genomes, archive_version);
        v_segments[no_segments]->add_raw(empty_ctg, nullptr, nullptr);		// To ensure that raw (special) segments are present in the archive
    }

    if (archive_version >= 3000)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->reset_prev_sample_name();

    return true;
}

// *******************************************************************************************
bool CAGCCompressor::Append(const string& _in_archive_fn, const string& _out_archive_fn, const uint32_t _verbosity, const bool _prefetch_archive, const bool _concatenated_genomes, const bool _adaptive_compression,
    const uint32_t no_threads)
{
    if (working_mode != working_mode_t::none)
        return false;

    in_archive_name = _in_archive_fn;
    out_archive_name = _out_archive_fn;
    prefetch_archive = _prefetch_archive;
    concatenated_genomes = _concatenated_genomes;
    adaptive_compression = _adaptive_compression;

    min_match_len = compression_params.min_match_len;
    uint32_t segment_size = compression_params.segment_size;
    verbosity = _verbosity;

    if (!load_file_type_info(in_archive_name))
        return false;

    working_mode = working_mode_t::pre_appending;

    if (!load_metadata())
        return false;
    working_mode = working_mode_t::appending;

    out_archive = make_shared<CArchive>(false, 32 << 20);

    if (!out_archive->Open(out_archive_name))
        return false;

    // !!! TODO (future): Add moving part of archive to the new one
    if (archive_version >= 3000 && archive_version < 4000)
        dynamic_pointer_cast<CCollection_V3>(collection_desc)->set_archives(in_archive, out_archive, no_threads, pack_cardinality, segment_size, kmer_length);

    no_samples_in_archive = collection_desc->get_no_samples();

    if (adaptive_compression)
        build_candidate_kmers_from_archive(no_threads);

    appending_init();

    working_mode = working_mode_t::appending;

    return true;
}

// *******************************************************************************************
void CAGCCompressor::AddCmdLine(const string& cmd_line)
{
    if (working_mode != working_mode_t::compression && working_mode != working_mode_t::appending)
        return;

    collection_desc->add_cmd_line(cmd_line);
}

// *******************************************************************************************
bool CAGCCompressor::Close(const uint32_t no_threads)
{
    bool r = true;

    if (working_mode == working_mode_t::none)
        r = false;
    else if (working_mode == working_mode_t::compression)
        r = close_compression(no_threads);
    else if (working_mode == working_mode_t::appending)
        r = close_compression(no_threads);

    working_mode = working_mode_t::none;

    return r;
}

// EOF
