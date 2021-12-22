// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include <iostream>
#include "../core/agc_basic.h"
#include <atomic>
#include <chrono>

using namespace std::chrono;

// *******************************************************************************************
CAGCBasic::CAGCBasic()
{
    compression_params.kmer_length = 0;
    compression_params.min_match_len = 0;
    compression_params.pack_cardinality = 1;

    pack_cardinality = 1;

    working_mode = working_mode_t::none;

    is_app_mode = true;
    kmer_length = 0;
    min_match_len = 0;

    verbosity = 0;
}

// *******************************************************************************************
CAGCBasic::~CAGCBasic()
{    
}

// *******************************************************************************************
bool CAGCBasic::load_metadata(const string &archive_name)
{
    if (prefetch_archive)
        in_archive = make_shared<CArchive>(true, ~0ull);           // ~0ull - special value - buffers whole archive
    else
        in_archive = make_shared<CArchive>(true, 32 << 10);        

    if (!in_archive->Open(archive_name))
    {
        if (is_app_mode)
            cerr << "Cannot open archive: " << archive_name << endl;
        return false;
    }
    
    int desc_sid = in_archive->GetStreamId("collection-desc");
    if (desc_sid < 0)
    {
        in_archive->Close();
        if (is_app_mode)
            cerr << "No collection description in archive" << endl;
        return false;
    }

    vector<uint8_t> v_desc_zstd;
    uint64_t tmp;

    if (!in_archive->GetPart(desc_sid, v_desc_zstd, tmp))
    {
        in_archive->Close();
        if (is_app_mode)
            cerr << "Problem with archive\n";
        return false;
    }

    vector<uint8_t> v_desc;
    v_desc.resize(tmp);
    ZSTD_decompress(v_desc.data(), v_desc.size(), v_desc_zstd.data(), v_desc_zstd.size());
    v_desc_zstd.clear();
    v_desc_zstd.shrink_to_fit();

    if (!collection_desc.deserialize(v_desc))
    {
        in_archive->Close();
        if (is_app_mode)
            cerr << "Cannot deserialize\n";
        return false;
    }

    vector<uint8_t> v_params;
    
    if (!in_archive->GetPart(in_archive->GetStreamId("params"), v_params, tmp))
    {
        in_archive->Close();
        if (is_app_mode)
            cerr << "Archive does not contain parameters section\n";
        return false;
    }

    auto p = v_params.begin();
    read(p, compression_params.kmer_length);
    read(p, compression_params.min_match_len);
    read(p, compression_params.pack_cardinality);

    kmer_length = compression_params.kmer_length;
    pack_cardinality = compression_params.pack_cardinality;
    min_match_len = compression_params.min_match_len;

    return true;
}

// *******************************************************************************************
void CAGCBasic::join_threads(vector<thread> &v_threads)
{
    for (auto& t : v_threads)
        t.join();

    v_threads.clear();
}

// *******************************************************************************************
void CAGCBasic::reverse_complement(contig_t& contig)
{
    reverse(contig.begin(), contig.end());
    for (auto& x : contig)
        if (x < 4)
            x = 3 - x;
}

// *******************************************************************************************
void CAGCBasic::reverse_complement_copy(const contig_t& src_contig, contig_t& dest_contig)
{
    dest_contig.clear();
    dest_contig.reserve(src_contig.size());

    for (auto p = src_contig.rbegin(); p != src_contig.rend(); ++p)
        dest_contig.emplace_back((*p < 4) ? 3 - *p : *p);
}

// EOL