// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
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
    compression_params.segment_size = 0;

    pack_cardinality = 1;

    working_mode = working_mode_t::none;

    is_app_mode = true;
    kmer_length = 0;
    min_match_len = 0;

    archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;

    if (archive_version < 2000)
        collection_desc = static_pointer_cast<CCollection>(make_shared<CCollection_V1>());
    else if(archive_version < 3000)
        collection_desc = static_pointer_cast<CCollection>(make_shared<CCollection_V2>());
    else if(archive_version < 4000)
        collection_desc = static_pointer_cast<CCollection>(make_shared<CCollection_V3>());

    verbosity = 0;
}

// *******************************************************************************************
CAGCBasic::~CAGCBasic()
{    
}

// *******************************************************************************************
bool CAGCBasic::load_file_type_info(const string& archive_name)
{
    vector<uint8_t> v_data;

    if (prefetch_archive)
        in_archive = make_shared<CArchive>(true, ~0ull, ss_prefix(archive_version));           // ~0ull - special value - buffers whole archive
    else
        in_archive = make_shared<CArchive>(true, 32 << 10, ss_prefix(archive_version));

    if (!in_archive->Open(archive_name))
    {
        if (is_app_mode)
            cerr << "Cannot open archive: " << archive_name << endl;
        return false;
    }

    m_file_type_info.clear();

    auto s_id = in_archive->GetStreamId("file_type_info");
    if (s_id < 0)
        return false;

    uint64_t n_items;

    if (!in_archive->GetPart(s_id, v_data, n_items))
        return false;

    auto p = v_data.begin();
    string key, val;

    for (size_t i = 0; i < n_items; ++i)
    {
        read(p, key);
        read(p, val);

        m_file_type_info.emplace(key, val);
    }

    archive_version = stoi(m_file_type_info["file_version_major"]) * 1000 + stoi(m_file_type_info["file_version_minor"]);

    if (archive_version < 2000)
        collection_desc = static_pointer_cast<CCollection>(make_shared<CCollection_V1>());
    else if (archive_version < 3000)
        collection_desc = static_pointer_cast<CCollection>(make_shared<CCollection_V2>());
    else if (archive_version < 4000)
        collection_desc = static_pointer_cast<CCollection>(make_shared<CCollection_V3>());

    return true;
}

// *******************************************************************************************
bool CAGCBasic::load_metadata_impl_v1()
{
    int desc_sid_v1 = in_archive->GetStreamId("collection-desc");

    vector<uint8_t> v_desc_zstd;
    uint64_t tmp;

    if (!in_archive->GetPart(desc_sid_v1, v_desc_zstd, tmp))
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

    auto collection_desc_v1 = dynamic_pointer_cast<CCollection_V1>(collection_desc);

    if (!collection_desc_v1->deserialize(v_desc))
    {
        in_archive->Close();
        if (is_app_mode)
            cerr << "Cannot deserialize\n";
        return false;
    }

    return true;
}

// *******************************************************************************************
bool CAGCBasic::load_metadata_impl_v2()
{
    int desc_main_sid_v2 = in_archive->GetStreamId("collection-main");
    int desc_details_sid_v2 = in_archive->GetStreamId("collection-details");

    vector<uint8_t> v_desc_zstd;
    uint64_t tmp;

    if (!in_archive->GetPart(desc_main_sid_v2, v_desc_zstd, tmp))
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

    bool expensive_collection_structures;
    if (is_app_mode)
    {
        if (working_mode == working_mode_t::decompression || working_mode == working_mode_t::none)
            expensive_collection_structures = false;
        else
            expensive_collection_structures = true;
    }
    else
        expensive_collection_structures = prefetch_archive;

    auto collection_desc_v2 = dynamic_pointer_cast<CCollection_V2>(collection_desc);

    if (!collection_desc_v2->deserialize_main(v_desc, expensive_collection_structures))
    {
        in_archive->Close();
        if (is_app_mode)
            cerr << "Cannot deserialize\n";
        return false;
    }

    while (in_archive->GetPart(desc_details_sid_v2, v_desc_zstd, tmp))
        collection_desc_v2->deserialize_details(v_desc_zstd, tmp,
            working_mode == working_mode_t::appending || working_mode == working_mode_t::pre_appending);

    return true;
}

// *******************************************************************************************
bool CAGCBasic::load_metadata_impl_v3()
{
    return true;

    // Do nothing here - loading is made in other place
}

// *******************************************************************************************
bool CAGCBasic::load_metadata()
{    
    if (archive_version >= 4000)
    {
        archive_version = 0;        // Invalid archive

        in_archive->Close();
        if (is_app_mode)
            cerr << "Unsupported archive version. Please use the most recent AGC application" << endl;
        return false;
    }

    if (archive_version < 2000)             // v1
        load_metadata_impl_v1();
    else if (archive_version < 3000)        // v2
        load_metadata_impl_v2();
    else if (archive_version < 4000)        // v3
        load_metadata_impl_v3();

    uint64_t tmp;
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

    if (archive_version >= 2000)
        read(p, compression_params.segment_size);
    else
        compression_params.segment_size = 0;

    kmer_length = compression_params.kmer_length;
    pack_cardinality = compression_params.pack_cardinality;
    min_match_len = compression_params.min_match_len;
    segment_size = compression_params.segment_size;

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
    if (contig.empty())
        return;

    int size = (int) contig.size();

    uint8_t* p = contig.data() + size - 1;
    uint8_t* q = contig.data();

    for (int i = 0; i < size / 2; ++i)
    {
        uint8_t x = (*p < 4) ? 3 - *p : *p;
        uint8_t y = (*q < 4) ? 3 - *q : *q;

        *q++ = x;
        *p-- = y;
    }

    if (size % 2)
        if (*p < 4)
            *p = 3 - *p;
}

// *******************************************************************************************
void CAGCBasic::reverse_complement_copy(contig_t& src_contig, contig_t& dest_contig)
{
    int size = (int) src_contig.size();
    
    dest_contig.resize(size);

    if (src_contig.empty())
        return;

    int i = 0;
    uint8_t* p = src_contig.data() + size - 1;
    uint8_t* q = dest_contig.data();

    switch (size % 4)
    {
    case 3:
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
        ++i;
    case 2:
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
        ++i;
    case 1:
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
        ++i;
    }

    for (; i < size; i += 4)
    {
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
        *q++ = (*p < 4) ? 3 - *p : *p;        --p;
    }
}

// EOL