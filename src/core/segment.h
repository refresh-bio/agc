#ifndef _SEGMENT_H
#define _SEGMENT_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include <vector>
#include <string>
#include <cinttypes>
#include <mutex>
#include <memory>
#include "../../libs/zstd.h"
#include "../core/lz_diff.h"
#include "../core/archive.h"
#include "../core/defs.h"

using namespace std;

class CSegment
{
    enum class internal_state_t {none, normal, packed};

    const uint8_t contig_separator = 0xffu;

    string name;
    shared_ptr<CArchive> in_archive;
    shared_ptr<CArchive> out_archive;
    uint32_t contigs_in_pack;
    uint32_t min_match_len;
    bool concatenated_genomes;

    int stream_id_ref;
    int stream_id_delta;

//    const int ref_compression_level = 15;
//    const int delta_compression_level = 19;

    internal_state_t internal_state;

    packed_block_t packed_ref_seq;
    packed_block_t packed_delta;
    uint64_t raw_ref_seq_size = 0;
    uint64_t raw_delta_size = 0;
    bool ref_transferred = false;

    CLZDiff lz_diff;

    uint32_t no_seqs;
    vector<contig_t> v_lzp;
    vector<contig_t> v_raw;

    uint64_t ref_size;
    uint64_t seq_size;
    uint64_t packed_size;
    mutex mtx;

    // *******************************************************************************************
    void bytes2tuples(const vector<uint8_t>& v_bytes, vector<uint8_t>& v_tuples)
    {
        uint8_t me = 0;
        
        if(!v_bytes.empty())
            me = *max_element(v_bytes.begin(), v_bytes.end());

        if (me < 4)
            bytes2tuples_impl<4, 4>(v_bytes, v_tuples);
        else if (me < 6)
            bytes2tuples_impl<3, 6>(v_bytes, v_tuples);
        else if (me < 16)
            bytes2tuples_impl<2, 16>(v_bytes, v_tuples);
        else
        {
            v_tuples = v_bytes;
            v_tuples.emplace_back(0x10u);
        }
    }

    // *******************************************************************************************
    void tuples2bytes(const vector<uint8_t>& v_tuples, vector<uint8_t>& v_bytes)
    {
        uint8_t marker = v_tuples.back();
        uint8_t no_bytes = marker >> 4;
        uint8_t trailing_bytes = marker & 0xf;
        uint32_t output_size = (uint32_t)((v_tuples.size() - 2) * no_bytes + trailing_bytes);

        v_bytes.reserve(output_size);

        if (no_bytes == 4)
            tuples2bytes_impl<4, 4>(v_tuples, v_bytes, output_size);
        else if (no_bytes == 3)
            tuples2bytes_impl<3, 6>(v_tuples, v_bytes, output_size);
        else if (no_bytes == 2)
            tuples2bytes_impl<2, 16>(v_tuples, v_bytes, output_size);
        else
            v_bytes.assign(v_tuples.begin(), v_tuples.begin() + v_tuples.size() - 1u);        
    }

    // *******************************************************************************************
    template<unsigned NO_BYTES, unsigned MULT>
    void bytes2tuples_impl(const vector<uint8_t>& v_bytes, vector<uint8_t>& v_tuples)
    {
        v_tuples.reserve((v_bytes.size() - NO_BYTES - 1u) / NO_BYTES + 1);

        size_t i;
        uint8_t c;

        for (i = 0; i + NO_BYTES <= v_bytes.size(); i += NO_BYTES)
        {
            c = 0;

            for (uint32_t j = 0; j < NO_BYTES; ++j)
                c = c * MULT + v_bytes[i + j];

            v_tuples.emplace_back(c);
        }
                
        for (c = 0; i < v_bytes.size(); ++i)
            c = c * MULT + v_bytes[i];

        v_tuples.emplace_back(c);

        v_tuples.emplace_back((NO_BYTES << 4) + (v_bytes.size() % NO_BYTES));        // marker for decoding
    }

    // *******************************************************************************************
    template<unsigned NO_BYTES, unsigned MULT>
    void tuples2bytes_impl(const vector<uint8_t>& v_tuples, vector<uint8_t>& v_bytes, const uint32_t output_size)
    {
        uint32_t i, j;

        v_bytes.resize(output_size);

        for (i = j = 0; j + NO_BYTES <= output_size; ++i, j += NO_BYTES)
        {
            uint8_t c = v_tuples[i];

            for (int k = NO_BYTES - 1u; k >= 0; --k)
            {
                v_bytes[(size_t) j + k] = c % MULT;
                c /= MULT;
            }
        }

        uint8_t c = v_tuples[i];

        uint32_t n = output_size % NO_BYTES;

        if(n)
            for (int k = n-1u; k >= 0; --k)
            {
                v_bytes[(size_t) j + k] = c % MULT;
                c /= MULT;
            }
    }

    // *******************************************************************************************
    void add_to_archive(const int stream_id, const contig_t& data, const int compression_level, bool buffered, ZSTD_CCtx* zstd_ctx)
    {
        size_t a_size = (uint32_t) ZSTD_compressBound(data.size());
        uint8_t *packed = new uint8_t[a_size+1u];
        uint32_t packed_size = (uint32_t) ZSTD_compressCCtx(zstd_ctx, (void *) packed, a_size, data.data(), data.size(), compression_level);
        packed[packed_size] = 0;      // ZSTD compression marker - plain (0)

        if(packed_size + 1u < (uint32_t) data.size())
        {
            vector<uint8_t> v_packed(packed, packed + packed_size + 1);
            if(buffered)
                out_archive->AddPartBuffered(stream_id, v_packed, data.size());
            else
                out_archive->AddPart(stream_id, v_packed, data.size());
        }
        else
        {
            if(buffered)
                out_archive->AddPartBuffered(stream_id, data, 0);
            else
                out_archive->AddPart(stream_id, data, 0);
        }

        delete[] packed;
    }

    // *******************************************************************************************
    void add_to_archive_tuples(const int stream_id, const contig_t& data, const int compression_level, bool buffered, ZSTD_CCtx* zstd_ctx)
    {
        vector<uint8_t> v_tuples;

        bytes2tuples(data, v_tuples);

        size_t a_size = (uint32_t) ZSTD_compressBound(v_tuples.size());
        uint8_t *packed = new uint8_t[a_size+1u];
        uint32_t packed_size = (uint32_t) ZSTD_compressCCtx(zstd_ctx, (void *) packed, a_size, v_tuples.data(), v_tuples.size(), compression_level);
        packed[packed_size] = 1;      // ZSTD compression marker - tuples (1)

        if(packed_size + 1u < (uint32_t) data.size())
        {
            vector<uint8_t> v_packed(packed, packed + packed_size + 1);
            if(buffered)
                out_archive->AddPartBuffered(stream_id, v_packed, data.size());
            else
                out_archive->AddPart(stream_id, v_packed, data.size());
        }
        else
        {
            if(buffered)
                out_archive->AddPartBuffered(stream_id, data, 0);
            else
                out_archive->AddPart(stream_id, data, 0);
        }

        delete[] packed;
    }

    // *******************************************************************************************
    void store_in_archive(const contig_t& data, bool buffered, ZSTD_CCtx* zstd_ctx)
    { 
        string stream_name = name + "-ref";

        stream_id_ref = out_archive->RegisterStream(stream_name);

//        uint32_t best_len = 0;
//        uint32_t best_cnt = 0;
        double best_frac = 0.0;
        double frac_limit = 0.5;

        for (uint32_t i = 4; i < 32; ++i)
        {
            uint32_t cnt = 0;
            uint32_t cur_size = 0;

            for (uint32_t j = 0; (size_t) j + i < data.size(); ++j)
            {
                cnt += data[j] == data[(size_t) j + i];
                cur_size += data[j] < 4;            // exclude non-ACGT from counting
            }

            double frac = 0.0;
            if (cur_size)
                frac = (double)cnt / cur_size;

            if (frac > best_frac)
            {
//                best_len = i;
//                best_cnt = cnt;
                best_frac = frac;

                if (best_frac >= frac_limit)
                    break;
            }
        }

        if (best_frac < 0.5)
            add_to_archive_tuples(stream_id_ref, data, 13, buffered, zstd_ctx);
        else
            add_to_archive(stream_id_ref, data, 19, buffered, zstd_ctx);
    }

    // *******************************************************************************************
    void store_in_archive(const vector<contig_t>& v_data, bool buffered, ZSTD_CCtx* zstd_ctx)
    {
        contig_t pack;

        size_t res_size = v_data.size();
        for (const auto& x : v_data)
            res_size += x.size();

        pack.reserve(res_size);

        for (auto& x : v_data)
        {
            pack.insert(pack.end(), x.begin(), x.end());
            pack.push_back(contig_separator);
        }

        if (stream_id_delta < 0)
            stream_id_delta = out_archive->RegisterStream(name + "-delta");

//        add_to_archive(stream_id_delta, pack, delta_compression_level, buffered, zstd_ctx);
        add_to_archive(stream_id_delta, pack, 17, buffered, zstd_ctx);
    }

    // *******************************************************************************************
    void store_compressed_delta_in_archive(bool buffered)
    {
        if (stream_id_delta < 0)
        {
            string stream_name = name + "-delta";
            stream_id_delta = out_archive->RegisterStream(stream_name);
        }

        if(buffered)
            out_archive->AddPartBuffered(stream_id_delta, packed_delta, raw_delta_size);
        else
            out_archive->AddPart(stream_id_delta, packed_delta, raw_delta_size);
    }

    void unpack(ZSTD_DCtx* zstd_ctx);

public:
    // *******************************************************************************************
    CSegment(const string &_name, shared_ptr<CArchive> _in_archive, shared_ptr<CArchive> _out_archive,
        const uint32_t _contigs_in_pack, const uint32_t _min_match_len, const bool _concatenated_genomes) :
        name(_name), in_archive(_in_archive), out_archive(_out_archive), 
        contigs_in_pack(_contigs_in_pack), min_match_len(_min_match_len), concatenated_genomes(_concatenated_genomes),
        no_seqs(0), ref_size(0), seq_size(0), packed_size(0)
    {
        stream_id_ref = -1;
        stream_id_delta = -1;
        internal_state = internal_state_t::none;

        lz_diff.SetMinMatchLen(min_match_len);
    };

    ~CSegment()
    {
    }

    uint32_t add_raw(const contig_t& s, bool buffered, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx);
    uint32_t add(const contig_t& s, bool buffered, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx);
    uint64_t estimate(const contig_t& s, ZSTD_DCtx* zstd_dctx);

    void get_coding_cost(const contig_t& s, vector<uint32_t> &v_costs, const bool prefix_costs, ZSTD_DCtx* zstd_dctx);

    void finish(bool buffered, ZSTD_CCtx* zstd_ctx);
    bool get_raw(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx);
    bool get(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx);
    void clear();
    uint64_t get_no_seqs();

    void appending_init();
};

// EOF
#endif