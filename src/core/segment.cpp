// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-03-16
// *******************************************************************************************

#include "../core/segment.h"

// *******************************************************************************************
uint32_t CSegment::add_raw(const contig_t& s, bool buffered, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);

    if (internal_state == internal_state_t::packed)
        unpack(zstd_dctx);

    if (v_raw.size() == contigs_in_pack)
    {
        store_in_archive(v_raw, buffered, zstd_cctx);
        v_raw.clear();
    }

    ++no_seqs;
    v_raw.emplace_back(s);

    return (uint32_t) (no_seqs - 1u);
}

// *******************************************************************************************
uint32_t CSegment::add(const contig_t& s, bool buffered, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);

    if (internal_state == internal_state_t::packed)
        unpack(zstd_dctx);

    if (no_seqs == 0)
    {
        lz_diff->Prepare(s);

        store_in_archive(s, buffered, zstd_cctx);

        ref_size = s.size() + 1;
    }
    else
    {
        if (v_lzp.size() == contigs_in_pack)
        {
            store_in_archive(v_lzp, buffered, zstd_cctx);
            v_lzp.clear();
        }

        contig_t delta;

        lz_diff->Encode(s, delta);

#ifdef IMPROVED_LZ_ENCODING
        if (delta.empty())       // same sequence as reference
            return 0;
#endif

        auto p = find(v_lzp.begin(), v_lzp.end(), delta);
        if (p != v_lzp.end())
            return no_seqs - distance(p, v_lzp.end());

        v_lzp.push_back(delta);

        seq_size += s.size() + 1;
        packed_size += delta.size() + 1;
    }

    ++no_seqs;

    return no_seqs - 1u;
}

// *******************************************************************************************
uint64_t CSegment::estimate(const contig_t& s, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);

    if (internal_state == internal_state_t::packed)
        unpack(zstd_dctx);

    if (ref_size == 0)
        return 0;
    else
    {
        contig_t delta;
        lz_diff->Encode(s, delta);
        return delta.size();
    }
}

// *******************************************************************************************
void CSegment::get_coding_cost(const contig_t& s, vector<uint32_t>& v_costs, const bool prefix_costs, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);

    if (internal_state == internal_state_t::packed)
        unpack(zstd_dctx);

    if (ref_size == 0)
        return;
    else
        lz_diff->GetCodingCostVector(s, v_costs, prefix_costs);
}

// *******************************************************************************************
void CSegment::finish(bool buffered, ZSTD_CCtx* zstd_ctx)
{
    if (!v_lzp.empty())
        store_in_archive(v_lzp, buffered, zstd_ctx);
    if (!v_raw.empty())
        store_in_archive(v_raw, buffered, zstd_ctx);
    if (!packed_delta.empty())
        store_compressed_delta_in_archive(buffered);
}

// *******************************************************************************************
bool CSegment::get_raw(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
    // Retrive pack of delta-coded contigs
    vector<uint8_t> raw_seq;
    vector<uint8_t> pack_raw_seq;
    vector<uint8_t> zstd_raw_seq;
    uint64_t raw_seq_size;

    stream_id_delta = in_archive->GetStreamId(name + "-delta");
    int part_id = id_seq / contigs_in_pack;
    int seq_in_part_id = id_seq % contigs_in_pack;

    in_archive->GetPart(stream_id_delta, part_id, zstd_raw_seq, raw_seq_size);

    if (raw_seq_size == 0)
        pack_raw_seq = move(zstd_raw_seq);
    else
    {
        pack_raw_seq.resize(raw_seq_size);
        ZSTD_decompressDCtx(zstd_ctx, pack_raw_seq.data(), pack_raw_seq.size(), zstd_raw_seq.data(), zstd_raw_seq.size());
    }

    // Retrive the requested delta-coded contig
    uint32_t b_pos = 0;
    uint32_t e_pos = 0;
    int cnt = 0;

    for (uint32_t i = 0; i < pack_raw_seq.size(); ++i)
    {
        if (pack_raw_seq[i] == contig_separator)
        {
            ++cnt;
            if (cnt == seq_in_part_id)
                b_pos = i + 1;
            else if (cnt == seq_in_part_id + 1)
            {
                e_pos = i;
                break;
            }
        }
    }

    ctg.assign(pack_raw_seq.begin() + b_pos, pack_raw_seq.begin() + e_pos);

    return true;
}

// *******************************************************************************************
bool CSegment::get(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
    // Retrive reference contig
    vector<uint8_t> ref_seq;
    vector<uint8_t> zstd_ref_seq;
    uint64_t ref_seq_size;

    stream_id_ref = in_archive->GetStreamId(name + "-ref");
    in_archive->GetPart(stream_id_ref, 0, zstd_ref_seq, ref_seq_size);

    if (ref_seq_size == 0)
        ref_seq = move(zstd_ref_seq);       // No compression
    else
    {
        ref_seq.resize(ref_seq_size);

        if (zstd_ref_seq.back() == 0)
            ZSTD_decompressDCtx(zstd_ctx, ref_seq.data(), ref_seq.size(), zstd_ref_seq.data(), zstd_ref_seq.size() - 1u);
        else
        {
            vector<uint8_t> v_tuples;
            v_tuples.resize(ref_seq_size+1);

            auto output_size = ZSTD_decompressDCtx(zstd_ctx, v_tuples.data(), v_tuples.size(), zstd_ref_seq.data(), zstd_ref_seq.size() - 1u);
            v_tuples.resize(output_size);
            tuples2bytes(v_tuples, ref_seq);
        }
    }

    if (id_seq == 0)
    {
        ctg = move(ref_seq);
        return true;
    }

    // Retrive pack of delta-coded contigs
    vector<uint8_t> delta_seq;
    vector<uint8_t> pack_delta_seq;
    vector<uint8_t> zstd_delta_seq;
    uint64_t delta_seq_size;

    stream_id_delta = in_archive->GetStreamId(name + "-delta");
    int part_id = (id_seq - 1) / contigs_in_pack;
    int seq_in_part_id = (id_seq - 1) % contigs_in_pack;

    in_archive->GetPart(stream_id_delta, part_id, zstd_delta_seq, delta_seq_size);

    if (delta_seq_size == 0)
        pack_delta_seq = move(zstd_delta_seq);
    else
    {
        pack_delta_seq.resize(delta_seq_size);
        ZSTD_decompressDCtx(zstd_ctx, pack_delta_seq.data(), pack_delta_seq.size(), zstd_delta_seq.data(), zstd_delta_seq.size());
    }

    if (contigs_in_pack > 1)
    {
        // Retrive the requested delta-coded contig
        uint32_t b_pos = 0;
        uint32_t e_pos = 0;
        int cnt = 0;

        for (uint32_t i = 0; i < pack_delta_seq.size(); ++i)
        {
            if (pack_delta_seq[i] == contig_separator)
            {
                ++cnt;
                if (cnt == seq_in_part_id)
                    b_pos = i + 1;
                else if (cnt == seq_in_part_id + 1)
                {
                    e_pos = i;
                    break;
                }
            }
        }

        delta_seq.assign(pack_delta_seq.begin() + b_pos, pack_delta_seq.begin() + e_pos);
    }
    else
        delta_seq.assign(pack_delta_seq.begin(), pack_delta_seq.end() - 1);
    
    // LZ decode delta-encoded contig
    ctg.clear();
    lz_diff->Decode(ref_seq, delta_seq, ctg);

    return true;
}

// *******************************************************************************************
void CSegment::appending_init()
{
    if (internal_state != internal_state_t::none)
        return;

    // Retrive reference contig
    contig_t ref_seq;
    
    int in_stream_id_ref = in_archive->GetStreamId(name + "-ref");
    int in_stream_id_delta = in_archive->GetStreamId(name + "-delta");

    int out_stream_id_ref = -1;
    int out_stream_id_delta = -1;

    if(in_stream_id_ref >= 0)
        out_stream_id_ref = out_archive->RegisterStream(name + "-ref");
    if(in_stream_id_delta >= 0)
        out_stream_id_delta = out_archive->RegisterStream(name + "-delta");

    // Copy of all parts except last one
    if (in_stream_id_ref >= 0)
    {
        in_archive->GetPart(in_stream_id_ref, packed_ref_seq, raw_ref_seq_size);
        out_archive->AddPart(out_stream_id_ref, packed_ref_seq, raw_ref_seq_size);
        ref_transferred = true;
        no_seqs = 1;
    }
    else
    {
        packed_ref_seq.clear();
        no_seqs = 0;
    }

    if (in_stream_id_delta >= 0)
    {
        vector<uint8_t> tmp_data;
        uint64_t tmp_meta;

        uint32_t no_parts = (uint32_t) in_archive->GetNoParts(in_stream_id_delta);
        for (uint32_t i = 0; i + 1 < no_parts; ++i)
        {
            in_archive->GetPart(in_stream_id_delta, tmp_data, tmp_meta);
            out_archive->AddPart(out_stream_id_delta, tmp_data, tmp_meta);
            no_seqs += contigs_in_pack;
        }

        in_archive->GetPart(in_stream_id_delta, packed_delta, raw_delta_size);
    }

    internal_state = internal_state_t::packed;

    stream_id_ref = out_stream_id_ref;
    stream_id_delta = out_stream_id_delta;
}

// *******************************************************************************************
void CSegment::clear()
{
    lock_guard<mutex> lck(mtx);

    no_seqs = 0;
    ref_size = 0;
    seq_size = packed_size = 0;
    v_raw.clear();
    v_lzp.clear();

    internal_state = internal_state_t::normal;
}

// *******************************************************************************************
uint64_t CSegment::get_no_seqs()
{
    lock_guard<mutex> lck(mtx);

    return no_seqs;
}

// *******************************************************************************************
void CSegment::unpack(ZSTD_DCtx* zstd_ctx)
{
    if (!packed_ref_seq.empty())
    {
        contig_t ref_seq;

        if (raw_ref_seq_size == 0)
            ref_seq = move(packed_ref_seq);       // No compression
        else
        {
            vector<uint8_t> v_tuples;
            v_tuples.resize(raw_ref_seq_size + 1);

            auto output_size = ZSTD_decompressDCtx(zstd_ctx, v_tuples.data(), v_tuples.size(), packed_ref_seq.data(), packed_ref_seq.size() - 1u);

            v_tuples.resize(output_size);

            if (packed_ref_seq.back() == 1)      // marker (tuples to bytes conversion needed)
                tuples2bytes(v_tuples, ref_seq);
            else
                ref_seq.swap(v_tuples);
        }

        packed_ref_seq.clear();
        packed_ref_seq.shrink_to_fit();

        lz_diff->Prepare(ref_seq);
        ref_size = ref_seq.size() + 1;
    }

    if (!packed_delta.empty())
    {
        contig_t delta_seq;

        if (raw_delta_size == 0)
            delta_seq = move(packed_delta);
        else
        {
            delta_seq.resize(raw_delta_size);
            ZSTD_decompressDCtx(zstd_ctx, delta_seq.data(), delta_seq.size(), packed_delta.data(), packed_delta.size());
        }

        packed_delta.clear();
        packed_delta.shrink_to_fit();

        v_lzp.clear();

        // Retrive the requested delta-coded contig
        uint32_t b_pos = 0;

        if (contigs_in_pack > 1)
        {
            for (uint32_t i = 0; i < delta_seq.size(); ++i)
                if (delta_seq[i] == contig_separator)
                {
                    v_lzp.emplace_back(delta_seq.begin() + b_pos, delta_seq.begin() + i);
                    b_pos = i + 1;
                }
        }
        else
            v_lzp.emplace_back(delta_seq.begin(), delta_seq.begin() + (delta_seq.size() - 1));

        no_seqs += (uint32_t) v_lzp.size();

        if (ref_size == 0)          // There is no reference sequence so the deltas are in fact raw sequences
            swap(v_raw, v_lzp);
    }

    internal_state = internal_state_t::normal;
}

// EOF
