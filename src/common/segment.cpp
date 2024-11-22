// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "segment.h"

// *******************************************************************************************
uint32_t CSegment::add_raw(const contig_t& s, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);

    if (internal_state == internal_state_t::packed)
        unpack(zstd_dctx);

    if (v_raw.size() == contigs_in_pack)
    {
        store_in_archive(v_raw, zstd_cctx);
        v_raw.clear();
    }

    ++no_seqs;
    v_raw.emplace_back(s);

    return (uint32_t) (no_seqs - 1u);
}

// *******************************************************************************************
uint32_t CSegment::add(const contig_t& s, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);

    if (internal_state == internal_state_t::packed)
        unpack(zstd_dctx);

    if (no_seqs == 0)
    {
        lz_diff->Prepare(s);

        store_in_archive(s, zstd_cctx);

        ref_size = s.size() + 1;
    }
    else
    {
        if (v_lzp.size() == contigs_in_pack)
        {
            store_in_archive(v_lzp, zstd_cctx);
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

        seq_size += s.size() + 1;
        packed_size += delta.size() + 1;

        v_lzp.emplace_back(move(delta));
    }

    ++no_seqs;

    return no_seqs - 1u;
}

// *******************************************************************************************
uint64_t CSegment::estimate(const contig_t& s, uint32_t bound, ZSTD_DCtx* zstd_dctx)
{
    if (ref_size == 0)
        return 0;

    {
        lock_guard<mutex> lck(mtx);

        if (internal_state == internal_state_t::packed)
            unpack(zstd_dctx);

        lz_diff->AssureIndex();
    }

    return lz_diff->Estimate(s, bound);
}

// *******************************************************************************************
void CSegment::get_coding_cost(const contig_t& s, vector<uint32_t>& v_costs, const bool prefix_costs, ZSTD_DCtx* zstd_dctx)
{
    if (ref_size == 0)
        return;

    {
        //    lock_guard<mutex> lck(mtx);
        lock_guard<mutex> lck(mtx);

        if (internal_state == internal_state_t::packed)
            unpack(zstd_dctx);
        lz_diff->AssureIndex();
    }

    lz_diff->GetCodingCostVector(s, v_costs, prefix_costs);
}

// *******************************************************************************************
size_t CSegment::get_ref_size() const
{
    return ref_size;
}

// *******************************************************************************************
void CSegment::finish(ZSTD_CCtx* zstd_ctx)
{
    if (!v_lzp.empty())
        store_in_archive(v_lzp, zstd_ctx);
    if (!v_raw.empty())
        store_in_archive(v_raw, zstd_ctx);
    if (!packed_delta.empty())
        store_compressed_delta_in_archive();
}

// *******************************************************************************************
bool CSegment::get_raw(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
    // Retrive pack of raw contigs
//    vector<uint8_t> pack_raw_seq;

    int part_id = id_seq / contigs_in_pack;
    int seq_in_part_id = id_seq % contigs_in_pack;

    uint8_t* pack_raw_seq = nullptr;
    contig_t buf;
    size_t pack_raw_seq_size = 0;

    vector<uint8_t> zstd_raw_seq;
    uint64_t raw_seq_size;

    if (!fast)
    {
        tie(stream_id_delta, ignore) = in_archive->GetPart(name + ss_delta_ext(archive_version), part_id, zstd_raw_seq, raw_seq_size);

        if (raw_seq_size == 0)
        {
            pack_raw_seq = zstd_raw_seq.data();
            pack_raw_seq_size = zstd_raw_seq.size();
        }
        else
        {
            buf.resize(raw_seq_size);
            ZSTD_decompressDCtx(zstd_ctx, buf.data(), buf.size(), zstd_raw_seq.data(), zstd_raw_seq.size());
            pack_raw_seq = buf.data();
            pack_raw_seq_size = buf.size();
        }
    }
    else
    {
        auto p_raw = pf_packed_raw_seq.find(part_id);
        
        if (p_raw == pf_packed_raw_seq.end())
        {
            tie(stream_id_delta, ignore) = in_archive->GetPart(name + ss_delta_ext(archive_version), part_id, zstd_raw_seq, raw_seq_size);

            if (pf_packed_raw_seq.size() >= pf_max_size)
                pf_packed_raw_seq.erase(pf_packed_raw_seq.begin());

            tie(p_raw, ignore) = pf_packed_raw_seq.insert(make_pair(part_id, vector<uint8_t>()));

            if (raw_seq_size == 0)
                p_raw->second = move(zstd_raw_seq);
            else
            {
                p_raw->second.resize(raw_seq_size);
                ZSTD_decompressDCtx(zstd_ctx, p_raw->second.data(), p_raw->second.size(), zstd_raw_seq.data(), zstd_raw_seq.size());
            }
        }

        pack_raw_seq = p_raw->second.data();
        pack_raw_seq_size = p_raw->second.size();
    }

    // Retrive the requested delta-coded contig
    uint32_t b_pos = 0;
    uint32_t e_pos = 0;
    int cnt = 0;

    for (uint32_t i = 0; i < pack_raw_seq_size; ++i)
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

    ctg.assign(pack_raw_seq + b_pos, pack_raw_seq + e_pos);

    return true;
}

// *******************************************************************************************
bool CSegment::get(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
    // Retrive reference contig
//    contig_t ref_seq;
    vector<uint8_t> zstd_ref_seq;

    vector<uint8_t> zstd_delta_seq;
    uint64_t delta_seq_size = 0;

    uint64_t ref_seq_size = 0;
    int part_id = (id_seq - 1) / contigs_in_pack;

    if (!fast)
    {
        tie(stream_id_ref, ignore, stream_id_delta, ignore) = in_archive->GetParts(
            name + ss_ref_ext(archive_version), 0, zstd_ref_seq, ref_seq_size,
            name + ss_delta_ext(archive_version), part_id, zstd_delta_seq, delta_seq_size);
    }
    else
    {
        if(ref_seq.empty())
            tie(stream_id_ref, ignore) = in_archive->GetPart(name + ss_ref_ext(archive_version), 0, zstd_ref_seq, ref_seq_size);

        auto p_delta = pf_packed_delta_seq.find(part_id);
        if (p_delta == pf_packed_delta_seq.end())
        {
            tie(stream_id_delta, ignore) = in_archive->GetPart(name + ss_delta_ext(archive_version), part_id, zstd_delta_seq, delta_seq_size);
            if (pf_packed_delta_seq.size() >= pf_max_size)
                pf_packed_delta_seq.erase(pf_packed_delta_seq.begin());
        }
    }

    if (ref_seq.empty())
    {
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
                v_tuples.resize(ref_seq_size + 1);

                auto output_size = ZSTD_decompressDCtx(zstd_ctx, v_tuples.data(), v_tuples.size(), zstd_ref_seq.data(), zstd_ref_seq.size() - 1u);

                v_tuples.resize(output_size);
                tuples2bytes(v_tuples, ref_seq);
            }
        }
    }

    if (id_seq == 0)
    {
        if (!fast)
        {
            ctg = move(ref_seq);
            ref_seq.clear();
            ref_seq.shrink_to_fit();
        }
        else
            ctg = ref_seq;

        return true;
    }

    uint8_t* pack_delta_seq;
    bool need_deallocate_pack_delta_seq = false;
    auto p_delta = pf_packed_delta_seq.begin();

    if (!fast)
    {
        if (delta_seq_size == 0)
        {
            pack_delta_seq = zstd_delta_seq.data();
            delta_seq_size = zstd_delta_seq.size();
        }
        else
        {
            pack_delta_seq = new uint8_t[delta_seq_size];
            need_deallocate_pack_delta_seq = true;
            ZSTD_decompressDCtx(zstd_ctx, pack_delta_seq, delta_seq_size, zstd_delta_seq.data(), zstd_delta_seq.size());
        }
    }
    else
    {
        p_delta = pf_packed_delta_seq.find(part_id);

        if (p_delta == pf_packed_delta_seq.end())
        {
            tie(p_delta, ignore) = pf_packed_delta_seq.insert(make_pair(part_id, make_pair(vector<uint8_t>(), vector<uint32_t>())));

            if (delta_seq_size == 0)
                p_delta->second.first = zstd_delta_seq;
            else
            {
                p_delta->second.first.resize(delta_seq_size);
                ZSTD_decompressDCtx(zstd_ctx, p_delta->second.first.data(), delta_seq_size, zstd_delta_seq.data(), zstd_delta_seq.size());
            }

            auto& sep_pos = p_delta->second.second;

            pack_delta_seq = p_delta->second.first.data();
            delta_seq_size = p_delta->second.first.size();

            if (contigs_in_pack == 1)
            {
                sep_pos.emplace_back(0);
                sep_pos.emplace_back(delta_seq_size);
            }
            else
            {
                sep_pos.emplace_back(0);

                for (uint32_t i = 0; i < delta_seq_size; ++i)
                    if (pack_delta_seq[i] == contig_separator)
                        sep_pos.emplace_back(i + 1);
            }
        }
        else
        {
            pack_delta_seq = p_delta->second.first.data();
            delta_seq_size = p_delta->second.first.size();
        }
    }

    // Retrive pack of delta-coded contigs
    contig_t delta_seq;
    int seq_in_part_id = (id_seq - 1) % contigs_in_pack;

    if (!fast)
    {
        if (contigs_in_pack > 1)
        {
            // Retrive the requested delta-coded contig
            uint32_t b_pos = 0;
            uint32_t e_pos = 0;
            int cnt = 0;

            for (uint32_t i = 0; i < delta_seq_size; ++i)
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

            delta_seq.assign(pack_delta_seq + b_pos, pack_delta_seq + e_pos);
        }
        else
            delta_seq.assign(pack_delta_seq, pack_delta_seq + delta_seq_size - 1);
    }
    else
        delta_seq.assign(pack_delta_seq + p_delta->second.second[seq_in_part_id], pack_delta_seq + p_delta->second.second[seq_in_part_id + 1] - 1);

    // LZ decode delta-encoded contig
    ctg.clear();
    lz_diff->Decode(ref_seq, delta_seq, ctg);

    if (need_deallocate_pack_delta_seq)
        delete[] pack_delta_seq;

    if (!fast)
    {
        ref_seq.clear();
        ref_seq.shrink_to_fit();
    }

    return true;
}

// *******************************************************************************************
bool CSegment::get_raw_locked(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
    lock_guard<mutex> lck(mtx);

    return get_raw(id_seq, ctg, zstd_ctx);
}

// *******************************************************************************************
bool CSegment::get_locked(const uint32_t id_seq, contig_t& ctg, ZSTD_DCtx* zstd_ctx)
{
    lock_guard<mutex> lck(mtx);

    return get(id_seq, ctg, zstd_ctx);
}

// *******************************************************************************************
void CSegment::appending_init()
{
    if (internal_state != internal_state_t::none)
        return;

    // Retrive reference contig
    contig_t ref_seq;
    
    int in_stream_id_ref = in_archive->GetStreamId(name + ss_ref_ext(archive_version));
    int in_stream_id_delta = in_archive->GetStreamId(name + ss_delta_ext(archive_version));

    int out_stream_id_ref = -1;
    int out_stream_id_delta = -1;

    if(in_stream_id_ref >= 0)
        out_stream_id_ref = out_archive->RegisterStream(name + ss_ref_ext(archive_version));
    if(in_stream_id_delta >= 0)
        out_stream_id_delta = out_archive->RegisterStream(name + ss_delta_ext(archive_version));

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
    bool empty_ctx = zstd_ctx == nullptr;

    if (!packed_ref_seq.empty())
    {
        contig_t ref_seq;

        if (raw_ref_seq_size == 0)
            ref_seq = move(packed_ref_seq);       // No compression
        else
        {
            vector<uint8_t> v_tuples;

            v_tuples.resize(raw_ref_seq_size + 1);

            if (zstd_ctx == nullptr)
                zstd_ctx = ZSTD_createDCtx();

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
            if (zstd_ctx == nullptr)
                zstd_ctx = ZSTD_createDCtx();

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

    if (empty_ctx && zstd_ctx != nullptr)
        ZSTD_freeDCtx(zstd_ctx);

    internal_state = internal_state_t::normal;
}

// EOF
