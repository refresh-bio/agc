// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include "../core/collection_v3.h"
#include <cassert>
#include <future>

// *******************************************************************************************
bool CCollection_V3::set_archives(shared_ptr<CArchive> _in_archive, shared_ptr<CArchive> _out_archive,
	uint32_t _no_threads, size_t _batch_size, uint32_t _segment_size, uint32_t _kmer_length)
{
	lock_guard<mutex> lck(mtx);

	in_archive = _in_archive;
	out_archive = _out_archive;

	batch_size = _batch_size;
	segment_size = _segment_size;
	no_threads = _no_threads;
	kmer_length = _kmer_length;

	if (in_archive == nullptr)
		return prepare_for_compression();
	else if (out_archive == nullptr)
		return prepare_for_decompression();
	else
		return prepare_for_appending_copy();
}

// *******************************************************************************************
bool CCollection_V3::prepare_for_compression()
{
	collection_samples_id = out_archive->RegisterStream("collection-samples");
	collection_contig_id = out_archive->RegisterStream("collection-contigs");
	collection_details_id = out_archive->RegisterStream("collection-details");	

	return true;
}

// *******************************************************************************************
bool CCollection_V3::prepare_for_appending_copy()
{
//	auto in_collection_samples_id = in_archive->GetStreamId("collection-samples");
	auto in_collection_contig_id = in_archive->GetStreamId("collection-contigs");
	auto in_collection_details_id = in_archive->GetStreamId("collection-details");

	collection_samples_id = out_archive->RegisterStream("collection-samples");
	collection_contig_id = out_archive->RegisterStream("collection-contigs");
	collection_details_id = out_archive->RegisterStream("collection-details");

	load_batch_sample_names();

	// in and out ids for collection-* must be the same!

	auto no_contig_batches = in_archive->GetNoParts(in_collection_contig_id);

	// Transfer all but the last one batch from in to out archive
	vector<uint8_t> data;
	uint64_t meta;

	for (size_t i = 0; i < no_contig_batches - 1; ++i)
	{
		in_archive->GetPart(in_collection_contig_id, i, data, meta);
		out_archive->AddPart(collection_contig_id, data, meta);
		in_archive->GetPart(in_collection_details_id, i, data, meta);
		out_archive->AddPart(collection_details_id, data, meta);
	}

	return true;
}

// *******************************************************************************************
bool CCollection_V3::prepare_for_appending_load_last_batch()
{
	lock_guard<mutex> lck(mtx);
	
	auto in_collection_contig_id = in_archive->GetStreamId("collection-contigs");
	auto in_collection_details_id = in_archive->GetStreamId("collection-details");

	auto no_contig_batches = in_archive->GetNoParts(in_collection_contig_id);

	// Transfer all but the last one batch from in to out archive
	vector<uint8_t> data;
	uint64_t meta;

	// Load last batch
	load_batch_contig_names(no_contig_batches - 1);
	load_batch_contig_details(no_contig_batches - 1);

	if (no_samples_in_last_batch == batch_size)
	{
		in_archive->GetPart(in_collection_contig_id, no_contig_batches - 1, data, meta);
		out_archive->AddPart(collection_contig_id, data, meta);
		in_archive->GetPart(in_collection_details_id, no_contig_batches - 1, data, meta);
		out_archive->AddPart(collection_details_id, data, meta);

		clear_batch_contig(no_contig_batches - 1);
	}

	return true;
}

// *******************************************************************************************
bool CCollection_V3::prepare_for_decompression()
{
	collection_samples_id = in_archive->GetStreamId("collection-samples");
	collection_contig_id = in_archive->GetStreamId("collection-contigs");
	collection_details_id = in_archive->GetStreamId("collection-details");

	load_batch_sample_names();

	return true;
}

// *******************************************************************************************
void CCollection_V3::complete_serialization()
{
	lock_guard<mutex> lck(mtx);

	store_batch_sample_names();
}

// *******************************************************************************************
void CCollection_V3::zstd_compress(ZSTD_CCtx*& cctx, vector<uint8_t>& v_input, vector<uint8_t>& v_output, int level)
{
	if (cctx == nullptr)
	{
		cctx = ZSTD_createCCtx();
	}

	v_output.resize(ZSTD_compressBound(v_input.size()));
	auto c_size = ZSTD_compressCCtx(cctx, v_output.data(), v_output.size(), v_input.data(), v_input.size(), level);

	v_output.resize(c_size);
}

// *******************************************************************************************
void CCollection_V3::zstd_decompress(ZSTD_DCtx*& dctx, vector<uint8_t>& v_input, vector<uint8_t>& v_output, size_t raw_size)
{
	if (dctx == nullptr)
		dctx = ZSTD_createDCtx();

	v_output.resize(raw_size);
	ZSTD_decompressDCtx(dctx, v_output.data(), v_output.size(), v_input.data(), v_input.size());
}

// *******************************************************************************************
void CCollection_V3::store_batch_sample_names()
{
	vector<uint8_t> v_data, v_tmp;

	determine_collection_samples_id();

	serialize_sample_names(v_tmp);

	zstd_compress(zstd_cctx_samples, v_tmp, v_data, 19);

	out_archive->AddPartBuffered(collection_samples_id, v_data, v_tmp.size());
}

// *******************************************************************************************
void CCollection_V3::load_batch_sample_names()
{
	vector<uint8_t> v_data, v_tmp;
	uint64_t raw_size;

	determine_collection_samples_id();

	in_archive->GetPart(collection_samples_id, v_tmp, raw_size);

	zstd_decompress(zstd_dctx_samples, v_tmp, v_data, raw_size);

	deserialize_sample_names(v_data);
}

// *******************************************************************************************
void CCollection_V3::store_batch_contig_names(uint32_t id_from, uint32_t id_to)
{
	vector<uint8_t> v_data, v_tmp;

	determine_collection_contig_id();

	serialize_contig_names(v_tmp, id_from, id_to);

	zstd_compress(zstd_cctx_contigs, v_tmp, v_data, 18);

	out_archive->AddPartBuffered(collection_contig_id, v_data, v_tmp.size());
}

// *******************************************************************************************
void CCollection_V3::load_batch_contig_names(size_t id_batch)
{
	vector<uint8_t> v_data, v_tmp;
	uint64_t raw_size;

	if (unpacked_contig_data_batch_id >= 0 && unpacked_contig_data_batch_id != (int) id_batch)
		clear_batch_contig(unpacked_contig_data_batch_id);

	determine_collection_contig_id();

	in_archive->GetPart(collection_contig_id, id_batch, v_tmp, raw_size);

	zstd_decompress(zstd_dctx_contigs, v_tmp, v_data, raw_size);

	deserialize_contig_names(v_data, id_batch * batch_size);

	unpacked_contig_data_batch_id = id_batch;
}

// *******************************************************************************************
void CCollection_V3::clear_batch_contig(size_t id_batch)
{
	size_t to_batch_id = min(sample_desc.size(), (id_batch + 1) * batch_size);

	for (size_t i = id_batch * batch_size; i < to_batch_id; ++i)
	{
		sample_desc[i].contigs.clear();
		sample_desc[i].contigs.shrink_to_fit();
	}
}

// *******************************************************************************************
void CCollection_V3::store_batch_contig_details(uint32_t id_from, uint32_t id_to)
{
	array<vector<uint8_t>, 5> v_data;
	array<vector<uint8_t>, 5> v_packed;
	vector<uint8_t> v_stream;

	determnine_collection_details_id();

	serialize_contig_details(v_data, id_from, id_to);

	if (no_threads >= 4)
	{
		vector<future<void>> v_fut;
		v_fut.reserve(5);

		for (int i = 0; i < 5; ++i)
			v_fut.emplace_back(async([&, i]() {zstd_compress(zstd_cctx_details[i], v_data[i], v_packed[i], 19); }));

		for (int i = 0; i < 5; ++i)
			v_fut[i].wait();
	}
	else
	{
		for (int i = 0; i < 5; ++i)
			zstd_compress(zstd_cctx_details[i], v_data[i], v_packed[i], 19);
	}

	for (int i = 0; i < 5; ++i)
	{
		append(v_stream, v_data[i].size());
		append(v_stream, v_packed[i].size());
	}

	for (int i = 0; i < 5; ++i)
		v_stream.insert(v_stream.end(), v_packed[i].begin(), v_packed[i].end());

	out_archive->AddPartBuffered(collection_details_id, v_stream, 0);
}

// *******************************************************************************************
void CCollection_V3::load_batch_contig_details(size_t id_batch)
{
	array<vector<uint8_t>, 5> v_data;
	array<vector<uint8_t>, 5> v_packed;
	vector<uint8_t> v_stream;

	uint64_t aux;

	if (unpacked_contig_data_batch_id >= 0 && unpacked_contig_data_batch_id != (int) id_batch)
		clear_batch_contig(unpacked_contig_data_batch_id);

	determnine_collection_details_id();

	in_archive->GetPart(collection_details_id, id_batch, v_stream, aux);

	uint8_t* ptr = v_stream.data();
	array<pair<uint32_t, uint32_t>, 5> a_sizes;

	for (int i = 0; i < 5; ++i)
	{
		read(ptr, a_sizes[i].first);
		read(ptr, a_sizes[i].second);
	}

	for (int i = 0; i < 5; ++i)
	{
		v_packed[i].assign(ptr, ptr + a_sizes[i].second);
		ptr += a_sizes[i].second;
	}

	if (no_threads >= 4)
	{
		vector<future<void>> v_fut;
		v_fut.reserve(5);

		for (int i = 0; i < 5; ++i)
			v_fut.emplace_back(async([&, i]() {zstd_decompress(zstd_dctx_details[i], v_packed[i], v_data[i], a_sizes[i].first); }));

		for (int i = 0; i < 5; ++i)
			v_fut[i].wait();
	}
	else
	{
		for (int i = 0; i < 5; ++i)
			zstd_decompress(zstd_dctx_details[i], v_packed[i], v_data[i], a_sizes[i].first);
	}

	deserialize_contig_details(v_data, id_batch * batch_size);

	unpacked_contig_data_batch_id = id_batch;
}

// *******************************************************************************************
void CCollection_V3::serialize_sample_names(vector<uint8_t>& v_data)
{
	append(v_data, (uint32_t) sample_desc.size());

	for (auto& x: sample_desc)
		append(v_data, x.name);
}

// *******************************************************************************************
void CCollection_V3::deserialize_sample_names(vector<uint8_t>& v_data)
{
	uint8_t* p = v_data.data();

	uint32_t no_samples;

	read(p, no_samples);

	sample_desc.resize(no_samples);

	for (size_t i = 0; i < no_samples; ++i)
	{
		read(p, sample_desc[i].name);
		sample_ids[sample_desc[i].name] = i;
	}
}

// *******************************************************************************************
vector<string> CCollection_V3::split_string(const string& s)
{
	auto p = s.begin();
	vector<string> components;

	while (true) 
	{
		auto q = find(p, s.end(), ' ');
		components.push_back(string(p, q));
		if (q == s.end())
			break;

		p = q + 1;
	}

	return components;
}

// *******************************************************************************************
string CCollection_V3::encode_split(vector<string>& prev_split, vector<string>& curr_split)
{
	string enc;

	for (size_t i = 0; i < curr_split.size(); ++i)
	{
		if (prev_split[i] == curr_split[i])
			enc.push_back(-127);						// same component marker
		else if (prev_split[i].size() != curr_split[i].size())
			enc.append(curr_split[i]);			
		else
		{
			char cnt = 0;
			uint32_t cmp_len = curr_split[i].size();

			auto p_ptr = prev_split[i].data();
			auto c_ptr = curr_split[i].data();

			for (uint32_t j = 0; j < cmp_len; ++j)
			{
				if (p_ptr[j] == c_ptr[j])
				{
					if (cnt == 100)
					{
						enc.push_back(-cnt);			// repetition marker
						cnt = 1;
					}
					else
						++cnt;
				}
				else
				{
					if (cnt)
					{
						enc.push_back(-cnt);			// repetition marker
						cnt = 0;
					}

					enc.push_back(c_ptr[j]);
				}
			}

			if (cnt)
				enc.push_back(-cnt);			// repetition marker
		}

		enc.push_back(' ');
	}

	enc.pop_back();			// remove final space

	return enc;
}

// *******************************************************************************************
string CCollection_V3::decode_split(vector<string>& prev_split, vector<string>& curr_split)
{
	string dec;
	string cmp;

	for (size_t i = 0; i < curr_split.size(); ++i)
	{
		if (curr_split[i].size() == 1 && curr_split[i].front() == -127)		// same component marker
		{
			dec.append(prev_split[i]);
			curr_split[i] = prev_split[i];
		}
		else
		{
//			uint32_t prev_idx = 0;
			cmp.clear();
			auto p_ptr = prev_split[i].data();

			for (auto c : curr_split[i])
			{
				if (c >= 0)
				{
					cmp.push_back(c);
					++p_ptr;
				}
				else
				{
					cmp.append(p_ptr, -c);
					p_ptr += -c;
				}
			}

			dec.append(cmp);
			curr_split[i] = move(cmp);
		}

		dec.push_back(' ');
	}

	dec.pop_back();			// remove final space

	return dec;
}

// *******************************************************************************************
void CCollection_V3::serialize_contig_names(vector<uint8_t>& v_data, uint32_t id_from, uint32_t id_to)
{
	append(v_data, id_to - id_from);

	string p_name;

	vector<string> sp_prev, sp_curr;

	for (auto p = sample_desc.begin() + id_from; p != sample_desc.begin() + id_to; ++p)
	{
		append(v_data, p->contigs.size());

		vector<string> prev_split;
		vector<string> curr_split;

		for (auto& x : p->contigs)
		{
			curr_split = split_string(x.name);

			if (curr_split.size() != prev_split.size())
				append(v_data, x.name);
			else
				append(v_data, encode_split(prev_split, curr_split));

			prev_split = move(curr_split);
		}
	}
}

// *******************************************************************************************
void CCollection_V3::deserialize_contig_names(vector<uint8_t>& v_data, size_t i_sample)
{
	uint8_t* p = v_data.data();

	uint32_t no_samples_in_curr_batch;
	uint32_t no_contigs_in_curr_sample;

	read(p, no_samples_in_curr_batch);

	for (size_t i = 0; i < no_samples_in_curr_batch; ++i)
	{
		read(p, no_contigs_in_curr_sample);

		auto& curr_sample = sample_desc[i_sample + i];

		curr_sample.contigs.resize(no_contigs_in_curr_sample);

		vector<string> prev_split;
		vector<string> curr_split;
		string enc;

		for (size_t j = 0; j < no_contigs_in_curr_sample; ++j)
		{
			read(p, enc);

			curr_split = split_string(enc);

			if (curr_split.size() != prev_split.size())
				curr_sample.contigs[j].name = enc;
			else
				curr_sample.contigs[j].name = decode_split(prev_split, curr_split);

			prev_split = move(curr_split);
		}
	}

	// important only for appending mode
	no_samples_in_last_batch = no_samples_in_curr_batch;	
}

// *******************************************************************************************
void CCollection_V3::serialize_contig_details(array<vector<uint8_t>, 5>& v_data, uint32_t id_from, uint32_t id_to)
{
	append(v_data[0], id_to - id_from);

	clear_in_group_ids();

	for (auto p = sample_desc.begin() + id_from; p != sample_desc.begin() + id_to; ++p)
	{
		append(v_data[0], p->contigs.size());

		uint32_t pred_raw_length = segment_size + kmer_length;

		for (auto& x : p->contigs)
		{
			append(v_data[0], x.segments.size());

			for (auto& seg : x.segments)
			{
				int prev_in_group_id = get_in_group_id(seg.group_id);

				uint32_t e_group_id = seg.group_id;
				uint32_t e_in_group_id;

				if (prev_in_group_id == -1)
					e_in_group_id = (int) seg.in_group_id;
				else
				{
					if (seg.in_group_id == 0)
						e_in_group_id = 0;
					else if ((int)seg.in_group_id == prev_in_group_id + 1)
						e_in_group_id = 1;
					else
						e_in_group_id = (uint32_t)zigzag_encode(seg.in_group_id, prev_in_group_id + 1) + 1u;
				}
				
				uint32_t e_raw_length = (uint32_t)zigzag_encode(seg.raw_length, pred_raw_length);

				append(v_data[1], e_group_id);
				append(v_data[2], e_in_group_id);
				append(v_data[3], e_raw_length);
				append(v_data[4], (uint32_t)seg.is_rev_comp);

				if ((int) seg.in_group_id > prev_in_group_id && seg.in_group_id > 0)
					set_in_group_id(seg.group_id, seg.in_group_id);
			}
		}
	}
}

// *******************************************************************************************
void CCollection_V3::deserialize_contig_details(array<vector<uint8_t>, 5>& v_data, size_t i_sample)
{
	array<vector<uint32_t>, 5> v_det;
	
	uint8_t* p = v_data[0].data();

	uint32_t no_samples_in_curr_batch;
	uint32_t no_contigs_in_curr_sample;
	uint32_t no_segments_in_curr_contig;
	size_t no_items = 0;

	read(p, no_samples_in_curr_batch);

	for (size_t i = 0; i < no_samples_in_curr_batch; ++i)
	{
		read(p, no_contigs_in_curr_sample);

		auto& curr_sample = sample_desc[i_sample + i];

		curr_sample.contigs.resize(no_contigs_in_curr_sample);

		for (size_t j = 0; j < no_contigs_in_curr_sample; ++j)
		{
			read(p, no_segments_in_curr_contig);

			curr_sample.contigs[j].segments.resize(no_segments_in_curr_contig);

			no_items += no_segments_in_curr_contig;
		}
	}

	for(int i = 1; i < 5; ++i)
	{
		v_det[i].resize(no_items);
		
		p = v_data[i].data();

		for (size_t j = 0; j < no_items; ++j)
			read(p, v_det[i][j]);
	}

	no_items = 0;

	clear_in_group_ids();

	uint32_t pred_raw_length = segment_size + kmer_length;

	for (size_t i = 0; i < no_samples_in_curr_batch; ++i)
	{
		auto& curr_sample = sample_desc[i_sample + i];

		for (size_t j = 0; j < curr_sample.contigs.size(); ++j)
		{
			auto& curr_contig = curr_sample.contigs[j];

			for (size_t k = 0; k < curr_contig.segments.size(); ++k, ++no_items)
			{
				uint32_t c_group_id = v_det[1][no_items];

				curr_contig.segments[k].group_id = c_group_id;
				int prev_in_group_id = get_in_group_id(c_group_id);

				uint32_t e_in_group_id = v_det[2][no_items];
				uint32_t c_in_group_id;

				if (prev_in_group_id == -1)
					c_in_group_id = e_in_group_id;
				else
				{
					if (e_in_group_id == 0)
						c_in_group_id = 0;
					else if (e_in_group_id == 1)
						c_in_group_id = prev_in_group_id + 1;
					else
						c_in_group_id = (uint32_t)zigzag_decode(e_in_group_id - 1u, prev_in_group_id + 1);

				}

				curr_contig.segments[k].in_group_id = c_in_group_id;

				uint32_t c_raw_length = (uint32_t)zigzag_decode(v_det[3][no_items], pred_raw_length);
				curr_contig.segments[k].raw_length = c_raw_length;

				curr_contig.segments[k].is_rev_comp = (bool)v_det[4][no_items];

				if ((int)c_in_group_id > prev_in_group_id && c_in_group_id > 0)
					set_in_group_id(c_group_id, c_in_group_id);
			}
		}
	}
}

// *******************************************************************************************
void CCollection_V3::store_contig_batch(uint32_t id_from, uint32_t id_to)
{
	lock_guard<mutex> lck(mtx);

	if (no_threads > 1)
	{
		future<void> fut_contigs = async([&]() {this->store_batch_contig_names(id_from, id_to); });
		store_batch_contig_details(id_from, id_to);
		fut_contigs.wait();
	}
	else
	{
		store_batch_contig_names(id_from, id_to);
		store_batch_contig_details(id_from, id_to);
	}

	for (auto p = sample_desc.begin() + id_from; p != sample_desc.begin() + id_to; ++p)
	{
		p->contigs.clear();
		p->contigs.shrink_to_fit();
	}
}

// *******************************************************************************************
bool CCollection_V3::register_sample_contig(const string& sample_name, const string& contig_name)
{
	string short_contig_name = extract_contig_name(contig_name);
	string stored_sample_name = sample_name;

	lock_guard<mutex> lck(mtx);

	if (sample_name.empty())
		stored_sample_name = short_contig_name;

	if (stored_sample_name != prev_sample_name)
	{
		auto q = sample_ids.find(stored_sample_name);
		if (q != sample_ids.end())
			return false;		// sample of the same name was already registered (prior to previous sample_name)

		uint32_t sample_id = (uint32_t)sample_ids.size();
		sample_ids[stored_sample_name] = sample_id;
		sample_desc.emplace_back(stored_sample_name);

		prev_sample_name = stored_sample_name;
	}

	sample_desc.back().contigs.emplace_back(contig_desc_t(contig_name));

	return true;
}

// *******************************************************************************************
void CCollection_V3::reset_prev_sample_name()
{
	lock_guard<mutex> lck(mtx);

	prev_sample_name.clear();
}

// *******************************************************************************************
void CCollection_V3::add_segment_placed(const string& sample_name, const string& contig_name, const uint32_t place, const uint32_t group_id, const uint32_t in_group_id, const bool is_rev_comp, const uint32_t raw_length)
{
	lock_guard<mutex> lck(mtx);

	string stored_sample_name = sample_name;

	if (sample_name.empty())
		stored_sample_name = extract_contig_name(contig_name);

	if (placing_sample_name != stored_sample_name)
	{
		placing_sample_name = stored_sample_name;
		placing_sample_id = sample_ids.find(stored_sample_name)->second;
	}

	for (auto& x : sample_desc[placing_sample_id].contigs)
	{
		if (x.name == contig_name)
		{
			if (place >= x.segments.size())
				x.segments.resize(place + 1);

			x.segments[place] = segment_desc_t(group_id, in_group_id, is_rev_comp, raw_length);

			return;
		}
	}
}

// *******************************************************************************************
void CCollection_V3::add_segments_placed(vector<segments_to_place_t>& segments_to_place)
{
	lock_guard<mutex> lck(mtx);

	for (const auto& desc : segments_to_place)
	{
		auto p = sample_ids.find(desc.sample_name);

		if (p == sample_ids.end())
		{
			assert("Wrong sample name\n");
			return;
		}

		for (auto& x : sample_desc[p->second].contigs)
		{
			if (x.name == desc.contig_name)
			{
				if (desc.seg_part_no >= x.segments.size())
					x.segments.resize(desc.seg_part_no + 1);

				x.segments[desc.seg_part_no] = segment_desc_t(desc.group_id, desc.in_group_id, desc.is_rev_comp, desc.data_size);

				break;
			}
		}
	}
}
// *******************************************************************************************
bool CCollection_V3::get_reference_name(string& reference_name)
{
	lock_guard<mutex> lck(mtx);

	if (sample_desc.empty())
		return false;

	reference_name = sample_desc.front().name;

	return true;
}

// *******************************************************************************************
bool CCollection_V3::get_samples_list(vector<string>& v_samples)
{
	lock_guard<mutex> lck(mtx);

	v_samples.clear();
	v_samples.reserve(sample_desc.size());

	for (auto& x : sample_desc)
		v_samples.emplace_back(x.name);

	sort(v_samples.begin(), v_samples.end());

	return true;
}

// *******************************************************************************************
bool CCollection_V3::get_contig_list_in_sample(const string& sample_name, vector<string>& v_contig_names)
{
	lock_guard<mutex> lck(mtx);

	auto p = sample_ids.find(sample_name);

	if (p == sample_ids.end())
		return false;		// Error: no such a sample

	if (sample_desc[p->second].contigs.empty())
		load_batch_contig_names(p->second / batch_size);

	v_contig_names.clear();
	v_contig_names.reserve(sample_desc[p->second].contigs.size());

	for (auto& x : sample_desc[p->second].contigs)
		v_contig_names.emplace_back(x.name);

	return true;
}

// *******************************************************************************************
bool CCollection_V3::get_sample_desc(const string& sample_name, vector<pair<string, vector<segment_desc_t>>>& sample_desc_)
{
	lock_guard<mutex> lck(mtx);

	sample_desc_.clear();

	auto p = sample_ids.find(sample_name);

	if (p == sample_ids.end())
		return false;		// Error: no such a sample

	if (sample_desc[p->second].contigs.empty())
	{
		load_batch_contig_names(p->second / batch_size);

		load_batch_contig_details(p->second / batch_size);
	}

	sample_desc_.reserve(sample_desc[p->second].contigs.size());

	for (auto& x : sample_desc[p->second].contigs)
		sample_desc_.emplace_back(x.name, x.segments);

	return true;
}

// *******************************************************************************************
bool CCollection_V3::get_contig_desc(const string& sample_name, string& contig_name, vector<segment_desc_t>& contig_desc)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);

	contig_desc.clear();

	auto p = sample_ids.find(sample_name);

	if (p == sample_ids.end())
		return false;		// Error: no such a sample

	if (sample_desc[p->second].contigs.empty())
		load_batch_contig_names(p->second / batch_size);

	if (sample_desc[p->second].contigs.empty() || sample_desc[p->second].contigs.front().segments.empty())
		load_batch_contig_details(p->second / batch_size);
	
	for (auto& x : sample_desc[p->second].contigs)
	{
		if (extract_contig_name(x.name) == short_contig_name)
		{
			contig_desc = x.segments;
			contig_name = x.name;
			return true;
		}
	}

	return false;
}

// *******************************************************************************************
bool CCollection_V3::is_contig_desc(const string& sample_name, const string& contig_name)
{
	lock_guard<mutex> lck(mtx);

	string short_contig_name = extract_contig_name(contig_name);

	auto p = sample_ids.find(sample_name);

	if (p == sample_ids.end())
		return false;		// Error: no such a sample

	if (sample_desc[p->second].contigs.empty())
		load_batch_contig_names(p->second / batch_size);

	for (auto& x : sample_desc[p->second].contigs)
		if (extract_contig_name(x.name) == contig_name)
			return true;

	return false;
}

// *******************************************************************************************
vector<string> CCollection_V3::get_samples_for_contig(const string& contig_name)
{
	lock_guard<mutex> lck(mtx);

	vector<string> v_samples;

	string short_contig_name = extract_contig_name(contig_name);

	size_t no_batches = (sample_desc.size() + batch_size - 1) / batch_size;

	for (size_t i = 0; i < no_batches; ++i)
	{
		if (sample_desc[i * batch_size].contigs.empty())
			load_batch_contig_names(i);

		size_t to_batch_id = min(sample_desc.size(), (i + 1) * batch_size);

		for (size_t j = i * batch_size; j < to_batch_id; ++j)
		{
			for (auto& x : sample_desc[j].contigs)
				if(extract_contig_name(x.name) == short_contig_name)
					v_samples.emplace_back(sample_desc[j].name);
		}

		clear_batch_contig(i);
	}

	return v_samples;
}

// *******************************************************************************************
size_t CCollection_V3::get_no_samples()
{
	lock_guard<mutex> lck(mtx);

	return sample_desc.size();
}

// *******************************************************************************************
int32_t CCollection_V3::get_no_contigs(const string& sample_name)
{
	lock_guard<mutex> lck(mtx);

	auto p = sample_ids.find(sample_name);

	if (p == sample_ids.end())
		return -1;		// Error: no such a sample

	if (sample_desc[p->second].contigs.empty())
		load_batch_contig_names(p->second / batch_size);

	return (int32_t) sample_desc[p->second].contigs.size();
}

// EOF
