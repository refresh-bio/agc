// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.1
// Date   : 2024-03-12
// *******************************************************************************************

#include "../core/archive.h"
#include "defs.h"

#include <iostream>
#include <algorithm>

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#endif

// *******************************************************************************************
CArchive::CArchive(const bool _input_mode, const size_t _io_buffer_size, const string& _lazy_prefix)
{
	input_mode = _input_mode;
	io_buffer_size = _io_buffer_size;

	if(input_mode)			// Ignore lazy_prefix in output mode
		lazy_prefix = _lazy_prefix;
}

// *******************************************************************************************
CArchive::~CArchive()
{
	Close();
}

// *******************************************************************************************
bool CArchive::Open(const string &file_name)
{
	lock_guard<mutex> lck(mtx);

	if (f_in.IsOpened())
		f_in.Close();
	if (f_out.IsOpened())
		f_out.Close();

	if (input_mode)
		f_in.Open(file_name, io_buffer_size);
	else
		f_out.Open(file_name);

	if (!f_in.IsOpened() && !f_out.IsOpened())
		return false;

	if (input_mode)
		deserialize();

	f_offset = 0;

	return true;
}

// *******************************************************************************************
bool CArchive::Close()
{
	lock_guard<mutex> lck(mtx);

	if (!f_in.IsOpened() && !f_out.IsOpened())
		return false;

	if (input_mode)
		f_in.Close();
	else
	{
		flush_out_buffers();
		serialize();
		f_out.Close();
	}

	return true;
}

// *******************************************************************************************
/*size_t CArchive::write_fixed(const uint64_t x)
{
	f_out.WriteUInt(x, 8);

	return 8;
}

// *******************************************************************************************
size_t CArchive::write(const uint64_t _x)
{
	int no_bytes = 0;
	uint64_t x = _x;

	for (size_t tmp = x; tmp; tmp >>= 8)
		++no_bytes;
	
	f_out.Put(no_bytes);

	for (int i = no_bytes; i; --i)
		f_out.Put((x >> ((i - 1) * 8)) & 0xff);

	return no_bytes + 1;
}*/

// *******************************************************************************************
size_t CArchive::write(const string &s)
{
	f_out.Write(s);
	f_out.Put(0);

	return s.size() + 1;
}

// *******************************************************************************************
size_t CArchive::read(string& s)
{
	s.clear();

	while (true)
	{
		int c = f_in.Get();
		if (c == EOF)
			return 0;

		if (c == 0)
			return s.size() + 1;

		s.push_back((char)c);
	}

	return 0;
}

// *******************************************************************************************
bool CArchive::serialize()
{
	size_t footer_size = 0;

	// Store stram part offsets
	footer_size += write(v_streams.size());

	for (auto& stream : v_streams)
	{
		size_t p = footer_size;

		footer_size += write(stream.stream_name);
		footer_size += write(stream.parts.size());
		footer_size += write(stream.raw_size);

		for (auto& part : stream.parts)
		{
			footer_size += write(part.offset);
			footer_size += write(part.size);
		}

		stream.packed_size += footer_size - p;
	}

	write_fixed(footer_size);

	return true;
}

// *******************************************************************************************
bool CArchive::deserialize()
{
	size_t footer_size;
	size_t file_size = f_in.FileSize();

	f_in.Seek(file_size - 8ull);
	read_fixed(footer_size);

	f_in.Seek(file_size -(size_t)(8 + footer_size));

	// Read stream part offsets
	size_t n_streams;
	read(n_streams);

	v_streams.resize(n_streams, stream_t());

	rm_streams.reserve(2 * n_streams);

	for (size_t i = 0; i < n_streams; ++i)
	{
		auto& stream_second = v_streams[i];

		read(stream_second.stream_name);
		read(stream_second.cur_id);
		read(stream_second.raw_size);

		stream_second.parts.resize(stream_second.cur_id);
		for (size_t j = 0; j < stream_second.cur_id; ++j)
		{
			read(stream_second.parts[j].offset);
			read(stream_second.parts[j].size);
		}

		stream_second.cur_id = 0;

		if(!is_lazy_str(stream_second.stream_name))
			rm_streams[stream_second.stream_name] = i;
	}
	
	f_in.Seek(0);

	return true;
}

// *******************************************************************************************
int CArchive::RegisterStream(const string &stream_name)
{
	lock_guard<mutex> lck(mtx);

	// Before adding new stream check if stream_name is already registered
	auto p = rm_streams.find(stream_name);
	if (p != rm_streams.end())
		return (int)p->second;

	int id = (int) v_streams.size();

	v_streams.emplace_back(stream_t());

	v_streams[id].cur_id = 0;
	v_streams[id].stream_name = stream_name;
	v_streams[id].raw_size = 0;
	v_streams[id].packed_size = 0;
	v_streams[id].packed_data_size = 0;

	rm_streams[stream_name] = id;

	return id;
}

// *******************************************************************************************
int CArchive::get_stream_id(const string& stream_name)
{
	if (is_lazy_str(stream_name))
		de_lazy();

	auto p = rm_streams.find(stream_name);
	if (p != rm_streams.end())
		return (int)p->second;

	return -1;
}

// *******************************************************************************************
int CArchive::GetStreamId(const string &stream_name)
{
	lock_guard<mutex> lck(mtx);

	return get_stream_id(stream_name);
}

// *******************************************************************************************
bool CArchive::add_part(const int stream_id, const vector<uint8_t>& v_data, const uint64_t metadata)
{
	v_streams[stream_id].parts.push_back(part_t(f_offset, v_data.size()));

	f_offset += write(metadata);
	f_out.Write(v_data.data(), v_data.size());

	f_offset += v_data.size();

	v_streams[stream_id].packed_size += f_offset - v_streams[stream_id].parts.back().offset;
	v_streams[stream_id].packed_data_size += v_data.size();

	return true;
}

// *******************************************************************************************
bool CArchive::AddPart(const int stream_id, const vector<uint8_t> &v_data, const uint64_t metadata)
{
	lock_guard<mutex> lck(mtx);
	
	return add_part(stream_id, v_data, metadata);
}

// *******************************************************************************************
int CArchive::AddPartPrepare(const int stream_id)
{
	lock_guard<mutex> lck(mtx);
	
	v_streams[stream_id].parts.push_back(part_t(0, 0));

	return static_cast<int>(v_streams[stream_id].parts.size()) - 1;
}

// *******************************************************************************************
bool CArchive::AddPartComplete(const int stream_id, const int part_id, const vector<uint8_t>& v_data, const uint64_t metadata)
{
	lock_guard<mutex> lck(mtx);
	
	v_streams[stream_id].parts[part_id] = part_t(f_offset, v_data.size());

	f_offset += write(metadata);
	f_out.Write(v_data.data(), v_data.size());

	f_offset += v_data.size();

	v_streams[stream_id].packed_size += f_offset - v_streams[stream_id].parts[part_id].offset;
	v_streams[stream_id].packed_data_size += v_data.size();

	return true;
}

// *******************************************************************************************
bool CArchive::AddPartBuffered(const int stream_id, const vector<uint8_t>& v_data, const uint64_t metadata)
{
	lock_guard<mutex> lck(mtx);

	m_buffer[stream_id].emplace_back(v_data, metadata);

	return true;
}

// *******************************************************************************************
bool CArchive::flush_out_buffers()
{
	for (auto& x : m_buffer)
		for (auto& y : x.second)
			add_part(x.first, y.first, y.second);

	m_buffer.clear();

	return true;
}

// *******************************************************************************************
bool CArchive::FlushOutBuffers()
{
	lock_guard<mutex> lck(mtx);

	return flush_out_buffers();
}

// *******************************************************************************************
void CArchive::SetRawSize(const int stream_id, const size_t raw_size)
{
	lock_guard<mutex> lck(mtx);
	
	v_streams[stream_id].raw_size = raw_size;
}

// *******************************************************************************************
size_t CArchive::GetRawSize(const int stream_id)
{
	lock_guard<mutex> lck(mtx);
	
	return v_streams[stream_id].raw_size;
}

// *******************************************************************************************
bool CArchive::get_part(const int stream_id, vector<uint8_t>& v_data, uint64_t& metadata)
{
	auto& p = v_streams[stream_id];

	if (p.cur_id >= p.parts.size())
		return false;

	v_data.resize(p.parts[p.cur_id].size);

	f_in.Seek(p.parts[p.cur_id].offset);

	if (p.parts[p.cur_id].size != 0)
		read(metadata);
	else
	{
		metadata = 0;
		p.cur_id++;
		return true;
	}

	f_in.Read(v_data.data(), p.parts[p.cur_id].size);

	p.cur_id++;

	return true;
}

// *******************************************************************************************
bool CArchive::GetPart(const int stream_id, vector<uint8_t> &v_data, uint64_t &metadata)
{
	lock_guard<mutex> lck(mtx);
	
	return get_part(stream_id, v_data, metadata);
}

// *******************************************************************************************
pair<int, bool> CArchive::GetPart(const string& stream_name, vector<uint8_t> &v_data, uint64_t &metadata)
{
	lock_guard<mutex> lck(mtx);
	
	int stream_id = get_stream_id(stream_name);

	if (stream_id < 0)
		return make_pair(-1, false);

	return make_pair(stream_id, get_part(stream_id, v_data, metadata));
}

// *******************************************************************************************
tuple<int, bool, int, bool> CArchive::GetParts(
	const string& stream_name1, vector<uint8_t>& v_data1, uint64_t& metadata1,
	const string& stream_name2, vector<uint8_t>& v_data2, uint64_t& metadata2)
{
	lock_guard<mutex> lck(mtx);

	bool res1 = false;
	bool res2 = false;

	int stream_id1 = get_stream_id(stream_name1);
	int stream_id2 = get_stream_id(stream_name2);

	if (stream_id1 >= 0)
		res1 = get_part(stream_id1, v_data1, metadata1);

	if (stream_id2 >= 0)
		res2 = get_part(stream_id2, v_data2, metadata2);

	return make_tuple(stream_id1, res1, stream_id2, res2);
}

// *******************************************************************************************
bool CArchive::get_part(const int stream_id, const int part_id, vector<uint8_t>& v_data, uint64_t& metadata)
{
	auto& p = v_streams[stream_id];

	if ((size_t)part_id >= p.parts.size())
		return false;

	v_data.resize(p.parts[part_id].size);

	f_in.Seek(p.parts[part_id].offset);

	if (p.parts[part_id].size != 0)
		read(metadata);
	else
	{
		metadata = 0;
		return true;
	}

	f_in.Read(v_data.data(), p.parts[part_id].size);

	return true;
}

// *******************************************************************************************
bool CArchive::GetPart(const int stream_id, const int part_id, vector<uint8_t> &v_data, uint64_t &metadata)
{
	lock_guard<mutex> lck(mtx);

	return get_part(stream_id, part_id, v_data, metadata);
}

// *******************************************************************************************
pair<int, bool> CArchive::GetPart(const string& stream_name, const int part_id, vector<uint8_t> &v_data, uint64_t &metadata)
{
	lock_guard<mutex> lck(mtx);

	int stream_id = get_stream_id(stream_name);

	if (stream_id < 0)
		return make_pair(-1, false);

	return make_pair(stream_id, get_part(stream_id, part_id, v_data, metadata));
}

// *******************************************************************************************
tuple<int, bool, int, bool> CArchive::GetParts(
	const string& stream_name1, const int part_id1, vector<uint8_t>& v_data1, uint64_t& metadata1,
	const string& stream_name2, const int part_id2, vector<uint8_t>& v_data2, uint64_t& metadata2)
{
	lock_guard<mutex> lck(mtx);

	bool res1 = false;
	bool res2 = false;

	int stream_id1 = get_stream_id(stream_name1);
	int stream_id2 = get_stream_id(stream_name2);

	if (stream_id1 >= 0)
		res1 = get_part(stream_id1, part_id1, v_data1, metadata1);

	if (stream_id2 >= 0)
		res2 = get_part(stream_id2, part_id2, v_data2, metadata2);

	return make_tuple(stream_id1, res1, stream_id2, res2);
}

// *******************************************************************************************
size_t CArchive::GetNoStreams()
{
	lock_guard<mutex> lck(mtx);

	return v_streams.size();
}

// *******************************************************************************************
size_t CArchive::GetNoParts(const int stream_id)
{
	lock_guard<mutex> lck(mtx);

	if (stream_id < 0 || (size_t)stream_id >= v_streams.size())
		return 0;

	return v_streams[stream_id].parts.size();
}

// *******************************************************************************************
size_t CArchive::GetStreamPackedSize(const int stream_id)
{
	lock_guard<mutex> lck(mtx);

	if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
		return 0;

	return v_streams[stream_id].packed_size;
}

// *******************************************************************************************
size_t CArchive::GetStreamPackedDataSize(const int stream_id)
{
	lock_guard<mutex> lck(mtx);

	if (stream_id < 0 || stream_id >= static_cast<int>(v_streams.size()))
		return 0;

	return v_streams[stream_id].packed_data_size;
}

// EOF
