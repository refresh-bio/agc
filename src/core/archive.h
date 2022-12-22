#ifndef _ARCHIVE_H
#define _ARCHIVE_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.0
// Date   : 2022-12-22
// *******************************************************************************************

#include <cstdio>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>
#include <thread>
#include <mutex>
#include "../core/io.h"

using namespace std;

class CArchive
{
	bool input_mode;
	CInFile f_in;
	COutFile f_out;
	size_t io_buffer_size;

	size_t f_offset;

	struct part_t{
		size_t offset;
		size_t size;

		part_t() : offset(0), size(0)
		{};

		part_t(size_t _offset, size_t _size) : offset(_offset), size(_size)
		{};
	};

	typedef struct {
		string stream_name;
		size_t cur_id;
		size_t raw_size;
		size_t packed_size;
		size_t packed_data_size;
		vector<part_t> parts;
	} stream_t;

	map<int, vector<pair<vector<uint8_t>, uint64_t>>> m_buffer;

	vector<stream_t> v_streams;
	unordered_map<string, size_t> rm_streams;
	string lazy_prefix;

	mutex mtx;

	bool serialize();
	bool deserialize();

	// *******************************************************************************************
	inline bool is_lazy_str(const string& str)
	{
		if (!input_mode)
			return false;

		if (lazy_prefix.empty())
			return false;

		if (str.size() <= lazy_prefix.size())
			return false;

		const char* p = str.c_str();
		const char* q = lazy_prefix.c_str();

		while (*p == *q)
			++p, ++q;

		return *q == 0;
	}

	// *******************************************************************************************
	void de_lazy()
	{
		rm_streams.reserve(2 * v_streams.size());

		for(size_t i = 0; i < v_streams.size(); ++i)
			if(is_lazy_str(v_streams[i].stream_name))
				rm_streams[v_streams[i].stream_name] = i;

		lazy_prefix.clear();
	}

	// *******************************************************************************************
	template<typename T>
	size_t write_fixed(const T x)
	{
		f_out.WriteUInt(static_cast<uint64_t>(x), 8);

		return 8;
	}

	// *******************************************************************************************
	template<typename T>
	size_t write(const T _x)
	{
		int no_bytes = 0;
		uint64_t x = static_cast<uint64_t>(_x);

		for (size_t tmp = x; tmp; tmp >>= 8)
			++no_bytes;

		f_out.Put(no_bytes);

		for (int i = no_bytes; i; --i)
			f_out.Put((x >> ((i - 1) * 8)) & 0xff);

		return no_bytes + 1;
	}

	// *******************************************************************************************
	size_t write(const string &s);
	
	// *******************************************************************************************
	template<typename T>
	size_t read_fixed(T& x)
	{
		x = static_cast<T>(f_in.ReadUInt(8));

		return 8;
	}

	// *******************************************************************************************
	size_t read(string& s);

	// *******************************************************************************************
	template<typename T>
	size_t read(T& x)
	{
		int no_bytes = f_in.Get();

		x = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			x <<= 8;
			x += static_cast<T>(f_in.Get());
		}

		return no_bytes + 1;
	}

	// *******************************************************************************************
	bool add_part(const int stream_id, const vector<uint8_t>& v_data, const uint64_t metadata);
	bool flush_out_buffers();
	int get_stream_id(const string& stream_name);
	bool get_part(const int stream_id, vector<uint8_t>& v_data, uint64_t& metadata);
	bool get_part(const int stream_id, const int part_id, vector<uint8_t>& v_data, uint64_t& metadata);

public:
	CArchive(const bool _input_mode, const size_t _io_buffer_size = 64 << 20, const string& _lazy_prefix = "");
	~CArchive();

	bool Open(const string &file_name);
	bool Close();

	int RegisterStream(const string &stream_name);
	int GetStreamId(const string &stream_name);

	size_t GetStreamPackedSize(const int stream_id);
	size_t GetStreamPackedDataSize(const int stream_id);

	bool AddPart(const int stream_id, const vector<uint8_t>& v_data, const uint64_t metadata = 0);
	int AddPartPrepare(const int stream_id);
	bool AddPartComplete(const int stream_id, const int part_id, const vector<uint8_t>& v_data, const uint64_t metadata = 0);
	bool AddPartBuffered(const int stream_id, const vector<uint8_t>& v_data, const uint64_t metadata = 0);

	bool FlushOutBuffers();

	bool GetPart(const int stream_id, vector<uint8_t> &v_data, uint64_t &metadata);
	bool GetPart(const int stream_id, const int part_id, vector<uint8_t>& v_data, uint64_t& metadata);

	pair<int, bool> GetPart(const string &stream_name, vector<uint8_t> &v_data, uint64_t &metadata);
	pair<int, bool> GetPart(const string& stream_name, const int part_id, vector<uint8_t> &v_data, uint64_t &metadata);

	tuple<int, bool, int, bool> GetParts(
		const string &stream_name1, vector<uint8_t> &v_data1, uint64_t &metadata1,
		const string& stream_name2, vector<uint8_t>& v_data2, uint64_t& metadata2);
	tuple<int, bool, int, bool> GetParts(
		const string& stream_name1, const int part_id1, vector<uint8_t> &v_data1, uint64_t &metadata1,
		const string& stream_name2, const int part_id2, vector<uint8_t> &v_data2, uint64_t &metadata2);

	void SetRawSize(const int stream_id, const size_t raw_size);
	size_t GetRawSize(const int stream_id);

	size_t GetNoStreams();
	size_t GetNoParts(const int stream_id);
};

// EOF
#endif