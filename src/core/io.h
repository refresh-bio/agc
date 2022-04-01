#ifndef _IO_H
#define _IO_H

// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-03-16
// *******************************************************************************************

#include <algorithm>
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>

#include <iostream>

using namespace std;

#ifndef _WIN32
#define my_fseek	fseek
#define my_ftell	ftell
#else
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
#include <fcntl.h>
#include <io.h>
#endif

// *******************************************************************************************
// Buffered input file
class CInFile
{
	size_t BUFFER_SIZE = 0;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	size_t buffer_filled;

	size_t file_size;
	size_t before_buffer_bytes;

public:
	// *******************************************************************************************
	CInFile() : f(nullptr), buffer(nullptr), buffer_pos(0), buffer_filled(0), file_size(0), before_buffer_bytes(0)
	{};

	// *******************************************************************************************
	~CInFile()
	{
		if (f)
			fclose(f);
		if (buffer)
			delete[] buffer;
	}

	// *******************************************************************************************
	bool Open(const string &file_name, const size_t _BUFFER_SIZE = 128 << 20)
	{
		if (f)
			return false;

		f = fopen(file_name.c_str(), "rb");
		if (!f)
			return false;

		my_fseek(f, 0, SEEK_END);
		file_size = my_ftell(f);
		my_fseek(f, 0, SEEK_SET);
		before_buffer_bytes = 0;

		if (_BUFFER_SIZE == ~0ull)
			BUFFER_SIZE = file_size;
		else
			BUFFER_SIZE = _BUFFER_SIZE;

		buffer = new uint8_t[BUFFER_SIZE];

		buffer_pos = 0;

		buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);

		return true;
	}

	// *******************************************************************************************
	bool Close()
	{
		if (f)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return true;
	}

	// *******************************************************************************************
	bool IsOpened()
	{
		return f != nullptr;
	}

	// *******************************************************************************************
	int Get()
	{
		if (buffer_pos < buffer_filled)
			return buffer[buffer_pos++];

		if (feof(f))
			return EOF;

		before_buffer_bytes += buffer_filled;

		buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
		if (buffer_filled == 0)
			return EOF;

		buffer_pos = 0;

		return buffer[buffer_pos++];
	}

	// *******************************************************************************************
	bool UnGet()
	{
		if (buffer_pos)
		{
			--buffer_pos;
			return true;
		}

		return false;
	}

	// *******************************************************************************************
	uint64_t ReadUInt(const int no_bytes)
	{
		uint64_t x = 0;
		uint64_t shift = 0;

		for (int i = 0; i < no_bytes; ++i)
		{
			uint64_t c = Get();
			x += c << shift;
			shift += 8;
		}

		return x;
	}

	// *******************************************************************************************
	uint64_t ReadUIntVar()
	{
		uint64_t x = 0;
		uint64_t c = Get();

		if ((c >> 7) == 0)		// [0, 0x8000)
		{
			x = c << 8;
			x += (uint64_t)Get();
		}
		else if ((c >> 6) == 0b10)	// [0x8000, 0x400000)
		{
			x = (c & 0x3f) << 16;
			x += (uint64_t) Get() << 8;
			x += (uint64_t)Get();
		}
		else if ((c >> 6) == 0b11)	// [0x80000, 0x4000000)
		{
			x = (c & 0x3f) << 24;
			x += (uint64_t)Get() << 16;
			x += (uint64_t)Get() << 8;
			x += (uint64_t)Get();
		}

		return x;
	}
		
	// *******************************************************************************************
	void Read(uint8_t *ptr, size_t size)
	{
		if (before_buffer_bytes + buffer_pos + size > file_size)
			size = file_size - (before_buffer_bytes + buffer_pos);

		size_t to_read = size;

		while (buffer_pos + to_read > BUFFER_SIZE)
		{
			memcpy(ptr, buffer + buffer_pos, BUFFER_SIZE - buffer_pos);
			ptr += BUFFER_SIZE - buffer_pos;
			to_read -= BUFFER_SIZE - buffer_pos;

			before_buffer_bytes += buffer_filled;
			buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
			buffer_pos = 0;
		}

		memcpy(ptr, buffer + buffer_pos, to_read);
		buffer_pos += to_read;
	}

	// *******************************************************************************************
	bool Eof() const
	{
		return before_buffer_bytes + buffer_pos >= file_size;
	}

	// *******************************************************************************************
	bool Seek(const size_t requested_pos)
	{
		if (requested_pos >= before_buffer_bytes && requested_pos < before_buffer_bytes + buffer_filled)
			buffer_pos = requested_pos - before_buffer_bytes;
		else
		{
			before_buffer_bytes = requested_pos;
			my_fseek(f, requested_pos, SEEK_SET);
			buffer_filled = fread(buffer, 1, BUFFER_SIZE, f);
			buffer_pos = 0;
		}

		return true;
	}

	// *******************************************************************************************
	size_t FileSize() const
	{
		if (f)
			return file_size;
		else
			return 0;
	}

	// *******************************************************************************************
	size_t GetPos() const
	{
		return before_buffer_bytes + buffer_pos;
	}
};

// *******************************************************************************************
// Buffered output file
class COutFile
{
	size_t BUFFER_SIZE;

	FILE *f;
	uint8_t *buffer;
	size_t buffer_pos;
	bool success;
	bool use_stdout;

public:
	// *******************************************************************************************
	COutFile() : BUFFER_SIZE (8u << 20), f(nullptr), buffer(nullptr), buffer_pos(0), success(false), use_stdout(false)
	{};

	// *******************************************************************************************
	~COutFile()
	{
		if (f)
			Close();
		if (buffer)
			delete[] buffer;
	}

	// *******************************************************************************************
	bool Open(const string &file_name, const size_t _BUFFER_SIZE = 8 << 20)
	{
		if (f)
			return false;

		use_stdout = file_name.empty();

		if (use_stdout)
		{
			f = stdout;
#ifdef _WIN32
			_setmode(_fileno(f), _O_BINARY);
#endif
		}
		else
		{
			f = fopen(file_name.c_str(), "wb");
			if (!f)
				return false;
		}

		BUFFER_SIZE = _BUFFER_SIZE;
		buffer = new uint8_t[BUFFER_SIZE];
		buffer_pos = 0;
		success = true;

		return true;
	}

	// *******************************************************************************************
	bool Close()
	{
		if (buffer_pos)
		{
			success &= fwrite(buffer, 1, buffer_pos, f) == buffer_pos;
			buffer_pos = 0;
		}

		fflush(f);

		if (f && !use_stdout)
		{
			fclose(f);
			f = nullptr;
		}
		if (buffer)
		{
			delete[] buffer;
			buffer = nullptr;
		}

		return success;
	}

	// *******************************************************************************************
	bool IsOpened()
	{
		return f != nullptr;
	}

	// *******************************************************************************************
	void Put(const uint8_t c)
	{
		if (buffer_pos == BUFFER_SIZE)
		{
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;
			buffer_pos = 0;
		}

		buffer[buffer_pos++] = c;
	}

	// *******************************************************************************************
	void Write(const uint8_t *p, size_t n)
	{
		uint8_t *q = (uint8_t *)p;

		while (buffer_pos + n > BUFFER_SIZE)
		{
			size_t small_n = BUFFER_SIZE - buffer_pos;
			memcpy(buffer + buffer_pos, q, small_n);
			success &= fwrite(buffer, 1, BUFFER_SIZE, f) == BUFFER_SIZE;

			buffer_pos = 0;
			n -= small_n;
			q += small_n;
		}

		memcpy(buffer+buffer_pos, q, n);
		buffer_pos += n;
	}

	// *******************************************************************************************
	void WriteUInt(const uint64_t _x, const int no_bytes)
	{
		uint64_t x = _x;

		for (int i = 0; i < no_bytes; ++i)
		{
			Put(x & 0xff);
			x >>= 8;
		}
	}

	// [0, 0x8000) -> 0 [15-bit value]
	// [0x8000, 0x400000) -> 10 [22-bit value]
	// [0x400000, 0x40000000) -> 11 [30-bit value]
	// Larger values are not supported!
	void WriteUIntVar(const uint64_t x)
	{
		if (x < 0x8000ull)		
		{
			Put((uint8_t) (x >> 8));
			Put((uint8_t) (x & 0xff));
		}
		else if (x < 0x400000ull)
		{
			Put((uint8_t) (0x80 + (x >> 16)));
			Put((uint8_t) ((x >> 8) & 0xff));
			Put((uint8_t) (x & 0xff));
		}
		else if (x < 0x40000000ull)
		{
			Put((uint8_t) (0xc0 + (x >> 24)));
			Put((uint8_t) ((x >> 16) & 0xff));
			Put((uint8_t) ((x >> 8) & 0xff));
			Put((uint8_t) (x & 0xff));
		}
		else
			cerr << "Too large value\n";
	}

	// *******************************************************************************************
	void Write(const string &s)
	{
		Write((uint8_t*)s.c_str(), s.size());
	}

	// *******************************************************************************************
	void Write(const string &s, const size_t start_pos, const size_t len)
	{
		Write((uint8_t*)s.c_str() + start_pos, len);
	}
};

// EOF
#endif