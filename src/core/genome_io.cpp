// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2024, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.2
// Date   : 2024-11-21
// *******************************************************************************************

#include "genome_io.h"
#include <algorithm>
#include <numeric>
#include <cctype>
#include <cstring>

// *******************************************************************************************
CGenomeIO::CGenomeIO()
{
	writing = false;

	out = nullptr;

	is_gzipped = false;
	use_stdout = false;

	buffer = nullptr;
	buffer_pos = 0;
	buffer_filled = 0;
}

// *******************************************************************************************
CGenomeIO::~CGenomeIO()
{
	Close();

	if (sdf)
		delete sdf;
	if (sif)
		delete sif;
}

// *******************************************************************************************
bool CGenomeIO::Open(const string& _file_name, const bool _writing)
{
	if (out || sif)
		return false;

	use_stdout = _file_name.empty();
	writing = _writing;

	if (writing)
	{
		if (use_stdout)
		{
			out = stdout;
#ifdef _WIN32
			_setmode(_fileno(out), _O_BINARY);
#endif
		}
		else
		{
			out = fopen(_file_name.c_str(), "wb");
			if (!out)
				return false;
			setvbuf(out, nullptr, _IOFBF, write_buffer_size);

		}

		buffer = nullptr;
	}
	else
	{
		sif = new refresh::stream_in_file(_file_name);
		if (!sif->is_open())
			return false;

		sdf = new refresh::stream_decompression(sif);

		buffer = new uint8_t[read_buffer_size];
	}

	buffer_pos = 0;
	buffer_filled = 0;

	return true;
}

// *******************************************************************************************
bool CGenomeIO::Close()
{
	if (writing)
	{
		if (out)
		{
			fflush(out);
			if(!use_stdout)
				fclose(out);
			out = nullptr;
		}
	}
	else
	{
		if (sif)
		{
			delete sdf;
			delete sif;
			sif = nullptr;
			sdf = nullptr;
		}
	}

	if (buffer)
		delete[] buffer;
	buffer = nullptr;

	return true;
}

// *******************************************************************************************
// Can be used only just after opening the file, i.e., prior to any reads
size_t CGenomeIO::FileSize() 
{
	size_t s = 0;

	if (!writing)
	{
		const size_t loc_buf_size = 1 << 25;
		char* loc_buf = new char[loc_buf_size];
		size_t readed;

		while (true)
		{
			sdf->read(loc_buf, loc_buf_size, readed);
			if (!readed)
				break;

			s += readed;
		}

		delete[] loc_buf;

		sif->restart();
		sdf->restart(sif);
	}

	return s;
}

// *******************************************************************************************
bool CGenomeIO::ReadContig(string& id, contig_t& contig)
{
	return !writing && read_contig(id, contig, false);
}

// *******************************************************************************************
bool CGenomeIO::ReadContigConverted(string& id, contig_t& contig)
{
	return !writing && read_contig(id, contig, true);
}

#if 0
// *******************************************************************************************
bool CGenomeIO::SaveContigConverted(const string& id, const contig_t& contig, const uint32_t line_length)
{
	return save_contig(id, contig, line_length, true);
}

// *******************************************************************************************
bool CGenomeIO::SaveContig(const string& id, const contig_t& contig, const uint32_t line_length)
{
	return save_contig(id, contig, line_length, false);
}
#endif

// *******************************************************************************************
bool CGenomeIO::SaveContigDirectly(const string& id, const contig_t& contig, const uint32_t gzip_level)
{
	return save_contig_directly(id, contig, gzip_level);
}

// *******************************************************************************************
bool CGenomeIO::fill_buffer()
{
	if (buffer_pos < buffer_filled)
	{
		memmove(buffer, buffer + buffer_pos, buffer_filled - buffer_pos);
		buffer_filled -= buffer_pos;
		buffer_pos = 0;
	}
	else
	{
		buffer_pos = 0;
		buffer_filled = 0;
	}

	size_t to_read = read_buffer_size - buffer_filled;
	size_t readed;
	
	sdf->read((char*) buffer + buffer_filled, to_read, readed);
	
	buffer_filled += readed;

	return buffer_filled != 0;
}

// *******************************************************************************************
bool CGenomeIO::read_contig_raw(string& id, contig_t& contig)
{
	if (!sif)
		return false;

	id.clear();
	contig.clear();

	// Read id
	while (true)
	{
		if (eof())
			if (!fill_buffer())
				return false;
		
		uint8_t c = buffer[buffer_pos++];
		if (c == '\n' || c == '\r')
			break;

		id.push_back((char) c);
	}

	if (!id.empty())
		id.erase(id.begin());

	// Read contig
	while (true)
	{
		int next_id_pos = find_contig_end();

		if (next_id_pos >= 0)
		{
			contig.insert(contig.end(), buffer + buffer_pos, buffer + next_id_pos);
			buffer_pos = (size_t) next_id_pos;
			break;
		}

		contig.insert(contig.end(), buffer + buffer_pos, buffer + buffer_filled);
		buffer_pos = buffer_filled;
		if (!fill_buffer())
			break;
	}

	return !id.empty() && !contig.empty();
}

// *******************************************************************************************
bool CGenomeIO::ReadContigRaw(string& id, contig_t& contig)
{
	return read_contig_raw(id, contig);
}

// *******************************************************************************************
int CGenomeIO::find_contig_end()
{
	auto p = find(buffer + buffer_pos, buffer + buffer_filled, '>');

	if (p == buffer + buffer_filled)
		return -1;

	return (int) (p - buffer);
}

// *******************************************************************************************
bool CGenomeIO::read_contig(string& id, contig_t& contig, const bool converted)
{
	if (!read_contig_raw(id, contig))
		return false;

	size_t len = contig.size();
	size_t in_pos = 0;
	size_t out_pos = 0;

	if(converted)
		for (; in_pos < len; ++in_pos)
		{
			auto c = contig[in_pos];
			if (c > 64)
				contig[out_pos++] = cnv_num[c];
		}
	else
		for (; in_pos < len; ++in_pos)
		{
			auto c = contig[in_pos];
			if (c > 64)
				contig[out_pos++] = c;
		}

	contig.resize(out_pos);

	return true;
}

#if 0
// *******************************************************************************************
bool CGenomeIO::save_contig(const string& id, const contig_t& contig, const uint32_t line_length, const bool converted)
{
	if(!writing || (!out && !gz_out))
		return false;

	if (is_gzipped)
	{
		gzputc(gz_out, '>');
		gzputs(gz_out, id.c_str());
		gzputc(gz_out, '\n');
	}
	else
	{
		putc('>', out);
		fputs(id.c_str(), out);
		putc('\n', out);
	}

	if (converted)
		save_contig_imp_cnv(contig, line_length);
	else
		save_contig_imp(contig, line_length);

	return true;
}
#endif

// *******************************************************************************************
bool CGenomeIO::save_contig_directly(const string& id, const contig_t& contig, const uint32_t gzip_level)
{
	string tmp = ">" + id + "\n";

	if (gzip_level)
	{
		auto overhead = gzip_zero_compressor.get_overhead(tmp.size());
		gzip_zero_compressor_buffer.resize(tmp.size() + overhead);
		auto packed_size = gzip_zero_compressor.compress(tmp.data(), tmp.size(), gzip_zero_compressor_buffer.data(), gzip_zero_compressor_buffer.size());

		if (fwrite(gzip_zero_compressor_buffer.data(), 1, packed_size, out) != packed_size)
			return false;
	}
	else
	{
		if (fwrite(tmp.data(), 1, tmp.size(), out) != tmp.size())
			return false;
	}

	return fwrite(contig.data(), 1, contig.size(), out) == contig.size();
}

#if 0
// *******************************************************************************************
void CGenomeIO::save_contig_imp(const contig_t& contig, const uint32_t line_length)
{
	size_t to_save = contig.size();
	auto p = contig.data();

	if (is_gzipped)
	{
		for (; to_save > line_length; to_save -= line_length)
		{
			gzwrite(gz_out, p, line_length);
			gzputc(gz_out, '\n');

			p += line_length;
		}

		for (uint32_t i = 0; i < to_save; ++i)
			gzputc(gz_out, *p++);

		gzputc(gz_out, '\n');
	}
	else
	{
		for (; to_save > line_length; to_save -= line_length)
		{
			fwrite(p, 1, line_length, out);
			putc('\n', out);

			p += line_length;
		}

		for (uint32_t i = 0; i < to_save; ++i)
			fputc(*p++, out);

		putc('\n', out);
	}
}

// *******************************************************************************************
void CGenomeIO::save_contig_imp_cnv(const contig_t& contig, const uint32_t line_length)
{
	size_t to_save = contig.size();
	uint8_t* line = new uint8_t[line_length + 1u];
	line[line_length] = '\n';
	auto p = contig.data();

	for(; to_save > line_length; to_save -= line_length)
	{
		uint32_t i;
		uint8_t *q = line;

		switch (i = line_length % 8)
		{
		case 7:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 6:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 5:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 4:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 3:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 2:	*q++ = cnv_num[*p++]; [[fallthrough]];
		case 1:	*q++ = cnv_num[*p++]; 
		}
		
		for (; i < line_length; i += 8)
		{
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
			*q++ = cnv_num[*p++];
		}

		is_gzipped ? gzwrite(gz_out, line, line_length + 1) : fwrite(line, 1, line_length + 1, out);
	}

	for (uint32_t i = 0; i < to_save; ++i)
		is_gzipped ? gzputc(gz_out, cnv_num[*p++]) : fputc(cnv_num[*p++], out);

	is_gzipped ? gzputc(gz_out, '\n') : putc('\n', out);

	delete[] line;
}
#endif

// EOF
