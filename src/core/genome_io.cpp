// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.0
// Date   : 2022-03-16
// *******************************************************************************************

#include "../core/genome_io.h"
#include <algorithm>
#include <numeric>
#include <cctype>
#include <cstring>

// *******************************************************************************************
CGenomeIO::CGenomeIO()
{
	writing = false;

	in = nullptr;
	out = nullptr;
	gz_in = nullptr;
	gz_out = nullptr;

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
}

// *******************************************************************************************
bool CGenomeIO::Open(const string& _file_name, const bool _writing)
{
	if (in || out || gz_in || gz_out)
		return false;

	is_gzipped = _file_name.length() > 3 && _file_name.substr(_file_name.length() - 3, 3) == ".gz"s;
	use_stdout = _file_name.empty();
	writing = _writing;

	if (writing)
	{
		if (is_gzipped)
		{
			gz_out = gzopen(_file_name.c_str(), "w3");
			if (!gz_out)
				return false;
			gzbuffer(gz_out, (uint32_t) gz_buffer_size);
		}
		else
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
				setvbuf(out, nullptr, _IOFBF, buffer_size);

			}
		}

		buffer = nullptr;
	}
	else
	{
		if (is_gzipped)
		{
			gz_in = gzopen(_file_name.c_str(), "r");
			if (!gz_in)
				return false;
			gzbuffer(gz_in, (uint32_t) gz_buffer_size);
		}
		else
		{
			in = fopen(_file_name.c_str(), "rb");
			if (!in)
				return false;
		}

		buffer = new uint8_t[buffer_size];
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
		if (gz_out)
		{
			gzclose(gz_out);
			gz_out = nullptr;
		}
		else if (out)
		{
			fflush(out);
			if(!use_stdout)
				fclose(out);
			out = nullptr;
		}
	}
	else
	{
		if (gz_in)
		{
			gzclose(gz_in);
			gz_in = nullptr;
		}
		else if (in)
		{
			fclose(in);
			in = nullptr;
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
		if (is_gzipped)
		{
			while (fill_buffer())
			{
				s += buffer_filled;
				buffer_pos = buffer_filled;
			}
			gzseek(gz_in, 0, SEEK_SET);
		}
		else
		{
			fseek(in, 0, SEEK_END);
			s = my_ftell(in);
			fseek(in, 0, SEEK_SET);
		}
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

	size_t to_read = buffer_size - buffer_filled;
	size_t readed;
	
	if (is_gzipped)
		readed = gzread(gz_in, buffer + buffer_filled, (uint32_t) to_read);
	else
		readed = fread(buffer + buffer_filled, 1, to_read, in);
	
	buffer_filled += readed;

	return buffer_filled != 0;
}

// *******************************************************************************************
bool CGenomeIO::read_contig_raw(string& id, contig_t& contig)
{
	if (!in && !gz_in)
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

	return p - buffer;
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
		case 7:	*q++ = cnv_num[*p++];
		case 6:	*q++ = cnv_num[*p++];
		case 5:	*q++ = cnv_num[*p++];
		case 4:	*q++ = cnv_num[*p++];
		case 3:	*q++ = cnv_num[*p++];
		case 2:	*q++ = cnv_num[*p++];
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

// EOF
