// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 3.1
// Date   : 2024-03-12
// *******************************************************************************************

#include "../core/agc_decompressor_lib.h"
#include "agc-api.h"
#include <cstring>

// *******************************************************************************************
char** agc_internal_cnv_vec2list(vector<string>& vec);

// *******************************************************************************************
// C++ part
// *******************************************************************************************

// *******************************************************************************************
CAGCFile::CAGCFile()
{
	agc = std::make_unique<CAGCDecompressorLibrary>(false);
	is_opened = false;
}

// *******************************************************************************************
CAGCFile::~CAGCFile()
{
}

// *******************************************************************************************
bool CAGCFile::Open(const std::string& file_name, bool prefetching)
{
	if (agc->IsOpened())
		return false;

	is_opened = agc->Open(file_name, prefetching);

	return is_opened;
}

// *******************************************************************************************
bool CAGCFile::Close()
{
	if (!is_opened)
		return false;

	return agc->Close();
}

// *******************************************************************************************
int CAGCFile::GetCtgLen(const std::string& sample, const std::string& name)	const
{
	if (!is_opened)
		return -1;

	return (int) agc->GetContigLength(sample, name);
}

// *******************************************************************************************
int CAGCFile::GetCtgSeq(const std::string& sample, const std::string& name, int start, int end, std::string& buffer) const
{
	if (!is_opened)
		return -1;

	return agc->GetContigString(sample, name, start, end, buffer);
}

// *******************************************************************************************
int CAGCFile::NSample() const
{
	if (!is_opened)
		return -1;

	return agc->GetNoSamples();
}

// *******************************************************************************************
int CAGCFile::NCtg(const std::string& sample)	const
{
	if (!is_opened)
		return -1;

	return agc->GetNoContigs(sample);
}

// *******************************************************************************************
int CAGCFile::GetReferenceSample(std::string& sample) const
{
	if (!is_opened)
		return -1;

	agc->GetReferenceSample(sample);

	return 0;
}

// *******************************************************************************************
int CAGCFile::ListSample(std::vector<std::string>& samples) const
{
	if (!is_opened)
		return -1;

	agc->ListSamples(samples);

	return 0;
}

// *******************************************************************************************
int CAGCFile::ListCtg(const std::string& sample, std::vector<std::string>& names) const
{
	if (!is_opened)
		return -1;

	agc->ListContigs(sample, names);

	return 0;
}

// *******************************************************************************************
// C part
// *******************************************************************************************
agc_t* agc_open(char* fn, int prefetching)
{
	agc_t* agc = new CAGCFile();
	bool r = agc->Open(fn, (bool)prefetching);

	if (!r)
	{
		delete agc;
		agc = NULL;
	}

	return agc;
}

// *******************************************************************************************
int agc_close(agc_t* agc)
{
	if (!agc)
		return -1;

	return agc->Close() ? 0 : -1;
}

// *******************************************************************************************
int agc_n_sample(const agc_t* agc)
{
	if (!agc)
		return -1;

	return agc->NSample();
}

// *******************************************************************************************
int agc_get_ctg_seq(const agc_t* agc, const char* sample, const char* name, int start, int end, char* buf)
{
	if (!agc)
		return -1;

	string buffer;

	if (agc->GetCtgSeq(sample ? sample : "", name, start, end, buffer) != 0)
		return -1;

	strcpy(buf, buffer.data());

	return (int) buffer.size();
}

// *******************************************************************************************
int agc_get_ctg_len(const agc_t* agc, const char* sample, const char* name)
{
	if (!agc)
		return -1;

	return agc->GetCtgLen(sample ? sample : "", name);
}

// *******************************************************************************************
int agc_n_ctg(const agc_t* agc, const char* sample)
{
	if (!agc)
		return -1;

	return agc->NCtg(sample);
}

// *******************************************************************************************
char* agc_reference_sample(const agc_t* agc)
{
	if (!agc)
		return NULL;

	string sample;

	if (agc->GetReferenceSample(sample) < 0)
		return NULL;

	char* c_sample = (char*) malloc(sample.size() + 1);
	strcpy(c_sample, sample.c_str());

	return c_sample;
}

// *******************************************************************************************
char** agc_list_sample(const agc_t* agc, int* n_sample)
{
	if (!agc)
		return NULL;

	vector<string> v_samples;

	agc->ListSample(v_samples);

	*n_sample = (int)v_samples.size();

	return agc_internal_cnv_vec2list(v_samples);
}

// *******************************************************************************************
char** agc_internal_cnv_vec2list(vector<string>& vec)
{
	char** list = (char**)malloc(sizeof(char*) * (vec.size() + 1));

	if (!list)
		return NULL;

	char** p = list;

	for (auto& s : vec)
	{
		*p = (char*)malloc(s.size() + 1);
		if (!*p)
		{
			for (auto q = list; q != p; ++q)
				free(*q);
			free(list);

			return NULL;
		}

		strcpy(*p, s.data());
		++p;
	}

	*p = NULL;

	return list;
}

// *******************************************************************************************
char** agc_list_ctg(const agc_t* agc, const char* sample, int* n_ctg)
{
	if (!agc)
		return NULL;

	vector<string> v_contigs;

	agc->ListCtg(sample, v_contigs);

	*n_ctg = (int)v_contigs.size();

	return agc_internal_cnv_vec2list(v_contigs);
}

// *******************************************************************************************
int agc_list_destroy(char** list)
{
	for (char** p = list; p; ++p)
		free(*p);

	free(list);

	return 0;
}

// *******************************************************************************************
int agc_list_destroy(char* sample)
{
	free(sample);

	return 0;
}


// EOF