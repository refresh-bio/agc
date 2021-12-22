// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#ifndef AGC_API_H
#define AGC_API_H

#ifdef __cplusplus
// *******************************************************************************************
// C++ API (available only when using C++ compiler)
// *******************************************************************************************
#include <vector>
#include <string>
#include <memory>

// *******************************************************************************************
class CAGCFile
{
	std::unique_ptr<class CAGCDecompressorLibrary> agc;
	bool is_opened;

public:
	CAGCFile();
	~CAGCFile();

	/**
	 * @param file_name		file name
	 * @param prefetching	true to preload whole file into memory (faster if you plan series of sequence queries), false otherwise
	 *
	 * @return false for error
	 */
	bool Open(const std::string& file_name, bool prefetching = true);

	/**
	 * @return true for success and false for error
	 */
	bool Close();

	/**
	 * Get the length of a contig. Return an error if _name_ is not present, or if
	 * _name_ is not unique but _sample_ is an empty string,
	 *
	 * @param sample   sample name; can be an empty string
	 * @param name     contig name
	 *
	 * @return contig length, or <0 for errors
	 */
	int GetCtgLen(const std::string& sample, const std::string& name) const;

	/**
	 * @param sample   sample name; can be an empty string
	 * @param name     contig name
	 * @param start    start offset
	 * @param end      end offset
	 * @param buf      sequence to be written; user should allocate memory (returned value)
	 *
	 * @return contig length, or <0 for errors
	 */
	int GetCtgSeq(const std::string& sample, const std::string& name, int start, int end, std::string& buffer) const;

	/**
	 * @return the number of samples
	 */
	int NSample() const;

	/**
	 * @param sample   sample name
	 *
	 * @return the number of contigs in sample
	 */
	int NCtg(const std::string& sample) const;

	/**
	 * @param samples  vector of strings with sample names (returned value)
	 *
	 * @return number of samples to be written to
	 */
	int ListSample(std::vector<std::string>& samples) const;

	/**
	 * @param sample    sample name; can be an empty string
	 * @param names     vector of strings with contig names (returned value)
	 *
	 * @return number of contigs in the sample
	 */
	int ListCtg(const std::string& sample, std::vector<std::string>& names) const;
};

typedef CAGCFile agc_t;
#define EXTERNC extern "C"
#else
typedef struct agc_t agc_t;
#define EXTERNC
#endif

// *******************************************************************************************
// C version of the API (can be used in C or C++ code)
// *******************************************************************************************

/**
 * @param fn			file name
 * @param prefetching	1 to preload whole file into memory (faster if you plan series of sequence queries), 0 otherwise
 *
 * @return NULL for error
 */
EXTERNC agc_t* agc_open(char* fn, int prefetching);

/**
 * @param fp   agc handle
 *
 * @return 0 for success and -1 for error
 */
EXTERNC int agc_close(agc_t* agc);

/**
 * Get the length of a contig. Return an error if _name_ is not present, or if
 * _name_ is not unique but _sample_ is NULL, 
 *
 * @param agc      agc handle
 * @param sample   sample name; can be NULL
 * @param name     contig name
 *
 * @return contig length, or <0 for errors
 */
EXTERNC int agc_get_ctg_len(const agc_t *agc, const char *sample, const char *name);

/**
 * @param agc      agc handle
 * @param sample   sample name; can be NULL
 * @param name     contig name
 * @param start    start offset
 * @param end      end offset
 * @param buf      sequence to be written; user should allocate memory (returned value)
 *
 * @return contig length, or <0 for errors
 */
EXTERNC int agc_get_ctg_seq(const agc_t *agc, const char *sample, const char *name, int start, int end, char *buf);

/**
 * @param agc      agc handle
 *
 * @return the number of samples
 */
EXTERNC int agc_n_sample(const agc_t* agc);

/**
 * @param agc      agc handle
 * @param sample   sample name
 *
 * @return the number of contigs in sample
 */
EXTERNC int agc_n_ctg(const agc_t *agc, const char *sample);

/**
 * @param agc       agc handle
 * @param n_sample  number of samples to be written to (returned value)
 *
 * @return array of NULL-terminated strings. Use agc_list_destroy() to deallocate.
 */
EXTERNC char **agc_list_sample(const agc_t *agc, int *n_sample);

/**
 * @param agc       agc handle
 * @param sample    sample name; can be NULL
 * @param n_ctg     number of contigs in the sample (returned value)
 *
 * @return array of NULL-terminated strings. Use agc_list_destroy() to deallocate.
 */
EXTERNC char **agc_list_ctg(const agc_t *agc, const char *sample, int *n_ctg);

/**
 * Deallocate an array of strings returned by agc_list_samples or agc_list_ctg
 *
 * @param list      array to deallocate
 */
EXTERNC int agc_list_destroy(char **list);

#endif

// EOF
