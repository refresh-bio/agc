// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021, S.Deorowicz, A.Danek, H.Li
//
// Version: 1.0
// Date   : 2021-12-17
// *******************************************************************************************

#include "application.h"
#include <algorithm>
#include <fstream>
#include <iterator>
#include "../core/utils.h"
#include "../../libs/ketopt.h"

// *******************************************************************************************
bool CApplication::parse_params(const int argc, const char** argv)
{
    if (argc == 1)
    {
        usage();
        return false;
    }

    execution_params.mode = string(argv[1]);

    if (argc == 2)
    {
        if (execution_params.mode == "create")
            usage_create();
        else if (execution_params.mode == "append")
            usage_append();
        else if (execution_params.mode == "getset")
            usage_getset();
        else if (execution_params.mode == "getctg")
            usage_getctg();
        else if (execution_params.mode == "listset")
            usage_listset();
        else if (execution_params.mode == "listctg")
            usage_listctg();
        else if (execution_params.mode == "info")
            usage_info();
        else
        {
            cerr << "Unknown mode: " << execution_params.mode << endl;
            return false;
        }
    }
    else
    {
        if (execution_params.mode == "create")
            return parse_params_create(argc - 1, argv + 1);
        else if (execution_params.mode == "append")
            return parse_params_append(argc - 1, argv + 1);
        else if (execution_params.mode == "getset")
            return parse_params_getset(argc - 1, argv + 1);
        else if (execution_params.mode == "getctg")
            return parse_params_getctg(argc - 1, argv + 1);
        else if (execution_params.mode == "listset")
            return parse_params_listset(argc - 1, argv + 1);
        else if (execution_params.mode == "listctg")
            return parse_params_listctg(argc - 1, argv + 1);
        else if (execution_params.mode == "info")
            return parse_params_info(argc - 1, argv + 1);
        else
        {
            cerr << "Unknown mode: " << execution_params.mode << endl;
            return false;
        }
    }

    return false;
}

// *******************************************************************************************
void CApplication::usage() const
{
	cerr << AGC_VERSION << endl;
    cerr << "Usage: agc <command> [options]\n";
    cerr << "Command:\n";
    cerr << "   create   - create archive from FASTA files\n";
    cerr << "   append   - add FASTA files to existing archive\n";
    cerr << "   getset   - extract sample from archive\n";
    cerr << "   getctg   - extract contig from archive\n";
    cerr << "   listset  - list sample names in archive\n";
    cerr << "   listctg  - list sample and contig names in archive\n";
    cerr << "   info     - show some statistics of the compressed data\n";
    cerr << "Note: run agc <command> to see command-specific options\n";
}

// *******************************************************************************************
void CApplication::usage_create() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc create [options] <ref.fa> [<in1.fa> ...] > <out.agc>\n";
    cerr << "Options:\n";
    cerr << "   -b <int>       - batch size " << execution_params.pack_cardinality.info() << "\n";
    cerr << "   -c             - concatenated genomes in a single file (default: " << bool2string(execution_params.concatenated_genomes) << ")\n";
    cerr << "   -d             - do not store cmd-line (default: " << bool2string(execution_params.store_cmd_line) << ")\n";
	cerr << "   -f             - use fast mode (default: " << bool2string(!execution_params.reproducibility_mode) << ")\n";
	cerr << "   -i <file_name> - file with FASTA file names (alterantive to listing file names explicitely in command line)\n";
    cerr << "   -k <int>       - k-mer length" << execution_params.k.info() << "\n";
    cerr << "   -l <int>       - min. match length " << execution_params.min_match_length.info() << "\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
	cerr << "   -s <int>       - expected segment size " << execution_params.segment_size.info() << "\n";
    cerr << "   -t <int>       - no of threads " << execution_params.no_threads.info() << "\n";
    cerr << "   -v <int>       - verbosity level " << execution_params.verbosity.info() << "\n";
}

// *******************************************************************************************
bool CApplication::parse_params_create(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int i, c;

	while ((c = ketopt(&o, argc, argv, 1, "t:b:s:k:l:cdfi:o:v:", 0)) >= 0) {
		if (c == 't') {
			execution_params.no_threads.assign(atoi(o.arg));
		} else if (c == 'b') {
			execution_params.pack_cardinality.assign(atoi(o.arg));
		} else if (c == 's') {
			execution_params.segment_size.assign(atoi(o.arg));
		} else if (c == 'k') {
			execution_params.k.assign(atoi(o.arg));
		} else if (c == 'l') {
			execution_params.min_match_length.assign(atoi(o.arg));
		} else if (c == 'c') {
			execution_params.concatenated_genomes = true;
		} else if (c == 'f') {
			execution_params.reproducibility_mode = false;
		} else if (c == 'd') {
			execution_params.store_cmd_line = false;
		} else if (c == 'i') {
			if (!load_file_names(o.arg, execution_params.input_names))
				return false;
		} else if (c == 'o') {
			execution_params.out_archive_name = o.arg;
			execution_params.use_stdout = false;
		} else if (c == 'v') {
			execution_params.verbosity.assign(atoi(o.arg));
		}
	}

	if (o.ind >= argc) {
		cerr << "No reference file name\n";
		return false;
	}

	execution_params.input_names.insert(execution_params.input_names.begin(), string(argv[o.ind]));

	for (i = o.ind + 1; i < argc; ++i)
		execution_params.input_names.emplace_back(argv[i]);

	return true;
}

// *******************************************************************************************
void CApplication::usage_append() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc append [options] <in.agc> [<in1.fa> ...] > <out.agc>\n";
    cerr << "Options:\n";
    cerr << "   -c             - concatenated genomes in a single file (default: " << bool2string(execution_params.concatenated_genomes) << ")\n";
    cerr << "   -d             - do not store cmd-line (default: " << bool2string(execution_params.store_cmd_line) << ")\n";
	cerr << "   -f             - use fast mode (default: " << bool2string(!execution_params.reproducibility_mode) << ")\n";
	cerr << "   -i <file_name> - file with FASTA file names (alterantive to listing file names explicitely in command line)\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
	cerr << "   -t <int>       - no of threads " << execution_params.no_threads.info() << "\n";
    cerr << "   -v <int>       - verbosity level " << execution_params.verbosity.info() << "\n";
}

// *******************************************************************************************
bool CApplication::parse_params_append(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int i, c;

	while ((c = ketopt(&o, argc, argv, 1, "t:cdfi:o:v:", 0)) >= 0) {
		if (c == 't') {
			execution_params.no_threads.assign(atoi(o.arg));
		} else if (c == 'c') {
			execution_params.concatenated_genomes = true;
		} else if (c == 'f') {
			execution_params.reproducibility_mode = false;
		} else if (c == 'd') {
			execution_params.store_cmd_line = false;
		} else if (c == 'i') {
			if (!load_file_names(o.arg, execution_params.input_names))
				return false;
		} else if (c == 'o') {
			execution_params.out_archive_name = o.arg;
			execution_params.use_stdout = false;
		} else if (c == 'v') {
			execution_params.verbosity.assign(atoi(o.arg));
		}
	}

	if (o.ind >= argc) {
		cerr << "No archive name\n";
		return false;
	}

	execution_params.in_archive_name = argv[o.ind];

	for (i = o.ind + 1; i < argc; ++i)
		execution_params.input_names.emplace_back(argv[i]);

	return true;
}

// *******************************************************************************************
void CApplication::usage_getset() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc getset [options] <in.agc> <sample_name1> [<sample_name2> ...] > <out.fa>\n";
    cerr << "Options:\n";
    cerr << "   -l <int>       - line length " << execution_params.line_length.info() << "\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
    cerr << "   -t <int>       - no of threads " << execution_params.no_threads.info() << "\n";
    cerr << "   -v <int>       - verbosity level " << execution_params.verbosity.info() << "\n";
}

// *******************************************************************************************
bool CApplication::parse_params_getset(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int c, i;

	execution_params.prefetch = true;

	while ((c = ketopt(&o, argc, argv, 1, "t:l:o:v:", 0)) >= 0) {
		if (c == 't') {
			execution_params.no_threads.assign(atoi(o.arg));
		} else if (c == 'l') {
			execution_params.line_length.assign(atoi(o.arg));
		} else if (c == 'o') {
			execution_params.output_name = o.arg;
			execution_params.use_stdout = false;
		} else if (c == 'v') {
			execution_params.verbosity.assign(atoi(o.arg));
		}
	}

	if (o.ind >= argc) {
		cerr << "No archive name\n";
		return false;
	}

	execution_params.in_archive_name = argv[o.ind];

	if (o.ind + 1 >= argc) {
		cerr << "No sample name\n";
		return false;
	}

	for (i = o.ind + 1; i < argc; ++i)
		execution_params.sample_names.emplace_back(argv[i]);

	return true;
}

// *******************************************************************************************
void CApplication::usage_getctg() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc getctg [options] <in.agc> <contig1> [<contig2> ...] > <out.fa>\n";
    cerr << "       agc getctg [options] <in.agc> <contig1@sample1> [<contig2@sample2> ...] > <out.fa>\n";
    cerr << "       agc getctg [options] <in.agc> <contig1:from-to>[<contig2:from-to> ...] > <out.fa>\n";
    cerr << "       agc getctg [options] <in.agc> <contig1@sample1:from-to> [<contig2@sample2:from-to> ...] > <out.fa>\n";
    cerr << "Options:\n";
    cerr << "   -l <int>       - line length " << execution_params.line_length.info() << "\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
    cerr << "   -t <int>       - no of threads " << execution_params.no_threads.info() << "\n";
    cerr << "   -p             - disable file prefetching (useful for short queries)" << "\n";
    cerr << "   -v <int>       - verbosity level " << execution_params.verbosity.info() << "\n";
}

// *******************************************************************************************
bool CApplication::parse_params_getctg(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int c, i;

	execution_params.prefetch = true;

	while ((c = ketopt(&o, argc, argv, 1, "t:l:o:pv:", 0)) >= 0) {
		if (c == 't') {
			execution_params.no_threads.assign(atoi(o.arg));
		} else if (c == 'l') {
			execution_params.line_length.assign(atoi(o.arg));
		} else if (c == 'o') {
			execution_params.output_name = o.arg;
			execution_params.use_stdout = false;
		} else if (c == 'p') {
			execution_params.prefetch = false;
		} else if (c == 'v') {
			execution_params.verbosity.assign(atoi(o.arg));
		}
	}

	if (o.ind >= argc) {
		cerr << "No archive name\n";
		return false;
	}

	execution_params.in_archive_name = argv[o.ind];

	if (o.ind + 1 >= argc) {
		cerr << "No contig name\n";
		return false;
	}

	for (i = o.ind + 1; i < argc; ++i)
		execution_params.contig_names.emplace_back(argv[i]);

	return true;
}

// *******************************************************************************************
void CApplication::usage_listset() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc listset [options] <in.agc> > <out.txt>\n";
    cerr << "Options:\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
}

// *******************************************************************************************
bool CApplication::parse_params_listset(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int c;

	execution_params.prefetch = false;

	while ((c = ketopt(&o, argc, argv, 1, "o:", 0)) >= 0) {
		if (c == 'o') {
			execution_params.output_name = o.arg;
			execution_params.use_stdout = false;
		}
	}

	if (o.ind >= argc) {
		cerr << "No archive name\n";
		return false;
	}

	execution_params.in_archive_name = argv[o.ind];

	return true;
}

// *******************************************************************************************
void CApplication::usage_listctg() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc listctg [options] <in.agc> <sample1> [<sample2> ...] > <out.txt>\n";
    cerr << "Options:\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
}

// *******************************************************************************************
bool CApplication::parse_params_listctg(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int c, i;

	execution_params.prefetch = false;

	while ((c = ketopt(&o, argc, argv, 1, "o:", 0)) >= 0) {
		if (c == 'o') {
			execution_params.output_name = o.arg;
			execution_params.use_stdout = false;
		}
	}

	if (o.ind >= argc) {
		cerr << "No archive name\n";
		return false;
	}

	execution_params.in_archive_name = argv[o.ind];

	if (o.ind + 1 >= argc) {
		cerr << "No sample name\n";
		return false;
	}

	for (i = o.ind + 1; i < argc; ++i)
		execution_params.sample_names.emplace_back(argv[i]);

	return true;
}

// *******************************************************************************************
void CApplication::usage_info() const
{
	cerr << AGC_VERSION << endl;
	cerr << "Usage: agc info [options] <in.agc> > <out.txt>\n";
    cerr << "Options:\n";
    cerr << "   -o <file_name> - output to file (default: output is sent to stdout)\n";
//    cerr << "   -v <int>       - verbosity level " << execution_params.verbosity.info() << "\n";      // Valid but hidden option
}

// *******************************************************************************************
bool CApplication::parse_params_info(const int argc, const char** argv)
{
	ketopt_t o = KETOPT_INIT;
	int c;

	execution_params.prefetch = false;
	execution_params.verbosity.assign(0);

	while ((c = ketopt(&o, argc, argv, 1, "o:v:", 0)) >= 0) {
		if (c == 'o') {
			execution_params.output_name = o.arg;
			execution_params.use_stdout = false;
		} else if (c == 'v') {
			execution_params.verbosity.assign(atoi(o.arg));
		}
	}

	if (o.ind >= argc) {
		cerr << "No archive name\n";
		return false;
	}

	execution_params.in_archive_name = argv[o.ind];

	return true;
}

// *******************************************************************************************
bool CApplication::load_file_names(const string &fn, vector<string>& v_file_names)
{
    fstream inf(fn);

    if (!inf.is_open())
    {
        cerr << "Cannot open file: " << fn << endl;
        return false;
    }

    v_file_names.assign(istream_iterator<string>(inf), istream_iterator<string>());

    return true;
}

// *******************************************************************************************
string CApplication::bool2string(bool x) const
{
    return x ? "true"s : "false"s;
}

// EOF
