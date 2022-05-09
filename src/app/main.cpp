// *******************************************************************************************
// This file is a part of AGC software distributed under MIT license.
// The homepage of the AGC project is https://github.com/refresh-bio/agc
//
// Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
//
// Version: 2.1
// Date   : 2022-05-06
// *******************************************************************************************

#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>

#include "../app/application.h"
#include "../core/agc_compressor.h"
#include "../core/agc_decompressor.h"

using namespace std;
using namespace std::chrono;

// *******************************************************************************************
int CApplication::Run(const int argc, const char** argv)
{
    if (!parse_params(argc, argv))
        return 0;

    for (int i = 0; i < argc; ++i)
        cmd_line += string(argv[i]) + " ";

    cmd_line.pop_back();

    auto t1 = chrono::high_resolution_clock::now();

    if (execution_params.mode == "create")
        create();
    else if (execution_params.mode == "append")
        append();
    else if (execution_params.mode == "getcol")
        getcol();
    else if (execution_params.mode == "getset")
        getset();
    else if (execution_params.mode == "getctg")
        getctg();
    else if (execution_params.mode == "listset")
        listset();
    else if (execution_params.mode == "listctg")
        listctg();
    else if (execution_params.mode == "info")
        info();
    else
    {
        cerr << "Unknown mode: " << execution_params.mode << endl;
        return 0;
    }

    auto t2 = chrono::high_resolution_clock::now();

    if(execution_params.verbosity() > 0)
        cerr << "***\nCompleted in           : " << duration_cast<duration<double>>(t2 - t1).count() << " s" << endl;

    return 0;
}

// *******************************************************************************************
bool CApplication::create()
{
    CAGCCompressor agc_c;

    bool r = agc_c.Create(
        execution_params.out_archive_name,
        execution_params.pack_cardinality(),
        execution_params.k(),
        execution_params.input_names.front(),
        execution_params.segment_size(),
        execution_params.min_match_length(),
        execution_params.concatenated_genomes,
        execution_params.adaptive_compression,
        execution_params.verbosity(),
        execution_params.no_threads());

    if (!r)
    {
        cerr << "Cannot create archive " << execution_params.out_archive_name << endl;
        return false;
    }

    if (execution_params.verbosity() > 0)
        cerr << "Start of compression\n";

    vector<pair<string, string>> v_sample_file_names;

    for (auto& fn : execution_params.input_names)
    {
        string sample_name = std::filesystem::path(fn).stem().string();
        v_sample_file_names.emplace_back(sample_name, fn);
    }

    if(r)
        r &= execution_params.reproducibility_mode ?
            agc_c.AddSampleFilesRep(v_sample_file_names, execution_params.no_threads()) :
            agc_c.AddSampleFiles(v_sample_file_names, execution_params.no_threads());

    if (r && execution_params.store_cmd_line)
        agc_c.AddCmdLine(cmd_line);

    r &= agc_c.Close(execution_params.no_threads());

    return r;
}

// *******************************************************************************************
bool CApplication::append()
{
    CAGCCompressor agc_c;

    bool r = agc_c.Append(execution_params.in_archive_name, execution_params.out_archive_name, execution_params.verbosity(), true, execution_params.concatenated_genomes, execution_params.adaptive_compression,
        execution_params.no_threads());

    if (!r)
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << " or create archive " << execution_params.out_archive_name << endl;
        return false;
    }

    vector<pair<string, string>> v_sample_file_names;

    for (auto& fn : execution_params.input_names)
    {
        string sample_name = std::filesystem::path(fn).stem().string();
        v_sample_file_names.emplace_back(sample_name, fn);
    }

    if (execution_params.verbosity() > 0)
        cerr << "Start of compression\n";

    if(r)
        r &= execution_params.reproducibility_mode ?
            agc_c.AddSampleFilesRep(v_sample_file_names, execution_params.no_threads()) :
            agc_c.AddSampleFiles(v_sample_file_names, execution_params.no_threads());

    if (r && execution_params.store_cmd_line)
        agc_c.AddCmdLine(cmd_line);

    r &= agc_c.Close(execution_params.no_threads());

    return r;
}

// *******************************************************************************************
bool CApplication::getcol()
{
    CAGCDecompressor agc_d(true);

    bool r = agc_d.Open(execution_params.in_archive_name, execution_params.prefetch);
        
    if (!r)
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << endl;
        return false;
    }

    r &= agc_d.GetCollectionFiles(
        execution_params.output_name,
        execution_params.line_length(), 
        execution_params.no_threads());
        
    r &= agc_d.Close();

    return r;
}

// *******************************************************************************************
bool CApplication::getset()
{
    CAGCDecompressor agc_d(true);

    bool r = agc_d.Open(execution_params.in_archive_name, execution_params.prefetch);

    if (!r)
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << endl;
        return false;
    }

    r &= agc_d.GetSampleFile(
        execution_params.output_name,
        execution_params.sample_names, 
        execution_params.line_length(), 
        execution_params.no_threads());
        
    r &= agc_d.Close();

    return r;
}

// *******************************************************************************************
bool CApplication::getctg()
{
    CAGCDecompressor agc_d(true);

    bool r = agc_d.Open(execution_params.in_archive_name, execution_params.prefetch);

    if (!r)
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << endl;
        return false;
    }

    r &= agc_d.GetContigFile(
        execution_params.output_name, 
        execution_params.contig_names, 
        execution_params.line_length(), 
        execution_params.no_threads());

    r &= agc_d.Close();

    return r;
}

// *******************************************************************************************
bool CApplication::listset()
{
    CAGCDecompressor agc_d(true);

    bool r = agc_d.Open(execution_params.in_archive_name, execution_params.prefetch);

    if (!r)
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << endl;
        return false;
    }

    vector<string> v_samples;

    agc_d.ListSamples(v_samples);

    COutFile outf;
    r &= outf.Open(execution_params.output_name);

    if (!r)
    {
        cerr << "Cannot open output file " << execution_params.output_name << endl;
        return false;
    }

    for (auto& s : v_samples)
        outf.Write(s + "\n");

    outf.Close();

    r &= agc_d.Close();

    return r;
}

// *******************************************************************************************
bool CApplication::listctg()
{
    CAGCDecompressor agc_d(true);

    bool r = agc_d.Open(execution_params.in_archive_name, execution_params.prefetch);

    if (!r)
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << endl;
        return false;
    }

    COutFile outf;
    r &= outf.Open(execution_params.output_name);

    if (!r)
    {
        cerr << "Cannot open output file " << execution_params.output_name << endl;
        return false;
    }

    for (auto& sn : execution_params.sample_names)
    {
        outf.Write(sn + "\n");

        vector<string> v_contigs;

        agc_d.ListContigs(sn, v_contigs);

        for (auto& s : v_contigs)
            outf.Write("   " + s + "\n");
    }

    outf.Close();

    r &= agc_d.Close();

    return r;
}

// *******************************************************************************************
bool CApplication::info()
{
    CAGCDecompressor agc_d(true);

    if (!agc_d.Open(execution_params.in_archive_name, execution_params.prefetch))
    {
        cerr << "Cannot open archive " << execution_params.in_archive_name << endl;
        return false;
    }

    vector<string> v_sample_names;
    vector<pair<string, string>> cmd_lines;
    uint32_t kmer_length;
    uint32_t min_match_len;
    uint32_t pack_cardinality;

    agc_d.ListSamples(v_sample_names);
    agc_d.GetCmdLines(cmd_lines);
    agc_d.GetParams(kmer_length, min_match_len, pack_cardinality);

    cerr << "No. samples      : " << v_sample_names.size() << endl;
    cerr << "k-mer length     : " << kmer_length << endl;
    cerr << "Min. match length: " << min_match_len << endl;
    cerr << "Batch size       : " << pack_cardinality << endl;
    cerr << "Command lines:" << endl;

    for (auto& cmd : cmd_lines)
        cout << cmd.second << " : " << cmd.first << endl;

    if (execution_params.verbosity() > 0)
    {
        map<string, string> m_file_type_info;

        agc_d.GetFileTypeInfo(m_file_type_info);

        cerr << "File type info:\n";
        for (auto& x : m_file_type_info)
            cout << "  " << x.first << " : " << x.second << endl;
    }

    agc_d.Close();

    return true;
}

// *******************************************************************************************
int main(int argc, char** argv)
{
    CApplication app;

    return app.Run(argc, const_cast<const char **>(argv));
}

// EOF
