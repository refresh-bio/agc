#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <random>
#include <agc-api.h>

void single_thread(const char *fn)
{
	// *** Open agc file
	CAGCFile agc;
	
	if(!agc.Open(fn, false))
	{
		std::cerr << "Cannot open " << fn << " file" << std::endl;
		return;
	}
	
	// *** Check no. of samples 
	std::cout << "No. samples: " << agc.NSample() << std::endl;

	// *** List samples
	std::vector<std::string> samples;
	int n_sample = agc.ListSample(samples);
	
	if(n_sample < 0)
	{
		std::cerr << "Error: cannot read list of samples" << std::endl;
		return;
	}
	
	std::cout << "Samples in file with no. of contigs" << std::endl;
	for(const auto &s : samples)
		std::cout << s << " : " << agc.NCtg(s) << std::endl;

	// *** List contigs in 0th sample, together with 
	std::cout << "\nContents of sample: " << samples.front() << std::endl;
	
	std::vector<std::string> contigs;
	int n_ctg = agc. ListCtg(samples.front(), contigs);
	
	if(n_ctg < 0)
	{
		std::cerr << "Error: cannot read list of contigs" << std::endl;
		return;
	}

	for(const auto &ctg : contigs)
		std::cout << ctg << " : " << agc.GetCtgLen(samples.front(), ctg) << std::endl;;

	// *** Print part of 0th contig of 0th sample
	int part_len = 2000;
	int from = 1000000;
	int to = from + part_len - 1;
	
	std::string seq;
	
	int ctg_len = agc.GetCtgLen(samples.front(), contigs.front());
	
	if(ctg_len <= to)
	{
		to = ctg_len - 1;
		if(ctg_len < part_len)
			from = 0;
		else
			from = to - (part_len - 1);
	}
	
	int seq_len = agc.GetCtgSeq(samples.front(), contigs.front(), from, to, seq);
	
	std::cout << from << " " << to << " " << seq_len << " " << seq << std::endl;
	
	// Query without sample name
	seq_len = agc.GetCtgSeq("", contigs.front(), from, to, seq);
	
	std::cout << from << " " << to << " " << seq_len << " " << seq << std::endl;	
	
	agc.Close();
}

void many_threads(const char *fn)
{
	CAGCFile agc;
	
	if(!agc.Open(fn, true))		// in case of many queries it is better to open agc file in prefetching mode
	{
		std::cerr << "Cannot open " << fn << " file" << std::endl;
		return;
	}

	std::vector<std::thread> v_thr;
	auto n_thr = std::thread::hardware_concurrency();
	
	std::cout << "Running " + std::to_string(n_thr) + " threads\n";
	
	// Multithreaded reading of contigs
	for(size_t i = 0; i < n_thr; ++i)
		v_thr.emplace_back([&agc, i]{
			std::vector<std::string> samples, contigs;
					
			agc.ListSample(samples);
			int n_sample = samples.size();

			std::mt19937 mt(i);
			int my_sample_id = mt() % n_sample;
				
			std::string sample_name = samples[my_sample_id];
			
			agc.ListCtg(sample_name, contigs);
			int n_contig = contigs.size();
		
			for(int i = 0; i < 10; ++i)
			{
				uint64_t sum_ctg = 0;
				
				std::string seq;
				std::string ctg_name = contigs[mt() % n_contig];
				agc.GetCtgSeq(sample_name, ctg_name, -1, -1, seq); 
				
				// Do anything - here we calculate sum of symbols
				for(auto c : seq)
					sum_ctg += static_cast<uint64_t>(c);
				
				std::cout << "Thread " + std::to_string(i) + " : " + sample_name + " " + ctg_name + " " + std::to_string(sum_ctg) + "\n";
			}
				
		});

	for(auto &t : v_thr)
		t.join();

	agc.Close();	
}

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		printf("Usage: example-agc-lib-c <agc_file>\n");
		return 0;
	}
	
	single_thread(argv[1]);
	
	many_threads(argv[1]);
	
	return 0;
}