#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <agc-api.h>

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		printf("Usage: example-agc-lib-c <agc_file>\n");
		return 0;
	}
	
	// *** Open agc file
	agc_t *agc = agc_open(argv[1], 0);
	if(!agc)
	{
		printf("Cannot open %s file\n", argv[1]);
		return 0;
	}
	
	// *** Check no. of samples 
	printf("No. samples: %d\n", agc_n_sample(agc));

	// *** List samples
	int n_sample;
	
	char **list_samples = agc_list_sample(agc, &n_sample);
	
	printf("Samples in file with no. of contigs\n");
	for(char **p = list_samples; *p; ++p)
		printf("%s : %d\n", *p, agc_n_ctg(agc, *p));

	// *** List contigs in 0th sample, together with 
	printf("\nContents of sample: %s\n", list_samples[0]);
	
	int n_ctg;
	char **list_contigs = agc_list_ctg(agc, list_samples[0], &n_ctg);
	
	for(char **p = list_contigs; *p; ++p)
		printf("%s : %d\n", *p, agc_get_ctg_len(agc, list_samples[0], *p));

	// *** Print part of 0th contig of 0th sample
	int part_len = 2000;
	int from = 1000000;
	int to = from + part_len - 1;
	
	char *seq = (char*) malloc(part_len + 1);
	
	int ctg_len = agc_get_ctg_len(agc, list_samples[0], list_contigs[0]);
	
	if(ctg_len <= to)
	{
		to = ctg_len - 1;
		if(ctg_len < part_len)
			from = 0;
		else
			from = to - (part_len - 1);
	}
	
	int seq_len = agc_get_ctg_seq(agc, list_samples[0], list_contigs[0], from, to, seq);
	
	printf("%d %d %d : %s\n", from, to, seq_len, seq);

	// Query without sample name
	seq_len = agc_get_ctg_seq(agc, NULL, list_contigs[0], from, to, seq);
	
	printf("%d %d %d : %s\n", from, to, seq_len, seq);

	free(seq);
	
	agc_list_destroy(list_samples);
	agc_list_destroy(list_contigs);
	
	agc_close(agc);
	
	return 0;
}