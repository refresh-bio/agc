#!/usr/bin/env python3
import sys
import py_agc_api as agc
import textwrap

        
# open AGC archive ("toy_ex/toy_ex.agc")
agc_arch = agc.CAGCFile()
if not agc_arch.Open("toy_ex/toy_ex.agc", True):
    print("Error: cannot open agc archive")
    sys.exit(1)

#Get number of samples in the archive
n = agc_arch.NSample();
print("No. samples: ", n)

#Get reference name
reference = agc_arch.GetReferenceSample()
print("Reference sample: ", reference)

#Get list of samples in the archive
set_list = agc.StringVector()


agc_arch.ListSample(set_list);
print("Samples in file with no. of contigs")
    
iterator = iter(set_list)
for i in range(n): #for i in range(len(set_list)):
    sample = set_list[i];
    no_ctg = agc_arch.NCtg(sample) #Get number of contigs for sample
    print(sample, ":", no_ctg)

    
#Get contents of the sample
sample=set_list[0]
print("\nContents of sample:", sample)
#Get number of contigs for the 0th sample
no_ctg = agc_arch.NCtg(sample)
#Get name of contigs of the 0th sample
ctg_list = agc.StringVector()
agc_arch.ListCtg(sample, ctg_list);
#Print name and length of each contig of the 0th sample
for i in range(no_ctg):
    ctg_len = agc_arch.GetCtgLen(sample, ctg_list[i]) #Get length of the of the contig
    print("length of contig", ctg_list[i],":",ctg_len) #print contig name and length
#Print part of contig in the sample (of length 5 if possible, from position 8, if possible)
    start=8 #start position of contig part
    length=5  #part length
    end=start+length-1  #end position of contig part
    ctg_len = agc_arch.GetCtgLen(sample, ctg_list[i]) #Get length of the of the contig
    if end >= ctg_len:
        end = ctg_len - 1
        if length > ctg_len:
            start = 0 #print whole contig
        else:
            start = end - (length - 1)
    print("\tsubsequence start:", start, ", end:", end, ", length:", length)
    #Get and print part of 0th contig in the 0th sample, from start to end
    seq = agc_arch.GetCtgSeq(sample, ctg_list[i], start, end)
    print("\tsubsequence length:", len(seq))
    print("\tpos:", start, "-", end, "len:", end-start+1, "seq:",seq);
    #Get (query without sample name) and print part of 0th contig in the 0th sample, from start to end (work only if contig name is unique among samples)
    seq = agc_arch.GetCtgSeq(ctg_list[i], start, end)
    print("\tsubsequence start:", start, ", end:", end, ", length:", length)
    print("\tsubsequence length:", len(seq))
    print("\tpos:", start, "-", end, "len:", end-start+1, "seq:",seq);
