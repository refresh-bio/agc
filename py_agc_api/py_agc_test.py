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

#Get list of samples in the archive
set_list = agc.StringVector()
agc_arch.ListSample(set_list);

print("Samples in file with no. of contigs")
for s in set_list:
    no_ctg = agc_arch.NCtg(s) #Get number of contigs for sample
    print(s, ":", no_ctg)

#Get contents of the first sample
print("\nContents of sample:", set_list[0])
#Get number of contigs for the 0th sample
no_ctg = agc_arch.NCtg(set_list[0])
#Get name of contigs of the 0th sample
ctg_list = agc.StringVector()
agc_arch.ListCtg(set_list[0], ctg_list);

#Print name and length of each contig of the 0th sample
for i in range(no_ctg):
    ctg_len = agc_arch.GetCtgLen(set_list[0], ctg_list[i]) #Get length of the of the contig
    print(ctg_list[i],":",ctg_len) #print contig name and length
    #Print part of contig in the sample
    start=8 #start position of contig part
    len=5  #part length
    end=start+len-1  #end position of contig part
    ctg_len = agc_arch.GetCtgLen(set_list[0], ctg_list[i]) #Get length of the of the contig
    if end >= ctg_len:
        end = ctg_len - 1
        if len > ctg_len:
            start = 0 #print whole contig
        else:
            start = end - (len - 1)
            #Get and print part of 0th contig in the 0th sample, from start to end
    seq = agc_arch.GetCtgSeq(set_list[0], ctg_list[i], start, end)
    print("\tpos:", start, "-", end, "len:", end-start+1, "seq:",seq);
    #Get (query without sample name) and print part of 0th contig in the 0th sample, from start to end (work only if contig name is unique among samples)
    seq = agc_arch.GetCtgSeq(ctg_list[i], start, end)
    print("\tpos:", start, "-", end, "len:", end-start+1, "seq:",seq);

