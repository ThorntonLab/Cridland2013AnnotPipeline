This README discuss the scripts needed to annotate TEs detected with the TE detection pipeline.

These scripts require the following output files 
1) fasta files generated by the detection pipeline have been assembled with Phrap (-forcelevel 10
-minscore 10 -minmatch 10)
2) these contigs have been aligned back to the reference genome with BLAST () #NOTE - the way the annotation script recognizes filenames these output files need to have the same name as the contig file with the addition of a .chromosome.output at the end of the file name 
3) these contigs have been aligned to a database containing all of the reference TE sequence only

annotate_TEs.pl

The input for this script is as follows

1) the directory where the contig files and the blast output files are located - these need to be in the same directory
2) the name of the lefthand contig file - this is the file that is created for upstream of the breakpoint - the script will parse this filename and then use it to reconstruct the file names of the right contig file and the blast output files associated with these
3) the estimate positions file generated by the TE detection pipeline - the line*_te_pos_estimate file - the pipeline checks both that there is a TE insertion and that the contig reconstructed aligns to the region of the genome estimated
4) the chrom_num file - this is the file with the numeric conversion for chromosome arms
5) the output file

The output file of this script has the following fields - this file can generate duplicate TE identifications for TEs shared with the reference - this is often the case when there are several small TEs very near each other in the reference and multiple contigs are created in this region - therefore this file needs to be edited for redundancy after it is generated

1) chromosome
2) upstream estimate
3) downstream estimate
4) a string with an id number, then a :, then a string with a preliminary description of what TE is located there - this was a first attempt at TE annotation, but was not what was used in the final paper
5) the word shared or novel to designate a TE shared with the refernce sequence (same TE at same location) or not
6) if novel - a number - this is for the upstream contig
   0 = complete reconstruction - one contig that contains the TE breakpoint
   1 = 2 contigs - one that is TE and one that aligns to the estimated genomic region, but do not have the exact breakpoint
   2 = did not resolve a contig
7) if novel - same as 6 except for the downstream contig

trim_repetition.pl

this script takes the output of the annotate_TEs.pl script and removes repetitive entries

the input for this script is

1) the output of annotate_TEs.pl
2) the output file - same format as input

###In the paper the data is then restricted to certain chromosome regions - this was done at this point in the process #########

make_best_te.pl 

This script determines the best match for each TE for each contig based on length of match and e-value

the input for this script is

1) the blast output file for the alignment to the TE database (this processes one file - I write a perl wrapper to go through each file in a directory)
2) the output file - this file is set up to append so you can run many instances and write to a single output file

the output file for this script is
    1) contig (contig name field can be parsed to chromosome number . ID number . left/right . fasta . Contig#) 
    2) TE
    3) length of TE match

make_best_te_contig.pl

This script determines which contig has the resolved TE breakpoing - this is only applicable for situations where a breakpoint is fully reconstructed

the input for this file is

1) the blast output file for the alignment to the reference sequence
2) the output file
3) the output file from annotate_TEs.pl following removal of redundant records

the output of this file is 
    1) which contig contains the breakpoint
    2) the total length of this contig

With the output files from make_best_te.pl and make_best_te_contig.pl the user can determine which TE has been annotated - In the paper we compare which TE is annotated with the upstream and downstream estimate for each event (left/right)

evaluate_shared.pl

This script identifies presence/absence calls for TEs shared with the reference

the input for this script is

1) the blast output file from the alignment back to the reference genome
2) the output file
3) the te_pos_estimate file
4) the file with the positions of the TEs in the reference 
5) distance for absence breakpoints 

