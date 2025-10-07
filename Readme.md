This workflow is for isolating RNA dark matter circles from a dataset that was provided by Rayan, in hopes that we can find new viruses or viroids kicking around in there. This document outlines the entire workflow up until the end of my rotation, so that future students can utilize the data and scripts for their own projects. This workflow starts with the 90% ORF clusters provided by Rayan, and ends with a set of sample origin annotated, reference genome and proteome filtered 35% identity ORF clusters. 


----- Step 1: Clustering -----

Clustering is based on the All vs All alignment function in usearch. A few things are needed for this:

1. Rayans 90% ORF clusters (can be found in the slack)
2. Artems copy of Usearch because the 32 bit version sucks
3. An EC2 instance with a ton of RAM

The resulting Uclust .UC file is the core of this workflow

==== Execution ====

```
usearch -cluster_fast orfs.cluster_size_over_10.fasta -id 0.35 -centroids circ_centroids.fasta -uc circ_clusters.uc
```
The resulting UC file had 11417413 instances, down from 14543248 (not that a good chunk of instances in the UC file are denoted by "C" and not a cluster member)

U clust output format: 

Field		Description
1		Record type S, H, C or N (see table below). C = cluster record
2		Cluster number (0-based).
3		Sequence length (S, N and H) or cluster size (C).
4		For H records, percent identity with target.
5		For H records, the strand: + or - for nucleotides, . for proteins.
6		Not used, parsers should ignore this field. Included for backwards compatibility.
7		Not used, parsers should ignore this field. Included for backwards compatibility.
8		Compressed alignment or the symbol '=' (equals sign). The = indicates that the query is 100% identical to the target sequence (field 10).
9		Label of query sequence (always present).
10		Label of target sequence (H records only).


----- Step 2: Sanity check -----

To check the distribution of obelisks and deltaviruses (which we know should be in the dataset), AA sequences are extracted from the centroids, which while not present in the UC file, can be done by 
extracting the SRA identifiers from column 9 of the UC file. Then the SRA identifiers can be run against Rayans 90% cluster file to pull the AA sequences. By turning the output FASTA file into a Diamond
DB, Diamond V2 can then be used to identify Obelisks and deltaviruses. 

What is needed:

1. The Usearch .UC file
2. Rayans 90% ORF clusters
3. an instance with Diamond V2
4. The RDRP1 RDRP5 data (this is the delta virus antigens)
5. The LOGAN obelisk search data (this is for oblin 1)

==== Execution ====

```
diamond blastp \
  --db cenDB.dmnd \
  --query Obelisk_DB_v1.1_Logan_Run_2_Oblin1.fa \
  --out Oblin1_results.m8 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --max-target-seqs 5 \
  --evalue 1e-5 \
  --threads 72
```
  plug and play with your target inputs...

===== Oblin Results =====

Total time = 6.253s
Reported 212678 pairwise alignments, 212678 HSPs.
60047 queries aligned.

===== DAG results =====

Total time = 5.777s
Reported 16007 pairwise alignments, 16007 HSPs.
9676 queries aligned

----- Step 3: Cluster Annotation (future users will likely start here) -----

This step is where the .UC file is annotated with both the species from which the sample is logged in the SRA as originating from. An AWK script is used to extract the SRA, KA, and contig number, then that data is 
compared against the SRA using a second script, and the awk output annotated with the species of origin. This annotated output is then concatenated with the .UC file to annotate the species of origin for each memeber of each cluster
(not just the centroids), and the resulting file split based on species parameters. 

What's needed:

1. sra_taxid.csv.zst file (see repo)
2. circle_contigs.fa.zst (or whatever Rayans original contig file is named)
3. AWK script to pull the SRA
4. AWK script to annotate the pulled SRA file
5. concatenation script (shotgun.sh - will need to be modified if you want to add more organisms)
6. Script to convert .UC to .tsv (polymorph.py)
7. Script to split the .tsv into designated bins based on species (Fly.grabber.R)


===== Execution of Fly.grabber.R on SRR.linked.triple.circles.tsv (annotated tsv) =====

Detected 13 columns and 11417413 rows
Processed 4145789 groups out of 4145789. 100% done. Time elapsed: 108s. ETA: 0s.
Output file summary (matches actual output files):
              output_file num_clusters total_sequences
              output file       clus num       seq num
1:          Saccharomyces        17076           56813
2: Caenorhabditis_elegans        15792           42257
3:             Drosophila        55754          121630
4:                no_hits      4057560         7073522

---- Step 4: Addition of contigs to .tsv file -----

This step adds the contigs of each SRA to the centroid of each cluster (trust me, trying to add the contig to each member made the file far too large). First, column 9, which contains the query, is split and the SRA extracted into a seperate file. 
The contigs are then pulled from the SRA. The resulting file is then merged with the Usearch output. 

What's needed:

1. gawk script to split the SRAs for your target file off (in this case, my saccharomyces file)
2. bash script to grep the contigs
3. python script to merge the files

==== Execution ====

==== gawk script ====

```
gawk -F'\t' 'NR==1 {
    print $0"\tSRR\tID\tKa"; next
}
{
    run=""; contig=""; ka="";
    # Extract SRR (SRR, DRR, ERR)
    if (match($9, /(SRR[0-9]+)/, m)) run=m[1];
    else if (match($9, /(DRR[0-9]+)/, m)) run=m[1];
    else if (match($9, /(ERR[0-9]+)/, m)) run=m[1];
    # Extract ID
    if (match($9, /_([0-9]+)_/, m)) ID=m[1];
    # Extract Ka
    if (match($9, /ka:f:([0-9.]+)/, m)) ka=m[1];
    print $0"\t"run"\t"ID"\t"ka;
}' output_clusters_Caenorhabditis_elegans_only.tsv > processed.C.elegans.circles.tsv
```
==== Contig extraction ====

```time zstd -dc tester.target.fa.zst | parallel --pipe --block 50M --jobs 20 grep -A1 -f tester.txt > extracted_sequences.fasta```

==== Merging with merger.py ====

```python3 merger.py c.elegans.X2.processed.circles.tsv c.elegans.circle.contigs.fa -o merged.output.tsv -v```

---- Step 5: Preparing files for filtering ----

Prior to filtering with Diamond2 and Bowtie2, the clusters for each species should be split into fasta files based on cluster size because you want to be able to trace an evolutionary history. Files also need to be checked because 
a malformed FASTA line will brick Diamond2. 

What's needed:

1. script to split species seggregated files based on cluster size (subsplitter.sh - this script also has a hit cutoff built in)
2. script to check the outputted fasta files.

==== Execution ====

==== subsplitter.sh ====

```
FILTER_SPECIES=Saccharomyces FILTER_MODE=count FILTER_CUTOFF=2 ./Subsplitter.sh Saccharomyces.X3.circles.w.contigs.tsv Saccharomyces.X3.0-3.fa 

 ^                                ^                      ^                                      ^                                ^ 
 
Filters based on the         Filters based on      Number of hits or                      Input                            Small output
header of the species        the number of hits    % of cluster that is
column. Defaults to          in the species        from the target species 
Saccharomyces unless         column. Defaults      for the filter modes
otherwise specified          to % unless 
                             specified


Saccharomyces.X3.3-5.fa     Saccharomyces.X3.5.fa
      ^                            ^  
      
   Medium output              Large clusters
```
==== subslitter.sh results ====

Filtered out 16833 clusters with <=2 saccharomyces hits
Splitting complete. Checking output files...
Small clusters (<3): 0 sequences
Medium clusters (3-5): 121 sequences
Large clusters (>5): 122 sequences

==== File checking ====

```
awk '/^>/ {
    if(header && !seq) print "Missing sequence for: " header
    header=$0; seq=""
} 
!/^>/ && NF > 0 {seq=seq $0} 
END {
    if(header && !seq) print "Missing sequence for: " header
}' Saccharomyces.X3.final.fasta
```

If a line is bad, cut it and move on. 

---- Step 6: Diamond and BT2 filtering ----

The target set of clusters from step 5 need to be filtered against the reference genome and proteome of the target creature to rule out endogenous RNA loops that can slip in. This is a pretty simple step and only involves running your data through Diamond,
then removing hits, then running the single filtered file through bowtie2 and removing hits from there to give a double negative list. This can all be done on the same machine. The cutoffs listed here are very aggressive, you may want to consider adjusting 
them on this step for something like human genetic material or just in general. 

What you need: 

1. EC2 instance with Diamond
2. script to remove Diamond hits (Diamond.breaker.py)
3. EC2 instance with bowtie 2
4. script to remove BT2 hits (Bowtie.breaker.py)

==== Execution ====

=== Diamond + breaker ====

```
diamond blastx -d yeastRP.dmnd -q Saccharomyces.X3.5.fa -o Saccharomyces.X3.5.dmnd.hits.m8 -f 6 --evalue 1e-5   --threads 64
```
Total time = 0.391s
Reported 350 pairwise alignments, 350 HSPs.
75 queries aligned.
```
python3 Diamond.breaker.py Saccharomyces.X3.5.dmnd.hits.m8 Saccharomyces.X3.5.fa Saccharomyces.X3.5.dmnd.null.fa
```
Filtering FASTA file...
Filtering complete:
  Sequences kept: 47
  Sequences removed: 75
  Total processed: 122

==== Bowtie2 + breaker ====

```
bowtie2 -p 8 -f -x reference_index -U Saccharomyces.X3.5.dmnd.null.fa -S Saccharomyces.X3.5.BT2.sam
```
47 reads; of these:
  47 (100.00%) were unpaired; of these:
    25 (53.19%) aligned 0 times
    15 (31.91%) aligned exactly 1 time
    7 (14.89%) aligned >1 times
46.81% overall alignment rate
```
Bowtie.breaker.py Saccharomyces.X3.5.BT2.sam Saccharomyces.X3.5.dmnd.null.fa Saccharomyces.X3.5.double.null.fa
```
Filtering FASTA file...
Filtering complete:
  Sequences kept: 25
  Sequences removed: 22
  Total processed: 47

---- Step 7: Searching in circles (Scanning the remaining loops for virus and viroid signatures) ----

Unfortunately, my search in yeast did not turn up any new viruses... 

However, there is now a FASTA file full of RNA circles with no clear origin. Whether they are viroids, viruses, transposons or just junk is completely unknown. To discern what they are, they are tested with multiple different methods. INFERNAL is a neat 
little program that can find ribozymes based on covariance models (CMs) - it uses runs inputted FASTAs against designated CMs to identify likely ribozymes and outputs a nice table. Publicly 
available CMs are not great, but Marcos de la Pena has been nice enough to gift his CMs for this, which are very powerful (Marcos.cm is the name of the CM file). Additionally, programs such as BLASTn can be used against the virus genome archive (for which a file can be found on this github), 
and Diamond can be deployed against a DB of viral RDRPs and other proteins (unfortunately I did not have time to do that, but the concept is the same as prior diamond runs). 

