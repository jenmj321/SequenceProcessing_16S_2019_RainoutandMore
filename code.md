# EvansLab_16S 2019 Sequence Processing for Jennifer, Tayler, Lukas, and Heather

Jennifer ran this code on the MSU HPCC between January 12 and January 18, 2020. 

Lukas's GLBRC code was used for the analysis. His protocol can be found [here.](https://github.com/GLBRC-TeamMicrobiome/TagAnalysis/blob/master/16S_USEARCH_Pipeline.md)

This Pipeline uses USEARCH through OTU and ZOTU mapping and then classification of OTUs and ZOTUs through silva (using USEARCH) and RDP (using QIIME2).

This pipeline assumes that you know how to login to the HPCC and use nano and that you have all of the necessary programs installed. 


## Quality Checking

Make sure you are in the folder with the extracted forward and reverse reads from Run 1

```
cd ./20191213_16S-V4_PE250
```

Files are stored as .gz so you need to unzip them. 

```
gunzip *.gz
```

### First look at the quality of raw unmerged seqs for run1

Make new directory for your quality files.
```
mkdir fastq_info_run1
```

Write the code to generate the quality files
```
nano fasta_info_fq.sh
!#/bin/bash
for fq in *.fastq
do
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_info $fq -output fastq_info_run1/$fq
done
```

make file executable and run the for loop create in your fasta_info_fq.sh file

```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```
move to the fastq_info directory
```
cd fastq_info_run1/
```

Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > Run1_fastq_lengths.txt
grep "^EE" * > Run1_fastq_EE.txt
```
look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

I saved these as csv files so that you can search them. 

### Look at the quality of raw unmerged seqs for run2 
Make sure you are in the folder with the extracted forward and reverse reads from Run 2
```
cd ./20200107_16S-V4_PE250/
```

Files are stored as .gz so you need to unzip them. 
```
gunzip *.gz
```

First look at the quality of raw unmerged seqs for run2
```
mkdir fastq_info_run2
```

Write the code to generate the quality files
```
nano fasta_info_fq.sh
!#/bin/bash
for fq in *.fastq
do
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_info $fq -output fastq_info_run2/$fq
done
```

make file executable and run the for loop create in your fasta_info_fq.sh file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```
move to the fastq_info directory
```
cd fastq_info_run2/
```
Now run the below code to summarize the fastq info for all of the forward and reverse reads
```
grep "^File" * > Run2_fastq_lengths.txt
grep "^EE" * > Run2_fastq_EE.txt
```
Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

I saved these files as .csv files


## Merge the forward and reverse sequences and trim adapters (for each run individually)

### Run 1

#### 2a) Merge Pairs for Run 1
Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html

-alnout gives you a human readable text file of the alignment and misalignments for each pair merged.

-tabbedout give you extensive information on the quality of the merge

This step takes approximately 1 minute
```
mkdir mergedfastq_run1
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout mergedfastq_run1/combined_merged_run1.fastq  -tabbedout mergedfastq_run1/combined_merged_run1_pair_report.txt -alnout mergedfastq_run1/combined_merged_run1_pair_aln.txt
```

If you are getting a low percentage of sequence pair merging (e.g. below 60%), consider trimming the reverse reads See bottom of script (“improving mergepairs”) for script on truncating reverse reads.
Tips on poor merging: https://drive5.com/usearch/manual/merge_badrev.html

For this run, pairs merged consistently around 70%. Here' the output summary. 


| Totals: |  |
| ------ | ------ |
| 11742936 | Pairs (11.7M) |
| 8287072 | Merged (8.3M, 70.57%) |
| 846082  | Alignments with zero diffs (7.21%) |
| 3428283 | Too many diffs (> 5) (29.19%) |
| 0  | Fwd tails Q <= 2 trimmed (0.00%) |
| 7 | Rev tails Q <= 2 trimmed (0.00%) |
| 27581 | No alignment found (0.23%) |
| 0 | Alignment too short (< 16) (0.00%) |
| 47251 | Staggered pairs (0.40%) merged & trimmed |
| 246.74 | Mean alignment length |
| 253.02 | Mean merged length |
| 0.26 | Mean fwd expected errors |
| 1.24 | Mean rev expected errors |
| 0.11 | Mean merged expected errors |

#### 2b) let's check sequence quality of the merged seqs using USEARCH's fastq_eestats2.

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 mergedfastq_run1/combined_merged_run1.fastq -output fastq_info_run1/combined_merged_run1_eestats2.txt
```
00:00 38Mb   CPU has 28 cores, defaulting to 10 threads

01:00 661Mb   100.0% Reading reads

8287072 reads, max len 479, avg 253.0

| Length | MaxEE 0.50 | MaxEE 1.00 | MaxEE 2.00 |
| ------  | ---------------- |  ----------------  | ---------------- |
|    50  |  8225396( 99.3%) |   8285708(100.0%) |   8286956(100.0%) |
|  100  |  8175961( 98.7%)  |  8282440( 99.9%)  |  8286880(100.0%) |
|   150  |  8141182( 98.2%)  |  8277908( 99.9%)  |  8286700(100.0%) |
 |  200  |  8110746( 97.9%)  |  8273450( 99.8%)  |  8285027(100.0%) |
|   250  |  8023131( 96.8%)  |  8231035( 99.3%)  |  8250200( 99.6%) |
|   300  |     5537(  0.1%)  |     6467(  0.1%)   |    6687(  0.1%) |
|   350    |   738(  0.0%)  |    1211(  0.0%)   |    1504(  0.0%) |
|   400   |     115(  0.0%)    |    324(  0.0%)   |     565(  0.0%) |
|   450     |     1(  0.0%)   |       6(  0.0%)     |     8(  0.0%) |


#### 2c) Optional Step: Remove any residual bases from adapter seqs using cut adapt
 I did not do this step. It wasn't working and Lukas said it was okay. 


### Run 2
#### 2a) Merge Pairs for Run 2
Make sure you are in the folder with the extracted forward and reverse reads from Run 2
https://www.drive5.com/usearch/manual/merge_options.html

-alnout gives you a human readable text file of the alignment and misalignments for each pair merged.

-tabbedout give you extensive information on the quality of the merge

This step takes approximately 5 minutes

```
mkdir mergedfastq_run2
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout mergedfastq_run2/combined_merged_run2.fastq  -tabbedout mergedfastq_run2/combined_merged_run2_pair_report.txt -alnout mergedfastq_run2/combined_merged_run2_pair_aln.txt
```

If you are getting a low percentage of sequence pair merging (e.g. below 60%), consider trimming the reverse reads See bottom of script (“improving mergepairs”) for script on truncating reverse reads.
Tips on poor merging: https://drive5.com/usearch/manual/merge_badrev.html

Pairs merged consistently around 70% 

| Totals: |  |
| --- | --- | 
|  10274096 | Pairs (10.3M) |
 |  7197875 | Merged (7.2M, 70.06%) |
 |  3086351 | Alignments with zero diffs (30.04%) |
|   3051209 | Too many diffs (> 5) (29.70%) |
|     25012 | No alignment found (0.24%) |
|         0 | Alignment too short (< 16) (0.00%) |
|     62342 | Staggered pairs (0.61%) merged & trimmed |
|    246.51 | Mean alignment length |
|    253.13 | Mean merged length |
|      0.36 | Mean fwd expected errors |
|      1.11 | Mean rev expected errors |
|      0.07 | Mean merged expected errors |


### 2b) let's check sequence quality of the merged seqs using USEARCH's fastq_eestats2.
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 mergedfastq_run2/combined_merged_run2.fastq -output fastq_info_run2/combined_merged_run2_eestats2.txt
```

00:00 38Mb   CPU has 28 cores, defaulting to 10 threads

01:13 661Mb   100.0% Reading reads

7197875 reads, max len 480, avg 253.1

| Length  |       MaxEE 0.50  |       MaxEE 1.00   |      MaxEE 2.00 |
| ------  | ---------------- |  ---------------- |  ---------------- |
|    50  |  7166676( 99.6%)  |  7196917(100.0%)  |  7197757(100.0%) |
 |  100  |  7110003( 98.8%) |   7192506( 99.9%) |   7197727(100.0%) |
|   150 |  7069774( 98.2%)  |  7185367( 99.8%)  |  7197645(100.0%) |
|   200  |  7042142( 97.8%)  |  7180142( 99.8%)  |  7196145(100.0%) |
|   250 |   6942755( 96.5%)  |  7120778( 98.9%)  |  7147706( 99.3%) |
|   300  |    15562(  0.2%)   |   17491(  0.2%)   |   18053(  0.3%) |
|   350   |    7180(  0.1%)   |    9014(  0.1%)   |    9805(  0.1%) |
 |  400   |     578(  0.0%)  |     1099(  0.0%)  |    1604(  0.0%) |
|   450   |       8(  0.0%)   |      19(  0.0%)   |     33(  0.0%) |

#### 2c) Optional Step: Remove any residual bases from adapter seqs using cut adapt
 I did not do this step. It wasn't working and Lukas said it was okay. 



## 3) Combine the two merged sequence files (or as many as you have Illumina runs)
Make sure you do not use replicate sample names between the two runs
This step is only necessary if more than one Illumina run are being analyzed together.

I made a new directory and copied the merged fastq files for each run. 
```
mkdir combined_runs_2020_01_13
cp ./20191213_16S-V4_PE250/mergedfastq_run1/combined_merged_run1.fastq ./combined_runs_2020_01_13/combined_merged_run1.fastq 
cp ./20200107_16S-V4_PE250/mergedfastq_run2/combined_merged_run2.fastq ./combined_runs_2020_01_13/combined_merged_run2.fastq 
```

concatenate the two separate fastq files from each run. 
```
cat combined_merged_run1.fastq combined_merged_run2.fastq > combined_merged_both_runs.fastq
```

check that the combined_merged_both_runs.fastq file has the same number of lines as the other two files combined

```
wc -l combined_merged_both_runs.fastq
```
61939788
```
wc -l combined_merged_run1.fastq
```
33148288
```
wc -l combined_merged_run2.fastq
```
28791500
they are the same number of lines! 

Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step. Any hyphens or underscores can be problematic and you need to remove these (use sed command and merged_cut files)

Search the file for any hyphens or underscores
```
grep -c "-" combined_merged_both_runs.fastq
```
116973
```
grep -c "_" combined_merged_both_runs.fastq
```


this prints out all of the lines where @ starts a line and - is in the same line. This will print a lot of sequence quality info but it's easy enough to sift through to double check that there are no sequence names printed which would mean that there's a - in a name. 
```
grep '^@.*-' combined_merged_both_runs.fastq
```

Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs using fastx_get_sample_names
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_get_sample_names combined_merged_both_runs.fastq -output combined_merged_both_runs_samples.txt
```


## 4) Filtering and Truncate the merged seqs to MaxEE and set length using fastq_filter
250 bp is the expected overlaps with 515F and 806R

-fastq_trunclen 250 trims the sequences to 250 bp

-fastq_maxee 1 make maxEE score 1 
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter combined_merged_both_runs.fastq -fastq_maxee 1 -fastq_trunclen 250 -fastaout combined_merged_both_runs_fil.fa
```

0:00 4.8Mb  FASTQ base 33 for file combined_merged_both_runs.fastq

00:00 39Mb   CPU has 28 cores, defaulting to 10 threads

02:16 628Mb   100.0% Filtering, 99.1% passed

15484947  Reads (15.5M)                   

86716  Discarded reads length < 250

46418  Discarded reads with expected errs > 1.00

15351813  Filtered reads (15.4M, 99.1%)


## 5) Filter so we only have unique sequences with fastx_uniques
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques combined_merged_both_runs_fil.fa  -fastaout uniques_combined_merged_both_runs_fil.fa -sizeout
```

00:42 5.0Gb   100.0% Reading combined_merged_both_runs_fil.fa

00:42 5.0Gb  CPU has 28 cores, defaulting to 10 threads      

00:56 8.2Gb   100.0% DF

00:57 8.3Gb  15351813 seqs, 3379201 uniques, 2282232 singletons (67.5%)

00:57 8.3Gb  Min size 1, median 1, max 148875, avg 4.54

01:32 5.9Gb   100.0% Writing uniques_combined_merged_both_runs_fil.fa


## 6) Cluster and map into OTUS and filter out singletons
There are two options here. (A) uses the traditional approach and clusters sequences into 0.97 identity cutoff OTUs. (B) uses unoise3 to identify ZOTUs.

I did both so you all can decide what you want to use! 

### Option 1: 97% OTUs 

#### Cluster into 0.97 OTUs using UPARSE and cluster_otus

This step will also denovo chimera check and filter out singletons.You can remove single sequences prior to clustering but singletons are also removed at the OTU clustering step (cluster_otus filters out OTUs <2 and unoise3 filters ZOTUs <8)

The clustering took a long time so I submitted it as a job. 

```
nano otu_cluster_97.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=05:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=10G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=EMAIL@msu.edu

########## Command Lines to Run ##########

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_otus uniques_combined_merged_both_runs_fil.fa -otus combined_merged_both_runs_otus.fa -uparseout combined_merged_both_runs_otus_uparse.txt -relabel OTU

###end of .sbatch
```

submitting the job 
```
sbatch otu_cluster_97.sbatch
```

check that the number of unique sequences are the same minus the chimeras 
```
grep -c '^>' uniques_combined_merged_both_runs_fil.fa
```
3379201
```
wc -l combined_merged_both_runs_otus_uparse.txt
```
1096969
It checks out! 

#### Mapping reads to traditional 0.97 OTUS
This took 3 hrs and 5 minutes to run

```
nano otu_mapping_97.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=05:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=10G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab combined_merged_both_runs.fastq -otus combined_merged_both_runs_otus.fa -uc combined_merged_both_runs_OTU_map.uc -otutabout combined_merged_both_runs_OTU_table.txt -biomout combined_merged_both_runs_OTU_jsn.biom -notmatchedfq combined_merged_run2_otu_unmapped.fq

###end of .sbatch
```

submitting the job 
```
sbatch otu_mapping_97.sbatch
```

### Option 2

#### Identify ZOTUs using unoise3

This step will also denovo chimera check and filter out low abundance ZOTUs. IMPORTANT: ZOTUs with less than 3 reads will be filtered out (i.e. -minsize 3) default setting filters out ZOTUs less than 8 reads
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -unoise3 uniques_combined_merged_both_runs_fil.fa -zotus combined_merged_both_runs_zotus.fa  -tabbedout combined_merged_both_runs_zotus_report.txt -minsize 3
```

I ran this as job because it took so long. 
```
nano otu_cluster_zotu.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -unoise3 uniques_combined_merged_both_runs_fil.fa -zotus combined_merged_both_runs_zotus.fa  -tabbedout combined_merged_both_runs_zotus_report.txt -minsize 3

###end of .sbatch
```

submitting the job 
```
sbatch otu_cluster_zotu.sbatch
```
Submitted batch job 53380671
this took 9 hours to run


#### Map ZOTUs back to samples
```
nano otu_cluster_step2_zotu.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

sed -i 's/Zotu/ZOTU/g' combined_merged_both_runs_zotus.fa

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab combined_merged_both_runs.fastq -zotus combined_merged_both_runs_zotus.fa -uc combined_merged_both_runs_ZOTU_map.uc -otutabout combined_merged_both_runs_ZOTU_table.txt -biomout combined_merged_both_runs_ZOTU_jsn.biom -notmatchedfq combined_merged_run2_ZOTU_unmapped.fq

###end of .sbatch
```

submitting the job 
```
sbatch otu_cluster_step2_zotu.sbatch
```

Submitted batch job 53445677

14947970 / 15484947 mapped to OTUs (96.5%)  


## Classifying taxonomy of OTUs and ZOTUs using QIIME2 against RDP and SILVA

### Using RDP database in QIIME2 

first, You need to capitalize the codons before you run the classifier
 
OTU rep set 
```
tr '[:lower:]' '[:upper:]'  < combined_merged_both_runs_otus.fa > combined_merged_both_runs_otus_CAP.fa
```

ZOU rep set
```
tr '[:lower:]' '[:upper:]'  < combined_merged_both_runs_zotus.fa > combined_merged_both_runs_zotus_CAP.fa
```
 
 
##### Classify the OTUs
``` 
 nano classify_OTU.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=50:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=30G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   classify_OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

source activate qiime2-2019.10 

#OTU
qiime tools import --type 'FeatureData[Sequence]' --input-path combined_merged_both_runs_otus_CAP.fa --output-path combined_merged_both_runs_otus.qza 

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads combined_merged_both_runs_otus.qza --o-classification combined_merged_both_runs_otus_RDP_taxonomy.qza 


###end of .sbatch
```

submitting the job 
```
sbatch classify_OTU.sbatch
```
Submitted batch job 53504191


##### Classify the ZOTUs
```
nano classify_ZOTU.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=30G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   classify_OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

source activate qiime2-2019.10

#OTU
qiime tools import --type 'FeatureData[Sequence]' --input-path combined_merged_both_runs_zotus_CAP.fa --output-path combined_merged_both_runs_zotus.qza

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads combined_merged_both_runs_zotus.qza --o-classification combined_merged_both_runs_zotus_RDP_taxonomy.qza 

###end of .sbatch
```
submitting the job 
```
sbatch classify_ZOTU.sbatch
```
Submitted batch job 53541077

I haven't figured out how to export the QIIME2 files yet. 


### Using Silva database and USEARCH

##### Classifying OTUs
```
nano classify_silva_OTU.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=30:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   taxa_class_OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax combined_merged_both_runs_otus.fa -db silva_16s_v123.fa -tabbedout combined_merged_both_runs_otus_silva_taxonomy.sintax -strand both

###end of .sbatch
```
submitting job
```
sbatch classify_silva_OTU.sbatch
```
Submitted batch job 53594826

##### Classifying OTUs
```
nano classify_silva_ZOTU.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=30:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   taxa_class_OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax combined_merged_both_runs_zotus.fa -db silva_16s_v123.fa -tabbedout combined_merged_both_runs_zotus_silva_taxonomy.sintax -strand both

###end of .sbatch
```
Submitting job to HPCC
```
sbatch classify_silva_ZOTU.sbatch
```
Submitted batch job 53619571



