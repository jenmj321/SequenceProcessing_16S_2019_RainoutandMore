# classifying taxonomy of OTUs and ZOTUs using QIIME2 against RDP

I am using RDP 16s database, downloaded on January 14, 2020. 
In order to do the classification, I first needed to reformat the database, because the RDP database is formatted very differently from what the input for Qiime is. 

## training a dataset in Qiime2 from RDP. 
 
This is a bit of a problem because RDP is formatted WAY differently than the QIIME classification files. 
 fortunately I found this website!!!  https://vdreis.com/rdp-database/
 
IMPORTANT NOTES

This code does NOT work on a MAC. I don't understand sed commands in MAC but it's way different. 

You can run this on the hpcc! 

also I changed all of the quotation marks and 's , because they didn't convert right from the website. 

## Making Taxonomy file from RDP download

call out all the lines starting with “>”, this way you get the reference ID and the taxonomic assignment
```
grep -e ">" current_Bacteria_unaligned.fa>RDP_tax.txt 
```

place taxonomic assignment on a new line
```
sed -i 's/Lineage=/\n/' RDP_tax.txt 
```

keep only  the ID – there is some extra info that is not needed
```
sed -i 's/\(>\w*\d*\).*/\1/' RDP_tax.txt 
```

all ‘newlines’ are deleted, so one long line is made
```
tr -d '\n' < RDP_tax.txt > RDP_tax_new.txt 
```

create a tab between the ID and the assignment
```
sed -i 's/Root;rootrank;/\t/g' RDP_tax_new.txt 
```

every ID and assignment is placed onto a new line
```
sed -i 's/>/\n/g' RDP_tax_new.txt 
```

first line is deleted as it is blank
```
sed -i '1{/^$/d}' RDP_tax_new.txt 
```



## Making the Fasta file:

This takes a surprisingly long time! 

copy database to become the fasta file only
```
scp current_Bacteria_unaligned.fa RDP.fa 
```

then you need to remove taxonomy assignment and any extra info, and converting all sequences to uppercase.

These steps take a surprisingly long time (like 7 hours), so I submitted them as a job. 

```
nano converting_classification_fasta.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=10G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

sed -i 's/\(>\w*\d*\).*/\1/' RDP.fa 
sed -i 's/\(.*\)/\U\1/' RDP.fa 

###end of .sbatch
```
submitting the job 
```
sbatch converting_classification_fasta.sbatch
```
Submitted batch job 53441447

### Step 2 Train the database for classification
training the database in Qiime

```
source activate qiime2-2019.10 

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path RDP.fa \
  --output-path RDP.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path RDP_tax_new.txt \
  --output-path ref-taxonomy.qza
```

This step took 8 hours 
```
nano classifier_extract.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=10G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

source activate qiime2-2019.10 

qiime feature-classifier extract-reads \
  --i-sequences RDP.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-min-length 300 \
  --p-max-length 450 \
  --o-reads ref-seqs.qza

###end of .sbatch
```

submitting the job 
```
sbatch classifier_extract.sbatch
```
Submitted batch job 53473568


### Then you need to train the classifier.
```
nano classifier_train.sbatch
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=10G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   OTU_97  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jones514@msu.edu

########## Command Lines to Run ##########

source activate qiime2-2019.10 

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza


###end of .sbatch
```

submitting the job 
```
sbatch classifier_train.sbatch
```
Submitted batch job 53502517



