#### THIS IS WHAT NEEDS TO HAPPEN NEXT
### Step 2 Train the database for classification
#### training the database in Qiime

source activate qiime2-2019.10 

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path RDP.fa \
  --output-path RDP.qza