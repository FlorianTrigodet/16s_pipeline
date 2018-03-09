# 16S rRNA pipeline
I used this pipeline for 16s analysis during my master and thesis

## Table of content

1. [IlluminaUtils](#preprocess)
  * Quality filtering (minoche et al.)
  * Merging
2. Swarm 
3. Chimera removal
4. Taxonomic assignation 
5. Complete table with OTUs and tax


**About raw data**
Raw data are in the `preprocess` directory. They are already demultiplexed 
and unziped.

## IlluminaUtils <a name="preprocess"></a>
IlluminaUtils python scripts by Meren were used to do the quality checking and 
paired end merging

```bash
cd preprocess

# Demultiplexed data are renamed according to their barcode/index
for line in `cat ../barcode_index.csv`
 do
 source=`echo $line | awk 'BEGIN{FS=";"}{print $1"_"$2}'`
 out=`echo $line | awk 'BEGIN{FS=";"}{print $3}'`
 cp `echo $source'_1_R1.fastq'` ./`echo $out'_R1.fastq'`
 cp `echo $source'_1_R2.fastq'` ./`echo $out'_R2.fastq'`
 done
```






