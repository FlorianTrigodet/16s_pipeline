# 16S rRNA pipeline
I used this pipeline for 16s analysis during my master and thesis


## Table of content

1. [Data description](#data)
2. [IlluminaUtils](#preprocess)
   * [Quality filtering (minoche et al.)](#qual)
   * [Merging](#merging)
3. [Swarm](#swarm)
4. [Chimera removal](#chimera)
5. [Taxonomic assignation](#taxa)
6. [Complete table with OTUs and tax](#otu_table)

## Data description <a name="data"></a>
We are interested in the bacterial communbity associated with the ennoblement of stainless steel (an spontaneous increase of the electrochemical potential). We previously observed a critical temperature set between 30°C and 40°C for the seawater in Brest (France) above which the ennoblement is no longer observed. Therfore we set up an experiment in which we exposed stainless steel to a range of temperature and collect the bacteria that have settled on the surface of the metal.

 
The temperature tested were 30, 33, 36, 38 and 40°C, with five replicates per conditions. 
Only three seawater tanks were available at the time of the expirement, so only three temperature tested simultanously, but each time with the 30°C conditions to ensure reproducibility
 

**About raw data**
Raw data are in the `preprocess` directory. They are already demultiplexed 
and unziped.

If the demultiplexed fastq were named according to the barcode and index. Samples were renamed using the barcode_index.csv and the following bash lines

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


## IlluminaUtils <a name="preprocess"></a>
IlluminaUtils python scripts by Meren were used to do the quality checking and paired end merging

Config files must be created for quality checking and merging

```bash
cd preprocess

# Create the first required file: qual-config.txt
ls *.fastq| awk 'BEGIN{FS="_R"}{print $1}' | uniq | awk 'BEGIN{print "sample\tr1\tr2"}{print $0 "\t" $0 "_R1.fastq\t" $0"_R2.fastq"}' > qual-config.txt


# Create the second required file: merge-config.txt 
ls *.fastq| awk 'BEGIN{FS="_R"}{print $1}' | uniq | awk 'BEGIN{print "sample\tr1\tr2"}{print $0 "\t" $0 "-QUALITY_PASSED_R1.fastq\t" $0"-QUALITY_PASSED_R2.fastq"}' > merge-config.txt
```

### Quality filtering (minoche et al.) <a name="qual"></a>
The quality filtering is made with IlluminaUtils using [Minoche _et al_](http://genomebiology.com/2011/12/11/R112) recommanded filtering parameters. 

```bash
cd preprocess

# Generate .ini files with barcode (.....) and primer associated for both R1 and R2.
# Exact match will be kept, and barcode + primers will be trimmed during filtering.
iu-gen-configs qual-config.txt --r1-prefix ^.....CCAGCAGC[C,T]GCGGTAA. --r2-prefix CCGTC[A,T]ATT[C,T].TTT[A,G]A.T

# loop for quality filtering
for i in *.ini
 do
 iu-filter-quality-minoche $i
 done
```

Output name: `*-QUALITY_PASSED_R1.fastq` and `*-QUALITY_PASSED_R2.fastq` 


### Merging <a name="merging"></a>

The merging function is similar to the quality filtering

```bash
iu-gen-configs merge-config.txt --r1-prefix ^.....CCAGCAGC[C,T]GCGGTAA. --r2-prefix CCGTC[A,T]ATT[C,T].TTT[A,G]A.T

# loop for quality filtering
for i in *.ini
 do
 iu-merge-pairs $i
 done
```

Once the quality filtering and merging done, I move the new fasta files into the `fasta_dir` and rename to only the sample name

```bash
for i in `ls preprocess/ | grep MERGED | sed 's/\(^.*\)_MERGED/\1/'`; do cp preprocess/$i"_MERGED" fasta_dir/$i".fasta"; done
```


## Swarm
OTUs clustering is performed by the swarm algorythm developped by [Mahé](https://peerj.com/articles/1420/). 
There is no need for a thereshold using this algorythm. 

All sequences are now in a fasta file per sample but the letters are sometime uppercase or lowercase. To avoid any confusion for the upcoming dereplication and clustering, I prefer to change all sequences to lowercase:

```bash
for i in fasta_dir/*.fasta; do awk '{if(!/>/){print tolower($0)}else{print $0}}' $i > $i.temp | mv $i.temp $i ; done
```

### Dereplication in each sample
Then we can dereplicate all sequences per sample. The idea is to have a fasta file with unique sequences and their frequences in the defline. 

```bash
for i in `ls fasta_dir/ | grep .fasta | sed 's/\(^.*\).fasta/\1/'`; do vsearch --derep_fulllength fasta_dir/$i.fasta --sizeout --relabel_sha1 --fasta_width 0 --output swarm/dereplicate/$i.derep.fasta ; done
```

Vsearch counting format is `>Seq;size=#;` and it needs to be changed to `>Seq_#` for swarm

```bash
cd swarm/dereplicate/
sed -i 's/size=/_/' *.derep.fasta
sed -i 's/;//g' *.derep.fasta
```

At this point we can create a contingency table of the unique sequence per sample. The resulting table will have all unique sequences as rows and samples as columns and filled with the sequence abundance. 
Special thanks to [Frederic Mahe](https://github.com/torognes/swarm/wiki/Working-with-several-samples) for that code.

```bash
cd swarm/dereplicate/

awk 'BEGIN {FS = "[>_]"}
           # Parse the sample files
           /^>/ {contingency[$2][FILENAME] = $3
           amplicons[$2] += $3
           if (FNR == 1) {
               samples[++i] = FILENAME
           }
          }
    END {# Create table header
          printf "amplicon"
          s = length(samples)
          for (i = 1; i <= s; i++) {
              printf "\t%s", samples[i]
          }
          printf "\t%s\n", "total"

          # Sort amplicons by decreasing total abundance (use a coprocess)
          command = "LC_ALL=C sort -k1,1nr -k2,2d"
          for (amplicon in amplicons) {
               printf "%d\t%s\n", amplicons[amplicon], amplicon |& command
          }
          close(command, "to")
          FS = "\t"
          while ((command |& getline) > 0) {
              amplicons_sorted[++j] = $2
          }
          close(command)

          # Print the amplicon occurrences in the different samples
          n = length(amplicons_sorted)
          for (i = 1; i <= n; i++) {
               amplicon = amplicons_sorted[i]
               printf "%s", amplicon
               for (j = 1; j <= s; j++) {
                   printf "\t%d", contingency[amplicon][samples[j]]
               }
               printf "\t%d\n", amplicons[amplicon]
          }}' *derep.fasta > ../amplicon_contingency_table.csv
```

### Dereplication for all sample
A dereplication at the sudy level is now required before the OTU clustering. It provides a unique fasta, with unique sequence and their abundance in the defline.


```bash
cd swarm/dereplicate/

export LC_ALL=C
cat *derep.fasta | \
awk 'BEGIN {RS = ">" ; FS = "[_\n]"}
     {if (NR != 1) {abundances[$1] += $2 ; sequences[$1] = $3}}
     END {for (amplicon in sequences) {
         print ">" amplicon "_" abundances[amplicon] "_" sequences[amplicon]}}' | \
sort --temporary-directory=$(pwd) -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > Temperature_effect.fasta
```

### Swarm
Swarm was used with default value -d 1, meaning only the difference of 1 base pair is taken to compare sequences. Other parameters include -t 8 the number of processors used, -s -w -l and -o the various outputs and -f for ...

```bash
cd swarm/

swarm -d 1 -f -t 8 -s Temperature_effect.stat -w OTUs_rep.fasta -o Temperature_effect.swarm -l Temperature_effect.log dereplicate/Temperature_effect.fasta
```

## Chimera removal <a name="chimera"></a>
In this case, the chimera detection is done on the OTUs representative sequence only. According to Mahé, Swarm should create OTUs out of Chimera sequences rather than having them included in an OTU. A lot of processor time can be spared when the chimera detection is not done on the global fasta file. 

The abundance format for each sequence needs to get from `>Seq_#` to `>Seq;size=#;` for Vsearch

```bash
cd swarm/

sed -i 's/_/;size=/' OTUs_rep.fasta
sed -i 's/^>.*/&;/' OTUs_rep.fasta
```

The chimera detection is done with Vsearch.

```bash
cd swarm/

vsearch --alignwidth 0 --uchime_denovo OTUs_rep.fasta --uchimeout Temperature_effect.uchimeout.txt
```


 






