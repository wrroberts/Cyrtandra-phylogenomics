# Hawaiian _Cyrtandra_ phylogenomics
## Scripts and code for Kleinkopf et al. Diversification of Hawaiian _Cyrtandra_ under the influence of incomplete lineage sorting and hybridization.

## Programs used for analyses:
##### Java
##### Python
##### R
##### Julia
##### [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
##### [Geneious](https://www.geneious.com/)
##### [MAFFT](https://mafft.cbrc.jp/alignment/software/)
##### [Phyutility](https://github.com/blackrim/phyutility)
##### [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
##### [Astral](https://github.com/smirarab/ASTRAL)
##### [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl)
##### [PhyloNet](https://bioinfocs.rice.edu/PhyloNet)
##### [HYBRIDCHECK](https://github.com/BenJWard/HybridCheck)

### Trimming the raw fastq files
##### Reads were trimmed separately for each read file for each sample using paired-end mode in Trimmomatic v.0.36.
```java
java -jar Trimmomatic-0.36 PE reads.fq.gz forward_paired.fq.gz forward_unpaired.fq.gz reverse_paired.fq.gz reverse_unpaired.fq.gz HEADCROP:13 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:36
```
##### Paired reads for each sample were merged in Geneious 8 and mapped to the probe set as reference. A single reference was constructed that contained each locus separated by 200 N's. This approach was used so that a single mapping run per species would be performed rather than separate mapping runs per locus per sample. Consensus sequences were extracted for each locus for each sample. In cases of multiple alleles, default settings in Geneious were allowed for ambiguity calls. Regions where coverage was not >5X were called as "?" and treated as missing data.

### Alignment and gene tree estimation
##### Individual loci were aligned using MAFFT v.7.2.71.
```bash
mafft --auto *.fasta --thread 12 > *.aln_cln.fasta
```
##### Alignments were trimmed using Phyutility v.2.7.1 and minimal column occupancy of 0.5
```java
java -jar phyutility.jar clean 0.5 -in *.fasta -out *.fasta.trimmed
```
##### Gene trees were estimated for each locus using RAxML v.8.2.9
```bash
raxml -T 12 -p 12345 -x 12345 -f a -# 100 -s *.fasta.trimmed -n *.fasta.raxml_bs -m GTRCAT
```
##### Use SumTrees to summarize bootstrap results onto the best-scoring maximum likelihood tree and collapse nodes with <30% bootstrap support
```python
python collapse_nodes_below_threshold.py DIR *.tre 30 > *.fasta.collapsed
```

### Species tree estimations
##### Concatenation and maximum likelihood tree with RAxML
##### Concatenate loci using phyutility
```java
java -jar phyutility.jar concat -in *.fasta.trimmed ... -out trim50-concat.fa
```
##### Maximum likelihood estimation with RAxML using the GTRGAMMA model and partitioning the alignment
```bash
raxml -T 12 -p 12345 -x 12345 -f a -# 200 -s trim50-concat.fa -n trim50-concat -m GTRGAMMA -q trim50-concat.model
```
##### Summary-coalescent species tree using ASTRAL v.5.6.1
##### First, getting branch support with default local posterior probabilities
```java
java -jar astral.5.6.1.jar -i trim50-collapse30.geneTrees.tre -o trim50-collapse30.geneTrees.astral
```
##### Second, getting branch support with multilocus bootstrapping
```java
java -jar astral.5.6.1.jar -i trim50-collapse30.geneTrees.tre -o trim50-collapse30.geneTrees.astral -b trim50-collapse30.geneTrees.boot -r 200
```

### Phylogenetic network estimation
### These analyses were performed for 3 island groups (Kaua'i, Oahu, and Maui-Nui). Gene trees were pruned to include only the taxa found on each of the 3 island groups prior to these analyses. The same approach detailed below was used for each island group.
##### Using SNaQ from the PhyloNetworks Julia package, we use the same input trees as for Astral
```julia
using PhyloNetworks
using PhyloPlots

# read the gene trees into julia
raxmlCF = readTrees2CF("trim50-collapse30.geneTrees.tre", writeTab=false, writeSummary=false)

# use the astral tree as the starting topology
tre = readTopology("trim50-collapse30.geneTrees.astral-mlbs.tre")

# estimate the best network for hmax number of hybridization
net0 = snaq!(tre, raxmlCF, hmax=0, filename="net0_raxml", seed=12345)
net1 = snaq!(net0, raxmlCF, hmax=1, filename="net1_raxml", seed=12345)
net2 = snaq!(net1, raxmlCF, hmax=2, filename="net2_raxml", seed=12345)
net3 = snaq!(net2, raxmlCF, hmax=3, filename="net3_raxml", seed=12345)
net4 = snaq!(net3, raxmlCF, hmax=4, filename="net4_raxml", seed=12345)
```
##### Using PhyloNet v.3.6.4 and maximum pseudo-likelihood
##### The same approach as above was used and analyses were performed over each of the 3 island groups.
##### The following commands were used in each nexus file: 
##### InferNetwork_MPL geneTreeList numReticulations -s startingTopology -po -x 10 -pl 12
```java
# estimate the best network using from 0 to 4 allowed reticulations
java -jar PhyloNet_3.6.4.jar net0.nex
java -jar PhyloNet_3.6.4.jar net1.nex
java -jar PhyloNet_3.6.4.jar net2.nex
java -jar PhyloNet_3.6.4.jar net3.nex
java -jar PhyloNet_3.6.4.jar net4.nex
```

### ABBA-BABA tests for gene flow
##### The R package HYBRIDCHECK was used to calculate statistics in predefined 4-population phylogenies to test for gene flow.
##### Analyses were performed as below using concatenated alignments for each island group. Any columns with gaps were removed.
```R
require(HybridCheck)

Analysis1 <- HC$new('Kauai-concat.NOGAPS.fa')

populations1 <- list(P1=c('Clong'), P2=c('ClXk'), P3=c('CKaua'), Others=c('CaffW','Cwawr','Cwawr2'))

Analysis1$setPopulations(populations1)

Analysis1$prepareFourTaxonTests()

popCombos <- list(c(P1="P3", P2="P2", P3="P1", A="Others"))

Analysis1$prepareFourTaxonTests(popCombos)

Analysis1$runFourTaxonTests(selections='ALL', blockLength=1000L)

fttResults1 <- Analysis1$tabulateFourTaxonTests(selections='TESTED', neat=T, global=T)

split_list <- split(fttResults1, fttResults1$P1)
set1 <- as.data.frame(split_list[[1]])
set2 <- as.data.frame(split_list[[2]])
```
