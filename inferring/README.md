# FACS-seq method to profile mutant library of regulatory element: NGS data processing package

## What is this?
This python script collection is one of the two the software subpackages of FACS-seq, used for the NGS data processing and calculation of response of each sensor mutant. The paper describing this program is currently under review. Please cite this paper if this program is useful to your work.

## General description of the algorithm and experiment
To understand the sequence-function relationship of regulatory genetic elements, we developed FACS-seq approach to profile the funcitonal consequences of thousands of mutants in parallel. Briefly, a synthetic library consisted of many mutants was prepared, whose response was perturbed by relevant mutations. Next, the library was treated with a variaty of ligand concentration and subjected to FACS. Here, the library was sorted into diffrent bins according to low to high response signals. The cells sorted into each bin were recovered, respectively and subjected to NGS analysis of the DNA payload carrying the mutations. This package is used in the final step to process the NGS data to infer the response of each mutant for each condition (ligand concentration). For details about the mathematical framework, see our paper. Basically, users only need to edit a configure file to pinpoint the corresponding bin of each NGS data, as well as some parameters used in bin setting during FACS experiment. The package can follow the instructions in the configure file to automatically process the NGS data, giving rise to the final output. [Here](./image/FACS-seq_schematics.png) is a schematic of FACS-seq experiment.


## How to use it?
### Step 1：Installation and dependency
1. Install Python version 2.7
2. Install Scipy version 0.19.1 or above
3. Install Matplotlib version 2.0.2 or above
4. Install Numpy version 1.13.1 or above
5. Install Pandas version 0.18.1 or above



### Step 2：Prepare the necessary files.
All these files (or subdirectories) should be organized under a common working directory together with the all .py scripts.
Please check the example files post at GitHub, which are described as below.


#### File 1: NGS files (.fastq or .fq extension) under one directory (see example_data/)
Note: Please try to keep the name of each file meaningful but as simple as possible. The files can also be compressed as .gz format.


#### File 2: mutant library file (see example_library.csv)
The mutant library file specifies the synthetic mutant libary used in FACS-seq experiment. This file is used to map the NGS read to each mutant for counting. It is at .csv formate **without header line**, in which there are two columns in order of id and sequence, respectively. **Use tab as delimiter**. It should be noted that **-, _ and ' '(space) should be eliminated from any id name. Avoid id like 'super-001', 'super_001' or 'super 001'**.

|1st column (id)|2nd column(nucleotide sequence)|
|---------------|-------------------------------|
|mutant1|ATGAATATCTTACATATATGTGTGACCTCAAAATGGTTCAATATTGACAACAAAATTGTCGATCACCGCCCTTGA|
|mutant2|ATGAAAATCTTACATATATGTGTGACCTCAAAATGGTTCAATATTGACAACAAAATTGTCGATCACCGCCCTTGA|
|.......|....................|


#### File 3: experiment design file (see example_experiment_configure.txt)
This file specifies the experiment design of FACS-seq. In particular, it defines the relation between each NGS raw data and its role in FACS-seq experiment. Usually, to profile the response of a genetic regulatory element of interest, we performed FACS-seq under M different conditions (M ligand concentrations). In each condition, the library is subjected to FACS and sorted into N bins. Hence, in the experiment design file, a N(row)-M(column) matrix is presented to define the role of each NGS raw data (one bin under one condition). **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=0uM|Ligand=100uM|Ligand=500uM|
|---|----------|------------|------------|
|P1|B0P1|B1P1|B2P1|
|P2|B0P2|B1P2|B2P2|
|P3|B0P3|B1P3|B2P3|
|P4|B0P4|B1P4|B2P4|
|P5|B0P5|B1P5|B2P5|
|P6|B0P6|B1P6|B2P6|
|..|....|....|....|

**The name of each NGS raw data should be the same as those defined in the configure file (see below)**


#### File 4: cell count file (see example_cell_count_configure.txt)
This file specifies the number of cells sorted into each bin during FACS-seq. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=0uM|Ligand=100uM|Ligand=500uM|
|---|----------|------------|------------|
|P1|10587|50564|50575|
|P2|100000|100000|100000|
|P3|100000|100000|100000|
|P4|100000|100000|100000|
|P5|100000|100000|100000|
|P6|51472|50923|51607|
|..|.....|.....|.....|


#### File 5: bin boundary file (see example_bin_boundary_configure.txt)
This file specifies the boundary of bins used in FACS. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=0uM|Ligand=100uM|Ligand=500uM|
|---|----------|------------|------------|
|P1|0.722,10.000|0.722,10.000|0.722,10.000|
|P2|0.316,0.722|0.316,0.722|0.316,0.722|
|P3|-0.076,0.316|-0.076,0.316|-0.076,0.316|
|P4|-0.460,-0.076|-0.460,-0.076|-0.460,-0.076|
|P5|-0.861,-0.460|-0.861,-0.460|-0.861,-0.460|
|P6|-10.000,-0.861|-10.000,-0.861|-10.000,-0.861|
|..|..............|..............|..............|

**In each cell, comma is used to seperate the upper and lower boundary of the relevant bin. Usually, fluorescence signal is used in FACS to define bins; hence, a log10(FLU) is expected here.**

Note that in the example file, we used log10(GFP/mCherry) to define bins, which is the reason here for negative bin boundary value. GFP is the reporter under control of the studied regulatory element, while mCherry is constitutively expressed used to normalize cell-to-cell variability. -10.000 and 10.000 are used as numerical lower and upper boundary of negative and positive infinite, respectively. Usually, the gating strategy is determined prior to the cell sorting, hence identical gating parameters for all conditions (sorting experiment) are expected. Even though, note that different gating parameters are also supported. 


#### File 6: bin occupation file (see example_bin_occupation_configure.txt)
This file specifies the proportion of cells sorted into each bin in the whole pooled library. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=0uM|Ligand=100uM|Ligand=500uM|
|---|----------|------------|------------|
|P1|0.01649|0.04676|0.05377|
|P2|0.20371|0.35796|0.39821|
|P3|0.28313|0.19516|0.17723|
|P4|0.23385|0.18492|0.16726|
|P5|0.18695|0.14078|0.12048|
|P6|0.07587|0.07443|0.08305|
|..|.......|.......|.......|

**Note that parameters here are different from those in file 4!!!** One additional experiment (analysis cytometry using mutant library) is needed to obtain the parameters here. Hence, the sum of each column here should be exactly 1. 



### Step 3: Set up the configure file (see example_configure.txt)
The configure file is used to set all the necessary parameters and tell the program where to find necessary files. 
**This file is in a two-column format using tab as delimiter.** Each line starts with one word (name of one parameter) separated with the following (content of this parameter) by a tab delimiter. We describe each parameter as below.
 
**fastq**: data directory and all NGS raw data (.fastq or .fq extension). The first item of the content should be the data directory, while the following items refer to the raw NGS data under this directory.

**sample-label**: the label for each NGS raw data file. The order of the label should corresponds to the order of the file names specified by the 'fastq' parameter. For simplicity, it is fine to use raw data file name without extension as label. Note that the labels specified here should be the same as the library name defined in the experiment design file (see above, Step 2, File 3).

**iniLib**: the NGS data profiling the mutant library prior to cell sorting, used as a reference to eliminate those mutants that are underrepresented in the library.

**prefix_nucl**: several (4-10) upstream nucleotides flanking the variable region, used to specify and cut the variable region from the sequencing read. These nucleotides should be located in the the PCR product of NGS library.

**suffix_nucl**: several (4-10) downstream nucleotides flanking the variable region, used to specify and cut the variable region from the sequencing read. These nucleotides should be located in the the PCR product of NGS library.

**variable_region_len**: number of nucleotides of the variable region on the construct, where the synthetic mutations are introduced. It is determined in the library design. Hence, the length specified here should be consistent with that of the mutant library file (see above, Step 2, File 2).

**list-seq**: the name of the mutant library file (see above, Step 2, File 2).

**experiment_configure**: the name of the experiment design file (see above, Step 2, File 3).

**cell_count_configure**: the name of the cell count file (see above, Step 2, File 4)

**bin_boundary_configure**: the name of the bin boundary file (see above, Step 2, File 5)

**bin_occupation_configure**: the name of the bin occupation file (see above, Step 2, File 6)

**ReadsThreshold**: threshold of read count used to eliminate mutants detected less than this threshold in the library prior to cell sorting. Default: 20.

**log10u_sigma_range**: range of log10u (average response, see our paper for details) and sigma, within which the search for optimal value is executed. **format: " log10u_lower_limit,log10u_upper_limit, sigma_lower_limit,sigma_upper_limit", e.g." -1.5,0.9,0.05,0.4"**. Note that a double quotes is expected and the first character should be blank. We suggest using experimental data to determine the limit for log10u and sigma search range.

**search**: density of 2D grid plane used to search for optimal log10u and sigma value. 

**prefix**: prefix used for naming of all output files, keep it simple without any ‘-’, ‘_’ and ‘ ’. For example, ‘screen20171001’ is fine.

see **example_configure.txt**.

parameter|value
---------|-----
prefix|Ilovemicrobe
fastqpath|example_data
fastq|example_data|Lib1.fq.gz|B0P1|B0P2
forward_prefixseq|GCAC
forward_suffixseq|GTTT
sample-label|dCas9R1,dCas9R2,NCR1,NCR2,plasmid
sgrna-len|20
list-seq|example_library.csv
experiment_configure|example_experiment_configure.txt
name_configure|example_naming_configure.txt
control_setting|NC
FDR_threshold|0.05
ReadsThreshold|20
hit_gene_calling|position
gene_sgRNA_position|example_coding_region_position.txt
Operon_gene_List|example_operon.txt

After Step 2 and 3, check your working directory. It should looks like below:
[here](./image/files_prepared_before_data_processing.png)

### Step 4：Run the pipeline
Open the command line window (for example, terminal in Macbook), cd to the working directory and run the analysis pipeline.
cd path_to_your_working_directory
python CRISPRscreen_main.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 
python CRISPRscreen_main.py example_configure.txt

Check [here](./image/successful_running.png) for the output during a successful running of the abovementioned test.

The program will create an 'error.log' file under the working directory, open this file to check whether anything wrong happens. Generally, no content suggests successful running. Please post your 'error.log' if you cannot figure out the bugs when using this tool.

For a typical Macbook (for example, 2.6 GHz processor and 8 GB memory), the example test can be finalized within 30 minutes. The rate-limiting step is the mapping of the raw NGS data. For a typical Macbook, we expect a processing speed of 20 million reads per hour. For a genome-scale library with 50 k members (10 sgRNAs per gene assuming 5 k genes encoded by a genome), 100-fold coverage leads to 5 million reads per library, thus roughly 4 NGS libraries processed per hour.

## Output files
The output files will be organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory (prefiex_results). We term this subdirectory 'result directory' thereafter.

[your working directory should be like this after running the test](./image/wkd_after_example_running.png)

You can find many sub directories under the result directory.

[your result directory after running the test](./image/resultdir_after_example_running.png)

Below is the description. For the mathematical processing, see our paper. **All .csv flat files use tab as delimiter unless mentioned**

### NGS raw data profile
-------------------------------------------------------------
#### read count of each sgRNA in each library (prefix_count/)
**prefix.countsummary.txt**: basic statistics of the mapping ratio of each NGS library with a header line using tab as delimiter.

File|Label|Reads|Mapped|Synerror|Unknown|Percentage|Zerocounts|GiniIndex
----|-----|-----|------|--------|-------|----------|----------|---------
example_data/plasmid.fq.gz|plasmid|1000000|838484|90356|71160|0.8385|1734|0.2167
...|...|...|...|...|...|...|...|...

Reads denote the number of reads in the raw data. Mapped denotes number of reads mapping perfectly to one member of the synthetic sgRNA library. Synerror refers to those reads with one indel mutation or more mismatch mutations. Unknown refers to those reads where no forward_prefixseq or forward_suffixseq can be identified. Percentage is the mapping ratio. Zerocount refers to sgRNA number in the *in silico* library without any corresponding read detected. GiniIndex is a metric reflecting the member abundance uniformity in a library. Bigger Gini index indicates more biased distribution with over- represented or diluted members. Generally, more stringent the selective condition is, bigger Gini index we can expect. 

**prefix.count.txt**: raw read count for each sgRNA in the *in silico* library **before normalization**.

sgRNA|dCas9R1|dCas9R2|NCR1|NCR2|plasmid
-----|-------|-------|----|----|-------
gspKb3332_817|12|11|11|9|8
...|...|...|...|...|...

**prefix.normalizeCount.txt**: read count for each sgRNA in the *in silico* library **after normalization of sequencing depth** (for details, see our paper). **This dataset is used for following data processing.**

sgRNA|Gene|dCas9R1|dCas9R2|NCR1|NCR2|plasmid
-----|----|-------|-------|----|----|-------
gspKb3332_817|gspK|12.07|11.08|10.95|8.95|7.98
...|...|...|...|...|...|...

[**prefix_Libray_Gini_Score.png**](./image/all_Libray_Gini_Score.png): schematic of Gini index for each library.

### sgRNA level statistics
-------------------------------------------------------------
#### removed sgRNAs (removed.sgRNA/)
We remove the over diluted sgRNAs with read count less than one threshold ('ReadsThreshold' described in the configure file part).

 1. **prefix.removed.sgRNA.txt**: a simple list flat file with one sgRNA each line

============================================================
#### biological replicate agreement (replicate_consistence/)
Files under this dierectory is a minotoring panel for biological replicate agreement. The replicate information is encoded by the experiment design file (Step 2, File 5). Generally, for N experiments with 2 replicates each, the program produces N scatter plots and N flat files to describe the consistence between replicates for each experiment. One summarizing flat file for all conditions is also given.

 1. **prefix_replicates_reads_statistics.txt**: the summary file
 
 condition|pearson correlation coefficient|P value
 ---------|-------------------------------|-------
 stress1|0.868|0.0
 control1|0.832|0.0
 ...|...|...

 2. **prefix_oneexperiment_replicates.txt**: two-column flat file of read count of each sgRNA in two replicates. For example, NCR1 and NCR2 (in our test data) are two replicates for one experiment.
 
 sgRNA|NCR1_abundance|NCR1_reads|NCR1_abundance_vs_initial|NCR2_abundance|NCR2_reads|NCR2_abundance_vs_initial
 -----|--------------|----------|-------------------------|--------------|----------|-------------------------
 gspKb3332_817|-16.22|10.95|0.46|-16.51|8.95|0.16
 ...|...|...|...|...|...|...

 3. [**prefix_oneexperiment_replicates.png**](./image/all_control1_replicates.png): schematic of replicate agreement of one particular experiment.

============================================================
#### sgRNA read count, abundance change, fitness score, etc (prefix_sgRNA_statistics/)
This directory stores all dataset about sgRNA metrics. It is generally organized at three levels (three sub directories):
 1. **Library level** (information of one NGS library). N files, N = number of rows in experiment design file.
 
 sgRNA|gene|plasmid_Log2_abundnace|plasmid_reads|plasmid_Log2_abundnace_vs_initial
 -----|----|----------------------|-------------|---------------------------------
 gspKb3332_817|gspK|-16.68|7.98|0.0
 ...|...|...|...|...|

 2. **Condition level** (information of one experiement (average of two replicate NGS library)). N files, N = number of columns in experiment design file.
 
 sgRNA|gene|control1_Log2_abundnace|control1_reads|control1_Log2_abundnace_vs_initial|control1_relative_deviation
 -----|----|----------------------|---------------|----------------------------------|---------------------------
 gspKb3332_817|gspK|-16.37|9.90|0.31|0.10
 ...|...|...|...|...|...

 3. **Phenotype level (only need to focus on this level for simplicity)** (information of one phenotype (selective condition normalized by the control condition), under combined_condition_level directory). N files, N = number of 'stress'(selective) conditions in experiment design file.
 
 relative_abundnace_change = Log<sub>2</sub> (read count selective condition / read count control condition)
 
 normalized_change (**it is used as sgRNA fitness score**) = relative_abundnace_change - median relative_abundnace_change of NC sgRNAs
 
 Zscore = normalized_change / sigma of NC sgRNA normalized_change normal distribution
 
 Quality: it is tagged as 'Good' if the averaged read count in control condition is above the threshold ('ReadsThreshold' described in the configure file). **Only 'Good' sgRNAs are used in gene level calculation**.
  
 sgRNA|gene|relative_abundnace_change|normalized_change|Zscore|Quality
 -----|----|-------------------------|-----------------|------|-------
 gspKb3332_817|gspK|0.22|-0.03|-0.04|Good
 ...|...|...|...|...

============================================================
#### NC sgRNA distribution (NCsgRNA_ND/)
Theoretically, fitness socre (log2 abundance change) of NC sgRNA should follow a normal distribution. We hence use a normal distribution to fit NC sgRNA fitness score data.

 1. **prefix_NCsgRNA_ND.txt**: normal distribution of NC sgRNA relative abundance changes (before normalization by median of NC sgRNA relative abundance change, referring to 'relative_abundnace_change' column in phenotype level sgRNA statistics)
 
 condition|median|mean|stdev
 ---------|------|----|-----
 phenotype1|0.25|0.16|0.71
 phenotype2|0.25|0.16|0.71
 ...|...|...|...

 2. **prefix_NCsgRNA_normalized_ND.txt**: normal distribution of NC sgRNA fitness scores (after normalization by median of NC sgRNA relative abundance change, referring to 'normalized_change' column in phenotype level sgRNA statistics)
 
 condition|median|mean|stdev
 ---------|------|----|-----
 phenotype1|0.0|-0.09|0.71
 phenotype2|0.0|-0.09|0.71
 ...|...|...|...

 3. [**prefix_phenotype_NCsgRNAND.png**](./image/all_essential_NCsgRNAND.png): schematic of NC sgRNA fitness score distribution.


### gene level statistics: 
-------------------------------------------------------------
#### FPR-score curve (prefix_quasigeneFPR/)
We use a NC sgRNA derived 'quasi' gene simulation approach (score approach thereafter, we use this method in our paper) (score = |gene fitness| * -Log<sub>10</sub>Pvalue_MWUtest) to determine the false positive rate (*FPR*) for each gene-phenotype association. Hence, for each studied phenotype, the program give 15 simulated *FPR*-score curves with 1 ~ 15 sgRNAs per quasi gene, respectively. Thus, 15 files describing these curves and [one figure file](./image/all_essential_quasigeneFPR.png) are in this sub directory.

============================================================
#### P value-Q value curve (prefix_Pvalue_Qvalue/)
We use a Storey-Tibshirani approach (PNAS 2003) to convert *FPR* into *Q* values. We also use a simple student t test method (sgRNA for one gene vs NC sgRNAs) (t test approach thereafter) to calculate another *P* value and convert it into *Q* values.

 1. **Qvalue_scoreFPR.txt**: score approach derived *P-Q* value curve for all phenotypes.
 
 phenotype1|phenotype2|...|Qvalue
 ----------|----------|---|------
 0.00015|...|...|0.001
 0.00084|...|...|0.005
 ...|...|...|...
 0.03180|...|...|0.1
 
  2. **Qvalue_TtestPvalue.txt**: t test approach derived *P-Q* value curve for all phenotypes, similar to above.
  
  3. Distributions of *FPR* (or t test *P*) for all genes (2N [figures](./image/all_essential_Ttest_pValue.png), N = number of phenotypes).
  
  4. Comparison of *P-Q* curves from two approaches (N [figures](./image/Ilovemicrobe_example_Pvalue_Qvalue.png), N = number of phenotypes).
 
============================================================
#### gene fitness, statistical significance, etc (prefix_gene_statistics/)
This directory stores all dataset about gene metrics. N files are generated, corresponding to N studied phenotypes.

sgRNAnumber: number of sgRNAs used to produce the metrics for this gene (see 'position' approach in our paper)

MedianRAC: gene fintess

MedianZ: gene fitness normalized by the sigma of NC sgRNA fitness normal distribution

-Log10Pvalue_MWUtest: MWU test of sgRNA (of one gene) fitness vs. NC sgRNA fitness

FPRvalue: score approach derived *FPR*; FDRvalue: *Q* value derived from *FPR*;

-Log10Pvalue_Ttest: student t test of sgRNA (of one gene) fitness vs. NC sgRNA fitness

Qvalue_Ttest: *Q* value derived from t test *P* value

gene|sgRNAnumber|MedianRAC|MedianZ|negative Log10Pvalue_MWUtest|FDRvalue|FPRvalue|negative Log10Pvalue_Ttest|Qvalue_Ttest
----|-----------|---------|-------|--------------------|--------|--------|------------------|------------
gspK|11|-0.19|-0.27|0.76|0.30|0.18|0.74|0.36
...|...|...|...|...|...|...|...|...

============================================================
#### FDR-score curve (prefix_quasigeneFDR/)
A [figure](./image/all_essential_quasigeneFDR.png) and a file describing FDR vs. score relation considering sgRNA number per gene profile in the library using a 'quasi' gene simulation approach. It is used to produce the volcano plot.

 phenotype1|phenotype2|...|FDR
 ----------|----------|---|------
 21.91|...|...|0.001
 4.36|...|...|0.005
 ...|...|...|...
 0.72|...|...|0.1

### operon level statistics (reorganize gene data according to the operon file)
-------------------------------------------------------------
#### reorganize gene fitness and statistical significance according to operon structures (prefix_operon_statistics/)
3N files corresponding to N phenotypes. 

**prefix_phenotype_RemoveGene_statistics.txt**: flat file describing removed genes due to the availability of sgRNAs

**prefix_phenotype_RemoveOperon_statistics.txt**: flat file describing removed operons due to the availability of genes

**prefix_phenotype_operon_statistics.txt**: reorganized gene fitness and statistical signifcance dataset. Each line refers to one operon. Different types of metrics are seperated by tab. The same type of metrics corresponding to multiple genes in one operon are seperated by comma. For simplicity, only part of gene metrics are included here.

one phenotype

gene|sgRNAnumber|MedianRAC|MedianZ|negative Log10Pvalue_MWUtest|FDRvalue|FPRvalue
----|-----------|---------|-------|--------------------|--------|--------|
mtlD,mtlR|8,7|-0.110,-0.248|-0.154,-0.347|0.510,0.724|0.472,0.300|0.410,0.184
...|...|...|...|...|...|...|
