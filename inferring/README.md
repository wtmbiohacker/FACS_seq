# FACS-seq method to profile mutant library of regulatory element: NGS data processing package

## What is this?
This python script collection is one of the two the software subpackages of FACS-seq, used for the NGS data processing and calculation of response of each sensor mutant.

## General description of the algorithm and experiment
To understand the sequence-function relationship of regulatory genetic elements, we developed FACS-seq approach to profile the funcitonal consequences of thousands of mutants in parallel. Briefly, a synthetic library consisted of many mutants was prepared, whose response was perturbed by relevant mutations. Next, the library was treated with a variaty of ligand concentration and subjected to FACS. Here, the library was sorted into diffrent bins according to low to high response signals. The cells sorted into each bin were recovered, respectively and subjected to NGS analysis of the DNA payload carrying the mutations. This package is used in the final step to process the NGS data to infer the response of each mutant for each condition (ligand concentration). For details about the mathematical framework, see our paper. Basically, users only need to edit a configure file to pinpoint the corresponding bin of each NGS data, as well as some parameters used in bin setting during FACS experiment. The package can follow the instructions in the configure file to automatically process the NGS data, giving rise to the final output. [Here](./image/FACS-seq_schematics.png) is a schematic of FACS-seq experiment.


## How to use it?
MacOS or Linux is supported environment to run this tool; for Windows users, some system shell commands in the scripts ('os.system(...)') need to be revised.

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

|Bin|Ligand=100uM|Ligand=500uM|
|---|------------|------------|
|P1|B1P1|B2P1|
|P2|B1P2|B2P2|
|P3|B1P3|B2P3|
|P4|B1P4|B2P4|
|P5|B1P5|B2P5|
|P6|B1P6|B2P6|

**The name of each NGS raw data should be the same as those defined in the configure file (see below)**


#### File 4: cell count file (see example_cell_count_configure.txt)
This file specifies the number of cells sorted into each bin during FACS-seq. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|Ligand=500uM|
|---|------------|------------|
|P1|50564|50575|
|P2|100000|100000|
|P3|100000|100000|
|P4|100000|100000|
|P5|100000|100000|
|P6|50923|51607|


#### File 5: bin boundary file (see example_bin_boundary_configure.txt)
This file specifies the boundary of bins used in FACS. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|Ligand=500uM|
|---|------------|------------|
|P1|0.722,10.000|0.722,10.000|
|P2|0.316,0.722|0.316,0.722|
|P3|-0.076,0.316|-0.076,0.316|
|P4|-0.460,-0.076|-0.460,-0.076|
|P5|-0.861,-0.460|-0.861,-0.460|
|P6|-10.000,-0.861|-10.000,-0.861|

**In each cell, comma is used to seperate the upper and lower boundary of the relevant bin. Usually, fluorescence signal is used in FACS to define bins; hence, a log10(FLU) is expected here.**

Note that in the example file, we used log10(GFP/mCherry) to define bins, which is the reason here for negative bin boundary value. GFP is the reporter under control of the studied regulatory element, while mCherry is constitutively expressed used to normalize cell-to-cell variability. -10.000 and 10.000 are used as numerical lower and upper boundary of negative and positive infinite, respectively. Usually, the gating strategy is determined prior to the cell sorting, hence identical gating parameters for all conditions (sorting experiment) are expected. Even though, note that different gating parameters are also supported. 


#### File 6: bin occupation file (see example_bin_occupation_configure.txt)
This file specifies the proportion of cells sorted into each bin in the whole pooled library. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|Ligand=500uM|
|---|------------|------------|
|P1|0.04676|0.05377|
|P2|0.35796|0.39821|
|P3|0.19516|0.17723|
|P4|0.18492|0.16726|
|P5|0.14078|0.12048|
|P6|0.07443|0.08305|

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


### Step 4：Run the pipeline
Put all necessary files mentioned above, as well as the data directory under the working directory.
Open the command line window, cd to the working directory and run the analysis pipeline.

cd path_to_your_working_directory

python HTsensor_main.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 

python HTsensor_main.py example_configure.txt

The program will create an 'error.log' file under the working directory, open this file to check whether anything wrong happens. Generally, no content suggests successful running. Please post your 'error.log' if you cannot figure out the bugs when using it.

For a typical Macbook (for example, 2.6 GHz processor and 8 GB memory), the example test can be finalized within 5 minutes. The rate-limiting step is the mapping of the raw NGS read. For a typical Macbook, we expect a processing speed of 20 million reads per hour.

## Output files
The output files will be organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory (prefiex_results). We term this subdirectory 'result directory' thereafter.

Three subdirectories are located under this result directory, namely, prefix_rawcount, prefix_cleandataset and prefix_optimization. These directories store results about mapping of the raw NGS data to synthetic library, mutant-centered read count data in multiple conditions, calculated mutant response, respectively.

Below is the description for each of them.

### prefix_rawcount/ (mapping of the raw NGS data to synthetic library)
-------------------------------------------------------------
**prefix.countsummary.txt**: basic statistics of the mapping ratio of each NGS library with a header line using tab as delimiter.

File|Label|Reads|Mapped|Synerror|Unknown|Percentage|Zerocounts|GiniIndex
----|-----|-----|------|--------|-------|----------|----------|---------
example_data/plasmid.fq.gz|plasmid|1000000|838484|90356|71160|0.8385|1734|0.2167
...|...|...|...|...|...|...|...|...

Reads denote the number of reads in the raw data. Mapped denotes number of reads mapping perfectly to one member of the synthetic mutant library. Synerror refers to those reads with one indel mutation or more mismatch mutations. Unknown refers to those reads where no forward_prefixseq or forward_suffixseq can be identified. Percentage is the mapping ratio. Zerocount refers to sgRNA number in the *in silico* library without any corresponding read detected. GiniIndex is a metric reflecting the member abundance uniformity in a library. Bigger Gini index indicates more biased distribution with over- represented or diluted members. Generally, more stringent the selective condition is (e.g. the bin with highest fluorescence), bigger Gini index we can expect.

**prefix.count.txt**: raw read count for each mutant in the *in silico* library.

sensor|Lib1|B1P1|B1P2|B1P3|B1P4|B1P5|B1P6|B2P1|B2P2|B2P3|B2P4|B2P5|B2P6
------|----|----|----|----|----|----|----|----|----|----|----|----|----
TnaC|4|0|0|0|7|11|4|0|0|0|10|18|0
....|...|...|...|...|...|...|...|...|...|...|...|...|...

**This dataset is used for following data processing.**

**prefix.unmapped.txt**: unmapped read in NGS raw data, having a header line using tab as delimiter

unmapped read|Lib1|B1P1|B1P2|B1P3|B1P4|B1P5|B1P6|B2P1|B2P2|B2P3|B2P4|B2P5|B2P6
-------------|----|----|----|----|----|----|----|----|----|----|----|----|----
TGGTGATGGCTACAGAAGGGCAAATCAAGGGCGGGTGGATCGACAATTTTGTTGTCAATTTGGAACCATTTTGAGGTCACACATATATGTAAGATATTCATAATGCACTTATCCTCGCAAGACACAGCCATGGTC|0|0|0|0|0|0|0|0|0|1|0|0|0\
....|...|...|...|...|...|...|...|...|...|...|...|...|...

[**prefix_Libray_Gini_Score.png**](./image/myexample_Libray_Gini_Score.png): schematic of Gini index for each library.

### prefix_cleandataset/ (mutant-centered read count data in multiple conditions)
-------------------------------------------------------------

**prefix.eliminate.txt**: mutants that are eliminated from further analysis, due to the over-diluted representation in the initial library; a simple list flat file with one sgRNA each line. The threshold used here is defined in the configure file.

**initial_ab.csv**: the absolute abundances in the initial library (prior to sorting) of mutants passing quality control (used in further analysis).

**Under this directory, the program also creates another sensor_ctab/**: storing a series of .csv files. Each file corresponds to one mutant, specifying the read count in each condition and each bin. Each file has a header line and index column using tab as delimiter. Here is an example.

  |Ligand=100uM|Ligand=500uM
  |------------|------------
P1|1144|1644
P2|2650|2560
P3|497|486
P4|304|386
P5|402|316
P6|289|277

### prefix_optimization/: calculated mutant response
-------------------------------------------------------------

**sensor_Log10u.csv**: the average response of each mutant (row) in each condition (column). This file has a header line and index column using tab as delimiter. Note that we assume the response of each mutant follows a normal distribution (this file presents average of this normal distribution). For details, see our paper. Here is an example. 

Sensor|Ligand=100uM|Ligand=500uM
------|------------|------------
TnaC|0.0579|0.4862
TnaC_D21F_TTT|-0.2709|0.0251
TnaC_D21L_CTT|0.0265|0.1548
...|...|...
 
**sensor_sigma.csv**: the response noise of each mutant (row) in each condition (column). This file has a header line and index column using tab as delimiter. Note that we assume the response of each mutant follows a normal distribution (this file presents the deviation of this normal distribution). For details, see our paper. Here is an example.

Sensor|Ligand=100uM|Ligand=500uM
------|------------|------------
TnaC|0.3889|0.1466
TnaC_D21F_TTT|0.1761|0.2094
TnaC_D21L_CTT|0.2153|0.2113
...|...|...
 
**sensor_negLog10P.csv**: the -log10 formatted probability (thus smaller this value, bigger probability) for the calculated response of each mutant (row) in each condition (column). This file has a header line and index column using tab as delimiter. Note that we use a maximized probability algorithm to calculate the response (this file presents this probability). For details, see our paper. Here is an example.

Sensor|Ligand=100uM|Ligand=500uM
------|------------|------------
TnaC|1113.7|1078.8
TnaC_D21F_TTT|44.7|10.2
TnaC_D21L_CTT|12.8|51.7
...|...|...

**Under this directory, the program also creates another heatmap/**: storing a series of .png plots. Each file corresponds to one mutant, specifying the heatmap of probability for each searched average response and noise in each condition. In our algorithm, we generated a 2D grid to search for the optimal average response and noise, while for each of them, one probability can be calculated. Here is an [**example**](./image/TnaC_D21P_CCC_Ligand=500uM_opt.png).
