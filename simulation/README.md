# FACS-seq method to profile mutant library of regulatory element: simulation package

## What is this?
This python script collection is one of the two the software subpackages of FACS-seq, used for simulation prior to real experiment to determine some essential parameters, such as the bin numers used in real sorting experiment.

## General description of the algorithm and experiment
Briefly, the simulation create an de novo mutant library; and specify the response for each of them. Then, in silico sorting and NGS experiments are performed to generate a simulated NGS data. These datasets are processed by all the same scripts as described in "inferring" package to calculate the "experimentally determined" response for each mutant. Finally, comparison is executed to check the consistency between the original setting and "experimentally determined" response.

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

#### File 1: experiment design file (see example_experiment_configure.txt)
This file specifies the experiment design of FACS-seq. In particular, it defines the relation between each NGS raw data and its role in FACS-seq experiment. Usually, to profile the response of a genetic regulatory element of interest, we performed FACS-seq under M different conditions (M ligand concentrations). In each condition, the library is subjected to FACS and sorted into N bins. Hence, in the experiment design file, a N(row)-M(column) matrix is presented to define the role of each NGS raw data (one bin under one condition). **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|
|---|------------|
|P1|B1P1|
|P2|B1P2|
|P3|B1P3|
|P4|B1P4|
|P5|B1P5|
|P6|B1P6|

**The name of each NGS raw data should be the same as those defined in the configure file (see below)**


#### File 2: cell count file (see example_cell_count_configure.txt)
This file specifies the number of cells sorted into each bin during FACS-seq. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|
|---|------------|
|P1|50564|
|P2|100000|
|P3|100000|
|P4|100000|
|P5|100000|
|P6|50923|


#### File 3: bin boundary file (see example_bin_boundary_configure.txt)
This file specifies the boundary of bins used in FACS. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|
|---|------------|
|P1|0.722,10.000|
|P2|0.316,0.722|
|P3|-0.076,0.316|
|P4|-0.460,-0.076|
|P5|-0.861,-0.460|
|P6|-10.000,-0.861|

**In each cell, comma is used to seperate the upper and lower boundary of the relevant bin. Usually, fluorescence signal is used in FACS to define bins; hence, a log10(FLU) is expected here.**

Note that in the example file, we used log10(GFP/mCherry) to define bins, which is the reason here for negative bin boundary value. GFP is the reporter under control of the studied regulatory element, while mCherry is constitutively expressed used to normalize cell-to-cell variability. -10.000 and 10.000 are used as numerical lower and upper boundary of negative and positive infinite, respectively. Usually, the gating strategy is determined prior to the cell sorting, hence identical gating parameters for all conditions (sorting experiment) are expected. Even though, note that different gating parameters are also supported. 


#### File 4: bin occupation file (see example_bin_occupation_configure.txt)
This file specifies the proportion of cells sorted into each bin in the whole pooled library. Similar to experiment design file (file 3), a N-M matrix is expected here, corresponding to N bins and M conditions. **This file is at .csv format using tab as delimiter. It also contains header row and index column, as one example shown below (3 conditions, 6 bins).**

|Bin|Ligand=100uM|
|---|------------|
|P1|0.04676|
|P2|0.35796|
|P3|0.19516|
|P4|0.18492|
|P5|0.14078|
|P6|0.07443|

**Note that parameters here are different from those in file 4!!!** One additional experiment (analysis cytometry using mutant library) is needed to obtain the parameters here. Hence, the sum of each column here should be exactly 1.

### Step 3: Set up the configure file (see example_configure.txt)
The configure file is used to set all the necessary parameters and tell the program where to find necessary files. 
**This file is in a two-column format using tab as delimiter.** Each line starts with one word (name of one parameter) separated with the following (content of this parameter) by a tab delimiter. We describe each parameter as below.

**experiment_configure**: the name of the experiment design file (see above, Step 2, File 1).

**cell_count_configure**: the name of the cell count file (see above, Step 2, File 2)

**bin_boundary_configure**: the name of the bin boundary file (see above, Step 2, File 3)

**bin_occupation_configure**: the name of the bin occupation file (see above, Step 2, File 4)

**LibSize**: number of mutants in synthetic library that is created in silico

**ReadsThreshold**: threshold of read count used to eliminate mutants detected less than this threshold in the library prior to cell sorting. Default: 20.

**log10u_sigma_range**: range of log10u (average response, see our paper for details) and sigma, within which the search for optimal value is executed. **format: " log10u_lower_limit,log10u_upper_limit, sigma_lower_limit,sigma_upper_limit", e.g." -1.5,0.9,0.05,0.4"**. Note that a double quotes is expected and the first character should be blank. We suggest using experimental data to determine the limit for log10u and sigma search range.

**search**: density of 2D grid plane used to search for optimal log10u and sigma value. 

**prefix**: prefix used for naming of all output files, keep it simple without any ‘-’, ‘_’ and ‘ ’. For example, ‘screen20171001’ is fine.

see **example_configure.txt**.


### Step 4：Run the pipeline
Put all necessary files mentioned above, as well as the data directory under the working directory.
Open the command line window, cd to the working directory and run the analysis pipeline.

cd path_to_your_working_directory

python HTsensor_simulation_main.py configure

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 

python HTsensor_simulation_main.py example_configure_simulation.txt

The program will create an 'error.log' file under the working directory, open this file to check whether anything wrong happens. Generally, no content suggests successful running. Please post your 'error.log' if you cannot figure out the bugs when using it.


## Output files
The output files will be organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory (prefiex_results). We term this subdirectory 'result directory' thereafter.

Three subdirectories are located under this result directory, namely, prefix_library, prefix_optimization and prefix_comparison. These directories store results about de novo created mutant library, calculated mutant response and the consistency between the original setting and "experimentally determined" response, respectively.

Below is the description for each of them.

### prefix_library/ (de novo created mutant library)
-------------------------------------------------------------
**my_library.txt**: de novo created mutant library, namely, the initial abundance, average response, and response noise of each mutant (row). This file has a header line and index column using tab as delimiter. Here is an example. 

Sensor|Abundance|Log10u|sigma
------|---------|------|-----
mutant0|0.0233|0.2066|0.1592
mutant1|0.0082|-0.4356|0.2399
...|...|...|...|

### prefix_optimization/: calculated mutant response (all the same as that of inferring package)
-------------------------------------------------------------

**sensor_Log10u.csv**: the average response of each mutant (row) in each condition (column). This file has a header line and index column using tab as delimiter. Note that we assume the response of each mutant follows a normal distribution (this file presents average of this normal distribution). For details, see our paper. Here is an example. 

Sensor|Ligand=100uM|Ligand=500uM
------|------------|------------
mutant0|0.0579|0.4862
mutant1|-0.2709|0.0251
...|...|...
 
**sensor_sigma.csv**: the response noise of each mutant (row) in each condition (column). This file has a header line and index column using tab as delimiter. Note that we assume the response of each mutant follows a normal distribution (this file presents the deviation of this normal distribution). For details, see our paper. Here is an example.

Sensor|Ligand=100uM|Ligand=500uM
------|------------|------------
mutant0|0.3889|0.1466
mutant1|0.1761|0.2094
...|...|...
 
**sensor_negLog10P.csv**: the -log10 formatted probability (thus smaller this value, bigger probability) for the calculated response of each mutant (row) in each condition (column). This file has a header line and index column using tab as delimiter. Note that we use a maximized probability algorithm to calculate the response (this file presents this probability). For details, see our paper. Here is an example.

Sensor|Ligand=100uM|Ligand=500uM
------|------------|------------
mutant0|1113.7|1078.8
mutant1|44.7|10.2
...|...|...

**Under this directory, the program also creates another heatmap/**: storing a series of .png plots. Each file corresponds to one mutant, specifying the heatmap of probability for each searched average response and noise in each condition. In our algorithm, we generated a 2D grid to search for the optimal average response and noise, while for each of them, one probability can be calculated.


### prefix_comparison/: consistency between the original setting and "experimentally determined" response
-------------------------------------------------------------

**simulation_data.csv**: the originally set and "experimentally determined" response for each mutant. This file has a header line and index column using tab as delimiter. Here is an example. 

Sensor|Abundance|Log10u|sigma|Calculated Log10u|Calculated sigma
------|---------|------|-----|-----------------|----------------
mutant0|0.0233|0.2066|0.1592|0.1895|0.1346
mutant1|0.0082|-0.4356|0.2399|-0.4556|0.3039
...|...|...|...|...|...

**Log10u_comparison.png** scatter plot to check the consistency of average response. Here is an [example](./image/Log10u_comparison.png)

**sigma_comparison.png** scatter plot to check the consistency of response noise. Here is an [example](./image/sigma_comparison.png)
