
README for HiC-BART(1.0.3)

Introduction
============

HiC-BART (Binding Analysis for Regulation of Transcription with Hi-C) is a bioinformatics tool for predicting functional transcription factors (TFs) that is accociated with change in chromatin interaction in the human or mouse genomes, given a differential interaction profile calculated from HiC interaction matrices. HiC-BART leverages over 7,000 human TF binding profiles and over 5,000 mouse TF binding profiles from the public domain (collected in Cistrome Data Browser) to make the prediction.

HiC-BART is implemented in Python and distributed as an open-source package along with necessary data libraries.

**HiC-BART is developed and maintained by the Chongzhi Zang Lab at the University of Virginia.**



# Installation
#### Prerequisites

HiC-BART uses Python's distutils tools for source installation. Before installing HiC-BART, please make sure Python3 (Python 3.3 or higher is recommended) is installed in the system, and the following python packages are installed:

- setuptools
- numpy
- pandas
- scipy

#### Install the full package (All data included, requires at least 15GB hard drive storage in the installation directory)

To install a source distribution of HiC-BART, unpack the distribution tarball and open up a command terminal. Go to the directory where you unpacked BART, and simply run the install script an install BART globally or locally. For example, if you want to install the package HiC-BART-v1.0.0-py3-full.tar.gz:

```shell
$ tar zxf HiC-BART-v1.0.0-py3-full.tar.gz
$ cd HiC-BART-v1.0.0-py3-full
```

Install with root/administrator permission (by default, the script will install python library and executable codes globally):

```shell
$ python setup.py install
```

If you want to install everything under your own directory, for example, a directory as /path/to/bart/, use these commands:

```shell
$ mkdir -p /path/to/bart/lib/pythonX.Y/site-packages 
$ export PYTHONPATH=/path/to/bart/lib/pythonX.Y/site-packages/:$PYTHONPATH 
$ python setup.py install --prefix /path/to/bart 
$ export PATH=/path/to/bart/bin/:$PATH
```

In this value, X.Y stands for the major–minor version of Python you are using (such as 3.5 ; you can find this with sys.version[:3] from a Python command line).

Configure environment variables

You’ll need to add those two lines in your bash file (varies on each platform, usually is `~/.bashrc` or `~/.bash_profile`) so that you can use the BART command line directly:

```shell
$ export PYTHONPATH= "/path/to/bart/lib/pythonX.Y/site-packages/:$PYTHONPATH"
$ export PATH="/path/to/bart/bin/:$PATH"
```



#### Install from source package without data libraries (recommended)

You can download the Human or Mouse Data Library separately under your own directory. In this case, you have to edit the config file (e.g. BART1.0.1/BART/bart.conf) after you unpack the source package to provide the directory for the data. For example, if you download the hg38_library.tar.gz (or mm10_library.tar.gz) and unpack it under /path/to/library, then you can modify the bart.conf file as:

`hg38_library_dir = /path/to/library/`

Then you can run the install script and install BART source package globally or locally same as the full package described above.



# Tutorial
Positional arguments {geneset,profile}


#### Usage

Given a differential interaction HiC matrices, predict functional transcription factors that regulate these genes.

**Usage**:	`hic_bart [-h] -ci <index control> -cm <matrix control> -ti <index treatment> -tm <matrix treatment> -s <species> [-o <outdir>] [-r <view regions>]`

**Example**:	hic_bart -ci control_5000_abs.bed -cm control_5000.matrix -ti treatment_5000_abs.bed -tm treatment_5000.matrix -s hg38 -o bart_output

**Input arguments**:

`-ci --controlindex <index control>`

HiC-Pro index bed file for the control hic sample

`-tm --controlmatrix <matrix control>`

HiC-Pro interaction matrix file for the control hic sample

`-ti --treatindex <index treatment>`

HiC-Pro index bed file for the treatment hic sample

`-tm --treatmatrix <matrix treatment>`

HiC-Pro interaction matrix file for the treatment hic sample

`-s --species <species>`

Species, please choose from "hg38" or "mm10".

`-r --region <viewregion>`

Regions to expand when finding interactions


`-o --outdir <outdir>`

If specified, all output files will be written to that directory. Default: the current working directory





#### Output files

1. **name_auc.txt** contains the ROC-AUC scores for all TF datasets in human/mouse, we use this score to measure the similarity of TF dataset to cis-regulatory profile, and all TFs are ranked decreasingly by scores. There will be one file for up-regulated and down regulated interactions each (treatment against control). The file should be like this:

|   AR_56206   | AUC = 0.854 |
| :----------: | :---------: |
|   AR_56205   | AUC = 0.846 |
|   AR_69287   | AUC = 0.844 |
|   AR_68090   | AUC = 0.835 |
|  EZH2_41813  | AUC = 0.820 |
| SUZ12_74199  | AUC = 0.819 |
| JARID2_41807 | AUC = 0.819 |
|   AR_51837   | AUC = 0.818 |

2. **name_bart_results.txt** is a ranking list of all TFs, which includes the Wilcoxon statistic score, Wilcoxon p value, standard Wilcoxon statistic score (zscore), maximum ROC-AUC score, rank score (relative rank of z score, p value and max auc) and Irwin-Hall P-value for each TF. The most functional TFs of input data are ranked first. There will be one file for up-regulated and down regulated interactions each (treatment against control).  The file should be like this:

   |   TF   | Wilcoxon statistics | Wilcoxon P-value | Z-score | max_AUC | relative_rank | Irwin-Hall P-value |
   | :----: | :-----------------: | :--------------: | :-----: | :-----: | :-----------: | ------------------ |
   |   AR   |       22.917        |    1.572e-116    |  3.228  |  0.854  |     0.004     | 3.126e-07          |
   | POLR3D |        5.770        |    3.972e-09     |  3.038  |  0.582  |     0.021     | 4.125e-05          |
   |  CTCF  |       25.105        |    2.169e-139    |  3.026  |  0.560  |     0.023     | 5.332e-05          |
   |  BRDU  |        4.206        |    1.302e-05     |  2.952  |  0.592  |     0.025     | 7.065e-05          |
   |   PR   |        6.690        |    1.118e-11     |  2.647  |  0.640  |     0.030     | 1.158e-04          |
   | RAD21  |        8.222        |    1.003e-16     |  3.008  |  0.537  |     0.034     | 1.712e-04          |

3. **interactions.bed** is the differential interaction profile from the HiC interaction matrices.