# MOPC

MOPC is a python program for getting the Multi-Omics Periphery and Core of cancer.


## Getting Started

1. In order to download MOPC, you should clone this repository via the commands.

​       `git clone https://github.com/wangbingbo2019/MOPC.git`

​       `cd MOPC`

2. Prepare the environment for running the python program and install  the dependencies libraries of CCP which are listed in requirements.txt.

3. Once you have prepared the environment and download the repository, you can run in the command line

​        `python MOPC.py` 
		`python visualization.py`

## -------------------------------------------------

##### Files you need to prepare to run MOPC.py

###### input files

newHnet-2015.txt (The human interactome, an undirected network.)
cancer (The folder of the multi-omics data for 15 cancers,including transcriptome differential expression, DNA differential methylation, somatic mutation, and copy number variation.)

###### output files

(1) zscore_result.txt 			(The LCC z-score result list of 15 cancers in four omics.)
	fc_result.txt 				(The cutoff list of 15 cancers in four omics.)
(2) Transcriptome_CLine.png 	(The CLine result of 15 cancers in Transcriptome omics.)
	Methylation_CLine.png		(The CLine result of 15 cancers in Methylation omics.)
	Somatic_mutation_CLine.png	(The CLine result of 15 cancers in Somatic mutation omics.)
	CNV_CLine.png				(The CLine result of 15 cancers in CNV omics.)
	Transcriptome_UCurve.png	(The UCurve result of 15 cancers in Transcriptome omics.)
	Methylation_UCurve.png		(The UCurve result of 15 cancers in Methylation omics.)
	Somatic_mutation_UCurve.png	(The UCurve result of 15 cancers in Somatic mutation omics.)
	CNV_UCurve.png				(The UCurve result of 15 cancers in CNV omics.)

## Citation

If you use the software, please cite.

Multi-Omics Peripheral and Core Regions of cancer by Bingbo Wang, Xianan Dong, Jie Hu, Lin Gao.





