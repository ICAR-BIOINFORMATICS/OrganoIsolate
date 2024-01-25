# OrganoIsolate (Organelle genome extraction and assembly Tool)
A user-friendly tool designed to facilitate the assembly of chloroplast and mitochondria genomes from the raw sequencing data.  This manual provides a step-by-step guide on how to use the pipeline effectively. 

### Before starting please ensure that you have the following dependencies installed before proceeding:


Python (version 3 or higher):https://docs.python-guide.org/starting/install3/linux/

Under the python two modules are required to load namely "pandas" (pip3 install pandas) and "matplotlib" (pip3 install matplotlib)


Trimmomatic: https://git.launchpad.net/ubuntu/+source/trimmomatic/log/?h=ubuntu/focal

## Download the OrganoIsolate package from the github using 

		https://github.com/ICAR-BIOINFORMATICS/OrganoIsolate.git 


 
## After downloading the pipeline package, you need to follow the steps below one by one. 
 
Step 1: Downloaded python scripts to a desired location on your computer 

		oraganoIsolate.py and OI.py

Step 2: Preparing the Input Data Before running the pipeline, make sure you have the following input data ready:
		
1. Raw sequencing reads in farstq.gz format

2. Respective Reference genome of chloroplast or mitochondria
   
4. Python scripts that we have downloaded in the previous steps



Running the Pipeline Follow these steps to run the ChloroMiner:

1: Open a terminal or command prompt on your computer (Ctrl+Alt+T).

2: Navigate to the directory where the pipeline is installed.

3: Execute the pipeline script using the following command:

        python3 organoIsolate.py -adap <Trimomatic_adapter_path> -mref <reference_mitochondria_genome_file> -cref <reference_chloroplast_genome_file> -threads <no of threads >



Here 
-adap is the path of adapters which downloaded under the Trimmomatric folder

-mref is the reference mitochondrial genome (fasta file) which you would like to take as a reference 

-cref is the reference chloroplast genome (fasta file) which you would like to take as a reference 

-threads represent the no of threads can be adjusted as per the computer system’s specification 

Note: python script, adapter path and reference genomes are space separated 

If you want to run the script to get one organelle genome then keep one parameter of your interest (either -mref or -cref)

## Developed by:

Dr. Samarth Godara, Scientist, Division of Computer Applications, ICAR-Indian Agricultural Statistics Research Institute (IASRI), Library Avenue, Pusa, New Delhi-110012, contact: samarth.godara@icar.gov.in.

Dr. Shbana Begam, Scientist, ICAR-National Institute for Plant Biotechnology, New Delhi-110012, contact: shbana.begam@icar.gov.in
