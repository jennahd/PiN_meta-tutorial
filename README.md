# PiN_meta-tutorial
Metagenomics tutorial for the 2024 Protistology Nordics Meeting

## Installing Software

### 1. Install Conda

If you don't already have it, please install Conda or Miniconda. It will make installing everything else easier.

You can find Miniconda with installation instructions per system here:
[https://docs.anaconda.com/free/miniconda/](https://docs.anaconda.com/free/miniconda/)

### 2. Install BLAST+

**Option 1:** Can install in a conda env (doesn't work for newer macs)
conda create -n blast
conda activate blast
conda install bioconda::blast

**Option 2:** Install from a dmg file
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

For dmg file, if you get a warning that it's from the web and can't open, right click and click open with and then it will install after you press okay

**Option 3:** Download executables
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Download and add executable files or symlinks to /user/local/bin or wherever you like to keep them.

### 3. Install Whokaryote
You can find infomation about this tool for distinguishing prokaryotic and eukaryotic metagenomic contigs here: [https://github.com/LottePronk/whokaryote](https://github.com/LottePronk/whokaryote)

Make a conda environment and install whokaryote
```
conda create -c bioconda -n whokaryote whokaryote
```

### 4. Install Metabat2
You can find information about this tool for binning of metagenomic contigs here: https://bitbucket.org/berkeleylab/metabat

Make a conda environment and then install metabat2 (it didn't work for me doing it in one step)
```
conda create -n metabat2
conda install bioconda::metabat2
```

### 5. Install Anvi'0

Follow the system-specifc instructions for installing anvi'o found here: https://anvio.org/install/

If newly installing on MacOS with an M1 or M2 chip, running "conda update conda” as outlined in the installation tutorial caused issues for me. If you encounter this issue I recommend uninstalling and reinstalling conda and then running the installation tutorial without this step

You may also encounter an error "Failed building wheel for datrie" during the anvio installation with pip (https://github.com/merenlab/anvio/issues/2215). To resolve this run "mamba install datrie", and then the pip command again
#To resolve this run "mamba install datrie", and then the pip command again
