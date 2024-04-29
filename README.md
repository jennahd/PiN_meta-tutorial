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

## Searching metagenomes

In this part of the tutorial we will search for four species of Nucleariida (_Parvularia atlantis_, _Pompholyxophrys punicea_, _Nuclearia simplex_, and _Fonticula alba_) is assembled metagenomes available on NCBI using their 18S rRNA gene sequences, which can be found in the file: sequences/Nucleariida_18SrRNAgenes.fasta

First, we want to find the taxids of metagenomes of interest. Take a look at the NCBI taxonomy browser (https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) and find the taxid for "moss metagenome".

Now we want to extract all whole genome shotgun (WGS) sequences projects that have this taxid. To retrieve those accession, we need to use a specific script from NCBI.

The taxid2wgs.pl script can be downloaded from: https://ftp.ncbi.nlm.nih.gov/blast/WGS_TOOLS/
I've downloaded it for you, and you can find it in the taxid2wgs folder.

Run taxid2wgs and collect moss metagenome WGS accessions in a database file that can be used for virtual blast searches
```
perl taxid2wgs.pl \
  -title moss_metagenome \
  -alias_file moss_metagenome \
  1675540
```

If you get the folling error when running taxid2wgs.pl "500 Can't verify SSL peers without knowing which Certificate Authorities to trust", install the perl module cpan Mozilla::CA, use sudo if necessary and if you don't have sudo access try installing local::lib.

No worries if it doesn't work! We will be using versions of these databases with only select accessions included.

Now let's try searching metagenomes using blastn_vdb (there is also tblastn_vdb for searching protein sequences, with the database files I've prepared that include only select metagenome accessions found in the taxid2wgs folder (bioreactor, freshwater, hydrothermal vent, lake water, moss, and soil metagenomes).

```
blastn_vdb \                                                               ✔  took 6s    base    at 14:41:40 
    -query sequences/Nucleariida_18SrRNAgenes.fasta \
    -db taxid2wgs/bioreactor_metagenome  \
    -out taxid2wgs/Nucleariida_vs_bioreactor_metagenome.tsv \
    -outfmt "6 qacc qlen sacc slen length evalue pident qcovs" \
    -max_hsps 1 \
    -perc_identity 93 \
    -qcov_hsp_perc 30
```
I've taken all of the hits found across non-animal metagenomes and inferred a maximum likelihood phylogeny including diversity across Opisthokonta and an outgroup of other Obozoa. You can take a look at the resulting tree using iToL and add the dataset files to colour the different sequences and add environmental source information (ADDING THESE FILES).


