# PiN_meta-tutorial

**Mini-metagenomics tutorial "Metagenomics is for Protists too"**

![Option1](https://github.com/jennahd/PiN_meta-tutorial/assets/22906565/97bc807a-c744-4634-99f0-11a2272150ee)

## Installing Tools

If you are comfortable using the command-line, please install the following tools prior to the tutorial. But, please don't worry if you are not comfortable using the command-line. We will also provide information on web-based metagenome search tools. Also, in that case, we recommend asking a neighbour who is using the command-line to follow along during the tutorial :)

If you have any issues installing software, don't worry! The output files for each step are provided in the relevant folder ending in "_ready". So even if you are only able to install some of the tools, it will be useful. The most important tool to install is [Anvi'o](https://anvio.org/).

### 1. Install Conda

If you have it installed already, please install [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.anaconda.com/free/miniconda/index.html), the minimal installer for conda. Conda is a package and environmental management and installation tool that can help to solve a lot of installation and errors running programs introducted by the fact that tools require different dependencies and versions of dependencies. Generally, it will make installing everything easier!

You can find Miniconda with system-specific installation instructions here:
[https://docs.anaconda.com/free/miniconda/](https://docs.anaconda.com/free/miniconda/)

### 2. Install BLAST+

The [BLAST+ tools](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) allow you to search various kinds of data stored at [NCBI](https://www.ncbi.nlm.nih.gov/)/[ENA](https://www.ebi.ac.uk/ena/browser/home)/[DDBJ](https://www.ddbj.nig.ac.jp/index-e.html), with data uploaded to each mirrored across the others. We are specifically interested in using their tools to search data whole genome shotgun (WGS) sequence projects virtually on the command-line (blastn_vdb and tblastn_vdb), without having to download the data locally.

There are three different options for installing BLAST+ outlined below.

**Option 1:** Install in conda environment

Unfortunately, this option doesn't work for newer macs with M1 or M2 chips. But otherwise you can make a conda environment and install BLAST+ there.
```
conda create -n blast
conda activate blast
conda install bioconda::blast
```

The two specific tools of interest, "blastn_vdb" and "tblastn_vdb", aren't included in the version of blast automatically installed. To get those you need to run the following.
```
conda update blast
```

**Option 2:** Install the application

You can access system-specific installations here:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

If you see a warning that it is an application from the web and can't be opened, right click and open it from there and it will install once you press okay.

**Option 3:** Download executables

You can also directly download the executables here:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Then directly add the executables or symlinks to the executables to /user/local/bin or wherever you like to keep them.

### 3. Install Whokaryote

Whokaryote ([Pronk and Medema, 2022](https://pubmed.ncbi.nlm.nih.gov/35503723/)) is a tool for distinguishing prokaryotic and eukaryotic contigs from metagenomic data based on gene structure. You can find information about the tool here: [https://github.com/LottePronk/whokaryote](https://github.com/LottePronk/whokaryote).

As an option, you can also run the tool TIARA ([Karlicki, Antonowicz, and Karnkowska](https://academic.oup.com/bioinformatics/article/38/2/344/6375939))as part of Whokaryotes, which is a similar tool that identifies eukaryotic contigs in metagenomic datasets using a deep-learing apprach. You can find information about the tool here: https://github.com/ibe-uw/tiara.

Another option for distinguishing prokaryotic and eukaryotic contigs in metagenomic dataset is EukRep ([West et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29496730/)): https://github.com/patrickwest/EukRep

The only way to install whokaryote is by making a conda environment. In our tests, this didn't work properly on all systems, so don't worry if it doesn't work for you (though you are of course always welcome to troubleshoot any installation issues!).
```
conda create -c bioconda -n whokaryote whokaryote
```

### 4. Install MetaBAT2

MetaBAT2 ([Kang et al., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/)) is a metagenomic binning tool that can be used to reconstruct metagenome assembled genomes from metagenomic asssemblies by clustering contigs based on their statistical properties and read coverage across different metagenomic datasets. You can find information about the tool here: https://bitbucket.org/berkeleylab/metabat.

The easiest way to install the tool is to make a conda environment and then install metabat2. You can install it as one step, but in our tests that didn't work on every system.
```
conda create -n metabat2
conda install bioconda::metabat2
```

### 5. Install Anvi'o

[Anvi'o](https://anvio.org/) is an "An open-source, community-driven analysis and visualization platform for microbial 'omics." that focuses on interactive visualization. They have extensive tutorials available for visualizing and interacting with data generated from various omics tools. We will use this tool to look at how well contigs are distinguised as eukaryotic by Whokaryote and how well metaBAT2 performed in terms of binning using our dataset.

Please follow the system-specifc instructions for installing anvi'o found here: https://anvio.org/install/

**DO NOT RUN "conda update conda”** as outlined in the installation tutorial, as this broke our versions of conda on some systems. If you missed this and encounter the issue, I recommend uninstalling and reinstalling conda and then going through the installation tutorial again without this step.

You may also encounter an error "Failed building wheel for datrie" during the step to install anvi'o with pip (https://github.com/merenlab/anvio/issues/2215). To resolve this, run `mamba install datrie`, and then the pip command again.

**Below are the instructions for the mini-metagenomics tutorial, they will be updated again before the meeting, so you can stop reading here :)**

## Part 1: Searching metagenomes

In this part of the tutorial we will search for four species of Nucleariida (_Parvularia atlantis_, _Pompholyxophrys punicea_, _Nuclearia simplex_, and _Fonticula alba_) is assembled metagenomes available on NCBI using their 18S rRNA gene sequences, which can be found in the file: sequences/Nucleariida_18SrRNAgenes.fasta

First, we want to find the taxids of metagenomes of interest. Take a look at the NCBI taxonomy browser (https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) and find the taxid for "moss metagenome".

Now we want to extract all whole genome shotgun (WGS) sequences projects that have this taxid. To retrieve those accession, we need to use a specific script from NCBI.

The taxid2wgs.pl script can be downloaded from: https://ftp.ncbi.nlm.nih.gov/blast/WGS_TOOLS/
I've downloaded it for you, and you can find it in the [taxid2wgs](https://github.com/jennahd/PiN_meta-tutorial/blob/main/taxid2wgs/taxid2wgs.pl) folder.

Run taxid2wgs and collect moss metagenome WGS accessions in a database file that can be used for virtual blast searches
```
perl taxid2wgs/taxid2wgs.pl \
  -title moss_metagenome \
  -alias_file moss_metagenome \
  1675540
```

If you get the folling error when running taxid2wgs.pl "500 Can't verify SSL peers without knowing which Certificate Authorities to trust", install the perl module `cpan Mozilla::CA`, use sudo if necessary and if you don't have sudo access try installing local::lib.

No worries if it doesn't work! We will be using versions of these databases with only select accessions included.

Now let's try searching metagenomes using blastn_vdb (there is also tblastn_vdb for searching protein sequences, with the database files I've prepared that include only select metagenome accessions found in the taxid2wgs folder (bioreactor, freshwater, hydrothermal vent, lake water, moss, and soil metagenomes).

```
blastn_vdb \
    -query sequences/Nucleariida_18SrRNAgenes.fasta \
    -db taxid2wgs/bioreactor_metagenome  \
    -out taxid2wgs/Nucleariida_vs_bioreactor_metagenome.tsv \
    -outfmt "6 qacc qlen sacc slen length evalue pident qcovs" \
    -max_hsps 1 \
    -perc_identity 93 \
    -qcov_hsp_perc 30
```
I've taken all of the hits found across non-animal metagenomes and inferred a maximum likelihood phylogeny including diversity across Opisthokonta and an outgroup of other Obozoa. You can take a look at the resulting tree found in the folder "tree" using iToL and add the dataset files also found in the "tree" folder to colour the different sequences and add environmental source information (ADDING THESE FILES).

## Part 2: Binning metagenomes

The first step is to download the metagenome assembly of a mixed culture including a eukaryote of interest.

I've already done that and retrieved the assembled metagenome of _Parvularia atlantis_ from https://figshare.com/articles/dataset/Genomic_data_for_Ministeria_vibrans_Parvularia_atlantis_Pigoraptor_vietnamica_and_Pigoraptor_chileana/19895962/1. It can be found in the folder "sequences". The metagenomes is also available on ENA at the project accession PRJEB52884 

**DON'T DO THESE STEPS**
I've also downloaded the metagenomic reads, performed read mapping against the metagenomic contigs with bowtie2 and used the resulting bam file to generate a coverage depth profile of reads mapped against our metagenomic contigs (found in the folder "metabat2") and an anvi'o profile database (found in the folder "anvio"). In addition, I've prepared an anvi'o contigs database and called single-copy genes and taxonomy (found in the folder "anvio"). You can find the steps to prepare these files here (but don't run them, they are time and memory intensive).

```
#Retrieve files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/006/ERR9765196/ERR9765196_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR976/006/ERR9765196/ERR9765196_2.fastq.gz

wget https://figshare.com/ndownloader/files/35315719
unzip 35315719
mv Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta \
    P_atlantis_metagenome.fasta
rm 35315719
#rm -r Parvularia_atlantis - careful when using remove -r!

#Bowtie2 read mapping

bowtie2-build \
    -f Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta \
    Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta

bowtie2 \
    -q \
    --fr \
    -x Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta \
    -1 ERR9765196_1.fastq.gz \
    -2 ERR9765196_2.fastq.gz \
    -p 20 | \
    samtools view \
    -h \
    -@ 20 \
    -bS - > Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.bam

rm Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.sam

samtools \
    view \
    -b \
     -@ 20 \
    -F 4 Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.bam \
    -o Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.mapped.bam

rm Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.bam

samtools \
    sort \
    -@ 20 \
    Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.mapped.bam \
    -o Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.mapped.sorted.bam

rm "Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.mapped.bam

samtools \
    index \
    -@ 20 \
    Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta.mapped.sorted.bam

#Get contigs depth profile
conda activate metabat2

jgi_summarize_bam_contig_depths \
    --outputDepth bowtie/P_atlantis_metagenome_depth.txt \
    bowtie2/Patl_draftAssembly_includesContamination.fasta.mapped.sorted.bam

#Prepare anvio'o contigs and profile databases
conda activate anvio-8

anvi-gen-contigs-database -f Parvularia_atlantis/Patl_draftAssembly_includesContamination.fasta \
                          -n P_atlantis_metagenome \
                          -T 20 \
                          -o anvio_contigs/P_atlantis_metagenome.db

anvi-run-hmms \
    -c anvio_contigs/P_atlantis_metagenome.db \
    -T 20

anvi-run-scg-taxonomy \
    -c anvio_contigs/P_atlantis_metagenome.db \
    -T 4 \
    -P 5

anvi-profile -i bowtie2/Patl_draftAssembly_includesContamination.fasta.mapped.sorted.bam \
             -c anvio_contigs/P_atlantis_metagenome.db \
             -T 20 \
             -S P_atlantis_metagenome \
             --cluster-contigs \
             -o anvio_profile
```

The next step is to sort the metagenomic contings into eukaryotic and prokaryotic fractions. We will do this using the program Whokaryote. I've already run protein calling using prodigal (the default in whokaryote) since that step takes some time. You can find the file in the "whokaryote" folder.

Check how many threads you have on your computer. Using 4 should be okay for most laptops
```
#Check what each of the flags means
whokaryote.py -h

#Run whokaryote
gunzip sequences/P_atlantis_metagenome.fasta.gz

whokaryote.py \
    --contigs P_atlantis_metagenome.fasta \
    --outdir whokaryote \
    --f \
    --model T \
    --threads 4 \
    --gff whokaryote/contigs_genes.gff
```

Now to use the output for anvi'o later, we want to **remove the header from the resulting file whokaryote/whokaryote_predictions_T.tsv**. You can do this using your favourite text editor or on the command-line with nano.

Next, we will bin our metagenomic contigs using MetaBAT2. I've already generated an depth coverage profile which indicates the abundance of the different metagenomic contigs (see above) that we will use in this step.
```
metabat2 \
    -i P_atlantis_metagenome.fasta \
    -o metabat2/bin \
    -a metabat2/P_atlantis_metagenome_depth.txt \
    -t 4
```

From the output we can then generate a tsv file with contig to bin information that we can use for anvio'o
```
for i in metabat2/*.fa ; do 
    bin=$(echo $i | cut -d "/" -f2 | cut -d "." -f1-2 | sed 's/\./_/g')
    grep ">" $i | cut -d ">" -f2 | sed "s/$/\t$bin/g" \
    >> metabat2/contig_bins.tsv
done
```

Now that we have all of our binning information we can add it to our anvi'o profile as a collection!

```
conda activate anvio-8

gunzip anvio/PROFILE.db.gz
gunzip anvio/P_atlantis_metagenome.db.gz

anvi-import-collection  whokaryote/whokaryote_predictions_T.tsv \
                       -p anvio/PROFILE.db \
                       -c anvio/P_atlantis_metagenome.db \
                       --contigs-mode \
                       -C Whokaryote

anvi-import-collection metabat2/contig_bins.tsv \
                       -p anvio/PROFILE.db \
                       -c anvio/P_atlantis_metagenome.db \
                       --contigs-mode \
                       -C MetaBAT
```

Finally, we can now open our prepared anvi'o file in interactive mode!
```
anvi-interactive \
    -p anvio/PROFILE.db \
    -c anvio/P_atlantis_metagenome.db
```
