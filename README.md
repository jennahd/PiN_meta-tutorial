# PiN_meta-tutorial

**Mini-metagenomics tutorial "Metagenomics is for Protists too"**

![Option1](https://github.com/jennahd/PiN_meta-tutorial/assets/22906565/97bc807a-c744-4634-99f0-11a2272150ee)

## Installing Tools

If you are comfortable using the command-line, please install the following tools prior to the tutorial. But, please don't worry if you are not comfortable using the command-line. We will also provide information on web-based metagenome search tools. Also, in that case, we recommend asking a neighbour who is using the command-line to follow along during the tutorial :)

If you have any issues installing software, don't worry! The output files for each step are provided in the relevant folder ending in "_ready". So even if you are only able to install some of the tools, it will be useful. The most important tool to install is [Anvi'o](https://anvio.org/).

### 1. Install Conda

If you have it installed already, please install [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.anaconda.com/free/miniconda/index.html), the minimal installer for conda. Conda is a package and environmental management and installation tool that can help to solve a lot of installation and errors running programs introduced by the fact that tools require different dependencies and versions of dependencies. Generally, it will make installing everything easier!

You can find Miniconda with system-specific installation instructions here:
[https://docs.anaconda.com/free/miniconda/](https://docs.anaconda.com/free/miniconda/)

### 2. Install BLAST+

The [BLAST+ tools](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) allow you to search various kinds of data stored at [NCBI](https://www.ncbi.nlm.nih.gov/)/[ENA](https://www.ebi.ac.uk/ena/browser/home)/[DDBJ](https://www.ddbj.nig.ac.jp/index-e.html), with data uploaded to each mirrored across the others. We are specifically interested in using their tools to search whole genome shotgun (WGS) sequence projects virtually on the command-line (blastn_vdb and tblastn_vdb), without having to download the data locally.

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

As an option, you can also run the tool TIARA ([Karlicki, Antonowicz, and Karnkowska](https://academic.oup.com/bioinformatics/article/38/2/344/6375939))as part of Whokaryote, which is a similar tool that identifies eukaryotic contigs in metagenomic datasets using a deep-learning approach. You can find information about the tool here: https://github.com/ibe-uw/tiara.

Another option for distinguishing prokaryotic and eukaryotic contigs in metagenomic dataset is EukRep ([West et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29496730/)): https://github.com/patrickwest/EukRep

The only way to install whokaryote is by making a conda environment. In our tests, this didn't work properly on all systems, so don't worry if it doesn't work for you (though you are of course always welcome to troubleshoot any installation issues!).
```
conda create -c bioconda -n whokaryote whokaryote
```

If you get the error "UnsatisfiableError: The following specifications were found to be incompatible with each other:" a potential work-around is outlined here: [https://github.com/LottePronk/whokaryote/issues/8](https://github.com/LottePronk/whokaryote/issues/8)
```
conda create -n whokaryote -y && conda activate whokaryote
conda install -c conda-forge -c bioconda tiara numpy=1.19.4 whokaryote
```

### 4. Install MetaBAT2

MetaBAT2 ([Kang et al., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/)) is a metagenomic binning tool that can be used to reconstruct metagenome assembled genomes from metagenomic assemblies by clustering contigs based on their statistical properties and read coverage across different metagenomic datasets. You can find information about the tool here: https://bitbucket.org/berkeleylab/metabat.

The easiest way to install the tool is to make a conda environment and then install metabat2. You can install it as one step, but in our tests that didn't work on every system.
```
conda create -n metabat2
conda install bioconda::metabat2
```

### 5. Install Anvi'o

[Anvi'o](https://anvio.org/) is an "An open-source, community-driven analysis and visualization platform for microbial 'omics." that focuses on interactive visualization. They have extensive tutorials available for visualizing and interacting with data generated from various omics tools. We will use this tool to look at how well contigs are distinguished as eukaryotic by Whokaryote and how well metaBAT2 performed in terms of binning using our dataset.

Please follow the system-specific instructions for installing anvi'o found here: https://anvio.org/install/

**DO NOT RUN "conda update conda”** as outlined in the installation tutorial, as this broke our versions of conda on some systems. If you missed this and encounter the issue, I recommend uninstalling and reinstalling conda and then going through the installation tutorial again without this step.

You may also encounter an error "Failed building wheel for datrie" during the step to install anvi'o with pip (https://github.com/merenlab/anvio/issues/2215). To resolve this, run `mamba install datrie`, and then the pip command again.

**Below are the instructions for the mini-metagenomics tutorial, you can stop reading here :)**

## Tutorial Overview

Throughout the tutorial, please pay attention to which commands you should run today and which are provided only as an example for future reference. For each section, the files that should be output are available in the respective folder named "_ready" (in case that step doesn't work for you), while additional files are provided in the folders named "_extra". If one step doesn't work for you (e.g., you weren't able to install that software), you can add "_ready" to the relevant folder for subsequent steps. For example, if you are unable to install whokaryote, for future steps that require output files from this program you can change the folder name from [whokaryote](whokaryote) to [whokaryote_ready](whokaryote_ready).

For the tutorial, we will be searching for several species of Nucleariidae, namely (_Parvularia atlantis_, _Pompholyxophrys punicea_, _Nuclearia simplex_, and _Fonticula alba_) in assembled WGS metagenomic datasets, searching the available genomes (_Parvularia atlantis_ and _Fonticula alba_) against SRA datasets, and perform binning on the metagenome of the mixed culture that the _Parvularia atlantis_ genome was obtained from using both automated binning and interactive visualization methods. Nucleariidae is a group of protists that belongs to Opisthokonta (within Obazoa) that together with Fungi forms the Holomycota. You can find more information about Nucleariidae in this recent review: [Gabaldón, Völcker, and Torruella, 2022](https://doi.org/10.1016/j.protis.2022.125895) and information about how the _Parvularia atlantis_ genome was obtained through several assembly and manual curation steps of a mixed culture metagenome here: [Ocaña-Pallarès et al., 2022](https://doi.org/10.1038/s41586-022-05110-4)

The tutorial presentation including information on metagenomics for protists is available here: [PiN_meta-tutorial.pptx](PiN_meta-tutorial.pptx).

**As a reminder, please help your neighbours and take the opportunity to discuss your results and any issues with each other :)
**

## Tutorial Part 1: Searching Metagenomes

There is a lot of publically available metagenomic data hosted on NCBI/ENA/DDBJ in the form of sequence read archive (SRA) and whole genome shotgun (WGS) sequencing projects, but how to search these datasets isn't immediately obvious.

### Searching WGS metagenomes - web

In this part of the tutorial we will search for four species of Nucleariidae (_Parvularia atlantis_, _Pompholyxophrys punicea_, _Nuclearia simplex_, and _Fonticula alba_) is assembled metagenomes available on NCBI using their 18S rRNA gene sequences, which can be found in the file: [sequences/Nucleariidae_18SrRNAgenes.fasta](sequences/Nucleariidae_18SrRNAgenes.fasta)

We can find information on which kinds of metagenomic datasets are available by searching "metagenomes" in the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi), which will give you a list of metagenome taxonomies available.

Go ahead and take a look. Now, find and click on "lake water metagenome". Here you will see the Taxonomy ID associated with metagenomes from this environment "1647806". On the right-hand side of the webpage, we can also see the various Entrez records from each Database associated with this taxid. There are currently 18,820 SRA projects, which includes amplicon (metabarcoding), metatranscriptome, and metagenome sequencing data. Of these, 1931 have used the strategy "WGS", indicating that they are shotgun metagenomes. However, unfortunately, only 36 (1.9%) are provided as metagenome assemblies (see the Entrez "Assembly" records). **As a reminder, please upload your assemblies to public repositories and not only your raw sequence data!**

Let's take a look at a second metagenomic environment "moss metagenome" and make note of the associated taxid "1675540". We will now try searching for several species of Nucleariidae in this environment using the web-version of BLASTN. **Note: WGS projects aren't included in the standard BLAST databases nt and nr. To search them, we need to specifically search WGS projects and provide a relevant taxid (you can't search them all at once)**

1. Navigate to the NCBI Nucleotide BLAST (blastn) page: https://blast.ncbi.nlm.nih.gov/Blast.cgi
2. Upload the file [sequences/Nucleariidae_18SrRNAgenes.fasta](sequences/Nucleariidae_18SrRNAgenes.fasta)
3. Under "Database", select "Whole-genome shotgun contigs (wgs)" from the drop-down menu.
4. Under "Limit by" we will keep "Organism" selected and in the provided box type "moss metagenome" or its taxid "1675540".
5. Under "Program Selection" choose the "Somewhat similar sequences (blastn)" option and then press "BLAST".

Running the search will take a bit, so in the meantime you can move on to the next step.

_When it's done running, here are some things to consider:_
- Which species do we find closer hits to in moss metagenomes?
- How many would you say are "close" hits to each species? Are they the same species? The same genus?
- Do we find the same top hits for the different searched species?
- Do you expect that the top hits will cluster together with the searched species in a phylogenetic tree?

### Searching SRA metagenomes - web

It can be time and resource intensive to search SRA projects, and it is most common to search for specific marker genes in metagenomic sequence reads (e.g., mOTUs or mTAGs). But there is now a tool available online that we can use to search entire genomes against metagenomic reads found in the SRA archive! [Branchwater Metagenome Query](https://branchwater.jgi.doe.gov) can be used to search sequences (50kb or longer) against reads and within a few minutes will return SRA metagenomes with hits alongside available metadata.

Let's try it out with:
1. The [_Parvularia atlantis_ genome](sequences/GCA_943704415.1_Parvularia_atlantis_v1_genomic.fna).
2. The [_Fonticula alba_ genome](sequences/GCA_000388065.2_Font_alba_ATCC_38817_V2_genomic.fna).
3. Then find and download a genome (or long sequence fragment) from your favourite protist and try searching it.

When the searches are done running, here are some things to consider:
- What does the location data tell us about the global distribution of hits?
- Which kinds of environments are most hits from? Are there any patterns?
- Try adusting the cANI (containment average nucleotide identity) to 0.97 (species-level) and 0.94. Which specific environments are the closest relatives of each species found in?

### Searching WGS metagenomes - command-line

It can be time and data storage intensive to iteratively download metagenomes to search them. Thankfully, NCBI provides a script that can be used to collect WGS accessions associated with each taxid into a database file that we can then search virtually using "blastn_vdb" (nucleotide sequence search) or "tblastn_vdb" (protein sequence search).

Running the NCBI taxid2wgs.pl script can be problematic as it used Perl modules that can be tricky to install on some systems. I've therefore gone ahead and put curated database files (only select accessions included so that searches are faster) for several metagenome environments where I know Nucleariidae are found in the [taxid2wgs](taxid2wgs) folder.

Let's try searching metagenomes using blastn_vdb and the taxid2wgs database files. If needed, first activate your blast conda environment.

Here's an example for moss metagenomes:
```
blastn_vdb \
    -query sequences/Nucleariidae_18SrRNAgenes.fasta \
    -db taxid2wgs/moss_metagenome  \
    -out taxid2wgs/Nucleariidae_vs_moss_metagenome.tsv \
    -outfmt "6 qacc qlen sacc slen length evalue pident qcovs" \
    -max_hsps 1 \
    -perc_identity 93 \
    -qcov_hsp_perc 30
```
If you get an error, try running the command again. It is probably a random network error and it will work eventually (I have emailed NCBI about this and apparently there is no fix).

Take a look at the resulting "taxid2wgs/Nucleariida_vs_moss_metagenome.tsv" file. The results for the searches can also be found in [taxid2wgs_ready](taxid2wgs_ready).

**If at a later point you want to make your own database files and search many accessions, below are some suggestions on how to do that. But you don't need to do that today. So feel free to move on to the phylogenetic tree section**

The taxid2wgs.pl script can be downloaded from: https://ftp.ncbi.nlm.nih.gov/blast/WGS_TOOLS/
I've downloaded it for you, and you can find it in the [taxid2wgs_extra](taxid2wgs/taxid2wgs.pl) folder.

Here is an example of how you would make a database file for moss metagenomes:
```
perl taxid2wgs_extra/taxid2wgs.pl \
  -title moss_metagenome \
  -alias_file moss_metagenome \
  1675540
```
If you get the folling error when running taxid2wgs.pl "500 Can't verify SSL peers without knowing which Certificate Authorities to trust", install the perl module `cpan Mozilla::CA`, use sudo if necessary and if you don't have sudo access try installing local::lib. This fix worked for me, but it doesn't seem to for every system, so you might need to do some trouble-shooting.

I've also run the script on the metagenome taxids found in this file: [taxid2wgs_extra/taxids.tsv](taxid2wgs_extra/taxids.tsv) and put the resuting database files here: [taxid2wgs_extra/taxid2wgs_databases](taxid2wgs_extra/taxid2wgs_databases) for your reference.

When running database files with a large number of accessions, you can recieve a random network error that will kill the search. To resolve this, I tend to make a file with a list of all accessions and then iteratively search each accession (you can provide a specific accession to the -db flag) in a loop (or submit jobs in a loop for each accession to a cluster, such as UPPMAX). That way, if one search fails it doesn't kill the entire search and I can then run the specific accession searches that failed again.

### Resulting phylogenetic tree

I've taken all of the hits found across non-animal metagenomes and inferred a maximum likelihood phylogeny for you to take a look at.

More specifically, I included reference 18S rRNA gene sequence diversity from across Opisthokonta (representatives from taxa found in [EukProt v3](https://evocellbio.com/eukprot/) and unicellular Holozoa found in [EukRibo v1](https://zenodo.org/records/6327891)), with select additional sequences from NCBI, and an outgroup of other Obazoa (Apusomonadidae and Breviatea). In addition, I've added sequences annotated as "Rotosphaerida" (another formal name for nucleariids) from (PR2)[https://pr2-database.org/]. I then used the tool CD-HIT (https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki) to reduce redundancy and removed 100% identical sequences. I retrieved the metagenomic contigs with hits from NCBI using the Entrez E-utilities(https://www.ncbi.nlm.nih.gov/books/NBK25500/). I then extracted the 18S rRNA gene region using [Barrnap](https://github.com/tseemann/barrnap). However, Barrnap doesn't always work so well for protists, and using an hmm profile specific for your group of interest is probably the best option. Of note, about a third of metagenomic contigs also encoded the 28S rRNA gene on the same fragment. I then combined these metagenomic 18S rRNA gene sequences with the reference sequences and sequences annotated as "Rotosphaerida" from and a recent study of eukaryotic diversity in various envioronments using long-read metabarcoding: [Jamy et al., 2022](https://doi.org/10.1038/s41559-022-01838-4). I then alignment the sequences using [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html) (mafft-einsi), trimmed poorly aligned and sparse regions with [TrimAL](https://vicfero.github.io/trimal/) (-gt 0.2), and inferred a maximum-likelihood phylogeny using IQ-TREE(http://www.iqtree.org/) (with the GTR+FO+R7 model of evolution selected, and 1000 SH-aLRT and 1000 ufboot supports).

You can find the resulting phylogeny here: (tree/ALL_sequences.einsi.gt20perc.trim.treefile)[tree/ALL_sequences.einsi.gt20perc.trim.treefile].

1. Open the tree using [iTOL](https://itol.embl.de/).
2. Adjust the visualiztion of the tree as you would like and set the root to the outgroup.
3. Add colours for the different groups by dragging the file (tree/colors_styles_template.txt)[tree/colors_styles_template.txt] onto the tree.
4. Add environmental information for the metagenomic and long-read metabarcoding sequences by dragging the file (tree/dataset_binary_template.txt)[tree/dataset_binary_template.txt] onto the tree.

Some things to consider:
- daf
- adff

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
