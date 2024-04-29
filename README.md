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

## Binning metagenomes

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
whokaryote.py \
    --contigs P_atlantis_metagenome.fasta \
    --outdir whokaryote \
    --f \
    --model T \
    --threads 4 \
    --gff whokaryote/contigs_genes.gff
```

Now to use the output for anvi'o later, we want to remove the header from the resulting file whokaryote/whokaryote_predictions_T.tsv. You can do this using your favourite text editor or on the command-line with nano.

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
