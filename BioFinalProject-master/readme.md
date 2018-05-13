# Bioinformatic Final Project
Final Project on 1051.​Theory and Practice of Bioinformatics at NCCU CS

### Goal
Reproduce [Shiao, M.-S. S. et al. Expression Divergence of Chemosensory Genes between Drosophila sechellia and Its Sibling Species and Its Implications for Host Shift. Genome Biol Evol 7, 2843–58 (2015)](https://www.ncbi.nlm.nih.gov/pubmed/26430061)

### Members
1. 魏孝全: Drosophila simulans (GSE67862)
2. [Blueswen](https://github.com/Blueswen): Write Nextflow Script, R Script
3. 江易倫:
4. 李恭儀: Drosophila sechelia (GSE67861) 資料下載、執行、校對、產生 heat map
5. 楊宗翰: Drosophila simulans (GSE67862) 資料下載、執行、校對、cummeRbund 視覺化程式修改、heatmap
6. 楊子萲: Drosophila sechelia (GSE67861) 序列資料下載、執行、校對、產生 heat map

### Data
1. [Sample Data](https://drive.google.com/open?id=0BzrrKiA0aYbbNUljM0d0Unk2N1U): Dataset used at Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks
2. [Experiment Data and Results](https://drive.google.com/drive/folders/0B6ckDKQ9DXOjSzZPcmhRV0Nad28?usp=sharing)

### Dependencies
1. [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml): RNA-Seq reads mapping
2. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/): Transcriptome assembly and differential expression analysis for RNA-Seq
3. [CummeRbund](http://compbio.mit.edu/cummeRbund/): R package for data visualiztion
4. [proto](https://cran.r-project.org/web/packages/proto/index.html): R package for argparse
5. [argparse](https://cran.r-project.org/web/packages/argparse/index.html): R package for parse argv
6. [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html): R package for plot
6. [NextFlow](https://www.nextflow.io/)

### Usage
1. Prepare data: fastq files, fasta file, gtf file
2. Prepare Bowtie index

  ```
  $ mkdir sample/genome
  $ bowtie2-build sample/genome.fa sample/genome/genome
  ```
3. Prepare fqlist file and grouplist files
  * fqlist file example:

    ```
    ID001_1.fastq,ID001_2.fastq
    ID002_1.fastq,ID002_2.fastq
    ID003_1.fastq,ID003_2.fastq
    ID004_1.fastq,ID004_2.fastq
    ...
    ```
  * grouplist file example:

    ```
    ID001,ID002,...
    ID003,ID004,...
    ...
    ```
4. Execute Nextflow script

  ```
  $ nextflow main.nf --p 8 --fqlist sample/fqlist.txt\
                     --fqpath sample/fq/GSE67587/\
                     --grouplist sample/grouplist.txt\
                     --genome sample/genome/genome
                     --gtf sample/gene.gtf\
                     --fa sample/genome.fa\
                     --output results/\
  ```
5. Prepare gene name file and gene families file
  * gene name file: Get gene name from cuff gene short name
  * gene families file example:

    ```
    Ir,6642,6682
    Obp,7093,7142
    Or,7151,7204
    ```
6. Execute R Script

  ```
  $ Rscript plot_hmap.r -i results/diff_out -l gene_name_list -s target_species_name -r gene_families_file -g group_label1,group_label2 -o results/plot
  ```

### Q&A
1. Install CummeRbund Problems
  1. [Install  CummeRbund](https://bioconductor.org/packages/release/bioc/html/cummeRbund.html)
  2. [Upgrade Bioconductor to the latest version available for this version of R](https://rdrr.io/bioc/BiocInstaller/man/BiocUpgrade.html)
  3. [Question: errors when installing cummeRbund package](https://support.bioconductor.org/p/61555/)

### Reference
1. [Shiao, M.-S. S. et al. Expression Divergence of Chemosensory Genes between Drosophila sechellia and Its Sibling Species and Its Implications for Host Shift. Genome Biol Evol 7, 2843–58 (2015).](https://www.ncbi.nlm.nih.gov/pubmed/26430061)
2. [Trapnell, C. et al. Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. Nat Protoc 7, 562–578 (2014).](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html)
