#!/bin/bash
mkdir sample/genome
bowtie2-build sample/genome.fa sample/genome/genome
nextflow main.nf --p 12 --fqlist sample/fqlist.txt\
                    --fqpath sample/fq/\
                    --grouplist sample/grouplist.txt\
                    --genome sample/genome/\
                    --gtf sample/genes.gtf\
                    --fa sample/genome.fa\
                    --output results/
