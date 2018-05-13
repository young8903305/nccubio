#!/bin/bash
mkdir sample/genome
bowtie2-build sample/genome.fa sample/genome/genome
nextflow main.nf --p 12 --fqlist sample/full_fqlist.txt\
                    --fqpath sample/fq/\
                    --grouplist sample/full_grouplist.txt\
                    --genome sample/genome/\
                    --gtf sample/genes.gtf\
                    --fa sample/genome.fa\
                    --output full_results/
