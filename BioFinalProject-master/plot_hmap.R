#source("http://bioconductor.org/biocLite.R")
#biocLite("cummeRbund")
library('cummeRbund')
library('proto')
library('argparse')
library('ggplot2')
Rstudio <- FALSE

#####################
## Process options ##
#####################

# Rscript plot_hmap.R -i diff_out -l all_3spe-exist+noDup.csv -s Dsec -r range.txt -g Dsec_M_TW,Dsec_F_TW -o plot

args <- NULL
parser <- ArgumentParser(description='Plot heatmap of cuffdiff diff_out.')
parser$add_argument('-i', '--input', help='input diff_out folder', required=TRUE)
parser$add_argument('-l', '--gene_list', help='gene name list file', required=TRUE)
parser$add_argument('-s', '--sp_name', help='target species name', required=TRUE)
parser$add_argument('-r', '--range_file', help='gene families file', required=TRUE)
parser$add_argument('-g', '--group_name', help='cuffdiff group name', required=TRUE)
parser$add_argument('-o', '--output', help='output folder', required=TRUE)

if(Rstudio){
  m_input <- ''
  m_gene_list <- ''
  m_sp_name <- ''
  m_range_file <- ''
  m_group_name <- ''
  m_output <- ''
  args <- parser$parse_args(c('--input',m_input,'--gene_list',m_gene_list,'--sp_name',m_sp_name,'--range_file',m_range_file,'--group_name',m_group_name,'--output',m_output))
}else{
  args <- parser$parse_args()
}

input <- args$input
gene_list <- args$gene_list
sp_name <- args$sp_name
range_file <- args$range_file
group_name <- args$group_name
output <- args$output

if(!file.exists(input)){
  stop("input diff_out folder doesn't exist.")
}

if(!file.exists(gene_list)){
  stop("gene list file doesn't exist.")
}

if(!file.exists(range_file)){
  stop("gene famliies file doesn't exist.")
}

if(!file.exists(output)){
  dir.create(output)
}

#####################
##    Load Data    ##
#####################

Gid <- read.csv(gene_list,header=TRUE)

if( sp_name %in% colnames(Gid) == FALSE){
  stop(paste(sp_name,"doesn't exist in gene list file."))
}

gene_range <- read.csv(range_file, header=FALSE, col.names=c('fname','start','end'))
cuff <- readCufflinks(input, rebuild=TRUE)

#####################
##       Main      ##
#####################

Gdiff <- diffData(genes(cuff))
myGid <- Gdiff$gene_id
myG <- getGenes(cuff,myGid)
my_df <- myG@annotation

for( k in 1:nrow(gene_range) ){
  targetGid <- c()
  gene_name <- c()
  targetFB <- as.character(Gid[gene_range[k,'start']:gene_range[k,'end'],sp_name])
  for( i in 1:length(targetFB) ){
    for ( j in 1:nrow(my_df) ){
      if( grepl(targetFB[i], my_df[[j,'gene_short_name']]) ){
        targetGid <- c( targetGid, my_df[[j,'gene_id']])
        gene_name <- c( gene_name, as.character(Gid[ i+gene_range[k,'start']-1, 'X.gene']))
      }
    }
  }
  targetG <- getGenes(cuff,targetGid)

  group_label <- strsplit(group_name,',')[[1]]
  hmap <- csHeatmap(targetG,cluster="none",fullnames=FALSE,labRow=FALSE)
  hmap <- hmap + scale_x_continuous(breaks=c(0.25,1.25), labels = c(group_label[1],group_label[2])) + theme(axis.text.x = element_text(size = 15, vjust=1,angle = 0))
  hmap <- hmap + scale_y_continuous(breaks=1:length(gene_name), labels = gene_name) + theme(axis.text.y = element_text(size = 5,vjust=0.5,hjust=1))
  hmap <- hmap + theme( plot.title = element_text(hjust = 0.5, lineheight=0.8, face="bold")) + ggtitle(gene_range[k,'fname'])

  if(Rstudio){
    plot(hmap)
  }else{
    jpeg(filename=paste(output,'/',gene_range[k,'fname'],'_hmap.jpg',sep = ""), width = 8, height = 8, units = 'in', res = 300)
    plot(hmap)
    dev.off()
  }
}
