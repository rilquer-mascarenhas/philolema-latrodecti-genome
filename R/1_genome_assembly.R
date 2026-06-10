setwd('/home/rilquer/lgeo_plat/')
require(tidyverse)

# GENOME ASSEMBLY####
## Pre-assembly####
# kmer distribution and genoscope to estimate
# homozygosity and expected genome size
# Compare to post-assembly

## Post-assembly####
# Working on GFA file
require(tidyverse)
# Reading GFA and wrangling to get contig info
plat_ctg <- read_lines(file = 'data/plat_gfa/plat_f0.asm.bp.p_ctg.noseq.gfa') %>%
  as_tibble() %>%
  filter(str_sub(string = value, start = 1,end = 1) == "S") %>%
  mutate(contig = str_split_i(value,'\t',2),
         orig_LN = str_split_i(value,'\t',4),
         orig_cov = str_split_i(value,'\t',5)) %>%
  mutate(LN = as.numeric(str_split_i(orig_LN,':',3)),
         cov = as.numeric(str_split_i(orig_cov,':',3))) %>%
  select(-c(value,orig_LN,orig_cov))
plat_ctg %>%
  ggplot(aes(x=reorder(contig,-LN),y=LN))+geom_point()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_blank())

### Stats####
# Calculating N50, N90, L50 and L90

# N50 and N90
# N50 is calculated by sorting contigs according to their lengths, and then
# taking the halfway point of the total genome length. The size of the contig at that halfway point is the N50 value
N_calc <- function(vec,N){
  cutpoint <- (N/100) * sum(vec)
  vec <- sort(vec,decreasing=T)
  # Cumulative sum and getting the indeces of cumsum that are below the cutpoint
  # N50 is the index for vec right after that, cause that is the contig length that
  # goes above the cutpoint
  return(vec[(max(which(cumsum(vec) < cutpoint)))+1])
}
N_calc(plat_ctg$LN,50)
N_calc(plat_ctg$LN,90)

L_calc <- function(vec,N){
  cutpoint <- (N/100) * sum(vec)
  vec <- sort(vec,decreasing=T)
  # Cumulative sum and getting the indeces of cumsum that are below the cutpoint
  # N50 is the index for vec right after that, cause that is the contig length that
  # goes above the cutpoint
  return(max(which(cumsum(vec) < cutpoint))+1)
}
L_calc(plat_ctg$LN,50)
L_calc(plat_ctg$LN,90)

# Distribution of contig size
ggplot(plat_ctg,aes(x=LN))+
  geom_histogram(color='black',fill='lightblue')+
  theme_bw()

### Coverage####
# Plotting with ggcoverage
# Modify the GFA information of the reads in each sequence to be
# similar to a BedTrack file
# This will probably allow us to visualize edges and gaps
plat_bed <- read_lines(file = 'data/plat_gfa/plat_f0.asm.bp.p_ctg.noseq.gfa') %>% 
  as_tibble() %>%
  filter(str_sub(string = value, start = 1,end = 1) == "A") %>%
  mutate(contig = str_split_i(value,'\t',2),
         start = as.numeric(str_split_i(value,'\t',3)),
         read_start = as.numeric(str_split_i(value,'\t',6)),
         read_end = as.numeric(str_split_i(value,'\t',7)),
         orig_id = str_split_i(value,'\t',8)) %>%
  mutate(read_id = as.numeric(str_split_i(orig_id,':',3)),
         end = start + (read_end - read_start)) %>% 
  select(-c(value,orig_id)) %>%
  select(contig,start,end,read_start,read_end,read_id)

# Retrieving coverage per contig
coverage <- lapply(unique(plat_bed$contig),function(x){
  bed <- plat_bed %>% filter(contig==x) %>%
    select(contig,start,end)
  cov <- apply(bed,MARGIN=1,FUN=function(y){
    return(seq(as.numeric(y[2]),as.numeric(y[3])))
  }) %>% unlist() %>% table() %>% as_tibble() %>%
    rename(pos = '.',coverage='n')
  return(cov)
})
names(coverage) <- unique(plat_bed$contig)
saveRDS(coverage,'rds/plat_asm_coverage.rds')

coverage <- readRDS('rds/plat_asm_coverage.rds')
# Function to plot coverage
plot_cov <- function(cov,window = 10000,name=NULL) {
  require(scales)
  windowed_cov <- cov %>% 
    # Create a custom index of a window to group rows,
    # so I can get the average coverage later on
    mutate(window = ceiling(row_number() / window) * window) %>%
    group_by(window) %>% 
    summarize(coverage = mean(coverage))
  ggplot(data=windowed_cov,aes(x=window,y=coverage))+
    geom_bar(stat='identity',color=NA,fill='darkgrey')+
    scale_x_continuous(labels = scales::comma)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
    ggtitle(paste0('Contig name: ',name,' | Length: ',
                   format(as.numeric(plat_ctg$LN[which(plat_ctg$contig==name)]), big.mark = ",", scientific = FALSE),
                   ' | Plot window: ',format(window,big.mark = ",", scientific = FALSE),'bp'),
            subtitle = paste0('Mean coverage: ',format(round(as.numeric(summary(cov$coverage)[4]),2),nsmall=2),
                              '  |  Median coverage: ',as.numeric(summary(cov$coverage)[3]),
                              '  |  Min coverage: ',as.numeric(summary(cov$coverage)[1]),
                              '  |  Max coverage: ',as.numeric(summary(cov$coverage)[6])))
}
lapply(1:length(coverage),function(i){
  ln <- as.numeric(plat_ctg$LN[which(plat_ctg$contig==names(coverage)[i])])
  plot_cov(coverage[[i]],name = names(coverage)[i],window = ln/1000)
  ggsave(paste0('output/plat_asm/contig_coverage_plots/',names(coverage)[i],'.png'),
         width = 15,height = 6)
})

# Table of coverage info per contig
lapply(1:length(coverage),function(i){
  return(c(contig = names(coverage)[i],
           ln = nrow(coverage[[i]])-1,
           mean_cov = format(round(as.numeric(summary(coverage[[i]]$coverage)[4]),2),nsmall=2),
           median_cov = as.numeric(summary(coverage[[i]]$coverage)[3]),
           min_cov = as.numeric(summary(coverage[[i]]$coverage)[1]),
           max_cov = as.numeric(summary(coverage[[i]]$coverage)[6])))
}) %>% bind_rows() %>% write_csv('output/plat_asm/contig_coverage.csv')

### Busco####
# It was run on bash on the fasta file
# Results saved in ipyrad folder with reference file