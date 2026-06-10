setwd('/home/rilquer/lgeo_plat/')
require(tidyverse)

# PHILOLEMA WGS####

wgs_metadata <- read_csv('data/sample_data/york_master_samplesheet.csv') %>%
  filter(wgs=='yes')

## Alignment stats ####
### Linear alignment####
samples <- list.files('data/plat_wgs/linear_aln/cov_files/',pattern = '1.cov') %>% 
  str_split_i(pattern = '\\.',i = 1)

cov <- lapply(1:length(samples),function(i){
  cov <- read_tsv(paste0('data/plat_wgs/linear_aln/cov_files/',samples[i],'.cov')) %>% 
    mutate(sample = samples[i]) %>%
    rename(contig = '#rname',length = 'endpos') %>%
    select(-startpos)
  return(cov)
}) %>% bind_rows() %>%
  # Adding metadata
  left_join(wgs_metadata %>% select(code,binomial,locality,region1,region2),
            by = c(sample = 'code'))

# Average depth per sample
cov %>% group_by(sample) %>%
  summarize(mean_depth = mean(meandepth),
            median_depth = median(meandepth),
            max_depth = max(meandepth),
            min_depth = min(meandepth)) %>%
  # Adding metadata
  left_join(wgs_metadata %>% select(code,binomial,locality,region1,region2),
            by = c(sample = 'code')) %>% 
  print(n=Inf) %>%
  write_csv('output/wgs/cov_stats/linear_avg_depth_per_sample.csv')

# Average depth per contig
cov %>% group_by(contig) %>%
  summarize(mean_depth = mean(meandepth),
            median_depth = median(meandepth),
            max_depth = max(meandepth),
            min_depth = min(meandepth)) %>% 
  print(n=Inf) %>% 
  write_csv('output/wgs/cov_stats/linear_avg_depth_per_contig.csv')

# Plots
# Mean Depth
ggplot(data=cov,aes(x=contig,y=locality,fill=meandepth))+
  geom_tile()+
  scale_fill_gradient2(name = 'Mean Depth',
                       low = '#2b83ba', mid = '#fffebd', high = '#d7191c',
                       midpoint = median(cov$meandepth))+
  theme_classic()+
  theme(panel.spacing.y = unit(0, "mm"),
        strip.text.y.left = element_text(angle = 0,face='italic'),
        strip.background = element_rect(fill = "gray90", color = "white",
                                        linewidth = 1.5),
        strip.placement = "outside",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=8))+
  facet_grid(binomial~., scales = "free_y", space = "free_y", switch = "y") +
  ggtitle('Linear Reference Alignment - Depth',subtitle = paste0('Mean: ',mean(cov$meandepth),'; Median: ',median(cov$meandepth)))
ggsave('output/wgs/cov_stats/linear_mean_depth.png',width = 13, height = 6,dpi = 800)

# Coverage
ggplot(data=cov,aes(x=contig,y=locality,fill=coverage))+
  geom_tile()+
  scale_fill_gradient2(name = '% Coverage',
                       low = '#2b83ba', mid = '#fffebd', high = '#d7191c',
                       midpoint = 50)+
  theme_classic()+
  theme(panel.spacing.y = unit(0, "mm"),
        strip.text.y.left = element_text(angle = 0,face='italic'),
        strip.background = element_rect(fill = "gray90", color = "white",
                                        linewidth = 1.5),
        strip.placement = "outside",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=8))+
  facet_grid(binomial~., scales = "free_y", space = "free_y", switch = "y") +
  ggtitle('Linear Reference Alignment - Coverage',subtitle = paste0('Mean: ',mean(cov$coverage),'; Median: ',median(cov$coverage)))
ggsave('output/wgs/cov_stats/linear_mean_coverage.png',width = 13, height = 6,dpi = 800)

# MapQ
ggplot(data=cov,aes(x=contig,y=locality,fill=meanmapq))+
  geom_tile()+
  scale_fill_gradient2(name = 'Mean MapQ',
                       low = '#2b83ba', mid = '#fffebd', high = '#d7191c',
                       midpoint = median(cov$meanmapq))+
    theme_classic()+
  theme(panel.spacing.y = unit(0, "mm"),
        strip.text.y.left = element_text(angle = 0,face='italic'),
        strip.background = element_rect(fill = "gray90", color = "white",
                                        linewidth = 1.5),
        strip.placement = "outside",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=8))+
  facet_grid(binomial~., scales = "free_y", space = "free_y", switch = "y") +
  ggtitle('Linear Reference Alignment - Mapping quality',subtitle = paste0('Mean: ',mean(cov$meanmapq),'; Median: ',median(cov$meanmapq)))
ggsave('output/wgs/cov_stats/linear_mean_mapq.png',width = 13, height = 6,dpi = 800)

# Stats per contig across all samples
plat_cov <- read_tsv('data/plat_wgs/linear_aln/cov_files/plat_wgs_all.cov')
mean(plat_cov$meandepth)
median(plat_cov$meandepth)

# Length of contiguous against depth
ggplot(plat_cov,aes(x=endpos,y=meandepth))+
  geom_point(alpha=0.7)+
  geom_hline(yintercept = 20,linetype='dashed',color='red')

# Length of contig against mean mapping quality
ggplot(plat_cov,aes(x=endpos,y=meanmapq))+
  geom_point(alpha=0.7)
  geom_hline(yintercept = 15,linetype='dashed',color='red')

### Graph alignment####
samples <- list.files('data/plat_wgs/graph_aln/cov_files/',pattern = '1.q20.cov') %>% 
    str_split_i(pattern = '\\.',i = 1) %>%
    # Adding outgroup samples
    append(c('SRR37165114','SRR37165121','SRR37165125','SRR37165136'))

# Filtered for mapq20
cov <- lapply(1:length(samples),function(i){
  cov <- read_tsv(paste0('data/plat_wgs/graph_aln/cov_files/',samples[i],'.q20.cov'))%>% 
    mutate(sample = samples[i]) %>%
    rename(contig = '#rname',length = 'endpos') %>%
    select(-startpos)
  return(cov)
}) %>% bind_rows() %>% 
  # Adding metadata
  left_join(wgs_metadata %>% select(code,binomial,locality,region1,region2),
            by = c(sample = 'code'))

# Average depth per sample
cov %>% group_by(sample) %>%
  summarize(mean_depth = mean(meandepth),
            median_depth = median(meandepth),
            max_depth = max(meandepth),
            min_depth = min(meandepth)) %>%
  # Adding metadata
  left_join(wgs_metadata %>% select(code,binomial,locality,region1,region2),
            by = c(sample = 'code')) %>% 
  print(n=Inf) %>%
  write_csv('output/wgs/cov_stats/graph_avg_depth_per_sample.csv')

# Average depth per contig
cov %>% group_by(contig) %>%
  summarize(mean_depth = mean(meandepth),
            median_depth = median(meandepth),
            max_depth = max(meandepth),
            min_depth = min(meandepth)) %>% 
  print(n=Inf) %>% 
  write_csv('output/wgs/cov_stats/graph_avg_depth_per_contig.csv')

# Plots
# Mean Depth
ggplot(data=cov,aes(x=contig,y=locality,fill=meandepth))+
  geom_tile()+
  scale_fill_gradient2(name = 'Mean Depth',
                       low = '#2b83ba', mid = '#fffebd', high = '#d7191c',
                       midpoint = median(cov$meandepth))+
  theme_classic()+
  theme(panel.spacing.y = unit(0, "mm"),
        strip.text.y.left = element_text(angle = 0,face='italic'),
        strip.background = element_rect(fill = "gray90", color = "white",
                                        linewidth = 1.5),
        strip.placement = "outside",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=8))+
  facet_grid(binomial~., scales = "free_y", space = "free_y", switch = "y") +
  ggtitle('Graph Reference Alignment - Depth',subtitle = paste0('Mean: ',mean(cov$meandepth),'; Median: ',median(cov$meandepth)))
ggsave('output/wgs/cov_stats/graph_mean_depth.png',width = 13, height = 6,dpi = 800)

# Coverage
ggplot(data=cov,aes(x=contig,y=locality,fill=coverage))+
  geom_tile()+
  scale_fill_gradient2(name = '% Coverage',
                       low = '#2b83ba', mid = '#fffebd', high = '#d7191c',
                       midpoint = 50)+
  theme_classic()+
  theme(panel.spacing.y = unit(0, "mm"),
        strip.text.y.left = element_text(angle = 0,face='italic'),
        strip.background = element_rect(fill = "gray90", color = "white",
                                        linewidth = 1.5),
        strip.placement = "outside",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=8))+
  facet_grid(binomial~., scales = "free_y", space = "free_y", switch = "y") +
  ggtitle('Graph Reference Alignment - Coverage',subtitle = paste0('Mean: ',mean(cov$coverage),'; Median: ',median(cov$coverage)))
ggsave('output/wgs/cov_stats/graph_mean_coverage.png',width = 13, height = 6,dpi = 800)

# MapQ
ggplot(data=cov,aes(x=contig,y=locality,fill=meanmapq))+
  geom_tile()+
  scale_fill_gradient2(name = 'Mean MapQ',
                       low = '#2b83ba', mid = '#fffebd', high = '#d7191c',
                       midpoint = max(cov$meanmapq)/2)+
  theme_classic()+
  theme(panel.spacing.y = unit(0, "mm"),
        strip.text.y.left = element_text(angle = 0,face='italic'),
        strip.background = element_rect(fill = "gray90", color = "white",
                                        linewidth = 1.5),
        strip.placement = "outside",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=8))+
  facet_grid(binomial~., scales = "free_y", space = "free_y", switch = "y") +
  ggtitle('Graph Reference Alignment - Mapping quality',subtitle = paste0('Mean: ',mean(cov$meanmapq),'; Median: ',median(cov$meanmapq)))
ggsave('output/wgs/cov_stats/graph_mean_mapq.png',width = 13, height = 6,dpi = 800)

# Stats per contig across all samples
plat_cov <- read_tsv('data/plat_wgs/linear_aln/cov_files/plat_wgs_all.cov')
mean(plat_cov$meandepth)
median(plat_cov$meandepth)

# Length of contiguous against depth
ggplot(plat_cov,aes(x=endpos,y=meandepth))+
  geom_point(alpha=0.7)+
  geom_hline(yintercept = 20,linetype='dashed',color='red')

# Length of contig against mean mapping quality
ggplot(plat_cov,aes(x=endpos,y=meanmapq))+
  geom_point(alpha=0.7)
geom_hline(yintercept = 15,linetype='dashed',color='red')
## SNPs info from vcftools ####
library(tidyverse)
library(cowplot)

##Set ggplot theme
theme_set(theme_bw())

##VARIANT QUALITY
var_qual <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.lqual",
                       delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
summary(var_qual$qual)

#We can see that most sites have quite high quality scores [~1000]
#Let's set a minimum quality filter of 100 and filter more strongly on other features

##VARIANT DEPTH##
var_depth <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.ldepth.mean",
                        delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
summary(var_depth$mean_depth)

#Since we all took different subsets, these values will likely differ slightly but clearly in 
#this case most variants have a depth of 17-20x whereas there are some extreme outliers. 
#We will redraw our plot to exclude these and get a better idea of the distribution of mean depth.

b + theme_light() + xlim(0, 10)


#This gives a better idea of the distribution. We could set our minimum coverage at the 
#5 and 95% quantiles but we should keep in mind that the more reads that cover a site, 
#the higher confidence our basecall is. 10x is a good rule of thumb as a minimum cutoff 
#for read depth, although if we wanted to be conservative, we could go with 15x.

#What is more important here is that we set a good maximum depth cufoff. As the outliers show, 
#some regions clearly have extremely high coverage and this likely reflects mapping/assembly 
#errors and also paralogous or repetitive regions. We want to exclude these as they will bias 
#our analyses. Usually a good rule of thumb is something the mean depth x 2 - so in this case 
#we could set our maximum depth at 40x.

#So we will set our minimum depth to 10x and our maximum depth to 40x.

##VARIANT MISSINGNESS##

var_miss <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.lmiss",
                       delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#Our cichlid data has a very promising missingness profile - clearly most individuals have a 
#call at almost every site. Indeed if we look at the summary of the data we can see this even 
#more clearly.

summary(var_miss$fmiss)

#Most sites have almost no missing data. Although clearly, there are some (as the max value shows). 
#This means we can be quite conservative when we set our missing data threshold. We will remove 
#all sites where over 10% of individuals are missing a genotype. One thing to note here is that 
#vcftools inverts the direction of missigness, so our 10% threshold means we will tolerate 90% 
#missingness (yes this is confusing and counterintuitive… but that’s the way it is!). Typically 
#missingness of 75-95% is used.

##Minor Allele Frequency ####

#Last of all for our per variant analyses, we will take a look at the distribution of allele 
#frequencies. This will help inform our minor-allele frequency (MAF) thresholds. As previously, 
#we read in the data:

var_freq <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.frq",
                       delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

#However, this is simply the allele frequencies. To find the minor allele frequency at each site, 
#we need to use a bit of dplyr based code.

# find minor allele frequency
var_freq$maf <- var_freq |> select(a1, a2) |> apply(1, function(z) min(z))

#Here we used apply on our allele frequencies to return the lowest allele frequency at each variant. 
#We then added these to our dataframe as the variable maf. Next we will plot the distribution.

d <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
summary(var_freq$maf)
# Wasp:
# Set to 0.1

# Dan:
#The upper bound of the distribution is 0.5, which makes sense because if MAF was more than this, 
#it wouldn’t be the MAF! How do we interpret MAF? It is an important measure because low MAF alleles 
#may only occur in one or two individuals. It is possible that some of these low frequency alleles 
#are in fact unreliable base calls - i.e. a source of error.

#With 16 individuals, there are 28 alleles for a given site. Therefore MAF = 0.04 is equivalent to 
#a variant occurring as one allele in a single individual (i.e. 28 * 0.04 = 1.12). Alternatively, 
#an MAF of 0.1 would mean that any allele would need to occur at least twice (i.e. 28 * 0.1 = 2.8).

#Setting MAF cutoffs is actually not that easy or straightforward. Hard MAF filtering (i.e. setting 
#a high threshold) can severely bias estimation of the site frequency spectrum and cause problems 
#with demographic analyses. Similarly, an excesss of low frequency, ‘singleton’ SNPs (i.e. only 
#occurring in one individual) can mean you keep many uniformative loci in your dataset that make 
#it hard to model things like population structure.

#Usually then, it is best practice to produce one dataset with a good MAF threshold and keep 
#another without any MAF filtering at all. For now however, we will set our MAF to 0.1

#Let's use program cowplot function plot_grid to create a single summary figure for these 
#five site based summaries

plot_grid(a, b, c, d, rows = 2)

## Individual based statistics ####

#As well as a our per variant statistics we generated earlier, we also calculated some individual 
#metrics too. We can look at the distribution of these to get an idea whether some of our individuals 
#have not sequenced or mapped as well as others. This is good practice to do with a new dataset. 
#A lot of these statistics can be compared to other measures generated from the data (i.e. principal 
#components as a measure of population structure) to see if they drive any apparent patterns in the data.

## Mean depth per individual ####

ind_depth <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.idepth",
                        delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)

e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "#fb6a4a", colour = "black", alpha = 0.3)
summary(ind_depth$depth)
# Wasp interpretation:
# Depth is overall low, btween 3 and 5, with one individual on 16x

#Because we are only plotting data for 16 individuals, the plot looks a little disjointed. While 
#there is some evidence that some individuals were sequenced to a higher depth than others, 
#there are no extreme outliers. So this doesn’t suggest any issue with individual sequencing depth.

## Missing data per individual ####

ind_miss  <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.imiss",
                        delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "#fb6a4a", colour = "black", alpha = 0.3)

# Wasp interpretation:
# Some few individuals have lot of missing data in the variant sites

#Again this shows us, the proportion of missing data per individual is very small indeed. It ranges 
#from 0.01-0.16, so we can safely say our individuals sequenced well.


## Het and inbreeding per ind ####

ind_het <- read_delim("Rilquer/york_postdoc/platrodecti/wgs/snp_stats/plat_wgs.het",
                      delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "#fb6a4a", colour = "black", alpha = 0.3)

# Wasp interpretation:
# Extra high inbreeding. Could this be due to lack of depth and underestimation of heterozygosity?

#All individuals have a slightly negative inbreeding coefficient suggesting that we observed a bit 
#less heterozygote genotypes in these individuals than we would expect under Hardy-Weinberg 
#equilibrium. However, here we combined samples from four species and thus violate the assumption 
#of Hardy-Weinberg equilibrium. We would expect slightly negative inbreeding coefficients due to 
#the Wahlund-effect. Given that all individuals seem to show similar inbreeding coefficients, 
#we are happy to keep all of them. None of them shows high levels of allelic dropout (strongly 
#negative F) or DNA contamination (highly positive F).

#Let's use program cowplot function plot_grid to create a single summary figure for these 
#eight sites indivdiual based summaries

plot_grid(a, b, c, d, e, f, g, rows = 2)

## WGS tree ####
# RAXML trees were run in Unity
require(ggtree) # Issues with ggtree installation