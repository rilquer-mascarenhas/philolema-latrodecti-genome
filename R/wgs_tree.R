require(ggtree)
tree <- read.tree('data/wgs/iqtree/plat_vg_nuDNA_filtered.min14.phy.varsites.phy.contree')

# Format tip labels
# Change to metadata binomial + original_code
tree$tip.label <- tree$tip.label %>% as_tibble() %>% rename(code='value') %>% 
  mutate(code = str_replace(code,'\\.q20.rmd.bam','')) %>%
  left_join(wgs_metadata %>% select(code,original_code,binomial)) %>% 
  unite(col = new_label,binomial,original_code,sep=' ') %>%
  mutate(new_label = str_replace_all(new_label,' NA','')) %>%
  select(new_label) %>% unlist() %>% as.character()

# Root tree
# Root with Bruchophagus gibbus
tree <- ape::root(tree,outgroup = 'Bruchophagus gibbus')

ggtree(tree)+
  #geom_nodelab(hjust = -0.2,size=2)+
  #geom_nodelab(aes(label = node), color = "red", hjust = -0.5)+
  #geom_nodepoint()+
  #geom_tiplab()+xlim(0,35)
  # geom_cladelabel(node=34,
  #                 label="San Juan, Panama and Arizona (USA)",
  #                 offset=.2, align=TRUE)+
  # geom_cladelabel(node=40,
  #                 label="San Juan, Panama and Hawaii (USA)", 
  #                 offset=.2, align=TRUE) +
  geom_cladelabel(node=29,
                  label="P. latrodecti", 
                  offset=.2, align=TRUE) + 
  geom_cladelabel(node=46,
                  label="Israel", 
                  offset=.2, align=TRUE) + 
  geom_cladelabel(node=26,
                  label="South Africa", 
                  offset=.2, align=TRUE)+
  geom_tiplab()+
  theme_tree2() + 
  xlim(0, 0.5) + 
  theme_tree()