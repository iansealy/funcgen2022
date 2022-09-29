library(tidyverse)

# load data
deseq_results <- 
  read_tsv('Amp.counts.tsv',
           col_types = cols(Chr = col_character(), 
                            Strand = col_character()))

# filter to sig genes
sig_data <- 
  filter(deseq_results, adjp < 0.05)
# select normalised count columns
norm_counts <- 
  select(sig_data, contains("normalised"))

# calculate correlation matrix
# the cor function works on columns so we have to transpose the matrix with the t() function
cor_matrix <- cor(t(norm_counts))

# add gene ids and make a data frame
colnames(cor_matrix) <- sig_data$Gene
cor_df <- tibble(
  Gene = sig_data$Gene,
  as.data.frame(cor_matrix)
)

# pivot the matrix to make it tidy
cor_long <- cor_df %>% 
  pivot_longer(cols = -Gene, names_to = 'Gene2', values_to = 'cor') %>% 
  # sort by Gene
  arrange(Gene, Gene2) %>% 
  # remove where the gene is the same or higher
  # also only keep ones where cor > 0.8
  filter(!(Gene == Gene2 | Gene > Gene2), cor > 0.8) %>% 
  # add a type column that is always "geneExprCor"
  mutate(type = "geneExprCor") %>% 
  # reorder the columns
  select(Gene, type, Gene2, cor)
# this is now the edges in the source -> target format

# write out network file
# change this filename if necessary
network_file <- 'Amp.sif'
write_tsv(select(cor_long, Gene, type, Gene2), 
          file = network_file, col_names = FALSE)

# get node info
node_info <- select(sig_data, 
                    Name = Gene, 
                    GeneName = Name, 
                    adjp, log2fc)

# write out nodes file
# change this filename if necessary
nodes_file <- 'Amp.nodes.txt'
write_tsv(node_info, file = nodes_file)

# make edges df
unite(cor_long, 'edge', Gene, type, sep = ' (') %>% 
  unite(col = 'edge', edge, Gene2, sep = ') ') %>% 
  # write out edges file
  write_tsv(., file = 'Amp.edges.txt')

