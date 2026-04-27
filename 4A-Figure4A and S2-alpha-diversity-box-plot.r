library(phyloseq)
library(RColorBrewer)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(mia)
library(ggpubr)
library(knitr)
library(dplyr)



otu <- read.table("psmin5_ASV_table.txt", sep = "\t", header = TRUE, row.names = 1)
taxon <- read.table("psmin5_taxa_table.txt", sep = "\t", header = TRUE, row.names = 1)
samples <- read.table("metadata.txt", sep = "\t", header = TRUE, row.names = 1)

taxon <- as.matrix(taxon)
TAX <- tax_table(taxon)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  sample_data(samples),
  tax_table(TAX))
  
ps

to_remove <- c("LP196","eDNA1")

ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)

ps

#The function below creates a table with selected (or all) diversity indicators.

alpha_div <- microbiome::alpha(ps, index = "all")
alpha_div_table <- kable(alpha_div)
write.table(alpha_div_table, file = "alpha_diversity_results.txt")


# Get the metadata from the `phyloseq` object

ps.meta <- meta(ps)



#Add the diversity indices to metadata

ps.meta$Shannon <- alpha_div$diversity_shannon 
ps.meta$InverseSimpson <-  alpha_div$diversity_inverse_simpson
ps.meta$Richness <- alpha_div$observed 


#order your groups for plotting 

ps.meta$cat <- factor(ps.meta$cat,levels=c("L-H","T-H","L-D","T-D","L-W","T-W","L-S","T-S"))

#define colors 
mycolors <- c(
  "L-H" = "darkolivegreen1",
  "T-H" = "darkolivegreen3",
  "L-D" = "brown1",
  "T-D" = "brown3",
  "L-W" = "cadetblue1",
  "T-W" = "cadetblue3",
  "L-S" = "magenta2",
  "T-S" = "magenta4")


#plot shannon 

shannon = ggplot(ps.meta,aes(x = cat, y = Shannon, color = cat)) +
geom_boxplot(outlier.shape = NA)+
labs(x = "", y = "Shannon Diversity") +
theme_bw()+
theme(panel.grid = element_blank())+
theme(axis.text.x=element_text(size=10))+
theme(legend.text = element_text(size=10))+
theme(legend.title=element_text(size=12))+
scale_colour_manual(values=mycolors,name="Site-type",
                         breaks=c("L-H","T-H","L-D","T-D","L-W","T-W","L-S","T-S"),
                         labels=c("L-H:Luminao-healthy","T-H:Tumon-healthy","L-D:Luminao-diseased","T-D:Tumon-diseased","L-W:Luminao-seawater","T-W:Tumon-seawater","L-S:Luminao-sediment","T-S:Tumon-sediment"))


shannon
ggsave("Figure3A-shannon_diversity_boxplot.tiff", units="in", width=6, height=3, dpi=300, compression = 'lzw')


richness = ggplot(ps.meta,aes(x = cat, y = Richness, color = cat)) +
geom_boxplot(outlier.shape = NA)+
labs(x = "", y = "ASV Richness") +
theme_bw()+
theme(panel.grid = element_blank())+
theme(axis.text.x=element_text(size=10))+
theme(legend.text = element_text(size=10))+
theme(legend.title=element_text(size=12))+
scale_colour_manual(values=mycolors,name="Site-type",
                         breaks=c("L-H","T-H","L-D","T-D","L-W","T-W","L-S","T-S"),
                         labels=c("L-H:Luminao-healthy","T-H:Tumon-healthy","L-D:Luminao-diseased","T-D:Tumon-diseased","L-W:Luminao-seawater","T-W:Tumon-seawater","L-S:Luminao-sediment","T-S:Tumon-sediment"))


richness
ggsave("Figure2A-Richiness_boxplot.tiff", units="in", width=6, height=3, dpi=300, compression = 'lzw')



simpson = ggplot(ps.meta,aes(x = cat, y = InverseSimpson, color = cat)) +
geom_boxplot(outlier.shape = NA)+
labs(x = "", y = "Inverse Simpson Diversity") +
theme_bw()+
theme(panel.grid = element_blank())+
theme(axis.text.x=element_text(size=10))+
theme(legend.text = element_text(size=10))+
theme(legend.title=element_text(size=12))+
scale_colour_manual(values=mycolors,name="Site-type",
                         breaks=c("L-H","T-H","L-D","T-D","L-W","T-W","L-S","T-S"),
                         labels=c("L-H:Luminao-healthy","T-H:Tumon-healthy","L-D:Luminao-diseased","T-D:Tumon-diseased","L-W:Luminao-seawater","T-W:Tumon-seawater","L-S:Luminao-sediment","T-S:Tumon-sediment"))


simpson
ggsave("FigureS2B-in-Simpson_diversity_boxplot.tiff", units="in", width=6, height=3, dpi=300, compression = 'lzw')
