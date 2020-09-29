rm(list = ls())

## Let's make a plot of WZA results from Lodgepole Pine DD0 data

pine <- read.csv("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/DD_0_Pine.csv")

library(ggplot2)

hits_to_highlight <- pine[(-log10(pine$top_candidate_p) > 25)&(pine$Z<10),]$contig

max_Z <- pine[pine$Z == max(pine$Z),]$contig

### Read in the raw p-values...
#raw <- read.csv("~/UBC/Convergence/Yeaman_Data/doi_10.5061_dryad.0t407__v1/yeaman_pine_spruce_convergence_archive2/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_raw", sep = "\t")
raw_1 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/scaffold750644_27863_31476.tsv", header = T, comment.char = "&")
raw_2 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/tscaffold1962_233534_237733.tsv", header = T, comment.char = "&")
raw_3 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/scaffold878015.2_134352_141334.tsv", header = T, comment.char = "&")
raw_4 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/tscaffold3400_439096_449060.tsv", header = T, comment.char = "&")
x = 10


corre <-cor( -log10(pine$top_candidate_p), y = pine$Z, method = "spearman")

cor_text = bquote("Spearman's" ~ rho ~ "=" ~ .(corre))

ggplot( data = pine, aes(x = -log10(top_candidate_p), y = Z))+
  geom_point(alpha = 0.5)+
  scale_y_continuous(expression(Z[W]))+
  scale_x_continuous(expression(-log[10]*"(Top-Candidate Test "*italic("p")*"-value)"))+
  annotate("text", label = cor_text, x = 60, y = 0)+
  theme_bw()



ggplot( data = raw_1, aes(x = pos_gcontig, y = -log10(DD_0)))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.0004004730154822944)), lty = 2, alpha = 0.6)+
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous("Position in Contig")

ggplot( data = raw_2, aes(x = pos_gcontig, y = -log10(DD_0)))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.0004004730154822944)), lty = 2, alpha = 0.6)+
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous("Position in Contig")

ggplot( data = raw_3, aes(x = pos_gcontig, y = -log10(DD_0)))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.0004004730154822944)), lty = 2, alpha = 0.6)+
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous("Position in Contig")

ggplot( data = raw_4, aes(x = pos_gcontig, y = -log10(DD_0)))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.0004004730154822944)), lty = 2, alpha = 0.6)+
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous("Position in Contig")
