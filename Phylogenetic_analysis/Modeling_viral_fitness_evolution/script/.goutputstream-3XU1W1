#!/usr/bin/env R

library(tidyverse)
library("ape")
library("ggplot2")
library("ggtree")
library(castor)
library(scales)
library(data.table)
library(cmdstanr)
library(tictoc)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(makedummies)
library(ggnewscale)
library(phangorn)
library(viridis)


args = commandArgs(trailingOnly=T)

##########args##########

transition.node.info.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/tree/Omicron_trees.parsimony.transition_branch_info.ba5.txt'
transition.node.info <- read.table(transition.node.info.name,header=T)

mut.gain.info.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/mutation_regression/2022_11_07/test_2022_11_07.3.method3.growth_gain.txt'
mut.gain.info <- read.table(mut.gain.info.name,header=T)

mut.info.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/GISAID_data/2022_10_15/metadata.mut_long.Spike.tsv'
mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)
mut.info <- mut.info %>% mutate(mut = gsub("Spike_","",mut), mut.mod = gsub("[A-Z]$","",mut))

ref.hap.info.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/mutation_regression/2022_11_07/test_2022_11_07.method3.ref_hap.txt'
ref.hap.info <- read.table(ref.hap.info.name,header=T)
ref.hap.info <- ref.hap.info %>% mutate(mut = gsub("Spike_","",mut))


lineage.interest <- "ba5"

tree.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/tree/nhx_tree/ba5.filt_annot.nwk'
tree <- read.tree(tree.name)


tip.df <- data.frame(tip_Id = 1:length(tree$tip), Id = tree$tip)


tree.info.df <- ggtree(tree)$data



metadata.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/GISAID_data/2022_10_15/metadata.tsv'
metadata <- fread(metadata.name,header=T,sep="\t",check.names=T)
metadata <- metadata %>% filter(Accession.ID %in% as.character(tip.df$Id))

max.date <- max(as.Date(metadata$Collection.date))

mut.info.interest <- mut.info %>% filter(Id %in% as.character(tip.df$Id))


mut.group.df <- data.frame()

for(i in 1:nrow(mut.gain.info)){

mut.group <- mut.gain.info$mut[i]
mut.v <- strsplit(mut.group, ":")[[1]]

temp.df <- data.frame(mut.group = mut.group, mut = mut.v)


mut.group.df <- rbind(mut.group.df,temp.df)

}


mut.info.merged <- mut.info.interest %>% inner_join(mut.group.df,by="mut")

mut.info.merged <- mut.info.merged %>% distinct(mut.group,Id,.keep_all=T)

mut.info.mat <- mut.info.merged %>% mutate(value = 1) %>% select(-mut,-mut.mod) %>% spread(key=mut.group,value=value)
mut.info.mat[is.na(mut.info.mat)] <- 0

mut.info.mat <- mut.info.mat %>% as.data.frame()
rownames(mut.info.mat) <- mut.info.mat$Id

mut.info.mat <- mut.info.mat %>% select(-Id)

mut.info.mat <- mut.info.mat[as.character(tip.df$Id),]

mut.gain.info.filtered <- mut.gain.info %>% filter(mut %in% colnames(mut.info.mat))

#mut.gain.info.filtered <- mut.gain.info.filtered %>% inner_join(ref.hap.info,by="mut")
#mut.gain.info.filtered <- mut.gain.info.filtered %>% mutate(mean.log = log(mean),mean.log.directed = ifelse(hap == 0, as.numeric(mean.log),-as.numeric(mean.log)))


mut.info.mat <- mut.info.mat %>% select(all_of(mut.gain.info.filtered$mut))

sum.ancestral.l <- list()


for(i in 1:ncol(mut.info.mat)){

mut.group.interest <- colnames(mut.info.mat)[i]

state.v <- as.numeric(mut.info.mat[,i]) + 1


fit.asr <- asr_max_parsimony(tree, state.v, Nstates=NULL, 
             transition_costs="all_equal",
             edge_exponent=0, weight_by_scenarios=TRUE,
             check_input=TRUE)


#fit.asr <- asr_mk_model(tree, state.v, 
#              Nstates = NULL, 
#              tip_priors = NULL, 
#              rate_model = "ER", 
#              transition_matrix = NULL, 
#              include_ancestral_likelihoods = TRUE, 
#              reroot = TRUE,
#              root_prior = "empirical", 
#              Ntrials = 1, 
#              optim_algorithm = "nlminb",
#              optim_max_iterations = 200,
#              optim_rel_tol = 1e-8,
#              store_exponentials = TRUE,
#              check_input =TRUE, 
#              Nthreads = 1)

ancstats <- as.data.frame(fit.asr$ancestral_likelihoods)
colnames(ancstats) <- c('absence','presence')
ancstats$node <- 1:tree$Nnode + Ntip(tree)

sum.ancestral.l[[mut.group.interest]] <- c(as.numeric(mut.info.mat[,i]),ancstats$presence)

}

Unzip <- function(...) rbind(data.frame(), ...)
sum.ancestral.df <- as.data.frame(t(do.call(Unzip, sum.ancestral.l)))
colnames(sum.ancestral.df) <- colnames(mut.info.mat)
rownames(sum.ancestral.df) <- 1:nrow(sum.ancestral.df)

mut.gain.info.filtered <- mut.gain.info %>% filter(mut %in% colnames(sum.ancestral.df))

sum.ancestral.df <- sum.ancestral.df %>% select(all_of(mut.gain.info.filtered$mut))

ref.hap.info.filtered <- ref.hap.info %>% filter(mut %in% colnames(sum.ancestral.df))
rownames(ref.hap.info.filtered) <- ref.hap.info.filtered$mut

ref.hap.info.filtered <- ref.hap.info.filtered[colnames(sum.ancestral.df),]

sum.ancestral.df <- t(t(sum.ancestral.df) - ref.hap.info.filtered$hap)

mut.effect.df <- t(t(sum.ancestral.df) * log(mut.gain.info.filtered$mean)) %>% as.data.frame()



predict.df.Re <- data.frame(node = tree.info.df$node, Re = exp(apply(mut.effect.df,1,sum)))



mut.info.spread <- mut.info.interest %>% select(Id,mut.mod) %>% mutate(value = 1) %>% spread(key=mut.mod,value =value)



tip.df.merged <- tip.df %>% left_join(mut.info.spread,by="Id")
tip.df.merged[is.na(tip.df.merged)] <- 0

mut.count.df <- data.frame(Id = tip.df.merged$Id, mut.count = apply(tip.df.merged[,c("F486","K444","L452","N460","R346")],1,sum))

color.mut.v <- brewer.pal(5, "Set1")
names(color.mut.v) <- c("R346","K444","L452","N460","F486")

print("tree")

col.v <- viridis(100)



n.seq <- nrow(tip.df)




node.df <- transition.node.info %>% select(node,mut)


plot.df <- tip.df.merged %>% arrange(tip_Id) %>% select(-tip_Id,-Id)
rownames(plot.df) <- tip.df.merged$Id
plot.df <- ifelse(plot.df == 1,"+","-")
plot.df <- plot.df[,c("R346","K444","L452","N460","F486")]

#plot.df <- plot.df[tip.df2$Id,]

g <- ggtree(tree,size=0.5,mrsd = max.date) %>% flip(5502,5531)
g <- g %<+% predict.df.Re + aes(color=Re) 
g <- g + scale_colour_gradientn(colors=c("gray70","red","orange","yellow"),breaks=c(1.25,1.4,1.55),limits=c(1.25,1.55),oob=squish) 
g <- g + new_scale_color()
g <- g %<+% node.df + aes(color=mut) 
g <- g + geom_point2(size=2)
g <- g + scale_color_manual(values=color.mut.v,breaks=names(color.mut.v),na.value = NA) 
g <- g + theme_tree2()
g <- g + scale_x_continuous(limits=c(2021,2024.25),breaks=seq(2021,2023,0.25))
g <- g + ylim(0,n.seq + 700)
g1 <- gheatmap(g, plot.df, offset = 0.2, width = 0.3,
               colnames_position = "top",
               colnames_offset_y = 0,colnames_angle=90,hjust=0,color=NA)
g1 <- g1 + scale_fill_manual(values=c('gray90','navy'),breaks=c("-","+"))
g1 <- g1 + ggtitle(lineage.interest)


pdf.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/tree/BA5_tree.Re.pdf'
pdf(pdf.name,width=8,height=5)
g1
dev.off()



#subtree
subtree.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/tree/ba5sub4_filt_bq11.nwk'
subtree <- read.tree(subtree.name)

subtree.tip.df <- data.frame(tip_Id = 1:length(subtree$tip), Id = subtree$tip)

mrca.subtree <- getMRCA(tree, subtree.tip.df$Id)

mrca.subtree.parent <- Ancestors(tree, mrca.subtree, type = c("parent"))

#g1 <- g1 + geom_text(aes(label=node), hjust=-.3,color="black")
g2 <- viewClade(g1,node=mrca.subtree.parent)
g2


pdf.name <- '/media/jampei/C6700DA7700D9F75/data/Sato_analysis/convergent_evolution/tree/BQ11_tree.Re.pdf'
pdf(pdf.name,width=8,height=5)
g2
dev.off()










