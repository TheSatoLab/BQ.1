#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(castor)
library(phangorn)
library(patchwork)
library(RColorBrewer)
library(ggnewscale)
 
args = commandArgs(trailingOnly=T)


#parameters
min.num.decendants <- 3
min.prop.decendants <- 0.7


#inputs
mut.interest.info.name <- "script/mut.interest.info.txt"
mut.interest.info <- read.table(mut.interest.info.name,header=T)

metadata.name <- args[1]
metadata <- fread(metadata.name,header=T,sep="\t",check.names=T)

mut.info.name <- args[2]
mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)


#preprocessing
mut.info <- mut.info %>% mutate(mut = gsub("Spike_","",mut), mut.mod = gsub("[A-Z]$","",mut))

mut.info.spread <- mut.info %>% select(Id,mut.mod) %>% mutate(value = 1) %>% spread(key=mut.mod,value=value)
mut.info.spread <- mut.info.spread %>% select(Id,F486,K444,L452,N460,R346)


#detect transition branches
plot.l <- list()
for(lineage.interest in c("BA.1","BA.2","BA.4","BA.5")){

  #read tree and preprocessing
  tree.name <- paste('../Trees/',lineage.interest,'.nwk',sep="")
  tree <- read.tree(tree.name)

  tree.info.df <- ggtree(tree)$data

  count_descendants <- function(node){
    n.descendants <- length(Descendants(tree, node, type = "tips")[[1]])
    return(n.descendants)
  }

  tree.info.df$num.decendants <- as.numeric(map(tree.info.df$node,count_descendants))

  tip.df <- data.frame(tip_Id = 1:length(tree$tip), Id = tree$tip)
  tip.df.merged <- tip.df %>% left_join(mut.info.spread,by="Id")
  tip.df.merged[is.na(tip.df.merged)] <- 0

  metadata.interest <- metadata %>% filter(Accession.ID %in% as.character(tip.df$Id))

  max.date <- max(as.Date(metadata.interest$Collection.date))
  max.x <- max(tree.info.df$x)

  mut.interest.v <- mut.interest.info %>% filter(lineage == lineage.interest) %>% pull(mut)

  mut.info.mat <- tip.df.merged[,mut.interest.v]


  #detect transition branches
  node.transition.df <- data.frame()
  branch_state.l <- list()

  for(mut.name in mut.interest.v){
    #ancestral state reconstruction
    state.v <- mut.info.mat %>% pull(mut.name)

    fit.asr <- asr_max_parsimony(tree, state.v + 1, Nstates=NULL, 
               transition_costs="all_equal",
               edge_exponent=0, weight_by_scenarios=TRUE,
               check_input=TRUE)

    state.mat <- fit.asr$ancestral_likelihoods
    state.node.v <- state.mat[,1]
    state.node.v <- ifelse(state.node.v>0.5,0,1)

    state.color <- c(state.v,state.node.v)

    branch_state.l[[mut.name]] <- state.color

    #extract transition branch
    tree.info.df.interest <- tree.info.df %>% select(parent,node,num.decendants) %>% mutate(state = state.color)

    tree.info.df.interest.parent <- tree.info.df.interest %>% select(node,state) %>% rename(parent = node, state.parent = state)
    tree.info.df.interest.merged <- merge(tree.info.df.interest,tree.info.df.interest.parent,by="parent")
    tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged %>% filter(state == 1, state.parent == 0)


    transition.node.info.df <- data.frame()
    for(i in 1:nrow(tree.info.df.interest.merged.0to1)){

      node.interest <- tree.info.df.interest.merged.0to1$node[i]
      date.interest <- max.date - round(365 * (max.x - tree.info.df$x[node.interest]))

      tip.v <- Descendants(tree, node.interest, type = "tips")[[1]]

      tip.df.with_mut <- tip.df.merged %>% filter(get({{mut.name}}) == 1, tip_Id %in% tip.v)

      num.tip.with_mut <-  nrow(tip.df.with_mut)

      mut.info.interest <- tip.df.with_mut %>% left_join(mut.info %>% filter(mut.mod == mut.name),by="Id")

      mut_type.major <- mut.info.interest %>% group_by(mut) %>% summarize(count = n()) %>% slice_max(count,n=1,with_ties = F) %>% pull(mut)

      temp.df <- data.frame(date = date.interest, mut_type = mut_type.major,num.tip.with_mut)
      transition.node.info.df <- rbind(transition.node.info.df,temp.df)

    }

    tree.info.df.interest.merged.0to1 <- cbind(tree.info.df.interest.merged.0to1,transition.node.info.df)
    tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged.0to1 %>% mutate(prop.num.tip.with_mut = num.tip.with_mut / num.decendants)
    
    # filtering transition nodes
    tree.info.df.interest.merged.0to1.filtered <- tree.info.df.interest.merged.0to1 %>% filter(num.tip.with_mut >= min.num.decendants, prop.num.tip.with_mut >= min.prop.decendants) %>% mutate(mut = mut.name)

    node.transition.df <- rbind(node.transition.df,tree.info.df.interest.merged.0to1.filtered)

  }

  #save transition branch info
  out.name <- paste("output/transition_branch.",lineage.interest,".txt",sep="")
  write.table(node.transition.df,out.name,col.names=T,row.names=F,sep="\t",quote=F)


  branch_state.df <- as.data.frame(branch_state.l)
  branch_state.df <- branch_state.df %>% mutate(total = apply(branch_state.df,1,sum), node = 1:nrow(branch_state.df))
  node.df <- node.transition.df %>% select(node,mut)


  #make plot
  plot.df <- tip.df.merged %>% arrange(tip_Id) %>% select(-tip_Id,-Id)
  rownames(plot.df) <- tip.df.merged$Id
  plot.df <- ifelse(plot.df == 1,"+","-")
  plot.df <- plot.df[,c("R346","K444","L452","N460","F486")]


  n.seq <- nrow(plot.df)

  plot.df.node.count <- branch_state.df %>% select(node,total)
  plot.df.node.count[is.na(plot.df.node.count)] <- 0

  col.branch.v <- c("gray75",brewer.pal(9, "Greens")[c(4,6)],"black")

  color.mut.v <- brewer.pal(6, "Set1")[c(1:2,4:6)]
  names(color.mut.v) <- c("R346","K444","L452","N460","F486")

  print("tree")

  g <- ggtree(tree,size=0.5,mrsd = max.date)
  g <- g %<+% plot.df.node.count + aes(color=total) 
  g <- g + scale_color_gradientn(colors=col.branch.v,breaks=seq(0,3),limits=c(0,3))
  g <- g + new_scale_color()
  g <- g %<+% node.df + aes(color=mut) 
  g <- g + geom_point2(size=2)
  g <- g + scale_color_manual(values=color.mut.v,breaks=names(color.mut.v),na.value = NA) 
  g <- g + theme_tree2()
  g <- g + scale_x_continuous(limits=c(2021,2024.25),breaks=seq(2021,2023,0.5))
  g <- g + ylim(0,n.seq + 700)


  g1 <- gheatmap(g, plot.df, offset = 0.2, width = 0.3,
                 colnames_position = "top",
                 colnames_offset_y = 0,colnames_angle=90,hjust=0,color=NA)
  g1 <- g1 + scale_fill_manual(values=c('gray90','navy'),breaks=c("-","+"))
  g1 <- g1 + ggtitle(lineage.interest)

  plot.l[[lineage.interest]] <- g1

}

#output plots
pdf.name <- "output/Omicron_trees.parsimony.pdf"
pdf(pdf.name,width=8,height=15)

wrap_plots(plot.l,ncol=1)

dev.off()




