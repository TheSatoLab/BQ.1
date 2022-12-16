#!/usr/bin/env R

library(tidyverse)
library(cmdstanr)


#parameters
generation_time <- 2.1
bin.size <- 1
S_w <- 10

iter_sampling <- 2000
iter_warmup <- 500
seed <- 1234
parallel_chains <- 4
adapt_delta <- 0.99
max_treedepth <- 15
chains <- 4


#inputs
stan.f.name <- "script/multinomial_mut_regression.stan"
multi_nomial_model <- cmdstan_model(stan.f.name)

mut.mat.name <- "input/S_substitution_profile.txt"
mut.mat <- read.table(mut.mat.name,header=T)

data.name <- "input/lineage_frequency.txt"
data <- read.table(data.name,header=T,sep="\t")


#output names
out.mut_effect.name <- "output/mutation_effect.txt"
out.S_haplotype_Re.name <- "output/S_haplotype_Re.txt"



#preprocessing
TS <- as.matrix(data.frame(X0 = 1, X1 = data$date.bin.num))

Y <- data %>% select(- date.bin.num)
group.df <- data.frame(group_Id = 1:ncol(Y), hap = colnames(Y))

mut.mat <- mut.mat[match(colnames(Y), mut.mat$hap),]
COEF <- cbind(mut.mat %>% select(-hap))
ref.hap.df <- data.frame(mut = colnames(COEF), hap = as.numeric(COEF[1,]))


COEF.2 <- t(t(COEF) - as.numeric(COEF[1,]))

Y <- Y %>% as.matrix()
Y_sum.v <- apply(Y,1,sum)

num.coef <- ncol(COEF.2)


data.stan <- list(K = ncol(Y),
                  N = nrow(Y),
                  D = num.coef,
                  TS = TS,
                  COEF = COEF.2,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  S_w = S_w,
                  Y_sum = Y_sum.v)

#stan fitting
fit.stan <- multi_nomial_model$sample(
      data = data.stan,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      seed = seed,
      parallel_chains = parallel_chains,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      chains = chains)

stat.info.growth_gain <- fit.stan$summary("growth_gain") %>% as.data.frame()
stat.info.growth_rate <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info.growth_gain <- stat.info.growth_gain %>% mutate(mut = colnames(COEF.2))
stat.info.growth_rate <- stat.info.growth_rate %>% mutate(hap = group.df$hap)

stat.info.growth_gain.q <- fit.stan$summary("growth_gain", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
stat.info.growth_rate.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)

stat.info.growth_gain <- stat.info.growth_gain %>% inner_join(stat.info.growth_gain.q,by="variable")
stat.info.growth_rate <- stat.info.growth_rate %>% inner_join(stat.info.growth_rate.q,by="variable")

write.table(stat.info.growth_rate,out.mut_effect.name,col.names=T,row.names=F,sep="\t",quote=F)
write.table(stat.info.growth_gain,out.mut_effect.name,col.names=T,row.names=F,sep="\t",quote=F)

