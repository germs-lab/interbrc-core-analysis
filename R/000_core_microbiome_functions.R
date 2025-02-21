# Load required packages
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(hillR)
library(microViz)
library(magrittr)
library(lazyeval)
library(minpack.lm)
library(Hmisc)
library(stats4)


# The functions below use abundance-occupancy distributions fitted to a neutral model, described by Shade and Stopnisek, 2019. Core microbial taxa are selected based on their contributions to overall microbial beta-diversity. In the described function, a core microbial taxa must contribute at least ~2% of variation to Bray-Curtis dissimilarity to be considered a 'core' microbial taxa -- but this value can be manipulated within the function. 

# Core microbial taxa whose abundance and occupancy are above the fitted neutral model's confidence intervals indicates greater occupancy across samples given abundance, suggesting deterministic selection by the plant. Alternatively, core microbial taxa that are below the fitted neutral model indicates greater abundance given lower occupancy; these taxa may be  dispersal limited.

# More info on the functions can be found in Shade and Stopnisek, 2019. Original functions were developed in VanWallendael et al 2021. Code was adapted and updated by Nicco Benucci (GLBRC, Bonito Lab) and Brandon Kristy (GLBRC, Evans Lab). 

# EXTRACT CORE FUNCTION: This function extracts a core microbial community based on abundnace occupancy distributions and each taxa's contributions to BC-dissimilarity. The threshold defined for this analysis is 2% (1.02 in the below function). 
ExtractCore <- function(physeq, Var, method, Group=NULL, Level=NULL){
  set.seed(37920)
  # input dataset needs to be rarified and minimum depth included
  if (min(sample_sums(physeq)) == max(sample_sums(physeq))) {
    nReads = min(sample_sums(physeq))
    #nReads %T>% print()
    rare <- physeq
    #rare %T>% print()
    taxon <- as(tax_table(rare), "matrix")
    taxon <- as.data.frame(taxon)
    dim(taxon)  %T>% print()
  } else {
    nReads = min(sample_sums(physeq))
    #nReads %T>% print()
    rare <-
      rarefy_even_depth(
        physeq, sample.size = nReads,
        trimOTUs = TRUE,
        replace = TRUE,
        verbose = FALSE)
    taxon <- as(tax_table(rare), "matrix")
    taxon <- as.data.frame(taxon)
  }
  # choosing a subset or using the whole phyloseq object as is
  if (is.null(Group)) {
    otu <- rare@otu_table %>% as("matrix")
    map <- rare@sam_data %>% as("data.frame")
  } else{
    sub_group <- substitute(Group)
    sub_set <- subset(sample_data(rare), eval(parse(text=sub_group)) %in% Level)
    physeq1 <- merge_phyloseq(otu_table(rare),
                              tax_table(rare),
                              refseq(rare),
                              sub_set)
    otu_table(physeq1) <- otu_table(physeq1)[which(rowSums(otu_table(physeq1)) > 0), ]
    otu <- physeq1@otu_table %>% as("matrix")
    map <- physeq1@sam_data %>% as("data.frame")
    print("Grouping Factor")
    map[,Group] %T>% print()
    taxon <- as(tax_table(physeq1), "matrix")
    taxon <- as.data.frame(taxon)
  }
  map$SampleID <- rownames(map)
  #print("Check: dimension of datframe and metadata")
  #dim(otu) %T>% print() # funcitons form magrittr package
  #dim(map) %T>% print() 
  # calculating occupancy and abundance
  otu_PA <-
    1 * ((otu > 0) == 1)                                             # presence-absence data
  otu_occ <-
    rowSums(otu_PA) / ncol(otu_PA)                                   # occupancy calculation
  otu_rel <-
    apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # mean relative abundance
  occ_abun <-
    add_rownames(as.data.frame(cbind(otu_occ, otu_rel)), "otu")     # combining occupancy and abundance data frame
  # NOTE! add_rownames is deprecated and generates a warning, a bug of tidyverse, 
  # alternative you can use:
  # occ_abun <- tibble::rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),"otu")
  # Ranking OTUs based on their occupancy
  # For caluclating raking index we included following conditions:
  # - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
  # - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
  Var <- enquo(Var) # lazy evaluation
  PresenceSum <-
    data.frame(otu = as.factor(row.names(otu)), otu) %>%
    gather(SampleID, abun,-otu) %>%
    left_join(map, by = "SampleID") %>%
    group_by(otu, !!Var) %>%
    dplyr::summarise(
      time_freq = sum(abun > 0) / length(abun), # frequency of detection between time points
      coreTime = ifelse(time_freq == 1, 1, 0)) %>% # 1 only if occupancy 1 with specific time, 0 if not
    group_by(otu) %>%
    dplyr::summarise(
      sumF = sum(time_freq),
      sumG = sum(coreTime),
      nS = length(!!Var)* 2,
      Index = (sumF + sumG) / nS) # calculating weighting Index based on number of points detected
  # PresenceSum %T>% print()
  # ranking otus
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by = "otu") %>%
    transmute(otu = otu,
              rank = Index) %>%
    arrange(desc(rank))
  #otu_ranked %T>% print()
  # calculating BC dissimilarity based on the 1st ranked OTU
  otu_start = otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start, ])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x)
    sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
  x_names <-
    apply(combn(ncol(start_matrix), 2), 2, function(x)
      paste(colnames(start_matrix)[x], collapse = "-"))
  # creating a data.frame and adding first OTU name as 1
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  # Initialize your data structures:  calculating the contribution of ranked OTUs to the BC similarity
  BCaddition <- NULL
  BCaddition <- rbind(BCaddition,df_s)
  # calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. 
  # Can be set to the entire length of OTUs in the dataset.
  # it might take some time if more than 5000 OTUs are included.
  for(i in 2:411){ #nrow(otu_ranked)
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    y <-
      apply(combn(ncol(start_matrix), 2), 2, function(y)
        sum(abs(start_matrix[, y[1]] - start_matrix[, y[2]])) / (2 * nReads))
    df_a <- data.frame(x_names, y)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  # Calculating the BC dissimilarity of the whole dataset (not needed if the second loop 
  # is already including all OTUs)
  z <-
    apply(combn(ncol(otu), 2), 2, function(z)
      sum(abs(otu[, z[1]] - otu[, z[2]])) / (2 * nReads))
  # overwrite the names here
  x_names <-
    apply(combn(ncol(otu), 2), 2, function(x)
      paste(colnames(otu)[x], collapse = "-"))
  df_full <- data.frame(x_names, z)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition, df_full, by='x_names')
  # ranking the obtained BC 
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    # Calculate mean Bray-Curtis dissimilarity
    summarise(MeanBC=mean(BC)) %>%            
    arrange(desc(-MeanBC)) %>%
    # Calculate proportion of the dissimilarity explained by the n number of ranked OTUs 
    mutate(proportionBC=MeanBC/max(MeanBC))
  #BC_ranked %T>% print()
  # Calculating the increase BC
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  BC_ranked <- drop_na(BC_ranked) 
  if (method=="elbow"){
    fo_difference <- function(pos){
      left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
      right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
      return(left - right)
    }
    BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
    elbow <- which.max(BC_ranked$fo_diffs)
    core_otus <- otu_ranked$otu[1:elbow]
    #core_otus %T>% print()
  }
  # Creating threshold for core inclusion - last call method using a
  # final increase in BC similarity of equal or greater than 5%
  else{
    lastCall <-
      as.numeric(as.character(dplyr::last(
        subset(BC_ranked, IncreaseBC >= 1.02)$rank)))
    core_otus <- otu_ranked$otu[1:lastCall]
    #core_otus %T>% print()
  }
  # Adding Core otus for creating occupancy abundance plot
  occ_abun$fill <- 'no'
  occ_abun$fill[occ_abun$otu %in% core_otus] <- "core"
  return_list <-
    list(core_otus, BC_ranked, otu_ranked, occ_abun, otu, map, taxon)
  return(return_list)
}

## NOTE: This technique requires en even sampling depth to perform, which requires rarefaction. I tested this function /wo rarefaction (using the mean sampling depth instead), and there were no differences in core community composition or identification.

# NEUTRAL MODEL FUNCTION: This function fits the extracted core OTU table from ExtractCore to a neutral model, identifying taxa that are above or below the fitted model predictions. This provides insights into potential taxa that may be deterministically selected by the plant host. 
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}