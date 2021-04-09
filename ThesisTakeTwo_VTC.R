###This is the code for my thesis. We start by generating the true
###values. Next, we generate the parameters based on the correlations
###of the true values. The centrality measures turned out weird,
###so I reran them. Next I created a master dataframe and then 
###found the mean across replications. Then I found the raw bias
###and standardized bias. Then there's a bunch of code that 
###needed to be fixed. Finally, I graphed the standardized bias
###and reran the metamodels. 

#Conditions
nodes <- c(7, 12, 20, 36, 57)
percentages <- c(0.2, 0.4, 0.6, 0.8)
N.vals <- c(274, 802, 2654, 34653)

#Packages
library(GeneNet)
library(corpcor)
#install.packages("rockchalk")
library(rockchalk)
library(huge)
library(qgraph)
library(igraph)
library(miscTools)
library(tidyverse)
library(data.table)
library(sjmisc)
library(dplyr)
library(lme4)
library(gtools)
library(plyr)
#install.packages("MuMIn")
library(MuMIn)

#Generate true data

set.seed(128937)

#tic("parameter_set")
for(rep in 1:1){
  for(a in 1:length(nodes)){
    for(b in 1:length(percentages)){
      ###Veronica
    pcor_net <- ggm.simulate.pcor(nodes[a], 
                                etaA = percentages[b],
                                stdprec=FALSE)
    write.csv(pcor_net, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Partial\\", rep, "-", a, "-", b, "-","-partial.csv"), row.names = FALSE)
    cor_net <- round(pcor2cor(pcor_net), 6)
    #Start here -- simulate data from these correlation matrices
    write.csv(cor_net, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Corr\\", rep, "-", a, "-", b, "-corr.csv"), row.names = FALSE)
    norm_data <- mvrnorm(100000, mu = rep(0, nodes[a]), Sigma = cor_net)
    write.csv(norm_data, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Norm\\", rep, "-", a, "-", b, "-norm.csv"), row.names = FALSE)
    #######
    #Generate non-normal data
    nonnorm_data <- exp(norm_data)
    write.csv(nonnorm_data, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_NonNorm\\", rep, "-", a, "-", b, "-nonnorm.csv"), row.names = FALSE)
    #Let's just get the means and covariance matrix of the non-normal data to make sure we have it empirically (since the below function holds only in expectation)
    cov_nonnorm <- cov(nonnorm_data)
    write.csv(cov_nonnorm, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CovNonNorm\\", rep, "-", a, "-", b, "-CovNonNorm.csv"), row.names = FALSE)
    mean_nonnorm <- colMeans(nonnorm_data)
    write.csv(mean_nonnorm, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_MeanNonNorm\\", rep, "-", a, "-", b, "-MeanNonNorm.csv"), row.names = FALSE)
    #######
    #OK, so here's where it gets confusing
    #We're going to use the nonparanormal transformation, which is a powerful way to relax the assumption of non-normality
    npn.transformed <- huge.npn(nonnorm_data)
    #Let's write out yet another correlation matrix -- this time for the nonparanormal transformed data
    cor.npn <- cor(npn.transformed)
    write.csv(cor.npn, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CorNPN\\", rep, "-", a, "-", b, "-cor_npn.csv"), row.names = FALSE)
    #Now we'll run a network analysis on the nonparanormal transformed data
    #The centrality values we get from this are our "population values"
    #I'm taking this just from Sacha Epskamp's page -- if you know a more expedient way to do it, go ahead of course!
    #Source: http://sachaepskamp.com/files/Cookbook.html
    Graph_lasso <- qgraph(cor.npn, graph = "glasso", layout = "spring",
                          sampleSize = nrow(nonnorm_data))
    ###Lindley
    centralitymeasures <- centrality(Graph_lasso) 
    graph <- as.igraph(Graph_lasso)
    e_central <- eigen_centrality(graph, directed = FALSE, scale = TRUE,
                                  weights = NULL, options = arpack_defaults)
    write.csv(centralitymeasures$InDegree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_degree\\", rep, "-", a, "-", b, "-degree.csv"), row.names = FALSE)
    write.csv(centralitymeasures$Closeness, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_close\\", rep, "-", a, "-", b, "-close.csv"), row.names = FALSE)
    write.csv(centralitymeasures$Betweenness, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_between\\", rep, "-", a, "-", b, "-between.csv"), row.names = FALSE)
    write.csv(e_central, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_eigen\\", rep, "-", a, "-", b, "-eigen.csv"), row.names = FALSE)
    
    
    }
  }
}
#toc()
####Veronica
#####################################################################
#As a secondary check, let's make sure the correlations are what they should be
correct.covariance <- function(the.data, a, b) {
  e_exp.a <- exp(mean(the.data[,a]) + .5*var(the.data[,a]))
  e_exp.b <- exp(mean(the.data[,b]) + .5*var(the.data[,b]))
  third.term <- ((exp(cov(the.data[,a], the.data[,b])))-1)
  e_exp.a*e_exp.b*third.term
}

#######################################################################
###Lindley
#run_times <- data.frame("Nodes", "Percentage", "Sample_Size", "Runtime")
con_kappa <- data.frame("Nodes", "Percentage", "Sample_Size", "Condition_number")

##Run conditions
set.seed(298438238)
for(rep in 1:500){
  for(a in 1:length(nodes)){
   for(b in 1:length(percentages)){
     for(c in 1:length(N.vals)){
      #start <- Sys.time()
      con_corr <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Corr\\", 1, "-", a, "-", b, "-corr.csv"))
      con_norm_data <- mvrnorm(N.vals[c], mu = rep(0, nodes[a]), Sigma = con_corr)
      write.csv(con_norm_data, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_Norm\\", rep, "-", a, "-", b, "-", c, "-cnorm.csv"), row.names = FALSE)
      con_nonnorm_data <- exp(con_norm_data)
      write.csv(con_nonnorm_data, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_NonNorm\\", rep, "-", a, "-", b, "-", c, "-cnonnorm.csv"), row.names = FALSE)
      con_cov_nonnorm <- cov(con_nonnorm_data)
      write.csv(con_cov_nonnorm, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_CovNonNorm\\", rep, "-", a, "-", b, "-", c, "-CovNonNorm.csv"), row.names = FALSE)
      con_mean_nonnorm <- colMeans(con_nonnorm_data)
      write.csv(con_mean_nonnorm, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_MeanNonNorm\\", rep, "-", a, "-", b, "-", c, "-MeanNonNorm.csv"), row.names = FALSE)
      con_nonnormal_corr <- cor(con_nonnorm_data)
      write.csv(con_nonnormal_corr, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_CorNonnormal\\", rep, "-", a, "-", b, "-", c, "-cor_non.csv"), row.names = FALSE)
      condition <- kappa(con_nonnormal_corr)
     # con_kappa <-rbind(con_kappa, data.frame(
      #  X.Nodes. = nodes[a],
       # X.Percentage. = percentages[b],
        #X.Sample_Size. = N.vals[c],
        #X.Condition_number. = condition))
      con_Graph_lasso <- qgraph(con_nonnormal_corr, graph = "glasso", layout = "spring",
                          sampleSize = nrow(nonnorm_data))
      con_centralitymeasures <- centrality(con_Graph_lasso) #Use whichever settings you had before here
      con_graph <- as.igraph(con_Graph_lasso)
      con_e_central <- eigen_centrality(con_graph, directed = FALSE, scale = TRUE,
                                  weights = NULL, options = arpack_defaults)
      write.csv(con_centralitymeasures$InDegree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_degree\\", rep, "-", a, "-", b, "-", c,  "-degree.csv"), row.names = FALSE)
      write.csv(con_centralitymeasures$Closeness, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_close\\", rep, "-", a, "-", b, "-", c, "-close.csv"), row.names = FALSE)
      write.csv(con_centralitymeasures$Betweenness, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_between\\", rep, "-", a, "-", b, "-", c, "-between.csv"), row.names = FALSE)
      write.csv(con_e_central, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_eigen\\", rep, "-", a, "-", b, "-", c, "eigen.csv"), row.names = FALSE)
      #run_times <-rbind(run_times, data.frame(
       # X.Nodes. = nodes[a],
        #X.Percentage. = percentages[b],
        #X.Sample_Size. = N.vals[c],
        #X.Runtime. = Sys.time() - start))
      
     }
    }
  }
}

###Veronica
sym  <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CorNPN\\", 1, "-", a, "-", b, "-cor_npn.csv"), header = TRUE, stringsAsFactors = FALSE)
is.sym <- matrix(NA, nrow(sym), ncol(sym)) #Call the putatively symmetric matrix 'sym'; then 'is.sym' will be the one where we record whether each pair is symmetric

for (a in 1:nrow(sym)) {
  for (b in 1:nrow(sym)) {
    if (round(sym[a,b], 3) == round(sym[b,a], 3)) {
      is.sym[a,b] <- 1}
    if (round(sym[a,b],3) != round(sym[b,a],3)) {
      is.sym[a,b] <- 0}
  }
}

#The degrees came out wonky, so I'm going to try something else.
#I'm worried about the other measures, so I'm just going to rerun
#everything

###Ignore this

for(rep in 1:500){
  for(a in 1:length(nodes)){
    for(b in 1:length(percentages)){
      for(c in 1:length(N.vals)){
        norm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Norm\\", 1, "-", a, "-", b, "-norm.csv"), header = TRUE, stringsAsFactors = FALSE)
        nonnorm_data <- exp(norm_data)
        cor.npn  <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CorNPN\\", 1, "-", a, "-", b, "-cor_npn.csv"), header = TRUE, stringsAsFactors = FALSE)
        roundCorr <- round(cor.npn, 0)
        rownames(cor.npn) <- colnames(cor.npn)
        Graph_lasso <- qgraph(cor.npn, graph = "glasso", layout = "spring",
                              sampleSize = nrow(nonnorm_data))
        graph <- as.igraph(Graph_lasso)
        p_centrality <- centrality(Graph_lasso, 
                   weighted = FALSE, signed = TRUE)
        eigen_cent <- eigen_centrality(graph, directed = FALSE, scale = TRUE)
        p_degree <- p_centrality[["OutDegree"]]
        p_between <- p_centrality[["Betweenness"]]
        p_close <- p_centrality[["Closeness"]]
        p_eigen <- eigen_cent$vector
        write.csv(p_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_degree\\", rep, "-", a, "-", b, "-degree.csv"), row.names = FALSE)
        write.csv(p_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_close\\", rep, "-", a, "-", b, "-close.csv"), row.names = FALSE)
        write.csv(p_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_between\\", rep, "-", a, "-", b, "-between.csv"), row.names = FALSE)
        write.csv(p_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_eigen\\", rep, "-", a, "-", b, "-eigen.csv"), row.names = FALSE)
        #condition
        con_nonnorm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_NonNorm\\", rep, "-", a, "-", b, "-", c, "-cnonnorm.csv"), header = TRUE, stringsAsFactors = FALSE)
        con_nonnormal_corr <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_CorNonnormal\\", rep, "-", a, "-", b, "-", c, "-cor_non.csv"), header = TRUE, stringsAsFactors = FALSE)
        rownames(con_nonnormal_corr) <- colnames(con_nonnormal_corr)
        con_Graph_lasso <- qgraph(con_nonnormal_corr, graph = "glasso", layout = "spring",
                                  sampleSize = nrow(nonnorm_data))
        c_graph <- as.igraph(con_Graph_lasso)
        c_centrality <- centrality(con_Graph_lasso, 
                                   weighted = FALSE, signed = TRUE)
        c_eigen_cent <- eigen_centrality(c_graph, directed = FALSE, scale = TRUE)
        c_degree <- c_centrality[["OutDegree"]]
        c_between <- c_centrality[["Betweenness"]]
        c_close <- c_centrality[["Closeness"]]
        c_eigen <- c_eigen_cent$vector
        write.csv(c_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_degree\\", rep, "-", a, "-", b, "-", c,  "-degree.csv"), row.names = FALSE)
        write.csv(c_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_close\\", rep, "-", a, "-", b, "-", c, "-close.csv"), row.names = FALSE)
        write.csv(c_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_between\\", rep, "-", a, "-", b, "-", c, "-between.csv"), row.names = FALSE)
        write.csv(c_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_eigen\\", rep, "-", a, "-", b, "-", c, "eigen.csv"), row.names = FALSE)
        
      }
    }
  }
}

###And this

for(rep in 1:500){
  for(a in 1:length(nodes)){
    for(b in 1:length(percentages)){
      for(c in 1:length(N.vals)){
        norm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Norm\\", 1, "-", a, "-", b, "-norm.csv"), header = TRUE, stringsAsFactors = FALSE)
        nonnorm_data <- exp(norm_data)
        cor.npn  <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CorNPN\\", 1, "-", a, "-", b, "-cor_npn.csv"), header = TRUE, stringsAsFactors = FALSE)
        roundCorr <- round(cor.npn, 0)
        rownames(cor.npn) <- colnames(cor.npn)
        Graph_lasso <- qgraph(cor.npn, graph = "glasso", layout = "spring",
                              sampleSize = nrow(nonnorm_data))
        graph <- as.igraph(Graph_lasso)
        p_centrality <- centrality(Graph_lasso, 
                                   weighted = FALSE, signed = TRUE)
        eigen_cent <- eigen_centrality(graph, directed = FALSE, scale = TRUE)
        p_degree <- p_centrality[["OutDegree"]]
        p_between <- p_centrality[["Betweenness"]]
        p_close <- p_centrality[["Closeness"]]
        p_eigen <- eigen_cent$vector
        write.csv(p_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_degree\\", rep, "-", a, "-", b, "-degree.csv"), row.names = FALSE)
        write.csv(p_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_close\\", rep, "-", a, "-", b, "-close.csv"), row.names = FALSE)
        write.csv(p_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_between\\", rep, "-", a, "-", b, "-between.csv"), row.names = FALSE)
        write.csv(p_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_eigen\\", rep, "-", a, "-", b, "-eigen.csv"), row.names = FALSE)
        #condition
        con_nonnorm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_NonNorm\\", rep, "-", a, "-", b, "-", c, "-cnonnorm.csv"), header = TRUE, stringsAsFactors = FALSE)
        con_nonnormal_corr <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_CorNonnormal\\", rep, "-", a, "-", b, "-", c, "-cor_non.csv"), header = TRUE, stringsAsFactors = FALSE)
        rownames(con_nonnormal_corr) <- colnames(con_nonnormal_corr)
        con_Graph_lasso <- qgraph(con_nonnormal_corr, graph = "glasso", layout = "spring",
                                  sampleSize = nrow(nonnorm_data))
        c_graph <- as.igraph(con_Graph_lasso)
        c_centrality <- centrality(con_Graph_lasso, 
                                   weighted = FALSE, signed = TRUE)
        c_eigen_cent <- eigen_centrality(c_graph, directed = FALSE, scale = TRUE)
        c_degree <- c_centrality[["OutDegree"]]
        c_between <- c_centrality[["Betweenness"]]
        c_close <- c_centrality[["Closeness"]]
        c_eigen <- c_eigen_cent$vector
        write.csv(c_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_degree\\", rep, "-", a, "-", b, "-", c,  "-degree.csv"), row.names = FALSE)
        write.csv(c_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_close\\", rep, "-", a, "-", b, "-", c, "-close.csv"), row.names = FALSE)
        write.csv(c_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_between\\", rep, "-", a, "-", b, "-", c, "-between.csv"), row.names = FALSE)
        write.csv(c_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_eigen\\", rep, "-", a, "-", b, "-", c, "eigen.csv"), row.names = FALSE)
        
      }
    }
  }
}

####Here's the code to get the centrality measures
#Lindley
for(rep in 1:500){
  for(a in 1:length(nodes)){
    for(b in 1:length(percentages)){
      for(c in 1:length(N.vals)){
        #norm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Norm\\", 1, "-", a, "-", b, "-norm.csv"), header = TRUE, stringsAsFactors = FALSE)
        #nonnorm_data <- exp(norm_data)
        #cor.npn  <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CorNPN\\", 1, "-", a, "-", b, "-cor_npn.csv"), header = TRUE, stringsAsFactors = FALSE)
       # roundCorr <- round(cor.npn, 0)
        #rownames(cor.npn) <- colnames(cor.npn)
      #  Graph_lasso <- qgraph(cor.npn, graph = "glasso", layout = "spring",
                        #      sampleSize = nrow(nonnorm_data))
       # graph <- as.igraph(Graph_lasso)
      #  E(graph)$weight <- NULL
       # p_centrality <- centrality(Graph_lasso, 
      #                             weighted = FALSE, signed = TRUE)
       # eigen_cent <- eigen_centrality(graph, directed = FALSE, scale = TRUE)
      #  p_degree <- p_centrality[["OutDegree"]]
       # p_between <- p_centrality[["Betweenness"]]
      #  p_close <- p_centrality[["Closeness"]]
       # p_eigen <- eigen_cent$vector
        #write.csv(p_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_degree\\", rep, "-", a, "-", b, "-degree.csv"), row.names = FALSE)
        #write.csv(p_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_close\\", rep, "-", a, "-", b, "-close.csv"), row.names = FALSE)
        #write.csv(p_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_between\\", rep, "-", a, "-", b, "-between.csv"), row.names = FALSE)
        #write.csv(p_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_eigen\\", rep, "-", a, "-", b, "-eigen.csv"), row.names = FALSE)
        #condition
        con_nonnorm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_NonNorm\\", rep, "-", a, "-", b, "-", c, "-cnonnorm.csv"), header = TRUE, stringsAsFactors = FALSE)
        con_nonnormal_corr <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_CorNonnormal\\", rep, "-", a, "-", b, "-", c, "-cor_non.csv"), header = TRUE, stringsAsFactors = FALSE)
        rownames(con_nonnormal_corr) <- colnames(con_nonnormal_corr)
        con_Graph_lasso <- qgraph(con_nonnormal_corr, graph = "glasso", layout = "spring",
                                  sampleSize = nrow(con_nonnorm_data))
        c_graph <- as.igraph(con_Graph_lasso)
        c_centrality <- centrality(con_Graph_lasso, 
                                   weighted = FALSE, signed = TRUE)
        c_eigen_cent <- eigen_centrality(c_graph, directed = FALSE, scale = TRUE)
        c_degree <- c_centrality[["OutDegree"]]
        c_between <- c_centrality[["Betweenness"]]
        c_close <- c_centrality[["Closeness"]]
        c_eigen <- c_eigen_cent$vector
        write.csv(c_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_degree\\", rep, "-", a, "-", b, "-", c,  "-degree.csv"), row.names = FALSE)
        write.csv(c_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_close\\", rep, "-", a, "-", b, "-", c, "-close.csv"), row.names = FALSE)
        write.csv(c_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_between\\", rep, "-", a, "-", b, "-", c, "-between.csv"), row.names = FALSE)
        write.csv(c_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_eigen\\", rep, "-", a, "-", b, "-", c, "-eigen.csv"), row.names = FALSE)
      }
    }
  }
}

for(a in 1:length(nodes)){
  for(b in 1:length(percentages)){
      norm_data <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_Norm\\", 1, "-", a, "-", b, "-norm.csv"), header = TRUE, stringsAsFactors = FALSE)
      nonnorm_data <- exp(norm_data)
      cor.npn  <- read.csv(file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_CorNPN\\", 1, "-", a, "-", b, "-cor_npn.csv"), header = TRUE, stringsAsFactors = FALSE)
      roundCorr <- round(cor.npn, 0)
      rownames(cor.npn) <- colnames(cor.npn)
      Graph_lasso <- qgraph(cor.npn, graph = "glasso", layout = "spring",
            sampleSize = nrow(nonnorm_data))
      graph <- as.igraph(Graph_lasso)
      p_centrality <- centrality(Graph_lasso, 
                                 weighted = FALSE, signed = TRUE)
      eigen_cent <- eigen_centrality(graph, directed = FALSE, scale = TRUE)
      p_degree <- p_centrality[["OutDegree"]]
      p_between <- p_centrality[["Betweenness"]]
      p_close <- p_centrality[["Closeness"]]
      p_eigen <- eigen_cent$vector
      write.csv(p_degree, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_degree\\", rep, "-", a, "-", b, "-degree.csv"), row.names = FALSE)
      write.csv(p_close, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_close\\", rep, "-", a, "-", b, "-close.csv"), row.names = FALSE)
      write.csv(p_between, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_between\\", rep, "-", a, "-", b, "-between.csv"), row.names = FALSE)
      write.csv(p_eigen, file = paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_eigen\\", rep, "-", a, "-", b, "-eigen.csv"), row.names = FALSE)
  }
}
#write.csv(run_times, "G:\\My Drive\\Thesis\\Data\\runtimes.csv")
write.csv(con_kappa, "G:\\My Drive\\Thesis\\Data\\kappa.csv")

###List of file names
file_name <- c("degree", "close", 
                "between", "eigen")

rep <- seq(1, 500, by=1)
nds <- seq(1, 5, by = 1)
per <- seq(1, 4, by = 1)
n <- seq(1, 4, by = 1)

###Write a function to compare results-I'm going to redo this.

#diff_df <- function(file_name, rep, nds, per, n){                                                                                          
 #        parameter <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_", file_name, "\\", 1, "-", nds, "-", per, "-", file_name, ".csv"))
  #       condition <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_", file_name, "\\", rep, "-", nds, "-", per, "-", n, "-", file_name, ".csv"))
   #      diff <- as.data.frame(matrix(0, ncol = ncol(parameter), nrow = nrow(parameter)))
    #       for(k in 1:nrow(parameter)){
     #        for(m in 1:ncol(parameter)){
      #         diff_value <- parameter[k,m] - condition[k,m]
       #        diff[k,m] <- diff_value
       #    }
        #   }
         #write.csv(diff, file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Diff_", file_name, "\\", rep, "-", nds, "-", per, "-", n, "-", file_name, ".csv"))

#}
#grid <- expand.grid(file_name = file_name, rep = rep, nds = nds, 
 #                   per = per, n = n, 
  #                  stringsAsFactors = FALSE)
#pmap(grid, diff_df)

###Analysis code starts here


#cent_data <- as.data.frame(matrix(0, ncol = 57))
#for(e in 1:length(file_name)){
 # for(b in 1:length(nds)){
  #  for(c in 1:length(per)){
   #   for(d in 1:length(n)){
    #      for(f in length(rep)){
     #        condition <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_", file_name[e], "\\", f, "-", b, "-", c, "-", d, "-", file_name[e], ".csv"))
      #       if(file_name[e] == "eigen"){
       #        condition_data <- as.data.frame(condition$vector)
        #     } else {
         #      condition_data <- as.data.frame(condition)
          #   }
            # for(i in 1:nrow(condition_data)){
             #   cent <- as.numeric(condition_data[i,])
              #  central <- append(central, cent, after = length(central))
               #  } 
#             condition_df <- data.frame(file_name = file_name[e], rep = f, nds = b, per = c, n = d)
 #            central_df <- as.data.frame(matrix(data = central, nrow = 1))
  #           data_c <- cbind(condition_df, central_df)
   #          cent_data <- rbindFill(cent_data, data_c)
    #      }
    #    }
     # }
    #}
#  }
  

# write.csv(cent_data, file <- "D:\\Lindley Backup\\Thesis\\data\\c-cent.csv")

###Take two: this code creates a master dataframe of all of the centrality measures.

cent_data_2 <- as.data.frame(matrix(ncol = 7))
x <- c("file_name", "rep", "nds", "per", "n", "x", "V2")
colnames(cent_data_2) <- x
for(e in 1:length(file_name)){
  for(b in 1:length(nds)){
    for(c in 1:length(per)){
      for(d in 1:length(n)){
        for(f in 1:length(rep)){
          condition <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_", file_name[e], "\\", f, "-", b, "-", c, "-", d, "-", file_name[e], ".csv"))
          condition_data <- condition
          for(i in 1:nrow(condition_data)){
            cent <- as.numeric(condition_data[i,])
            condition_data[i,2] <- as.numeric(i)
           }
          condition_df <- data.frame(file_name = file_name[e], rep = as.numeric(f),
                                     nds = as.numeric(b), per = as.numeric(c), 
                                     n = as.numeric(d)) %>%
             slice(rep(1:n(), each = nrow(condition_data)))
          data_c <- cbind(condition_df, condition_data)
          data_c_names <- colnames(data_c)
          cent_names <- colnames(cent_data_2)
          common_names <- intersect(data_c_names, cent_names)
          cent_data_2 <- rbind(cent_data_2[common_names], data_c[common_names])
        }
      }
    }
  }
}

write.csv(cent_data_2, file <- "D:\\Lindley Backup\\Thesis\\data\\c-cent.csv")
cent_data_2 <- read.csv("G:\\My Drive\\Thesis\\Thesis\\c-cent.csv", header = TRUE,
                        stringsAsFactors = FALSE )

###Find the mean across replications for conditions
mean_con <- cent_data_2 %>% 
  group_by(nds, per, n, V2, file_name) %>%
  summarize(mean_c = mean(x, na.rm = TRUE))

###Add new columns to mean_con to prepare for adding values
#mean_con <- mean_con %>%
 # add_column(raw_bias = NA,
  #           relative_bias = NA)



all_parameter <- data.frame(file_name = character(), 
                            nds = numeric(), per = numeric(),
                            true_val = numeric(), V2 = numeric())
###Create a dataframe with corresponding true values.
for(e in 1:length(file_name)){
  for(b in 1:length(nds)){
    for(c in 1:length(per)){
      parameter <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_", file_name[e], "\\", 1, "-", b, "-", c, "-", file_name[e], ".csv"))
        
        for(i in 1:nrow(parameter)){
          parameter[i,2] <- as.numeric(i)
        }
      parameter <- parameter %>%
        rename_all(funs(c("true_val", "V2")))
      #  parameter <- parameter %>% 
       # dplyr::rename(
        #  true_val = x
        #)
       
       parameter_df <- data.frame(file_name = as.character(file_name[e]), 
                                   nds = as.numeric(b), per = as.numeric(c)) %>%
          slice(rep(1:n(), each = nrow(parameter)))
        parameter_df$file_name <- as.character(parameter_df$file_name)
        data_p <- cbind(parameter_df, parameter)
        all_parameter <- rbind(all_parameter, data_p)
      }
    }
  }

write.csv(all_parameter, "G:\\My Drive\\Thesis\\Thesis\\all_parameter.csv")

##Combine the mean_con and all_parameter dataframes
mean_con <- full_join(mean_con, all_parameter, by = c("file_name", "nds", "per", "V2"))
     
###Add in some missing values
mean_con$true_val[mean_con$file_name == "eigen" & mean_con$nds == 1 & mean_con$per == 1] <- 1
mean_con$true_val[mean_con$file_name == "eigen" & mean_con$nds == 2 & mean_con$per == 1] <- 5.500059e-01


###Find the raw bias
rb_frame <- data.frame("nds" = numeric(), "per" = numeric(),"raw_bias" = numeric(),
                       "n" = numeric(), "V2" = numeric(), "file_name" = character(),
                       "mean_c" = numeric(), "true_val" = numeric())
for(e in 1:length(file_name)){
  for(b in 1:length(nds)){
    for(c in 1:length(per)){
      for(d in 1:length(n)){
        row_select <- mean_con %>% 
          filter(file_name == file_name[e] &
                   nds == b &
                   per == c &
                   n == d)
        for(i in 1:nrow(row_select)){
        row_select <- ungroup(row_select)
        rb <- row_select$mean_c[i] - row_select$true_val[i]
        rb_item <- row_select[i,]
        rb_item$raw_bias <- rb
        rb_frame <- rbind(rb_frame, rb_item)
        
        }
      }
    }
  }
}

###Join dataframes
mean_con <- full_join(mean_con, rb_frame, by = c("file_name", "nds", "per", "V2", "n"))

###Fix up the columns
mean_con <- mean_con %>% 
  rename(mean_c = mean_c.x,
         true_val = true_val.x
)
mean_con <- subset(mean_con, select = -c(mean_c.y, true_val.y) )


###Save mean_con so I don't mess it up

write.csv(mean_con, "C:\\Users\\Owner\\Google Drive\\Thesis\\Thesis\\mean_con.csv")
mean_con <- read.csv("C:\\Users\\Owner\\Google Drive\\Thesis\\Thesis\\mean_con.csv", header = TRUE,
                     stringsAsFactors = FALSE)

###Find sd for standardized bias

mean_c <- cent_data_2 %>%
  group_by(nds, per, n, V2, file_name) %>%
  summarize(sd_c = sd(x)) %>%
  ungroup()


###Add to mean_con
mean_con$sd_c <- mean_c$sd_c

###Add standardized bias to mean_con
mean_con <- mean_con %>%
  mutate(stnd_bias = raw_bias/sd_c)

###Write file again
write.csv(mean_con, "C:\\Users\\Owner\\Google Drive\\Thesis\\Thesis\\mean_con.csv")

###Ignore this until....
###There's some duplicate rows in mean_con. Let's get rid of them.
un_mean_con <- distinct(mean_con)

###Handle the zero sds
###Veronica
#Function to calculate the percentage of zeroes in a given list
pct.zero <- function(x) {
  numerator <- length(subset(x, x == 0))
  denominator <- length(x)
  numerator/denominator
}

mean_con$unique.node <- mean_con$nds*100 + mean_con$per*10 + mean_con$n

degree <- subset(mean_con, mean_con$file_name == "degree")
by(degree$sd_c, degree$unique.node, pct.zero)

between <- subset(mean_con, mean_con$file_name == "between")
by(between$sd_c, between$unique.node, pct.zero)

close <- subset(mean_con, mean_con$file_name == "close")
by(close$sd_c, close$unique.node, pct.zero)

eigen <- subset(mean_con, mean_con$file_name == "eigen")
by(eigen$sd_c, eigen$unique.node, pct.zero)

###Lindley
###subset to those with less than 50% sd = 0
degree_0_sd <- NA

between_0_sd <- c(543, 544, 532, 533, 523, 442, 443, 433)

close_0_sd <- c(541, 542, 543, 544, 531, 532, 533, 521, 
                522, 523, 511, 512, 441, 442, 443, 431,
                432, 433, 421, 422, 423, 411, 341, 342,
                343, 331, 332, 321, 231)

eigen_0_sd <- c(244, 224)

###Subset to those with stnd_bias < abs(100)
degree_100 <- degree %>%
  subset(stnd_bias >= 100 | stnd_bias <= -100)
d_100 <- c(341, 342, 421, 422, 423, 431, 432, 433, 441, 442, 443,
           511, 512, 521, 522, 523, 531, 532, 533, 541, 542, 543,
           544)
between_100 <- between %>%
  subset(stnd_bias >= 100 | stnd_bias <= -100)
b_100 <- c(331, 341, 342, 343, 411, 421, 422, 423, 431, 432, 433,
           441, 442, 443, 511, 512, 521, 522, 523, 531, 532, 533,
           541, 542, 543, 544)
close_100 <- close %>%
  subset(stnd_bias >= 100 | stnd_bias <= -100)
c_100 <- c(231, 321, 331, 332, 341, 342, 343, 411, 421, 422, 423,
           431, 432, 433, 441, 442, 511, 512, 521, 522, 523, 531,
           532, 533, 541, 542, 543, 544)
eigen_100 <- eigen %>%
  subset(stnd_bias >= 100 | stnd_bias <= -100)
e_100 <- c(111, 112, 113, 114, 121, 122, 123, 124, 211, 212, 213,
           214, 544)

###Subset datasets
degree <- degree %>%
  subset(!(.$unique.node %in% d_100))

between <- between %>%
  subset(!(.$unique.node %in% between_0_sd)) %>%
  subset(!(.$unique.node %in% b_100))

close <- close %>%
  subset(!(.$unique.node %in% close_0_sd)) %>%
  subset(!(.$unique.node %in% c_100))

eigen <- eigen %>%
  subset(!(.$unique.node %in% eigen_0_sd)) %>%
  subset(!(.$unique.node %in% e_100))

###Sort datasets for stnd_bias
degree <- degree %>%
  arrange(desc(stnd_bias))

between <- between %>%
  arrange(desc(stnd_bias))

close <- close %>%
 arrange(desc(stnd_bias))

eigen <- eigen %>%
  arrange(desc(stnd_bias))


###Don't need
###Find relative bias
rlb_frame <- data.frame("relative_bias" = numeric())
for(e in 1:length(file_name)){
  for(b in 1:length(nds)){
    for(c in 1:length(per)){
      for(d in 1:length(n)){
        row_select <- mean_con %>% 
          filter(file_name == file_name[e] &
                   nds == b &
                   per == c &
                   n == d)
        for(i in 1:nrow(row_select)){
          row_select <- ungroup(row_select)
          rlb <- row_select$raw_bias[i]/row_select$true[i]
          rlb_item <- row_select[i,]
          rlb_item$raw_bias <- rlb
          rlb_frame <- full_join(rlb_frame, rlb_item, by = c("relative_bias"))
        } 
      }
    }
  }
}

###Join dataframes
mean_con <- full_join(mean_con, rlb_frame, by = c("file_name", "nds", "per", "V2", "n"))





###Add new columns to mean_con to prepare for adding values
cent_data_2 <- cent_data_2 %>%
  add_column(par_val = NA)

###Create a dataframe of differences

###Ignore This

for(e in 1:length(file_name)){
  for(b in 1:length(nds)){
    for(c in 1:length(per)){
      for(d in 1:length(n)){
        parameter <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_", file_name[e], "\\", 1, "-", b, "-", c, "-", d, "-", file_name[e], ".csv"))
        #parameter %>% 
         # rename(
           #  par_val = x
          #)
        for(i in 1:nrow(parameter)){
          parameter_data <- as.numeric(parameter[i,])
         # parameter_data[i,2] <- as.numeric(i)
          cent_data_2[cent_data_2$file_name == file_name[e] & cent_data_2$nds == b &
                        cent_data_2$per == c & cent_data_2$n == d, cent_data_2$V2 == i]$par_val <- parameter_data[i,1]
        }
      }
    }
  }
}

###Create dataframe with true values
cent_data_2 <- full_join(cent_data_2, all_parameter, by = c("file_name", "nds", "per", "V2"))

###We need to make a new dataframe with all data, but with the cells
###excluded.


###Now add a difference column and find the difference
diff_frame <- data.frame("nds" = numeric(), "per" = numeric(), "diff_val" = numeric(),
                       "n" = numeric(), "V2" = numeric(), "file_name" = character(),
                        "X" = numeric(), "true_val" = numeric(),
                       "rep" = numeric(), "x" = numeric(),
                       "con_nodes" = numeric(), "con_per" = numeric(),
                       "con_n" = numeric(), "unique.node" = numeric())
for(e in 1:length(file_name)){
  for(b in 1:length(nds)){
    for(c in 1:length(per)){
      for(d in 1:length(n)){
        c_row_select <- cent_data_2 %>% 
          filter(file_name == file_name[e] &
                   nds == b &
                   per == c &
                   n == d)
        p_row_select <- all_parameter %>%
          filter(file_name == file_name[e] &
                   nds == b &
                   per == c)
        for(i in 1:nrow(c_row_select)){
          diff <- c_row_select$x[i] - p_row_select$true_val[i]
          diff_item <- c_row_select[i,]
          diff_item$diff_val <- diff
          diff_item$true_val <- p_row_select$true_val[i]
          diff_frame <- rbind(as.data.frame(diff_frame), as.data.frame(diff_item))

          
        }
      }
    }
  }
}

write.csv(diff_frame, "G:\\My Drive\\Thesis\\Thesis\\diff_frame.csv")

cent_data_2 <- cent_data_2 %>%
  group_by(per, file_name, n, nds, V2, rep) %>%
  mutate(diff_val = cent_data_2$x - all_parameter$true_val)
cent_data_2$diff_val <- cent_data_2$x - all_parameter$true_val

###We need to get the numerical values for the condition variables, so we'll do that here.

cent_data_2 <- cent_data_2 %>%
  mutate(con_nodes = case_when(
    nds == '1' ~ 7,
    nds == '2' ~ 12,
    nds == '3' ~ 20,
    nds == '4' ~ 36,
    nds == '5' ~ 57
  ))

cent_data_2 <- cent_data_2 %>%
  mutate(con_per = case_when(
  per == 1 ~ 0.2,
  per == 2 ~ 0.4,
  per == 3 ~ 0.6,
  per == 4 ~ 0.8
  ))

cent_data_2 <- cent_data_2 %>%
  mutate(con_n = case_when(
    n == 1 ~ 274,
    n == 2 ~ 802,
    n == 3 ~ 2654,
    n == 4 ~ 34653
  ))

###Also there's an annoying NA first row. Let's get rid of that
cent_data_2 <- cent_data_2[-1,]

###Now we're going to add another column
cent_data_2$unique.node <- cent_data_2$nds*100 + cent_data_2$per*10 + cent_data_2$n


###Time for metamodels!

###Degree metamodel
degree_cen <- cent_data_2 %>%
  filter(file_name == "degree") %>%
  subset(!(.$unique.node %in% degree_0_sd)) %>%
  subset(!(.$unique.node %in% d_100))
degree.model <- lmer(diff_val ~ con_nodes*con_per*con_n + (1 | rep), data = degree)


r_degree <- r.squaredGLMM(degree.model)

###Betweenness metamodel
betweenness <- cent_data_2 %>%
  filter(file_name == "between") %>%
  subset(!(.$unique.node %in% between_0_sd)) %>%
  subset(!(.$unique.node %in% b_100))
Betweenness.model <- lmer(diff_val ~ con_nodes*con_per*con_n + (1 | rep), data = betweenness)

r_between <- r.squaredGLMM(Betweenness.model)

###Closeness metamodel
closeness <- cent_data_2 %>%
  filter(file_name == "close") %>%
  subset(!(.$unique.node %in% close_0_sd)) %>%
  subset(!(.$unique.node %in% c_100))
closeness.model <- lmer(diff_val ~ con_nodes*con_per*con_n + (1 | rep), data = closeness)

r_close <- r.squaredGLMM(closeness.model)

###Eigenvector metamodel
eigenvector <- cent_data_2 %>%
  filter(file_name == "eigen") %>%
  subset(!(.$unique.node %in% eigen_0_sd)) %>%
  subset(!(.$unique.node %in% e_100))
eigen.model <- lmer(diff_val ~ con_nodes*con_per*con_n + (1 | rep), data = eigenvector)

r_eigen <- r.squaredGLMM(eigen.model)

###Stop ignoring: Fixing the stnd_bias and models
###Veronica
read.cent <- read.table("G:\\My Drive\\Thesis\\Thesis\\c-cent.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
read.con <- read.table("G:\\My Drive\\Thesis\\Thesis\\mean_con.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
#Generate a unique identifier for each cell

mean_con <- read.con
mean_con$unique.cell <- mean_con$nds*100 + mean_con$per*10 + mean_con$n

#Data management - get duplicate rows out
mean_con.subset <- mean_con[,c(-1)]
mean_con.subset <- dplyr::distinct(mean_con.subset)

#Get rid of extremely big values and also negative infinity
mean_con.subset <- subset(mean_con.subset, (mean_con.subset$stnd_bias > -100) & (mean_con.subset$stnd_bias < 100))


#Get the mean, minimum, and maximum values for each combination of unique cell and index
true.vals <- mean_con.subset[,c("unique.cell", "file_name", "true_val")]
the.means <- aggregate(true.vals, by = list(true.vals$file_name, true.vals$unique.cell), mean)
the.min <- aggregate(true.vals, by = list(true.vals$file_name, true.vals$unique.cell), min)
the.max <- aggregate(true.vals, by = list(true.vals$file_name, true.vals$unique.cell), max)

#More data management... get means, minimum, maximum in a form where they can be merged back in
the.means <- the.means[,c("Group.1", "Group.2", "true_val")]
names(the.means) <- c("file_name", "unique.cell", "mean_true")
the.min <- the.min[,c("Group.1", "Group.2", "true_val")]
names(the.min) <- c("file_name", "unique.cell", "min_true")
the.max <- the.max[,c("Group.1", "Group.2", "true_val")]
names(the.max) <- c("file_name", "unique.cell", "max_true")

#Lots of merge steps. The end result will be a dataset we'll call "all.together" that has all the values we just created
mean_con.subset <- merge(mean_con.subset, the.means, by = c("unique.cell", "file_name"))
mean_con.subset <- merge(mean_con.subset, the.min, by = c("unique.cell", "file_name"))
mean_con.subset <- merge(mean_con.subset, the.max, by = c("unique.cell", "file_name"))

all.together <- mean_con.subset





#Calculate deviations of each node's true value to mean, min, and max
all.together$minus_mean <- all.together$true_val - all.together$mean_true
all.together$minus_min <- all.together$true_val - all.together$min_true
all.together$minus_max <- all.together$true_val - all.together$max_true

#Generate selection matrix
#This is going to create three columns that basically form a selection matrix. If a node has a 1 for a given centrality measure and a given cell, it is the chosen node for metamodels and tabulations 
all.together$use.this.mean <- NA
all.together$use.this.min <- NA
all.together$use.this.max <- NA

#Note: This is going to throw a warning because for some cells there aren't any valid values of a given index -- e.g., no value of closeness in the cells with larger sample sizes and numbesr of nodes
for (a in unique(all.together$unique.cell)) {
  for (b in unique(all.together$file_name)) {
    the.rows <- which((all.together$unique.cell == a) & (all.together$file_name == b))
    cell.subset <- subset(all.together, ((all.together$unique.cell == a) & (all.together$file_name == b)))
    
    #Get the value that's closest to the mean
    the.min.mean <- min(cell.subset$minus_mean)
    closest.mean.list <- which(cell.subset$minus_mean == the.min.mean)
    if (length(closest.mean.list) > 1) {
      closest.mean <- sample(closest.mean.list, 1)
    } else {
      closest.mean <- closest.mean.list
    }
    out.mean <- rep(0, nrow(cell.subset))
    out.mean[closest.mean] <- 1
    all.together$use.this.mean[the.rows] <- out.mean
    
    #Get the value that's closest to the minimum
    the.min.min <- min(cell.subset$minus_min)
    closest.min.list <- which(cell.subset$minus_min == the.min.min)
    if (length(closest.min.list) > 1) {
      closest.min <- sample(closest.min.list, 1)
    } else {
      closest.min <- closest.min.list
    }
    out.min <- rep(0, nrow(cell.subset))
    out.min[closest.min] <- 1
    all.together$use.this.min[the.rows] <- out.min
    
    #Get the value that's closest to the maximum
    
    #Get the value that's closest to the max
    the.min.max <- min(cell.subset$minus_max)
    closest.max.list <- which(cell.subset$minus_max == the.min.max)
    if (length(closest.max.list) > 1) {
      closest.max <- sample(closest.max.list, 1)
    } else {
      closest.max <- closest.max.list
    }
    out.max <- rep(0, nrow(cell.subset))
    out.max[closest.max] <- 1
    all.together$use.this.max[the.rows] <- out.max
    
  }
}



#OK, now just get the selection matrix and true values
selection <- all.together[,c("unique.cell", "file_name", "V2", "true_val", "use.this.mean", "use.this.min", "use.this.max", "stnd_bias", "nds", "per", "n")]


#Now merge the selection matrix back with the replication-level data file
read.cent$unique.cell <- read.cent$nds*100 + read.cent$per*10 + read.cent$n
rep.level <- merge(read.cent, selection, by = c("unique.cell", "file_name", "V2"))
rep.level <- rep.level[order(rep.level$rep) , ]

#Lindley
###Add columns with numerical variable names
selection <- selection  %>%
  mutate(con_nodes = case_when(
    nds == '1' ~ '7',
    nds == '2' ~ '12',
    nds == '3' ~ '20',
    nds == '4' ~ '36',
    nds == '5' ~ '57'
  ))

selection  <- selection %>%
  mutate(con_per = case_when(
    per == 1 ~ 0.2,
    per == 2 ~ 0.4,
    per == 3 ~ 0.6,
    per == 4 ~ 0.8
  ))

selection  <- selection  %>%
  mutate(con_n = case_when(
    n == 1 ~ 274,
    n == 2 ~ 802,
    n == 3 ~ 2654,
    n == 4 ~ 34653
  ))


###Make graph of standardized bias
###Label facets
n.labs <- c("Sample size 274", "Sample size 802",
            "Sample size 2654", "Sample size 34653")
names(n.labs) <- c(274, 802, 2654, 34653)
file.labs <- c("Degree", "Betweenness", "Closeness", "Eigenvector")
names(file.labs) <- c("degree", "between", "close", "eigen")

library(tidyverse)
library(wesanderson)

###Graph for mean
selection %>%
  filter(con_per == 0.4 & use.this.mean == 1) %>%
  ggplot(aes(x= factor(con_nodes, levels=c('7', '12', '20', '36', '57')), y = stnd_bias, fill = con_nodes)) +
  geom_col() +
  facet_wrap(~con_n + file_name, labeller = labeller(con_n = n.labs, file_name = file.labs)) +
  labs(x = "Number of nodes", y = "Standardized bias") +
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1")) +
  ylim(-35, 30)

###min graphs
selection %>%
  filter(con_per == 0.4 & use.this.min == 1) %>%
  ggplot(aes(x= factor(con_nodes, levels=c('7', '12', '20', '36', '57')), y = stnd_bias, fill = con_nodes)) +
  geom_col() +
  facet_wrap(~con_n + file_name, labeller = labeller(con_n = n.labs, file_name = file.labs)) +
  labs(x = "Number of nodes", y = "Standardized bias") +
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1"))

###max graphs
selection %>%
  filter(con_per == 0.4 & use.this.max == 1) %>%
  ggplot(aes(x= factor(con_nodes, levels=c('7', '12', '20', '36', '57')), y = stnd_bias, fill = con_nodes)) +
  geom_col() +
  facet_wrap(~con_n + file_name, labeller = labeller(con_n = n.labs, file_name = file.labs)) +
  labs(x = "Number of nodes", y = "Standardized bias") +
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1"))

###Veronica
#Get each person's deviation
rep.level$deviation <- rep.level$x - rep.level$true_val

#We're going to now get 12 datasets -- mean, min, and max (3) for each metric (of which there are 4)
mean.between <- subset(rep.level, (rep.level$file_name == "between" & rep.level$use.this.mean == 1))
min.between <- subset(rep.level, (rep.level$file_name == "between" & rep.level$use.this.min == 1))
max.between <- subset(rep.level, (rep.level$file_name == "between" & rep.level$use.this.max == 1))
mean.close <- subset(rep.level, (rep.level$file_name == "close" & rep.level$use.this.mean == 1))
min.close <- subset(rep.level, (rep.level$file_name == "close" & rep.level$use.this.min == 1))
max.close <- subset(rep.level, (rep.level$file_name == "close" & rep.level$use.this.max == 1))
mean.degree <- subset(rep.level, (rep.level$file_name == "degree" & rep.level$use.this.mean == 1))
min.degree <- subset(rep.level, (rep.level$file_name == "degree" & rep.level$use.this.min == 1))
max.degree <- subset(rep.level, (rep.level$file_name == "degree" & rep.level$use.this.max == 1))
mean.eigen <- subset(rep.level, (rep.level$file_name == "eigen" & rep.level$use.this.mean == 1))
min.eigen <- subset(rep.level, (rep.level$file_name == "eigen" & rep.level$use.this.min == 1))
max.eigen <- subset(rep.level, (rep.level$file_name == "eigen" & rep.level$use.this.max == 1))

#################
#Make sure to load in the lmerTest package so that you can get statistical significance tests
library(lmerTest)
#Lindley
lmer.mean.between <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = mean.between)
summary(lmer.mean.between)
lmer.min.between <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = min.between)
summary(lmer.min.between)

lmer.max.between <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = max.between)
summary(lmer.max.between)

lmer.min.degree <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = min.degree)
summary(lmer.min.degree)

lmer.max.degree <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = max.degree)
summary(lmer.max.degree)

lmer.mean.degree <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = mean.degree)
summary(lmer.mean.degree)

lmer.min.close <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = min.close)
summary(lmer.min.close)

lmer.mean.close <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = mean.close)
summary(lmer.mean.close)

lmer.mean.close <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = max.close)
summary(lmer.max.close)

lmer.min.eigen <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = min.eigen)
summary(lmer.min.eigen)

lmer.max.eigen <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = max.eigen)
summary(lmer.max.eigen)

lmer.mean.eigen <- lmer(deviation ~ nds*per*n + (1|unique.cell), data = mean.eigen)
summary(lmer.mean.eigen)