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

#Generate data
set.seed(128937)
#tic("parameter_set")
for(rep in 1:1){
  for(a in 1:length(nodes)){
    for(b in 1:length(percentages)){
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
    centralitymeasures <- centrality(Graph_lasso) #Use whichever settings you had before here
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

#####################################################################
#As a secondary check, let's make sure the correlations are what they should be
correct.covariance <- function(the.data, a, b) {
  e_exp.a <- exp(mean(the.data[,a]) + .5*var(the.data[,a]))
  e_exp.b <- exp(mean(the.data[,b]) + .5*var(the.data[,b]))
  third.term <- ((exp(cov(the.data[,a], the.data[,b])))-1)
  e_exp.a*e_exp.b*third.term
}

#######################################################################

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

#write.csv(run_times, "G:\\My Drive\\Thesis\\Data\\runtimes.csv")
write.csv(con_kappa, "G:\\My Drive\\Thesis\\Data\\kappa.csv")

###List of file names
file_name <- c("degree", "close", 
                "between", "eigen")

rep <- seq(1, 500, by=1)
nds <- seq(1, 5, by = 1)
per <- seq(1, 4, by = 1)
n <- seq(1, 4, by = 1)

###Write a function to compare results

diff_df <- function(file_name, rep, nds, per, n){                                                                                          
         parameter <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Parameter_", file_name, "\\", 1, "-", nds, "-", per, "-", file_name, ".csv"))
         condition <- read.csv(file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Condition_", file_name, "\\", rep, "-", nds, "-", per, "-", n, "-", file_name, ".csv"))
         diff <- as.data.frame(matrix(0, ncol = ncol(parameter), nrow = nrow(parameter)))
           for(k in 1:nrow(parameter)){
             for(m in 1:ncol(parameter)){
               diff_value <- parameter[k,m] - condition[k,m]
               diff[k,m] <- diff_value
           }
           }
         write.csv(diff, file <- paste0("D:\\Lindley Backup\\Thesis\\data\\Diff_", file_name, "\\", rep, "-", nds, "-", per, "-", n, "-", file_name, ".csv"))

}
grid <- expand.grid(file_name = file_name, rep = rep, nds = nds, 
                    per = per, n = n, 
                    stringsAsFactors = TRUE)
pmap(grid, diff_df)

