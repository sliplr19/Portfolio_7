###This is my thesis code. Veronica and I both worked on it, so 
###I've labeled the sections we each worked on. Everything below
###the name is done by that person until it gets to another person's 
###name.

###Let me walk you through what we did. First, we created true values
###(using nonparanormal transformation) and parameters (not using
###nonparanormal transformation). Then my centrality measures
###were messed up, so I did some fiddling and got them to work. 
###Then I created a master long form dataframe of every centrality
###value of the parameters and all of the centrality values for
###the true values. I also found the mean, raw bias, and 
###standardized bias for the parameters.Then Veronica found the 
###nodes that are the minimum distance from the mean, minimum,
###and maximum nodes. Then I created pretty graphs of 
###the standardized biases.Then we did multilevel metamodels
###and used stargazer to come up with pretty tables...
###and that's the project!

#Conditions-determined by finding the quartiles of numerous
#empirical studies
#Lindley
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
    #Veronica
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
    ###Lindley
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

####Here's the code to get the centrality measures
###Lindley
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


###This code creates a master dataframe of all of the centrality measures.

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
cent_data_2 <- read.csv("D:\\Lindley Backup\\Thesis\\data\\c-cent.csv", header = TRUE,
                        stringsAsFactors = FALSE )

###Find the mean across replications for conditions
mean_con <- cent_data_2 %>% 
  group_by(nds, per, n, V2, file_name) %>%
  summarize(mean_c = mean(x, na.rm = TRUE))

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

##Combine the mean_con and all_parameter dataframes
mean_con <- full_join(mean_con, all_parameter, by = c("file_name", "nds", "per", "V2"))
     

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

#Veronica
#Function to calculate the percentage of zeroes in a given list
pct.zero <- function(x) {
  numerator <- length(subset(x, x == 0))
  denominator <- length(x)
  numerator/denominator
}

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

##Getting the values across reps/nodes for stnd_bias
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

###Lindley
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
mean_sel <- selection %>%
  filter(con_per == 0.4 & use.this.mean == 1) %>%
  ggplot(aes(x= factor(con_nodes, levels=c('7', '12', '20', '36', '57')), y = stnd_bias, fill = con_nodes)) +
  geom_col() +
  facet_wrap(~con_n + file_name, labeller = labeller(con_n = n.labs, file_name = file.labs)) +
  labs(x = "Number of nodes", y = "Standardized bias") +
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1")) +
  ylim(-35, 30)

###min graphs
min_sel <- selection %>%
  filter(con_per == 0.4 & use.this.min == 1) %>%
  ggplot(aes(x= factor(con_nodes, levels=c('7', '12', '20', '36', '57')), y = stnd_bias, fill = con_nodes)) +
  geom_col() +
  facet_wrap(~con_n + file_name, labeller = labeller(con_n = n.labs, file_name = file.labs)) +
  labs(x = "Number of nodes", y = "Standardized bias") +
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1"))

###max graphs
max_sel <- selection %>%
  filter(con_per == 0.4 & use.this.max == 1) %>%
  ggplot(aes(x= factor(con_nodes, levels=c('7', '12', '20', '36', '57')), y = stnd_bias, fill = con_nodes)) +
  geom_col() +
  facet_wrap(~con_n + file_name, labeller = labeller(con_n = n.labs, file_name = file.labs)) +
  labs(x = "Number of nodes", y = "Standardized bias") +
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1"))


#Get each person's deviation
#Veronica
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
#Model for betweenness at values close to the mean
lmer.mean.between <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = mean.between)
summary(lmer.mean.between)
#Just to give a sense of how we'd interpret this, there are negative main effects of nds and n
#so as the number of nodes and people increase, deviation from the rue value (i.e., inaccuracy) goes down
#but we have a positive interaction -- so the effect of each variable attenuates the effect of the other
#Nice explanation of this kind of finding: https://stats.stackexchange.com/questions/80050/two-negative-main-effects-yet-positive-interaction-effect
###Lindley
lmer.min.between <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = min.between)
summary(lmer.min.between)

lmer.max.between <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = max.between)
summary(lmer.max.between)

lmer.min.degree <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = min.degree)
summary(lmer.min.degree)

lmer.max.degree <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = max.degree)
summary(lmer.max.degree)

lmer.mean.degree <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = mean.degree)
summary(lmer.mean.degree)

lmer.min.close <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = min.close)
summary(lmer.min.close)

lmer.mean.close <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = mean.close)
summary(lmer.mean.close)

lmer.max.close <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = max.close)
summary(lmer.max.close)

lmer.min.eigen <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = min.eigen)
summary(lmer.min.eigen)

lmer.max.eigen <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = max.eigen)
summary(lmer.max.eigen)

lmer.mean.eigen <- lmer(deviation ~ nds.x*per.x*n.x + (1|unique.cell), data = mean.eigen)
summary(lmer.mean.eigen)

###Make tables
library(stargazer)

###Betweenness
##Veronica
class(lmer.mean.between) <- "lmerMod"
class(lmer.min.between) <- "lmerMod"
class(lmer.max.between) <- "lmerMod"
stargazer(lmer.mean.between, lmer.min.between, lmer.max.between, title="Models for Betweenness Centrality",
          show.dep.var.labels = FALSE,
          column.labels = c("Deviation at Mean", "Deviation at Minimum", "Deviation at Maximum"),
          covariate.labels = c("Number Nodes (NN)", "Percentage Connections (PC)", "Sample Size (SS)", "NN x PC", "NN x SS", "PC x SS", "NN x PC x SS"),
          align=TRUE)

##Lindley
###Degree

class(lmer.mean.degree) <- "lmerMod"
class(lmer.min.degree) <- "lmerMod"
class(lmer.max.degree) <- "lmerMod"
stargazer(lmer.mean.degree, lmer.min.degree, lmer.max.degree, title="Models for Degree Centrality",
          show.dep.var.labels = FALSE,
          column.labels = c("Deviation at Mean", "Deviation at Minimum", "Deviation at Maximum"),
          covariate.labels = c("Number Nodes (NN)", "Percentage Connections (PC)", "Sample Size (SS)", "NN x PC", "NN x SS", "PC x SS", "NN x PC x SS"),
          align=TRUE)

###Closeness

class(lmer.mean.close) <- "lmerMod"
class(lmer.min.close) <- "lmerMod"
class(lmer.max.close) <- "lmerMod"
stargazer(lmer.mean.close, lmer.min.close, lmer.max.close, title="Models for Closeness Centrality",
          show.dep.var.labels = FALSE,
          column.labels = c("Deviation at Mean", "Deviation at Minimum", "Deviation at Maximum"),
          covariate.labels = c("Number Nodes (NN)", "Percentage Connections (PC)", "Sample Size (SS)", "NN x PC", "NN x SS", "PC x SS", "NN x PC x SS"),
          align=TRUE)

###Eigenvector

class(lmer.mean.eigen) <- "lmerMod"
class(lmer.min.eigen) <- "lmerMod"
class(lmer.max.eigen) <- "lmerMod"
stargazer(lmer.mean.eigen, lmer.min.eigen, lmer.max.eigen, title="Models for Eigenvector Centrality",
          show.dep.var.labels = FALSE,
          column.labels = c("Deviation at Mean", "Deviation at Minimum", "Deviation at Maximum"),
          covariate.labels = c("Number Nodes (NN)", "Percentage Connections (PC)", "Sample Size (SS)", "NN x PC", "NN x SS", "PC x SS", "NN x PC x SS"),
          align=TRUE)