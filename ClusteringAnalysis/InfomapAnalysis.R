#------ Script: random validation for infomap -------
#author: Geutg
#date: 20.1.2020
#randomize the infomap running to validate the observed
#
# the observed L with 500 tryouts is: 7.4675047349999994495

library(vegan)
library(bipartite)
library(tidyverse)
library(magrittr)
library(ggplot2)
#?commsim # Take a look here for a list of algorithms for shuffling.


#------------- randomizing functions --------------
create_randomized_binary_networks <- function(network, a_seed){
  A <- 1*(network>0) # Make the matrix binary
  
  # Create a set of randomized networks
  null <- vegan::nullmodel(A, method = 'curveball')
  shuffled_matrices <- simulate(null, nsim = SIMULATION_NUMBER, burnin = 5000, seed = a_seed)
  
  # The output is a 3D array (rows and columns as original matrix and another dimension for the 1000 simulations)  
  return (shuffled_matrices)
}
create_randomized_weighted_networks <- function(network, a_seed){
  # Tested hypothesis: proteins and chaperons interact with a given number of partners
  # assign the values after the shuffeling of a binary mat with curveball
  shuffled_matrices <- create_randomized_binary_networks(network, a_seed)
  
  # Get the the PORITIVE values of the edges
  link_values <- unlist(network)
  link_values <- subset(link_values, link_values>0)
  
  # Distribute the weights back to each randomization. Can also be done without a
  # for loop using apply but I don't have time to figure it out. It is fast
  # anyways as these are small matrices.
  for (i in 1:SIMULATION_NUMBER){
    x <- shuffled_matrices[,,i]
    x[x==1] <- sample(link_values)
    shuffled_matrices[,,i] <- x
  }
  return(shuffled_matrices)
}
create_randomized_weighted_non_constrained_networks <- function(network, a_seed){
  # Tested hypothesis: interactions between proteins and chaperons are non-random and are not constrained.
  # It preserves grand sum and the cells of the matrix are shuffled
  A <- network*(network>0)
  
  null <- vegan::nullmodel(A, method = 'r00_samp')
  shuffled_matrices <- simulate(null, nsim = SIMULATION_NUMBER, seed = a_seed)
  
  return (shuffled_matrices)
}


#------------- other functions--------
get_id_from_ensg <- function(ens_id){
  return(Prot_Attrib[ens_id,]$ID)
}

infomap_link_file_from_matrix <- function(network_matrix, output_path){
  # initdata string for writing 
  file_lines <- "#source target weight"
    
  for(i in rownames(network_matrix)){
    for(j in colnames(network_matrix)){
      r <- network_matrix[i,j]
      if (r > 0) {
        chap_id <- get_id_from_ensg(i)
        prot_id <- get_id_from_ensg(j)
        r_line <-paste(chap_id, prot_id, r, sep = " ")
        file_lines <- append(file_lines, r_line) 
      }
    }
  }
  
  #save the data to a file
  writeLines(file_lines, output_path)
}

infomap_l_value_from_file <- function(input_file_path, output_path, output_name){
  #exampled
  #input_file_path <- "~/InfomapTool/temps.txt"
  #output_name <- "testree1"
  #output_path <- "~/InfomapTool/WorkingFolder"
  
  # build command line to run infomap - included the -2 flag to avoid submodules
  infomap_command <- paste("~/InfomapTool/Infomap", input_file_path, output_path, "-i link-list -u -2 -N 500 --out-name", output_name, "--tree")
  
  #Run infomap
  system(infomap_command)
  
  # Get the map equation value, L from the output tree file
  output_file_path <- paste(output_path,"/",output_name,'.tree', sep = '')
  output <- read_lines(output_file_path)[1]
  pos_start <- str_locate_all(output, 'codelength')[[1]][2,2]+2
  pos_end <- str_locate_all(output,' in ')[[1]][3,1]-1
  L <- parse_number(str_sub(output, pos_start, pos_end))

  return(L)
}

IsObservedSignificant <- function(clusters_file_path, randomizing_method, run_id = "none"){
  # read network *metrix* file
  obs_data <- read.table(clusters_file_path, row.names = 1, header = TRUE)
  # run the infomap with the observed
  #infomap_link_file_from_matrix(test_mat2, '/home/geut/Desktop/test_link_file.txt')
  infomap_link_file_from_matrix(obs_data, "~/InfomapTool/input/Observed_link_file.txt")
  L_obs <- infomap_l_value_from_file(input_file_path = "~/InfomapTool/input/Observed_link_file.txt", 
                                     output_name = "Observed_correlations", 
                                     output_path = "~/InfomapTool/output")

  L_vector <- double(SIMULATION_NUMBER)
  all_rand_mats <- randomizing_method(obs_data, a_seed = 42)
  for (i in 1:SIMULATION_NUMBER){
    link_file_path <- paste("~/InfomapTool/WorkingFolder/simulation_no_", i,".txt", sep = '')
    # create the input file for infomap then run it
    infomap_link_file_from_matrix(all_rand_mats[,,i], link_file_path)
    L_vector[i] <- infomap_l_value_from_file(input_file_path = link_file_path, 
                                             output_name = paste("sim_no_", i, sep = ''), 
                                             output_path = "~/InfomapTool/output/simulations")
    print(paste(i, "/", SIMULATION_NUMBER, " simulations done.", sep = ''))
  }
  print (L_obs)
  hist(L_vector, col="lightblue")
  abline(v = L_obs, col="red", lwd=1, lty=5)
  
  # save the L's to a file for farther usage
  output_lines <- paste(lapply(L_vector, toString), collapse  = "\n")
  L_file_name <- paste("~/Desktop/L_list_",SIMULATION_NUMBER, "sim_", run_id,".txt", sep = '')
  writeLines(output_lines, L_file_name)
  
  return (t.test(L_vector, mu=L_obs, conf.level=0.95))
}

#TODO
#--------------- post-run-processing -----------
file_to_histogram <- function(file_name, obs_L) {
  string_vec <- readLines(file_name)
  nums <- as.numeric(string_vec)
  
  hist(nums, col="lightblue", xlim = c(7.44, 7.96))
  abline(v = obs_L, col="red", lwd=1, lty=5)
  
  return (nums))
}

file_to_ggplot_histogram <- function(file_name, obs_L) {
  string_vec <- readLines(file_name)
  nums <- as.numeric(string_vec)
  barlines <- "darkblue"
  
  p <- ggplot(as.data.frame(nums), aes(x = nums)) +
    geom_histogram(colour = barlines, aes(fill = ..count..), bins = 30) +
    scale_x_continuous(name = "Simulated Infomap L Value",
                       breaks = seq(7.782, 7.794, 0.002),
                       limits = c(7.783, 7.7935)) +
    scale_y_continuous(name = "Count") +
    ggtitle("Frequency histogram of simulated L values", subtitle = "observed L value: 7.47") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  p
  ggsave("/home/geut/Desktop/L_values_histogram100.pdf")
  return (nums)
}

significance_of_obsrvd_value <- function(L_simulation, L_obs) {
  # see if a value is significantly less the a population.

  cases_count <- sum(L_simulation<L_obs) # Sum all the cases
  p_value <- cases_count/SIMULATION_NUMBER #Calculate the p-value
  print(p_value)
  return(p_value)
}

#------------- const and globals -------------------------
SIMULATION_NUMBER = 1000
OBS_L_VALUE = 7.467505
Prot_Attrib <- read.csv("/home/geut/R/RandomizeInfomap/Infomap_Data/protein_attrib.csv", row.names = 1, header = TRUE)
INPUT_DATA_FILE <- "/home/geut/R/RandomizeInfomap/Infomap_Data/Pan_CANCER_analysis.tab"


#------------- running the code that validates the observed ------
# run with curveball
ttest_data <- IsObservedSignificant(INPUT_DATA_FILE,
                                    create_randomized_weighted_networks,
                                    "curveball_1000prod")
the_obs_L_curveball <- ttest_data$null.value

numbs = file_to_ggplot_histogram('Desktop/L_list_1000sim_curveball_1000prod.txt', OBS_L_VALUE)


  significance_of_obsrvd_value(numbs, OBS_L_VALUE)

#------------- running on test data ------------------
IsObservedSignificant("/home/geut/R/RandomizeInfomap/Infomap_Data/pan_test_cut.tab",
                      create_randomized_weighted_non_constrained_networks_through_binary, "htshfsh")

#------------- code testing-------
test_mat1 <- read.table("/home/geut/R/RandomizeInfomap/Infomap_Data/pan_test.tab", row.names = 1, header = TRUE)
test_mat2 <- test_mat1[1:5, 1:5]
A <- test_mat2*(test_mat2>0)



L_obs <- 2
L_sim <- rnorm(1000, mean=2.3, sd=0.2)


qplot(L_sim)+geom_vline(xintercept = L_obs, color='red') #L_sim<L_obs # Condition
cases_count <- sum(L_sim<L_obs) # Sum all the cases
p_value <- cases_count/SIMULATION_NUMBER #Calculate the p-value




#------------- desposable ---------------
create_randomized_weighted_non_constrained_networks_through_binary <- function(network, a_seed){
  # Tested hypothesis: interactions between proteins and chaperons are non-random and are not constrained.
  # It preserves grand sum and the cells of the matrix are shuffled
  A <- 1*(network>0) # Make the matrix binary
  
  # Create a set of randomized networks
  null <- vegan::nullmodel(A, method = 'r00_samp')
  shuffled_matrices <- simulate(null, nsim = SIMULATION_NUMBER, seed = a_seed)
  
  # Get the the PORITIVE values of the edges
  link_values <- unlist(network)
  link_values <- subset(link_values, link_values>0)
  
  # Distribute the weights back to each randomization. Can also be done without a
  # for loop using apply but I don't have time to figure it out. It is fast
  # anyways as these are small matrices.
  for (i in 1:SIMULATION_NUMBER){
    x <- shuffled_matrices[,,i]
    x[x==1] <- sample(link_values)
    shuffled_matrices[,,i] <- x
  }
  return(shuffled_matrices)
}


# run again with r00_samp
ttest_data <- IsObservedSignificant(INPUT_DATA_FILE,
                                    create_randomized_weighted_non_constrained_networks,
                                    "r00_samp")
the_obs_L_r00_samp <- ttest_data$null.value

# run again with r00_samp through binary transtering
ttest_data <- IsObservedSignificant(INPUT_DATA_FILE,
                                    create_randomized_weighted_non_constrained_networks_through_binary,
                                    "binary_r00_samp")
the_obs_L_r00_samp_fixed <- ttest_data$null.value