

# This code was written by Alexander Y. Mitrophanov, PhD. Finalized in July 2023.

# This script is the CI generator for pifcalc_supervised.R -- for the supervised-learning version of our algorithm. 

# This script is analogous to pifcalc_main_CIs.R, so see that script for more comments.


library(readxl)

# -------------------- CODE FOLDERS ----------------------------------

# It is assumed that the project scripts, including this one, will be run from this folder ("curr").
# It is also assumed that this folder has a subfolder DATA for all the (input and output) data files.
# !!!!! Users: change directory names and locations as needed.

curr <- getwd()
DATA_FOLDER <- paste(curr, "/DATA/", sep = "")


# -------------------- CODE PARAMS ----------------------------------

analysis_var <- 3

rel_tolerance <- 1e-8
num_reps <- 10000
set.seed(2)


# ---------------------------------- FUNCTIONS ---------------------------

  
loglikeli <- function(par, pathog_data, neut_data, experim_data, var_num) { # (-1)*log-likelihood, actually
  # This is for the supervised-learning version of the algorithm.
  
  # var_num : 3 is HAT, 4 is Cis, 5 is Ola
  
  # mixture_p <- par[1] # prob of pathogenic
  # TWO NORMALS
  loglike_value <- nrow(pathog_data)*log(par[1]) + sum(log(   dnorm(pathog_data[[var_num]],  mean = par[2], sd = par[3]    )   ))
  loglike_value <- loglike_value + nrow(neut_data)*log(1 - par[1]) + sum(log(   dnorm(neut_data[[var_num]], mean = par[4], sd = par[5]) ))
  
  return(-loglike_value)
  # take the negative because the optimizer minimizes by default
  
}



init_par <- function(dset1, var_number) {
  
  # initial-parameters vector
  par_vec <- c(0,0,0,0,0)
  par_vec[1] <- 0.5 # NEUTRAL (FLAT) PRIOR
  
  subset_neutral <- subset(dset1, classification == "Neutral")
  subset_pathogenic <- subset(dset1, classification == "Pathogenic")
  
  # initialization using the method of moments
  par_vec[2] <- mean( subset_pathogenic[[var_number]] )
  par_vec[3] <- sd( subset_pathogenic[[var_number]] )
  par_vec[4] <- mean( subset_neutral[[var_number]] )
  par_vec[5] <- sd( subset_neutral[[var_number]] )
  
  return(par_vec)
}




PIFs_semi_super_2sets_boot3 <- function(dset_predict,dset_labeled,exper_d,path_d,neut_d,rel_toler,meth) {
  
  PIF_list <- list() 
  
  for (I in 3:5) { # variables HAT, Cis, Ola
    
    init <- init_par(dset_labeled, I)
    result <- optim(par = init, fn = loglikeli, pathog_data = path_d, neut_data = neut_d, experim_data = exper_d, 
                    var_num = I, control = list(reltol = rel_toler), method = meth)
    params <- result$par

    
    PIF_list[[I - 2]] <- params[1]*dnorm(dset_predict[[I]], mean = params[2], sd = params[3])/
      (  params[1]*dnorm(dset_predict[[I]], mean = params[2], sd = params[3]) +
           (1 - params[1])*dnorm(dset_predict[[I]], mean = params[4], sd = params[5] )  )
  }
  
  PIF_list[[4]] <- PIF_list[[1]] + (1 - PIF_list[[1]])*PIF_list[[2]]*PIF_list[[3]]
  
  return(PIF_list)
  
}


PIF_boot <- function(data,rel_tol,num_reps,analysis_var) {

  meth <- "Nelder-Mead"
  output_PIFs <- list()
  
  data_pathog <- subset(data, classification == "Pathogenic")
  data_neut <- subset(data,  classification == "Neutral")
  data_unlabeled <- subset(data, classification == "Experimental")
  
  num_pathog <- nrow(data_pathog)
  num_neut <- nrow(data_neut)
  num_unlabeled <- nrow(data_unlabeled)
  
  for (I in 1:num_reps) {
    
    sample_pathog <- sample(num_pathog, replace = T)
    sample_neut <- sample(num_neut, replace = T)
    sample_unlabeled <- sample(num_unlabeled, replace = T)
    
    data_pathog_sample <- data_pathog[sample_pathog,]
    data_neut_sample <- data_neut[sample_neut,]
    data_unlabeled_sample <- data_unlabeled[sample_unlabeled,]
    
    data_labeled_sample <- rbind(data_pathog_sample,data_neut_sample)
    
    output_PIFs[[I]] <- PIFs_semi_super_2sets_boot3(data,data_labeled_sample,
                                                    data_unlabeled_sample,data_pathog_sample,
                                                    data_neut_sample,rel_tol,meth)[[analysis_var]]
    
    
  }
  
  return(output_PIFs)
  
  
}



quant_intervals <- function(boot_output) {
  
  lower_quant <- .025
  upper_quant <- .975
  
  num_variants <- length(boot_output[[1]])
  num_replicates <- length(boot_output)
  
  upper <- 1:num_variants
  lower <- 1:num_variants
  
  for (I in 1:num_variants){
    
    samples <- c()
    for (J in 1:num_replicates) {
      samples[J] <- boot_output[[J]][I]
    }
    
    upper[I] <- quantile(samples, probs = upper_quant)
    lower[I] <- quantile(samples, probs = lower_quant)
    
  }
  
  return(data.frame(lo = lower,up = upper))
  
}




# -----------------------------------------------------------------
# -------------------- DATA LOADING ------------------------------
# -----------------------------------------------------------------

file_in <- paste(DATA_FOLDER, "BRCA2_clinvar_variants_NGS_data.xlsx.xlsx", sep = "")
datafr <- read_excel(file_in) # main data frame
datafr <- datafr[c(1, 3:6)] # taking out the 'HGVS protein' column
colnames(datafr)[1] <- "variant"
colnames(datafr)[2] <- "classification"
colnames(datafr)[3] <- "HAT" 
colnames(datafr)[4] <- "Cis"
colnames(datafr)[5] <- "Ola" 

# accomodating the new terminology
datafr[datafr == "VUS"] <- "Experimental"
datafr[datafr == "Benign"] <- "Neutral"


# ---------------------------- BOOTSTRAP -------------------------------- 


boot_out <- PIF_boot(datafr,rel_tolerance,num_reps,analysis_var)
boot_out <- quant_intervals(boot_out)

file_out <- paste("pifcalc_supervised_CIs_", as.integer(analysis_var), sep = "")
file_out <- paste(file_out, ".RData", sep = "")
file_out <- paste(DATA_FOLDER, file_out, sep = "")
save(boot_out, file = file_out)









