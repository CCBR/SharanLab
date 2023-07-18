

# Code developed by Alexander Y. Mitrophanov, PhD. Finalized in July 2023.

# This is the script that calculates and plots PIFs for the supervised-learning version of our algorithm.
# This script is analogous to (a reduced version of) pifcalc_main.R, which implements the main (i.e., semi-supervised-learning) version of our algorithm.
# Supervised learning is implemented via a modified likelihood function, which does not use VUS (i.e., Experimental) data points.

# !!!!! Before running this script, run pifcalc_supervised_CIs.R, which generates and saves CI (confidence interval) information.
# Make sure to run pifcalc_supervised_CIs.R with the same value of the analysis_var parameter as assigned in this script (see below).

# !!!!! See pifcalc_main.R for more comments on the corresponding functions.


# ----------------------- LIBRARIES ---------------------------------

# MAKE SURE THE LIBRARIES LISTED BELOW ARE INSTALLED ON THE SYSTEM!

library(readxl)
library(dplyr)
library(writexl)


# -------------------- CODE FOLDERS ----------------------------------

# It is assumed that the project scripts, including this one, will be run from this folder ("curr").
# it is also assumed  that this folder has a subfolder DATA for all the (input and output) data files.
# !!!!!!! Users: change directory names and locations as needed.

curr <- getwd()
DATA_FOLDER <- paste(curr, "/DATA/", sep = "")


# -------------------- CODE PARAMS ----------------------------------


analysis_var <- 4
# 4 is default -- the main mode to run the algorithm.
# 4 means all variables (HAT, Cis, and Ola) are included in the final PIF calculation.
# Other options are: 1, 2, and 3 for only HAT, only Cis, and only Ola, respectively.
# CHANGE analysis_var (and rerun this script) to generate results for other possible values of analysis_var.

# NOTE: this script calculates (and saves in an Excel spreadsheet) all the four types of PIFs  
# (corresponding to the four possible values of analysis_var),
# and plots all PIFs in a 4-subplot figure (and saves CIs in Excel files), but uses only the analysis_var value indicated above 
# (and, therefore, the corresponding variables)
# in classification-accuracy calculation and cross-validation. 
# CI information is saved in numbered files and then integrating into the main-output PIF file MANUALLY.


plot_PIFs_legend_flag <- 1 # "1" means display the legend in the PIF plot, "0" otherwise.


pathog_thresh <- 0.99 # PIF greater than that is Pathogenic
neutral_thresh <- 0.05 # PIF less than or equal to that is Benign 
# the rest of the variants are Indeterminate
threshs <- c(pathog_thresh,neutral_thresh)

# target tolerance for the numerical optimizer
rel_tolerance <- 1e-8 # 1e-8 is typical value, see documentation for the function optim()
# numerical optimization methods for the log-likelihood:
opt_method <- "Nelder-Mead"
# Nelder-Mead in the default for the optim() function. 


# to "standardize" randomization results
set.seed(2)



# ------------------------------------------------------------------------------------
# ---------------------- FUNCTION DEFINITIONS ----------------------------------------
# ------------------------------------------------------------------------------------


pl <- function(x) {
  print(x, quote = FALSE)
}


loglikeli <- function(par, pathog_data, neut_data, experim_data, var_num) { # (-1)*log-likelihood, actually
  # This log-likelihood corresponds to supervised learning (only pathog and neut data are used)
  
  # var_num : 3 is HAT, 4 is Cis, 5 is Ola
  
  # mixture_p <- par[1] # prob of pathogenic
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



PIFs_semi_super_2sets <- function(dset_predict,dset,relt,meth) {
  # dset is for training, dset_predict is for prediction
  
  path_d <- subset(dset,classification == "Pathogenic")
  neut_d <- subset(dset,classification == "Neutral")
  exper_d <- subset(dset,classification == "Experimental")
  
  PIF_list <- list() 
  params_list_init <- c()
  params_list_fit <- c()
  
  for (I in 3:5) { # variables HAT, Cis, Ola
    
    init <- init_par(dset, I)
    params_list_init <- c(params_list_init,init)
    
    result <- optim(par = init, fn = loglikeli, pathog_data = path_d, neut_data = neut_d, experim_data = exper_d, 
                    var_num = I, control = list(reltol = relt), method = meth)
    params <- result$par
    params_list_fit <- c(params_list_fit,params)
    
    PIF_list[[I - 2]] <- params[1]*dnorm(dset_predict[[I]], mean = params[2], sd = params[3])/
      (  params[1]*dnorm(dset_predict[[I]], mean = params[2], sd = params[3]) +
           (1 - params[1])*dnorm(dset_predict[[I]], mean = params[4], sd = params[5] )  )
    
    
  }
  
  PIF_list[[4]] <- PIF_list[[1]] + (1 - PIF_list[[1]])*PIF_list[[2]]*PIF_list[[3]]
  PIF_list[[5]] <- params_list_init
  PIF_list[[6]] <- params_list_fit
  return(PIF_list)
  
}



accuracy_on_full_data <- function(data_set,rel_tol,meth,thr,analysis_var) {
  
  pathog_thresh <- thr[1]
  neutral_thresh <- thr[2]
  
  result <- PIFs_semi_super_2sets(data_set,data_set,rel_tol,meth)
  Combined_p <- result[[4]]
  HAT_p <- result[[1]]
  Cis_p <- result[[2]]
  Ola_p <- result[[3]]
  
  counter <- 0
  undetected_pathogenic <- c()
  undetected_neutral <- c()
  intermediate_list <- c()
  
  num_variants <- nrow(data_set)
  num_pathog <- nrow(subset(data_set, classification == "Pathogenic"))
  num_neut   <- nrow(subset(data_set, classification == "Neutral"))
  
  
  for (I in 1:num_variants) {
    if ( (data_set[[I,2]] == "Neutral" & result[[analysis_var]][I] <= neutral_thresh)
         | (data_set[[I,2]] == "Pathogenic" & result[[analysis_var]][I] > pathog_thresh)) {counter <- counter + 1}
    if (data_set[[I,2]] == "Neutral" & result[[analysis_var]][I] > neutral_thresh) {
      undetected_neutral <- c(undetected_neutral,I)
    }
    if (data_set[[I,2]] == "Pathogenic" & result[[analysis_var]][I] <= pathog_thresh) {
      undetected_pathogenic <- c(undetected_pathogenic,I)
    }
    if (result[[analysis_var]][I] > neutral_thresh & result[[analysis_var]][I] <= pathog_thresh) {
      intermediate_list <- c(intermediate_list,I)
    }
  }
  
  data_set_p <- cbind(data_set,HAT_p,Cis_p,Ola_p,Combined_p)
  undpath <- data_set_p[undetected_pathogenic,]
  undneut <- data_set_p[undetected_neutral,]
  intermediate_vars <- data_set_p[intermediate_list,]
  pl(" ")
  pl("Undetected Pathogenic:")
  print(undpath)
  pl(" ")
  pl("Undetected Benign:")
  print(undneut)
  pl(" ")
  pl("Number of indeterminate variants:")
  print(nrow(intermediate_vars))
  pl("Indeterminate variants:")
  print(intermediate_vars)
  
  
  out <- list()
  out[[1]] <- result[[analysis_var]] # Combined_p
  out[[2]] <- counter/(num_pathog + num_neut)
  
  # this is the main output that later becomes the output file
  out[[3]] <- data_set_p 
  
  
  out[[4]] <- result[[5]]
  out[[5]] <- result[[6]]
  return(out)
  
}



plot_PIFs <- function(data_set,PIFs, plot_title = "", thr, analysis_var, DATA_FOLDER, legend_flag) {
  # FIGURE: the plot generated using this function was saved as a PDF (using the GUI in the plot); 
  # subplot labels are supposed to be added manually later. 

  vus_col <- "gray50"
  
  file_in <- paste("pifcalc_supervised_CIs_", as.character(analysis_var), sep = "")
  file_in <- paste(file_in, ".RData", sep = "")
  file_in <- paste(DATA_FOLDER, file_in, sep = "")
  load(file_in) # loading boot_out
  
  boot_out_CIs <- boot_out
  data_set <- cbind(data_set,PIFs)
  data_set_CIs <- data_set
  
  
  # ---------- saving CI info in Excel ------
  # ----- renaming columns for output ---------
  if (analysis_var == 1) { colnames(boot_out_CIs) <- c( "HAT_PIF_CI_lower", "HAT_PIF_CI_upper" ) }
  else if (analysis_var == 2) { colnames(boot_out_CIs) <- c( "Cisplatin_PIF_CI_lower", "Cisplatin_PIF_CI_upper" ) }
  else if (analysis_var == 3) { colnames(boot_out_CIs) <- c( "Olaparib_PIF_CI_lower", "Olaparib_PIF_CI_upper" )  }
  else if (analysis_var == 4) { colnames(boot_out_CIs) <- c( "(HAT + Cisplatin + Olaparib)_PIF_CI_lower", 
                                                             "(HAT + Cisplatin + Olaparib)_PIF_CI_upper" ) }
  
  data_set_CIs <- cbind(data_set_CIs, boot_out_CIs) 
  file_out_CIs <- paste(DATA_FOLDER, "CI_data_", as.character(analysis_var), "_supervised.xlsx", sep = "")
  write_xlsx(data_set_CIs, path = file_out_CIs)
  # ----------------------------------------
  
  
  data_set <- cbind(data_set, boot_out) 
  data_set <- arrange(data_set, PIFs)
  
  #dev.new()
  par(cex.main = 1.4, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L", font.main = 1)
  plot(data_set$PIFs, main = "", xlab = "", ylab = "")
  title(ylab = "Probability of impact on function (PIF)", line=2.4) # the "line" param is horizontal positioning
  title(xlab = expression( paste( italic("BRCA2 "), "variants in order of increasing PIF" )  ), line=2.4) # vertical positioning here
  title(main = plot_title, line = 0.2)
  if (legend_flag == 1) {
    legend(7,.88, legend=c("Pathogenic", "Benign", "VUS"), 
           col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  }
  abline(h = thr[1], lty = "dashed")
  abline(h = thr[2], lty = "dashed")
  
  for (I in 1:nrow(data_set)) {
    points(I,data_set$PIFs[I], col = vus_col)
  }
  
  for (I in 1:nrow(data_set)) {
    
    if (data_set[I,2] == "Neutral") {
      colr <- "blue"
      points(I,data_set$PIFs[I], col = colr)
    }
    else if (data_set[I,2] == "Pathogenic") {
      colr <- "red"
      points(I,data_set$PIFs[I], col = colr)
    }
    else {colr <- vus_col}
    
    # if (data_set$PIFs[I] <= thr[1] & data_set$PIFs[I] > thr[2] & analysis_var == 4) {
    #   text(I - 27, data_set$PIFs[I] + .01, labels = data_set[I,1], cex = .7, col = colr)
    # }  # decided no labels
    
   lines(  c(I,I), c(data_set$lo[I], data_set$up[I]), col = colr)  
    
    
  }
  
}



fold_CV_dataset_gen <- function(data_set,num_folds) {
  # This function generates num_folds data sets for num_folds-fold CV (cross-validation). These data sets are then used in
  # the function that actually does the CV and spits out its accuracy (averaged over the folds).
  
  # Subsets of labeled data!
  # This involves randomization!
  train_set <- subset(data_set, classification == "Neutral" | classification == "Pathogenic")
  N_train <- nrow(train_set)
  rnd_order <- sample(1:N_train)
  
  dsets <- list()
  fold_size <- N_train/num_folds
  
  for (I in 1:num_folds) { # fold number
    
    fold_start <- (I - 1)*fold_size + 1
    fold_end <- fold_start + fold_size - 1
    
    dsets[[I]] <- train_set[rnd_order[fold_start:fold_end],]
    
  }
  
  return(dsets)
  
}



K_fold_CV <- function(data_set,num_folds,rel_tol,meth,thresh_vec,analysis_var) {
  
  exp_data <- subset(data_set, classification == "Experimental")
  
  fold_data <- fold_CV_dataset_gen(data_set,num_folds)
  fold_size <- nrow(fold_data[[1]])
  accuracies <- 1:num_folds
  
  for (I in 1:num_folds) {
    
    test_set <- fold_data[[I]]
    train_set <- exp_data
    ind <- 1:num_folds
    ind <- ind[-I]
    
    for (J in ind) {
      train_set <- rbind(train_set,fold_data[[J]])
    }
    
    PIFs <- PIFs_semi_super_2sets(test_set,train_set,rel_tol,meth)[[analysis_var]]
    
    counter <- 0
    for (J in 1:fold_size) {
      
      if ( PIFs[J] > thresh_vec[1] & test_set$classification[J] == "Pathogenic" ) {
        counter <- counter + 1
      }
      else if ( PIFs[J] <= thresh_vec[2] & test_set$classification[J] == "Neutral" ) {
        counter <- counter + 1
      }
      
    }
    
    accuracies[I] <- counter/fold_size
    
  }
  
  return(mean(accuracies))
  
}




# -------------------------------------------------------------------------
# -------------------------------- DATA LOADING ---------------------------
# -------------------------------------------------------------------------

# The main data set will be contained in the data frame "datafr"

file_in <- paste(DATA_FOLDER, "BRCA2_clinvar_variants_NGS_data.xlsx", sep = "")
datafr <- read_excel(file_in) # main data frame
datafr <- datafr[c(1, 3:6)] # taking out the 'HGVS protein' column
colnames(datafr)[1] <- "variant"
colnames(datafr)[2] <- "classification"
colnames(datafr)[3] <- "HAT" 
colnames(datafr)[4] <- "Cis"
colnames(datafr)[5] <- "Ola" 

# accommodating the new terminology
datafr[datafr == "VUS"] <- "Experimental"
datafr[datafr == "Benign"] <- "Neutral"


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------- DATA ANALYSIS ------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ----------------- Fitting the model and plotting PIFs ----------------

pl(" --------------------------------------------------------- ")
pl("Accuracy on the training data (fitting to full data set):")
pl(" --------------------------------------------------------- ")
result <- accuracy_on_full_data(datafr,rel_tolerance,opt_method,threshs,analysis_var)
init_params_final <- result[[4]]
fit_params_final <- result[[5]]
acc <- result[[2]]


# -------------------- plotting a figure with 4 PIF subplots ------
dev.new()
par(mfrow = c(2,2))
plot_title <- ""
PIF_vecs <- result[[3]]

for (I in 1:4) {
  if (I == 1) {
    plot_PIFs_legend_flag <- 1
  }
  else {
    plot_PIFs_legend_flag <- 0
  }
  
  # if (I == 4) {
  #   plot_title <- "HAT + Cisplatin + Olaparib"
  # } else if (I == 1) {
  #   plot_title <- "HAT"
  # } else if (I == 2) {
  #   plot_title <- "Cisplatin"
  # } else if (I == 3) {plot_title <- "Olaparib"}
  
plot_PIFs(datafr,PIF_vecs[[5 + I]], plot_title, threshs, I, DATA_FOLDER, plot_PIFs_legend_flag)
}
# ----------------------------------------------------------------


data_and_probs <- result[[3]] # not just PIFs, but data and all probabs
pl(" ")
pl("Fitting accuracy (%):")
print(acc*100)


# --------------------------- CROSS-VALIDATION ---------------------------


# At each cross-validation step, the accuracy is calculated as the fraction of correctly classified variants,
# averaged over folds.

pl(" ")
pl(" ")
pl(" ---------------------------------------------------------------------- ")
pl("Vector of K values and K-fold CV accuracies (%), including LOOCV (last):")
pl(" ---------------------------------------------------------------------- ")
K_vals = c(3,6,9,43)
print(K_vals)
for (K in K_vals) { # 43 is the total number of labeled variants
  
  acc <- K_fold_CV(datafr,K,rel_tolerance,opt_method,threshs, analysis_var)
  print(acc*100)
  
}
pl(" ")


# ----------------------------- saving FS and PIF table ------------------------

# changing column names a bit and accommodating new terminology
colnames(data_and_probs) <- c("HGVS nucleotide","Functional classification","HAT FS","Cisplatin FS", "Olaparib FS",
                              "HAT_PIF", "Cisplatin_PIF", "Olaparib_PIF", "(HAT + Cisplatin + Olaparib)_PIF")
data_and_probs[data_and_probs == "Experimental"] <- "VUS"
data_and_probs[data_and_probs == "Neutral"] <- "Benign"

out_file <- paste(DATA_FOLDER, "BRCA2_clinvar_variants_supervised_PIFs.xlsx", sep = "")
write_xlsx(data_and_probs, path = out_file)




