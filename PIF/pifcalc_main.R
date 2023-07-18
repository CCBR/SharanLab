

# Code developed by Alexander Y. Mitrophanov, PhD. Finalized in July 2023.

# This is the main script that calculates and plots PIFs.
# !!!!! Before running it, run pifcalc_main_CIs.R, which generates and saves CI (confidence interval) information
# Make sure to run pifcalc_main_CIs with the same value of the analysis_var parameter as assigned in this script 
# (see below).


# ----------------------- LIBRARIES ---------------------------------

# MAKE SURE THE LIBRARIES LISTED BELOW ARE INSTALLED ON THE SYSTEM!

library(readxl)
library(dplyr)
library(nortest) # for lilliefors test
library(writexl)


# -------------------- CODE FOLDERS ----------------------------------

# It is assumed that the project scripts, including this one, will be run from *this* folder ("curr").
# It is also assumed  that this folder has a subfolder DATA for all the (input and output) data files.
# !!! Users: change directory names and locations as needed.

curr <- getwd()
DATA_FOLDER <- paste(curr, "/DATA/", sep = "")


# -------------------- CODE PARAMETERS ----------------------------------


analysis_var <- 4 # 4 is default -- the main mode to run the algorithm
# 4 means all variables (HAT, Cis, and Ola)
# 4 means all variables (HAT, Cis, and Ola) are included in the final PIF calculation.
# Other options are: 1, 2, and 3 for only HAT, only Cis, and only Ola, respectively.
# CHANGE analysis_var (and rerun this script) to generate results for other possible values of analysis_var.

# NOTE: this script calculates (and saves in an Excel spreadsheet) all the four types of PIFs  
# (corresponding to the four possible values of analysis_var),
# but uses only the analysis_var value indicated above (and, therefore, the corresponding variables)
# in PIF visualization, CI calculation and saving, classification-accuracy calculation, and cross-validation. 
# CI information is saved in numbered Excel files and then integrating into the main-output PIF file MANUALLY.


plot_PIFs_legend_flag <- 0 # "1" means display legend in PIF plot, "0" otherwise
# Select as needed for generating a particular PIF plot (see analysis_var above to choose).


pathog_thresh <- 0.99 # PIF greater than that is Pathogenic
neutral_thresh <- 0.05 # PIF less than or equal to that is Benign (formerly called Neutral)
# the rest of the variants are Indeterminate.
threshs <- c(pathog_thresh,neutral_thresh)


# target tolerance for the numerical optimizer
rel_tolerance <- 1e-8 # 1e-8 is typical value, see documentation for the function optim()
# numerical optimization methods for the log-likelihood:
opt_method <- "Nelder-Mead"
# Nelder-Mead in the default for the optim() function. Other available methods might give a speed increase but might be less robust

# to "standardize" randomization results
set.seed(2)

# !!! NOTE: in the DATA ANALYSIS section (below in this script), comment individual code blocks, as needed -- 
# if only a subset of the script's outputs (plots) need to be generated



# -------------------- GRAPHICS PARAMS ------------------------------

par(cex.main = 1.4, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L", font.main = 1)
# These should be the same parameters as the ones inside functions that invoke dev.new() (see function definitions below). 
# Primary tuning of these parameters was done by running the plot_PIFs function.
# For axis-label distance tuning, use line = 2.4
# (deviations are possible in specific cases).

# Also, the colors for Path, Benign, and VUS are red, blue, and gray50, respectively.
# FIGURES: all individual plots were exported as PDF (using the user interface in the plot) 
# and saved to a folder under meaningful names; these plots are to be assembled into figures 
# manually.


# ------------------------------------------------------------------------------------
# ---------------------- FUNCTION DEFINITIONS ----------------------------------------
# ------------------------------------------------------------------------------------


pl <- function(x) {
  # convenient function for text output
  
  print(x, quote = FALSE)
}


loglikeli <- function(par, pathog_data, neut_data, experim_data, var_num) { 
  # (-1)*log-likelihood, actually -- what we will minimize
  # The likelihood definition used here is for semi-supervised learning.
  
  # var_num : 3 is HAT, 4 is Cis, 5 is Ola
  
  # mixture_p <- par[1] # prob of pathogenic
  # TWO NORMALS
  loglike_value <- nrow(pathog_data)*log(par[1]) + sum(log(   dnorm(pathog_data[[var_num]],  mean = par[2], sd = par[3]    )   ))
  loglike_value <- loglike_value + nrow(neut_data)*log(1 - par[1]) + sum(log(   dnorm(neut_data[[var_num]], mean = par[4], sd = par[5]) ))
  loglike_value <- loglike_value + sum(log( (1 - par[1])*dnorm(experim_data[[var_num]],mean = par[4], sd = par[5]) +
                                              par[1]*dnorm(experim_data[[var_num]], mean = par[2], sd = par[3] )  ) )
  
  return(-loglike_value)
  # take the negative because the optimizer minimizes by default
  
}



init_par <- function(dset1, var_number) {
  # mixture parameter initialization
  
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
  # this function does -log-likelihood minimization
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
  # this function fits the mixture model to the whole data set 
  
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
  out[[3]] <- data_set_p
  out[[4]] <- result[[5]]
  out[[5]] <- result[[6]]
  return(out)
  
}



plot_PIFs <- function(data_set,PIFs, plot_title = "", thr, analysis_var, DATA_FOLDER, legend_flag) {
  # main function for PIF-plot generation

  vus_col <- "gray50"
  
  file_in <- paste("pifcalc_main_CIs_", as.character(analysis_var), sep = "")
  file_in <- paste(file_in, ".RData", sep = "")
  file_in <- paste(DATA_FOLDER, file_in, sep = "")
  print(file_in)
  load(file_in) # load boot_out
  
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
  file_out_CIs <- paste(DATA_FOLDER, "CI_data_", as.character(analysis_var), "_main.xlsx", sep = "")
  write_xlsx(data_set_CIs, path = file_out_CIs)
  # ----------------------------------------
  
  
  data_set <- cbind(data_set, boot_out) 
  data_set <- arrange(data_set, PIFs)

  
  dev.new()
  par(cex.main = 1.4, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L", font.main = 1)
  plot(data_set$PIFs, main = "", xlab = "", ylab = "")
  title(ylab = "Probability of impact on function (PIF)", line=2.4) # the "line" param is horizontal positioning
  title(xlab = expression( paste( italic("BRCA2 "), "variants in order of increasing PIF" )  ), line=2.4) # vertical positioning here
  title(main = plot_title, line = 0.2)
  if (legend_flag == 1) {
    legend(7,.5, legend=c("Pathogenic", "Benign", "VUS"), 
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
    
    if (data_set$PIFs[I] <= thr[1] & data_set$PIFs[I] > thr[2] ) {
      
      # if (analysis_var == 4) { # decided no labels
      #   text(I - 15, data_set$PIFs[I], labels = data_set[I,1], cex = .7, col = colr)
      # }
      
    }  
    
   lines(  c(I,I), c(data_set$lo[I], data_set$up[I]), col = colr )  
    
    
  }
  
  
}





# --------------------------------------- PIF-PIF plotting functions ---------------------------


# PIF_PIF_plot_orig <- function(data_with_PIFs) { # this is the original version, which uses linear scale to show PIFs
#   # Plots PIF(all_vars) against PIF(HAT).
#   # Can be straightforwardly extended to other PIF combinations.
#   
#   vus_col <- "gray50"
#   PIF_difference_thresh = .05 
#   
#   dev.new()
#   par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
#   plot(data_with_PIFs$HAT_p, data_with_PIFs$Combined_p, col = vus_col, main = "", xlab = "", ylab = "")
#   lines(c(0,1), c(0,1), lty = 2)
#   title(xlab = "PIF(HAT assay only)", line=2.6) # vertical positioning here
#   title(ylab = "PIF(HAT + Cisplatin + Olaparib assay)", line=2.6) # the "line" param is horizontal positioning
#   legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
#          col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
#   
#   for (I in 1:nrow(data_with_PIFs)) {
#     
#     if (data_with_PIFs[I,2] == "Neutral") {
#       colr <- "blue"
#       points(data_with_PIFs$HAT_p[I], data_with_PIFs$Combined_p[I], col = colr)
#     }
#     else if (data_with_PIFs[I,2] == "Pathogenic") {
#       colr <- "red"
#       points(data_with_PIFs$HAT_p[I], data_with_PIFs$Combined_p[I], col = colr)
#     }
#     else {
#       colr <- vus_col 
#       }
#     
#     if(abs(data_with_PIFs$HAT_p[I] - data_with_PIFs$Combined_p[I]) > PIF_difference_thresh) {
#       text(data_with_PIFs$HAT_p[I] + .08, data_with_PIFs$Combined_p[I] + .008, labels = data_with_PIFs[I,1], cex = .7, col = colr)
#     }
#     
#   }
#   
#   # dev.new()
#   # par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
#   # plot(log(data_with_PIFs$HAT_p), log(data_with_PIFs$Combined_p), col = vus_col)
#   # lines(c(-10,0), c(-10,0), lty = 2)
#   
#   
# }




PIF_PIF_plot <- function(data_with_PIFs, threshs) { # this is the new version, which uses log scale
  # Plots PIF(all_vars) against PIF(HAT).
  # Can be straightforwardly extended to other PIF combinations.
  
  vus_col <- "gray50"
  PIF_difference_thresh = .05 
  
  dev.new()
  par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
  
  plot(log(data_with_PIFs$HAT_p), log(data_with_PIFs$Combined_p), col = vus_col, main = "", xlab = "", ylab = "",
       xlim = c(-10.5, .5), ylim = c(-10.5,.5))
  lines(c(-10.5, .5), c(-10.5, .5), lwd = .1)
  abline(h = log(threshs[1]), lty = 2)
  abline(h = log(threshs[2]), lty = 2)
  abline(v = log(threshs[1]), lty = 2)
  abline(v = log(threshs[2]), lty = 2)
  
  title(xlab = "log( PIF(HAT assay only) )", line=2.6) # vertical positioning here
  title(ylab = "log( PIF(HAT + Cisplatin + Olaparib assay) )", line=2.6) # the "line" param is horizontal positioning
  # legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
  #        col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(log(data_with_PIFs$HAT_p[I]), log(data_with_PIFs$Combined_p[I]), col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(log(data_with_PIFs$HAT_p[I]), log(data_with_PIFs$Combined_p[I]), col = colr)
    }
    else {
      colr <- vus_col 
    }
    
    if(abs(data_with_PIFs$HAT_p[I] - data_with_PIFs$Combined_p[I]) > PIF_difference_thresh) {
      text(log(data_with_PIFs$HAT_p[I]) - .18, log(data_with_PIFs$Combined_p[I]) - .28, labels = data_with_PIFs[I,1], cex = .7, col = colr)
    }
    
  }

}



PIF_PIF_plot_all_vs_Cis <- function(data_with_PIFs, threshs) { # this uses log scale
  # Plots PIF(all_vars) against PIF(Cis)
  
  vus_col <- "gray50"
  PIF_difference_thresh = .05 
  
  dev.new()
  par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
  
  plot(log(data_with_PIFs$Cis_p), log(data_with_PIFs$Combined_p), col = vus_col, main = "", xlab = "", ylab = "",
       xlim = c(-10.5, .5), ylim = c(-10.5,.5))
  lines(c(-10.5, .5), c(-10.5, .5), lwd = .1)
  abline(h = log(threshs[1]), lty = 2)
  abline(h = log(threshs[2]), lty = 2)
  abline(v = log(threshs[1]), lty = 2)
  abline(v = log(threshs[2]), lty = 2)
  
  title(xlab = "log( PIF(Cisplatin assay only) )", line=2.6) # vertical positioning here
  title(ylab = "log( PIF(HAT + Cisplatin + Olaparib assay) )", line=2.6) # the "line" param is horizontal positioning
  # legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
  #        col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(log(data_with_PIFs$Cis_p[I]), log(data_with_PIFs$Combined_p[I]), col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(log(data_with_PIFs$Cis_p[I]), log(data_with_PIFs$Combined_p[I]), col = colr)
    }
    else {
      colr <- vus_col 
    }
    
    if(abs(data_with_PIFs$Cis_p[I] - data_with_PIFs$Combined_p[I]) > PIF_difference_thresh) {
      text(log(data_with_PIFs$Cis_p[I]) - .18, log(data_with_PIFs$Combined_p[I]) - .25, labels = data_with_PIFs[I,1], cex = .7, col = colr)
    }
    
  }
  
}


PIF_PIF_plot_all_vs_Ola <- function(data_with_PIFs, threshs) { # this uses log scale
  # Plots PIF(all_vars) against PIF(Ola)
  
  vus_col <- "gray50"
  PIF_difference_thresh = .05 
  
  dev.new()
  par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
  
  plot(log(data_with_PIFs$Ola_p), log(data_with_PIFs$Combined_p), col = vus_col, main = "", xlab = "", ylab = "",
       xlim = c(-10.5, .5), ylim = c(-10.5,.5))
  lines(c(-10.5, .5), c(-10.5, .5), lwd = .1)
  abline(h = log(threshs[1]), lty = 2)
  abline(h = log(threshs[2]), lty = 2)
  abline(v = log(threshs[1]), lty = 2)
  abline(v = log(threshs[2]), lty = 2)
  
  title(xlab = "log( PIF(Olaparib assay only) )", line=2.6) # vertical positioning here
  title(ylab = "log( PIF(HAT + Cisplatin + Olaparib assay) )", line=2.6) # the "line" param is horizontal positioning
  # legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
  #        col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(log(data_with_PIFs$Ola_p[I]), log(data_with_PIFs$Combined_p[I]), col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(log(data_with_PIFs$Ola_p[I]), log(data_with_PIFs$Combined_p[I]), col = colr)
    }
    else {
      colr <- vus_col 
    }
    
    if(abs(data_with_PIFs$Ola_p[I] - data_with_PIFs$Combined_p[I]) > PIF_difference_thresh) {
      text(log(data_with_PIFs$Ola_p[I]) - .18, log(data_with_PIFs$Combined_p[I]) - .25, labels = data_with_PIFs[I,1], cex = .7, col = colr)
    }
    
  }
  
}



PIF_PIF_plot_HAT_vs_Ola <- function(data_with_PIFs, threshs) { # this uses log scale
  # Plots PIF(HAT) -- y-axis -- against PIF(Ola)
  
  vus_col <- "gray50"
  PIF_difference_thresh = .05 
  
  dev.new()
  par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
  
  plot(log(data_with_PIFs$Ola_p), log(data_with_PIFs$HAT_p), col = vus_col, main = "", xlab = "", ylab = "",
       xlim = c(-10.5, .5), ylim = c(-10.5,.5))
  lines(c(-10.5, .5), c(-10.5, .5), lwd = .1)
  abline(h = log(threshs[1]), lty = 2)
  abline(h = log(threshs[2]), lty = 2)
  abline(v = log(threshs[1]), lty = 2)
  abline(v = log(threshs[2]), lty = 2)
  
  title(xlab = "log( PIF(Olaparib assay only) )", line=2.6) # vertical positioning here
  title(ylab = "log( PIF(HAT assay only) )", line=2.6) # the "line" param is horizontal positioning
  # legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
  #        col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(log(data_with_PIFs$Ola_p[I]), log(data_with_PIFs$HAT_p[I]), col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(log(data_with_PIFs$Ola_p[I]), log(data_with_PIFs$HAT_p[I]), col = colr)
    }
    else {
      colr <- vus_col 
    }
    
    if(abs(data_with_PIFs$Ola_p[I] - data_with_PIFs$HAT_p[I]) > PIF_difference_thresh) {
      text(log(data_with_PIFs$Ola_p[I]) - .18, log(data_with_PIFs$HAT_p[I]) - .25, labels = data_with_PIFs[I,1], cex = .7, col = colr)
    }
    
  }
  
}



PIF_PIF_plot_HAT_vs_Cis <- function(data_with_PIFs, threshs) { # this uses log scale
  # Plots PIF(HAT) -- y-axis -- against PIF(Cis)
  
  vus_col <- "gray50"
  PIF_difference_thresh = .05 
  
  dev.new()
  par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
  
  plot(log(data_with_PIFs$Cis_p), log(data_with_PIFs$HAT_p), col = vus_col, main = "", xlab = "", ylab = "",
       xlim = c(-10.5, .5), ylim = c(-10.5,.5))
  lines(c(-10.5, .5), c(-10.5, .5), lwd = .1)
  abline(h = log(threshs[1]), lty = 2)
  abline(h = log(threshs[2]), lty = 2)
  abline(v = log(threshs[1]), lty = 2)
  abline(v = log(threshs[2]), lty = 2)
  
  title(xlab = "log( PIF(Cisplatin assay only) )", line=2.6) # vertical positioning here
  title(ylab = "log( PIF(HAT assay only) )", line=2.6) # the "line" param is horizontal positioning
  # legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
  #        col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(log(data_with_PIFs$Cis_p[I]), log(data_with_PIFs$HAT_p[I]), col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(log(data_with_PIFs$Cis_p[I]), log(data_with_PIFs$HAT_p[I]), col = colr)
    }
    else {
      colr <- vus_col 
    }
    
    if(abs(data_with_PIFs$Cis_p[I] - data_with_PIFs$HAT_p[I]) > PIF_difference_thresh) {
      text(log(data_with_PIFs$Cis_p[I]) - .18, log(data_with_PIFs$HAT_p[I]) - .25, labels = data_with_PIFs[I,1], cex = .7, col = colr)
    }
    
  }
  
}




PIF_PIF_plot_Cis_vs_Ola <- function(data_with_PIFs, threshs) { # this uses log scale
  # Plots PIF(Cis) -- y-axis -- against PIF(Ola)
  
  vus_col <- "gray50"
  PIF_difference_thresh = .05 
  
  dev.new()
  par(cex.main = 1.5, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L")
  
  plot(log(data_with_PIFs$Ola_p), log(data_with_PIFs$Cis_p), col = vus_col, main = "", xlab = "", ylab = "",
       xlim = c(-10.5, .5), ylim = c(-10.5,.5))
  lines(c(-10.5, .5), c(-10.5, .5), lwd = .1)
  abline(h = log(threshs[1]), lty = 2)
  abline(h = log(threshs[2]), lty = 2)
  abline(v = log(threshs[1]), lty = 2)
  abline(v = log(threshs[2]), lty = 2)
  
  title(xlab = "log( PIF(Olaparib assay only) )", line=2.6) # vertical positioning here
  title(ylab = "log( PIF(Cisplatin assay only) )", line=2.6) # the "line" param is horizontal positioning
  # legend(.6,.3, legend=c("Pathogenic", "Benign", "VUS"), 
  #        col=c("red", "blue", vus_col), lty=c(0,0,0), pch = c(1,1,1), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(log(data_with_PIFs$Ola_p[I]), log(data_with_PIFs$Cis_p[I]), col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(log(data_with_PIFs$Ola_p[I]), log(data_with_PIFs$Cis_p[I]), col = colr)
    }
    else {
      colr <- vus_col 
    }
    
    if(abs(data_with_PIFs$Ola_p[I] - data_with_PIFs$Cis_p[I]) > PIF_difference_thresh) {
      text(log(data_with_PIFs$Ola_p[I]) - .18, log(data_with_PIFs$Cis_p[I]) - .25, labels = data_with_PIFs[I,1], cex = .7, col = colr)
    }
    
  }
  
}



# (END PIF-PIF plotting functions) -----------------------------------------------------------------------





init_data_check <- function(data_fr) {
  # performs an initial data check

  total_num <- nrow(data_fr)
  num_experim <- nrow(subset(data_fr,classification == "Experimental"))
  num_path <- nrow(subset(data_fr,classification == "Pathogenic"))
  num_neut <- nrow(subset(data_fr,classification == "Neutral"))

  pl(" ")
  pl("Pathogenic, Benign, and VUS variant numbers:")
  print(num_path)
  print(num_neut)
  print(num_experim)
  pl("Total number of labeled variants (Pathogenic + Benign):")
  print(num_path + num_neut)
  if (num_experim + num_path + num_neut == total_num) {pl("")}
  else {pl("Initial data check failed...")}
  pl("Total number of variants in the data set:")
  print(total_num)
  pl(" ")
  pl(" ")

}




init_distrib_check <- function(data) {
  # checks if the data distributions for Pathogenic and Benign are normal. Does normality tests, Q-Q plots

  data_path <- subset(data, classification == "Pathogenic")
  data_neut <- subset(data, classification == "Neutral")
  
  Ps_neut <- 1:3
  Ps_path <- 1:3
  
  for (I in 3:5) {
    Ps_neut[I - 2] <- lillie.test(data_neut[[I]])$p.value
    Ps_path[I - 2] <- lillie.test(data_path[[I]])$p.value
  }
  
  pl(" ------------ Normality-test P-values (unadjusted) --------------------- ")
  pl("Neutral variants:")
  print(Ps_neut)
  pl("Pathogenic variants:")
  print(Ps_path)
  pl(" ")
  
  
  # ------------------ Q-Q plots ------------------------------------
  
  titles_path <- c("Pathogenic: HAT", "Pathogenic: Cisplatin", "Pathogenic: Olaparib")
  titles_ben <- c("Benign: HAT", "Benign: Cisplatin", "Benign: Olaparib")
  
  dev.new()
  par(cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2, font.lab = 1, font = 1, bty = "L",font.main = 1)
  par(mfcol = c(3,2), oma = c(0,1,0,0))
  
  for (I in 3:5) {
    qqnorm(data_neut[[I]], pch = 1, col = "blue", main = titles_ben[I - 2], ylim = c(-10,2.5) )
    qqline(data_neut[[I]], col = "blue", lwd = 1.5)
    
  }
  
  for (I in 3:5) {
    qqnorm(data_path[[I]], pch = 1, col = "red", main = titles_path[I - 2], ylim = c(-10,2.5))
    qqline(data_path[[I]], col = "red", lwd = 1.5)
  }
  
}




fold_CV_dataset_gen <- function(data_set,num_folds) {
  # This function generates num_folds data sets for num_folds-fold CV (cross-validation). 
  # These data sets are then used in K-fold_CV() --
  # the function that actually does the CV and spits out its accuracy (averaged over the folds).
  
  
  # Subsets of labelled data.
  # this involves randomization.
  train_set <- subset(data_set, classification == "Neutral" | classification == "Pathogenic")
  N_train <- nrow(train_set)
  rnd_order <- sample(1:N_train)
  
  dsets <- list()
  fold_size <- floor(N_train/num_folds)
  
  for (I in 1:num_folds) { # fold number

    fold_start <- (I - 1)*fold_size + 1
    fold_end <- fold_start + fold_size - 1

    dsets[[I]] <- train_set[rnd_order[fold_start:fold_end],]

  }
  

  return(dsets)
  
}



K_fold_CV <- function(data_set,num_folds,rel_tol,meth,thresh_vec,analysis_var) {
  # performs cross-validation
  
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



plot_var_PIFs <- function(datafr_total, list_intermed) {
  # This function plots one-variable PIF vs FS for HAT, Cis, and Ola.

  vus_col <- "gray50"
  
  # x-axis limits
  m1 <- min(datafr_total$HAT)
  m2 <- min(datafr_total$Cis)
  m3 <- min(datafr_total$Ola)
  m <- min(c(m1, m2, m3))
  M1 <- max(datafr_total$HAT)
  M2 <- max(datafr_total$Cis)
  M3 <- max(datafr_total$Ola)
  M <- min(c(M1, M2, M3))
  
  dev.new()
  par(cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.3, font.lab = 1, font = 1, bty = "L",font.main = 1)
  par(mfcol=c(3,1), oma = c(0,1,0,0))
  par(mar = c(7, 4.1, 4.1, 2.1)) 

  srt <- arrange(datafr_total, HAT)
  plot(srt$HAT, srt$HAT_p, xlab = "Functional score", ylab = "PIF", col = vus_col, main = "HAT assay only",
       ylim = c(0,1), xlim = c(m,M)  )
  legend(0,1, legend=c("Pathogenic", "Benign", "VUS", "Indeterminate"), 
         col=c("red", "blue", vus_col, "magenta"), lty=c(0,0,0,0), pch = c(1,1,1,1), cex=1)
  for (I in 1:nrow(srt)) {
    if (srt[I,2] == "Pathogenic") { points(srt$HAT[I], srt$HAT_p[I], col = "red") }
    if (srt[I,2] == "Neutral") { points(srt$HAT[I], srt$HAT_p[I], col = "blue") }
    if ( any(srt[[I,1]] == list_intermed) == TRUE ) { points(srt$HAT[I], srt$HAT_p[I], col = "magenta", cex = 2)  }
  }

  srt <- arrange(datafr_total, Cis)
  plot(srt$Cis, srt$Cis_p, xlab = "Functional score", ylab = "PIF", col = vus_col, main = "Cisplatin assay only",
       ylim = c(0,1), xlim = c(m,M))
  for (I in 1:nrow(srt)) {
    if (srt[I,2] == "Pathogenic") { points(srt$Cis[I], srt$Cis_p[I], col = "red") }
    if (srt[I,2] == "Neutral") { points(srt$Cis[I], srt$Cis_p[I], col = "blue") }
    if ( any(srt[I,1] == list_intermed) ) { points(srt$Cis[I], srt$Cis_p[I], col = "magenta", cex = 2)  }
  }

  srt <- arrange(datafr_total, Ola)
  plot(srt$Ola, srt$Ola_p, xlab = "Functional score", ylab = "PIF", col = vus_col, main = "Olaparib assay only",
       ylim = c(0,1), xlim = c(m,M))
  for (I in 1:nrow(srt)) {
    if (srt[I,2] == "Pathogenic") { points(srt$Ola[I], srt$Ola_p[I], col = "red") }
    if (srt[I,2] == "Neutral") { points(srt$Ola[I], srt$Ola_p[I], col = "blue") }
    if ( any(srt[I,1] == list_intermed) ) { points(srt$Ola[I], srt$Ola_p[I], col = "magenta", cex = 2)  }
  }

}




vis_Cis_Ola <- function(data_with_PIFs) {
  # This function compares Cis and Ola FS, does regression.
  
  vus_col <- "gray50"
  
  dev.new()
  par(cex.main = 1.4, cex.lab = 1.2, cex.axis = 1, font.lab = 1, font = 1, bty = "L", font.main = 1)
  plot(data_with_PIFs$Cis, data_with_PIFs$Ola, col = vus_col, xlab = "", ylab = "", xlim = c(-11,3), ylim = c(-11,3))
  
  title(xlab = "Cisplatin functional score", line=2.4) # vertical positioning here
  title(ylab = "Olaparib functional score", line=2.4) # the "line" param is horizontal positioning
  legend(-5,-7, legend=c("Pathogenic", "Benign", "VUS","Regression line (P < 0.00001)", "95% prediction band"), 
         col=c("red", "blue", vus_col,"black","black"), lty=c(0,0,0,1,2), pch = c(1,1,1,NA,NA), cex=1)
  
  
  for (I in 1:nrow(data_with_PIFs)) {
    
    if (data_with_PIFs[I,2] == "Neutral") {
      colr <- "blue"
      points(data_with_PIFs$Cis[I], data_with_PIFs$Ola[I], col = colr)
    }
    else if (data_with_PIFs[I,2] == "Pathogenic") {
      colr <- "red"
      points(data_with_PIFs$Cis[I], data_with_PIFs$Ola[I], col = colr)
    }
    else {
      colr <- vus_col 
    }
    
  }
  
  fit <- lm(formula = Ola ~ Cis, data = data_with_PIFs) 
  print(summary(fit))
  new_data <- data.frame(Cis = seq(from = -11, to = 3, by = .01))
  pred <- predict(object = fit, newdata = new_data,  interval = "predict", level = .95)
  lines(new_data[[1]],pred[,1]) # fit
  lines(new_data[[1]], pred[,2], lty=2) # lower pred interval
  lines(new_data[[1]], pred[,3], lty=2) # upper pred interval

}




plot_hist_distr <- function(dset1, init_pars, fit_pars) {
  # For HAT, Cis, and Ola, this function plots data histograms and the corresponding mixture p.d.f.
  
  dev.new()
  par(mfcol=c(3,1), oma = c(0,1,0,0))
  par(cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.3, font.lab = 1, font = 1, bty = "L",font.main = 1)
  par(mar = c(7, 4.1, 4.1, 2.1)) 

  
  num_bins <- 35
  c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
  # c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
  col_v <- c1


  right_boundary <-  6
  left_boundary <- -11
  upper_boundary <- 1

  
  # params reformatting 
  init_params <- list()
  init_params[[1]] <- init_pars[1:5] # HAT
  init_params[[2]] <- init_pars[6:10] # Cis
  init_params[[3]] <- init_pars[11:15] # Ola
  
  fit_params <- list()
  fit_params[[1]] <- fit_pars[1:5] # HAT
  fit_params[[2]] <- fit_pars[6:10] # Cis
  fit_params[[3]] <- fit_pars[11:15] # Ola
  
  main_title <- c("HAT", "Cisplatin", "Olaparib")
  # histograms and densities
  for (I in 1:3) {
    
    # histogram
    var_number <- I + 2
    M <- max(dset1[[var_number]])
    ax <- seq(from = left_boundary, to = right_boundary, length.out = num_bins + 1)
    H1 <- hist(dset1[[var_number]], breaks = ax, probability = TRUE,
               col = col_v, main = main_title[I],
               ylim = c( 0, upper_boundary), xlim = c(left_boundary, right_boundary),  
               ylab = "Probability density", xlab = "Functional score")
    
    
    x_variable <- seq(from = left_boundary, to = right_boundary, by = .01)
    
    # density before fitting
    pars <- init_params[[I]]
    prior <- pars[1]
    mix_density <- prior*dnorm(x_variable, mean = pars[2],sd =  pars[3]) +
      (1 - prior)*dnorm(x_variable, mean = pars[4],sd =  pars[5])
    lines(x_variable, mix_density, col = "black", lty = 2, lwd = 2)
    
    
    # density after fitting
    pars <- fit_params[[I]]
    prior <- pars[1]
    mix_density <- prior*dnorm(x_variable, mean = pars[2],sd =  pars[3]) +
      (1 - prior)*dnorm(x_variable, mean = pars[4],sd =  pars[5])
    lines(x_variable, mix_density, col = "black", lwd = 2)
    
    if (I == 1) {
      legend(-11,1, legend=c("Mixture model before fitting", "Mixture model after fitting", "Experimental-data histogram"), 
             col=c("black", "black","black"), lty=c(2,1,NA), pch = c(NA,NA,22), 
             cex=1.1, pt.bg = c(col_v,col_v,col_v), pt.cex = c(1,1,2))
    }
    
    
  }

  
}





# -------------------------------------------------------------------------
# -------------------------------- DATA LOADING ---------------------------
# -------------------------------------------------------------------------

# The main data set will be contained in the dataframe "datafr".

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



# ------------------ Initial data checking --------------------------------

init_data_check(datafr)

pl(" ")
readline(prompt = "Press a key...")
pl(" ")

init_distrib_check(datafr)



# ----------------- Fitting the model and plotting PIFs ----------------

pl(" --------------------------------------------------------- ")
pl("Accuracy on the training data (fitting to full data set):")
pl(" --------------------------------------------------------- ")
result <- accuracy_on_full_data(datafr,rel_tolerance,opt_method,threshs,analysis_var)
init_params_final <- result[[4]]
fit_params_final <- result[[5]]
acc <- result[[2]]

if (analysis_var == 4) {
  plot_title <- "HAT + Cisplatin + Olaparib"
} else if (analysis_var == 1) {
  plot_title <- "HAT"
} else if (analysis_var == 2) {
  plot_title <- "Cisplatin"
} else if (analysis_var == 3) {plot_title <- "Olaparib"}
plot_PIFs(datafr,result[[1]], plot_title, threshs, analysis_var, DATA_FOLDER, plot_PIFs_legend_flag) # result[[1]] is the PIFs for the variables in analysis_var
data_and_probs <- result[[3]] # data and all PIFs in this dataframe
pl(" ")
pl("Fitting accuracy (%):")
print(acc*100)


# ------------------- PIF-PIF plots -----------------------------------------


PIF_PIF_plot(data_and_probs, threshs)
PIF_PIF_plot_all_vs_Cis(data_and_probs, threshs)
PIF_PIF_plot_all_vs_Ola(data_and_probs, threshs)
PIF_PIF_plot_HAT_vs_Cis(data_and_probs, threshs)
PIF_PIF_plot_HAT_vs_Ola(data_and_probs, threshs)
PIF_PIF_plot_Cis_vs_Ola(data_and_probs, threshs)


# # # --------------------------- CROSS-VALIDATION ---------------------------


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


# # ------------------------ Parameter analysis ----------------------------------

pl(" ------------------------------------------------------------------ ")
pl("      Parameter analysis: initial (1st row) vs fitted (2nd row) mixture-density (two normals) params")
pl(" ------------------------------------------------------------------ ")

template <- matrix(1:30, nrow = 2, ncol = 15)
template[1,] <- init_params_final
template[2,] <- fit_params_final
print(template)
pl(" ")
pl(" ")

# # ------------------------- plot var PIFs -----------------------------------------
# !!!! Not included in the final set of figures!!!!
# INSERTED ARE THE REAL INTDETERM VARIANT NAMES FOR analysis_var = 4
# intermed <- c("c.C7142A", "c.C8710T","c.C8668A", "c.T7068G", "c.T8692G")
# plot_var_PIFs(data_and_probs,intermed)


# # ----------------- plot Cis-Ola regression ----------------------------------------

pl(" ------------------------------------------------------------------ ")
pl("                      Cis vs Ola comparison          ")
pl(" ------------------------------------------------------------------ ")
# !!!! Not included in the final set of figures!!!!
# vis_Cis_Ola(data_and_probs)
# pl(" ")
# pl(" ")


# # ---------------- plot histograms and distribution densities ---------------

plot_hist_distr(datafr,init_params_final,fit_params_final)



# # ----------------------------- saving variable-level and PIF table ------------------------

# changing column names a bit and accommodating new terminology
colnames(data_and_probs) <- c("HGVS nucleotide","Functional classification","HAT FS","Cisplatin FS", "Olaparib FS",
                              "HAT_PIF", "Cisplatin_PIF", "Olaparib_PIF", "(HAT + Cisplatin + Olaparib)_PIF")
data_and_probs[data_and_probs == "Experimental"] <- "VUS"
data_and_probs[data_and_probs == "Neutral"] <- "Benign"

out_file <- paste(DATA_FOLDER, "BRCA2_clinvar_variants_Main_PIFs.xlsx", sep = "")
write_xlsx(data_and_probs, path = out_file)









