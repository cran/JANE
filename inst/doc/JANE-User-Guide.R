## ----include = FALSE----------------------------------------------------------

knitr::opts_chunk$set(
  strip.white = FALSE,
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 4,
  fig.align = "center"
)

options(knitr.table.format = "html")
library(JANE) # load JANE
set.seed(1) # for reproducibility


## ----GS_install, eval = FALSE-------------------------------------------------
# # Current release from CRAN
# install.packages("JANE")
# 
# # Development version from GitHub
# # install.packages("devtools")
# devtools::install_github("a1arakkal/JANE")
# 
# # Development version with vignettes
# devtools::install_github("a1arakkal/JANE", build_vignettes = TRUE)

## ----eval = FALSE, message = FALSE--------------------------------------------
# # Load JANE
# library(JANE)
# 
# # Vignette
# RShowDoc("JANE-User-Guide", package = "JANE")

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 3,
#      model = "NDH",
#      noise_weights = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2:5,
#      K = 3,
#      model = "RS",
#      noise_weights = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2:5,
#      K = 2:10,
#      model = "RSR",
#      noise_weights = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 5,
#      model = "RS",
#      noise_weights = TRUE,
#      family = "poisson",
#      guess_noise_weights = 1L)

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 5,
#      model = "RS",
#      noise_weights = TRUE,
#      family = "lognormal",
#      guess_noise_weights = -3.5) # log-scale mean

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 5,
#      model = "RSR",
#      noise_weights = TRUE,
#      family = "bernoulli",
#      guess_noise_weights = 0.2) # expected noise edge proportion

## ----eval = FALSE-------------------------------------------------------------
# # Specify the 3 x 2 matrix containing the 1 x 2 mean vectors of the 3 bivariate normals
# mus <- matrix(c(-1,-1,  # Mean vector 1
#                  1,-1,  # Mean vector 2
#                  1, 1), # Mean vector 3
#               nrow = 3,
#               ncol = 2,
#               byrow = TRUE)
# 
# # Specify the 2 x 2 x 3 array containing the 2 x 2 precision matrices of the 3 bivariate normals
# omegas <- array(c(diag(rep(7,2)),  # Precision matrix 1
#                   diag(rep(7,2)),  # Precision matrix 2
#                   diag(rep(7,2))), # Precision matrix 3
#                   dim = c(2,2,3))
# 
# # Simulate a network
# sim_A(N = 100L,
#       model = "NDH",
#       mus = mus,
#       omegas = omegas,
#       p = rep(1/3, 3),
#       params_LR = list(beta0 = 1.0),
#       remove_isolates = TRUE)

## ----eval = !identical(Sys.getenv("IN_PKGDOWN"), "true"), message=FALSE-------
# Specify the 3 x 2 matrix containing the 1 x 2 mean vectors of the 3 bivariate normals
mus <- matrix(c(-1,-1,  # Mean vector 1
                 1,-1,  # Mean vector 2
                 1, 1), # Mean vector 3
              nrow = 3,
              ncol = 2, 
              byrow = TRUE)

# Specify the 2 x 2 x 3 array containing the 2 x 2 precision matrices of the 3 bivariate normals
omegas <- array(c(diag(rep(7,2)),  # Precision matrix 1
                  diag(rep(7,2)),  # Precision matrix 2
                  diag(rep(7,2))), # Precision matrix 3
                  dim = c(2,2,3))

desired_density <- 0.1 # Target network density
min_density <- desired_density * 0.99  # Lower bound for acceptable density
max_density <- desired_density * 1.01  # Upper bound for acceptable density
n_act <- 100L # Number of actors in the network

density <- Inf # Initialize density to enter while loop
beta0 <- 0.5 # Initial value for intercept parameter
n_while_loop <- 0 # Counter for outer loop iterations
max_its <- 100 # Maximum number of iterations
change_beta0 <- 0.1  # Amount to adjust beta0 by

# Adjust beta0 until simulated network has the approximate desired density
while(! (density >= min_density & density <= max_density) ){
  
  if(n_while_loop>max_its){
    break
  }
  
  n_retry_isolate <- 0
  retry_isolate <- T
  
  # Retry until a network with no isolates is generated (this while loop is optional)
  while(retry_isolate){
    
    sim_data <- sim_A(N = n_act, 
                      model = "NDH",
                      mus = mus, 
                      omegas = omegas, 
                      p = rep(1/3, 3), 
                      params_LR = list(beta0 = beta0),
                      remove_isolates = TRUE)
    
    n_retry_isolate <- n_retry_isolate + 1
    
    # Accept network if no isolates remain, or if retried more than 10 times at the same beta0
    if(nrow(sim_data$A) == n_act | n_retry_isolate>10){
      retry_isolate <- F
    }
    
  }
  
  # Compute network density
  density <- igraph::graph.density(
    igraph::graph_from_adjacency_matrix(sim_data$A, mode = "undirected")
    )

  # Adjust beta0 based on density feedback
  if (density > max_density)  {
    beta0 <- beta0 - change_beta0
  }
  
  if (density < min_density)  {
    beta0 <- beta0 + change_beta0
  }
  
  n_while_loop <- n_while_loop + 1
  
}    

A <- sim_data$A # Final simulated adjacency matrix
igraph::graph.density(igraph::graph_from_adjacency_matrix(A, mode = "undirected")) # Verify density

## ----eval = FALSE-------------------------------------------------------------
# # Specify the 3 x 2 matrix containing the 1 x 2 mean vectors of the 3 bivariate normals
# mus <- matrix(c(-1,-1,  # Mean vector 1
#                  1,-1,  # Mean vector 2
#                  1, 1), # Mean vector 3
#               nrow = 3,
#               ncol = 2,
#               byrow = TRUE)
# 
# # Specify the 2 x 2 x 3 array containing the 2 x 2 precision matrices of the 3 bivariate normals
# omegas <- array(c(diag(rep(7,2)),  # Precision matrix 1
#                   diag(rep(7,2)),  # Precision matrix 2
#                   diag(rep(7,2))), # Precision matrix 3
#                   dim = c(2,2,3))
# 
# # Simulate a network
# sim_A(N = 100L,
#       model = "RSR",
#       family = "poisson",
#       mus = mus,
#       omegas = omegas,
#       p = rep(1/3, 3),
#       params_LR = list(beta0 = 1),
#       params_weights = list(beta0 = 2),
#       noise_weights_prob = 0.1,
#       mean_noise_weights = 1,
#       remove_isolates = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 2:10,
#      model = "RSR",
#      noise_weights = FALSE,
#      control = list(IC_selection = "Total_BIC"))

## ----eval = FALSE-------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 3,
#      model = "RSR",
#      noise_weights = FALSE,
#      control = list(IC_selection = "Total_ICL",
#                     n_start = 10))

## ----eval = FALSE-------------------------------------------------------------
# # Specify starting values
# D <- 3
# K <- 5
# N <- nrow(sim_data$A)
# n_interior_knots <- 5L
# U <- matrix(stats::rnorm(N*D), nrow = N, ncol = D)
# omegas <- stats::rWishart(n = K, df = D+1, Sigma = diag(D))
# mus <- matrix(stats::rnorm(K*D), nrow = K, ncol = D)
# p <- extraDistr::rdirichlet(n = 1, rep(3,K))[1,]
# Z <-  extraDistr::rdirichlet(n = N, alpha = rep(1, K))
# beta <- stats::rnorm(n = 1 + 2*(1 + n_interior_knots))
# 
# my_starting_values <- specify_initial_values(A = network_adjacency_matrix,
#                                              D = D,
#                                              K = K,
#                                              model = "RSR",
#                                              n_interior_knots = n_interior_knots,
#                                              U = U,
#                                              omegas = omegas,
#                                              mus = mus,
#                                              p = p,
#                                              Z = Z,
#                                              beta = beta)
# 
# # Run JANE using my_starting_values (no need to specify D and K as function will
# # determine those values from my_starting_values)
# JANE(A = network_adjacency_matrix,
#      initialization = my_starting_values,
#      model = "RSR",
#      noise_weights = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# # Specify prior hyperparameters
# D <- 3
# K <- 5
# n_interior_knots <- 5L
# a <- rep(1, D)
# b <- 3
# c <- 4
# G <- 10*diag(D)
# nu <- rep(2, K)
# e <- rep(0.5, 1 + (n_interior_knots + 1))
# f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
# 
# my_prior_hyperparameters <- specify_priors(D = D,
#                                            K = K,
#                                            model = "RS",
#                                            n_interior_knots = n_interior_knots,
#                                            a = a,
#                                            b = b,
#                                            c = c,
#                                            G = G,
#                                            nu = nu,
#                                            e = e,
#                                            f = f)
# 
# # Run JANE using supplied prior hyperparameters (no need to specify D and K
# # as function will determine those values from my_prior_hyperparameters)
# JANE(A = network_adjacency_matrix,
#      initialization = "GNN",
#      model = "RS",
#      noise_weights = FALSE,
#      control = list(priors = my_prior_hyperparameters))

## ----eval = FALSE-------------------------------------------------------------
# # Specify first set of prior hyperparameters
# D <- 3
# K <- 5
# n_interior_knots <- 5L
# a <- rep(1, D)
# b <- 3
# c <- 4
# G <- 10*diag(D)
# nu <- rep(2, K)
# e <- rep(0.5, 1 + (n_interior_knots + 1))
# f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
# 
# my_prior_hyperparameters_1 <- specify_priors(D = D,
#                                              K = K,
#                                              model = "RS",
#                                              n_interior_knots = n_interior_knots,
#                                              a = a,
#                                              b = b,
#                                              c = c,
#                                              G = G,
#                                              nu = nu,
#                                              e = e,
#                                              f = f)
# 
# # Specify second set of prior hyperparameters
# D <- 2
# K <- 3
# n_interior_knots <- 5L
# a <- rep(1, D)
# b <- 3
# c <- 4
# G <- 10*diag(D)
# nu <- rep(2, K)
# e <- rep(0.5, 1 + (n_interior_knots + 1))
# f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
# 
# my_prior_hyperparameters_2 <- specify_priors(D = D,
#                                              K = K,
#                                              model = "RS",
#                                              n_interior_knots = n_interior_knots,
#                                              a = a,
#                                              b = b,
#                                              c = c,
#                                              G = G,
#                                              nu = nu,
#                                              e = e,
#                                              f = f)
# 
# # Create nested list
# my_prior_hyperparameters <- list(my_prior_hyperparameters_1,
#                                  my_prior_hyperparameters_2)
# 
# # Run JANE using supplied prior hyperparameters (no need to specify D and K
# # as function will determine those values from my_prior_hyperparameters)
# JANE(A = network_adjacency_matrix,
#      initialization = "GNN",
#      model = "RS",
#      noise_weights = FALSE,
#      control = list(priors = my_prior_hyperparameters))

## ----eval = FALSE-------------------------------------------------------------
# # Specify prior hyperparameters as unevaluated calls
# n_interior_knots <- 5L
# e <- rep(0.5, 1 + (n_interior_knots + 1))
# f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
# my_prior_hyperparameters <- specify_priors(model = "RS",
#                                            n_interior_knots = n_interior_knots,
#                                            a = quote(rep(1, D)),
#                                            b = b,
#                                            c = quote(D + 1),
#                                            G = quote(10*diag(D)),
#                                            nu = quote(rep(2, K)),
#                                            e = e,
#                                            f = f)
# 
# # Run JANE using supplied prior hyperparameters
# JANE(A = network_adjacency_matrix,
#      D = 2:5,
#      K = 2:10,
#      initialization = "GNN",
#      model = "RS",
#      noise_weights = FALSE,
#      control = list(priors = my_prior_hyperparameters))

## ----eval=FALSE---------------------------------------------------------------
# # Set up parallel plan with 5 workers (cores)
# future::plan(future::multisession, workers = 5)
# 
# # Run JANE in parallel
# res_parallel <- JANE(A = network_adjacency_matrix,
#                      D = 2:5,
#                      K = 3:10,
#                      initialization = "GNN",
#                      model = "RSR")
# 
# # Reset to sequential processing
# future::plan(future::sequential)

## ----eval=FALSE---------------------------------------------------------------
# JANE(A = network_adjacency_matrix,
#      D = 2,
#      K = 3,
#      initialization = "GNN",
#      model = "RSR",
#      case_control = TRUE,
#      control = list(n_control = 20))

## ----eval = FALSE-------------------------------------------------------------
# print(res)

## ----eval = FALSE-------------------------------------------------------------
# summary(res)

## ----eval = FALSE-------------------------------------------------------------
# summary(res, true_labels = true_labels_vec)
# summary(res, initial_values = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# plot(res)

## ----eval = FALSE-------------------------------------------------------------
# plot(res, type = "misclassified", true_labels = true_labels_vec)

## ----eval = FALSE-------------------------------------------------------------
# plot(res, type = "uncertainty")

## ----eval = FALSE-------------------------------------------------------------
# plot(res, type = "trace_plot")

