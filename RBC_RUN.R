library(tidyverse)
library(MASS)
library(Matrix)
library(Rcpp)
library(nleqslv)
library(gEcon)



rbc_object <- make_model('Desktop/RBC/rbc_test.gcn')



rbc_object <- set_free_par(
  model = rbc_object,
  list(
    rho = 0.99,
    #beta = 0.99,
    #alpha = ...,
    psi = 0#,
    #eta = ...
  )
)



rbc_object <- steady_state(rbc_object)
rbc_object <- solve_pert(rbc_object)



get_par_values(rbc_object)
get_ss_values(rbc_object)
re_solved(rbc_object)



path <- matrix(
  c(0.05),
  nrow = 1,
  ncol = 1
)

rbc_object_sim <- simulate_model(
  rbc_object,
  variables  = c('K_s'),
  shocks     = 'epsilon_Z',
  shock_path = path,
  sim_length = 100
)

plot_simulation(rbc_object_sim)



random_rbc_object_sim <- random_path(
  rbc_object,
  variables = c('Y','W', 'r', 'C', 'L_s'),
  sim_length = 100
)

plot_simulation(random_rbc_object_sim)


###MONTE CARLO????







