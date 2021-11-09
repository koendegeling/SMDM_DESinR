##-----------------------------------------------------------------------------------------##
##                                                                                         ##
##          Discrete Event Simulation in R to support healthcare decision making           ##
##                                                                                         ##
##                by Koen Degeling, Michiel van de Ven, Hendrik Koffijberg                 ##
##                                                                                         ##
##-----------------------------------------------------------------------------------------##
#
# This script has been prepared for the the "Discrete Event Simulation in R to support 
# healthcare decision making" short course presented to the Society for Medical Decision 
# Making (SMDM). It implements a discrete event simulation (DES) that models the pathway of 
# individuals who are referred for a potential knee replacement, which substantially improves
# their quality of life. These individuals have to go through several steps/checks before 
# they may potentially receive surgery for the knee replacement to be placed. Although the
# case study does not enforce a maximum number of surgeries (or other resources like tests),
# it is demonstrated how resources are implemented so that the example can be easily extended
# to consider resource constraints. More details about the hypothetical case study are 
# provided throughout the code and in the corresponding slide deck that is available from the
# GitHub repository (see below).
#
# Please reach out with any questions or suggestions (koen.degeling@unimelb.edu.au) and check
# the GitHub repository for potential updates and extensions:
# - https://github.com/koendegeling/SMDM_DESinR
#
# This 3_probabilistic_analysis.R script demonstrates how the health economic simulation 
# model that has been defined in the 2_health_economics.R script can be analyzed through a 
# probabilistic analysis.
#
# This script contains the following sections:
# 1. INITIALIZATION   sets up the environment for the analysis
# 2. FUNCTIONS        defines functions that are used to obtain parameter values and run the
#                     simulation model for a given set of parameter values
# 3. DETERMINISTIC    running the simulation model for the deterministic parameter values
# 4. PROBABILISTIC    runiing the probabilistic analysis of the simulation model
#
# This script was written using R Studio version 1.4.1106 and R version 4.0.3. The versions
# of the packages used are specified in where they are loaded in the script.
#
# Please note that time is modeled in years in this script. Furthermore, this script uses the
# 'snake_case' to define objects, such as parameters and functions, whereas the 'CamelCase'
# is used to refer to individuals' attributes in the simulation. This is done deliberately to
# clearly distinguish between the two.
#
#
#
#
## 1. INITIALIZATION ----

# This section sets up the analysis by installing and loading the required packages.

# Uncomment to clear the Global Environment and Console 
rm(list = ls()); cat("\14");

# Uncomment to install the used packages
#install.packages(pkgs = c("simmer", "flexsurv", "data.table", "MASS"))




## 2. FUNCTIONS ----

# The fn_pars_deterministic() function returns a list with all required parameters to run
# the simulation and their values in the deterministic base case. Using the same random
# number seed, these parameter should yield the same outcomes as those simulated in the 
# 2_health_economics.R script.
fn_pars_deterministic <- function() {
  
  pars <- list(
    
    d_death_shape = 0.1,
    d_death_rate  = 0.005,
    
    d_intake_shape = 2.0,
    d_intake_scale = 0.3,
    t_intake = 1/365.25,
    c_intake = 132,
    p_eligible_intake = 0.7,
    
    d_testing_shape = 2.5,
    d_testing_scale = 0.2,
    t_testing = 1/365.25,
    c_testing = 623,
    
    d_consult_shape = 3.0,
    d_consult_scale = 0.1,
    t_consult = 1/365.25,
    c_consult = 184,
    p_eligible_consult = 0.9,
    
    d_surgery_shape = 1.5,
    d_surgery_scale = 0.5,
    t_surgery = 1/52,
    c_surgery = 8132,
    
    d_recovery_shape = 2.0,
    d_recovery_scale = 2/52,
    t_recovery = 3/12,
    c_recovery = 4576,
    
    d_followup_shape  = 3.0,
    d_followup_scale  = 1.0,
    t_followup = 1/365.25,
    c_followup = 132,
    n_followup_rounds = 5,
    
    u_prior_surgery = 0.6,
    u_after_surgery = 0.9,
    
    p_eligible_intake_1  = 0.770,
    p_eligible_consult_1 = 0.945,
    c_surgery_1  = 12716,
    c_recovery_1 = 6329
    
  )
  
  return(pars)
  
}

# Like the fn_pars_deterministic() function, the fn_pars_probabilistic() function return a
# list of the parameters required to run the simulation and a value for each of them. The
# difference, however, is that for some parameters values in the fn_pars_probabilistic() 
# function are sampled from probability distributions to quantify the parameter uncertainty
# in the probabilistic analyses. Note that, like the deterministic parameter values, the 
# distributions used here are hypothetical.
fn_pars_probabilistic <- function() {
  
  # Library for multivariate normal distributions to sample correlated parameters that 
  # define the time-to-event distributions.
  library(MASS)   # v7.3-53
  
  # For the time-to-event distributions, parameter values are sampled from multivariate
  # normal distributions using the estimated coefficients and variance-covariance matrix
  # that can be obtained from the maximum likelihood estimation of such distributions. Note
  # that some parameters are represented by coefficients that are log-transformed, so those
  # need to be exponentiated to obtain the actual parameter values on real scale. The
  # transformations used here are in line with those of the flexsurv package.
  d_death_est  = c(shape = 0.1039571, rate = -5.3681496)
  d_death_vcov = matrix(data = c(2.060495e-05, -0.0005321442, -5.321442e-04, 0.0157431743),
                        nrow = 2, dimnames = list(c('shape', 'rate'), c('shape', 'rate')))
  
  d_intake_est  = c(shape = 0.7196607, scale = -1.2288401)
  d_intake_vcov = matrix(data = c(0.0012167280, 0.0002492552, 0.0002492552, 0.0005252385),
                         nrow = 2, dimnames = list(c('shape', 'scale'), c('shape', 'scale')))
  
  d_testing_est  = c(shape = 0.9709345, scale = -1.5933248)
  d_testing_vcov = matrix(data = c(0.0011847217, 0.0001934389, 0.0001934389, 0.0003184552),
                          nrow = 2, dimnames = list(c('shape', 'scale'), c('shape', 'scale')))
  
  d_consult_est  = c(shape = 1.076840, scale = -2.326732)
  d_consult_vcov = matrix(data = c(0.0012398980, 0.0001762389, 0.0001762389, 0.0002571628),
                          nrow = 2, dimnames = list(c('shape', 'scale'), c('shape', 'scale')))
  
  d_surgery_est  = c(shape = 0.4705519, scale = -0.7108770)
  d_surgery_vcov = matrix(data = c(0.001181326, 0.0003221130, 0.0003221130, 0.0008682243),
                          nrow = 2, dimnames = list(c('shape', 'scale'), c('shape', 'scale')))
  
  d_recovery_est  = c(shape = 0.6879554, scale = -3.2609574)
  d_recovery_vcov = matrix(data = c(0.0012262444, 0.0002597393, 0.0002597393, 0.0005602357),
                           nrow = 2, dimnames = list(c('shape', 'scale'), c('shape', 'scale')))
  
  d_followup_est  = c(shape = 0.7012124, scale = -3.2337920)
  d_followup_vcov = matrix(data = c(0.0012300962, 0.0002582606, 0.0002582606, 0.0005462214),
                           nrow = 2, dimnames = list(c('shape', 'scale'), c('shape', 'scale')))
  
  # Sampling the joint time-to-event distribution parameters
  d_death_coefs    <- mvrnorm(n = 1, mu = d_death_est, Sigma = d_death_vcov)
  d_intake_coefs   <- mvrnorm(n = 1, mu = d_intake_est, Sigma = d_intake_vcov)
  d_testing_coefs  <- mvrnorm(n = 1, mu = d_testing_est, Sigma = d_testing_vcov)
  d_consult_coefs  <- mvrnorm(n = 1, mu = d_consult_est, Sigma = d_consult_vcov)
  d_surgery_coefs  <- mvrnorm(n = 1, mu = d_surgery_est, Sigma = d_surgery_vcov)
  d_recovery_coefs <- mvrnorm(n = 1, mu = d_recovery_est, Sigma = d_recovery_vcov)
  d_followup_coefs <- mvrnorm(n = 1, mu = d_followup_est, Sigma = d_followup_vcov)
  
  # Probability of being eligible, with the probability sampled from a beta distribution, and
  # the probability for the experimental strategy determined based based on the log odds
  # compared to the control strategy, in which uncertainty is represented.
  p_eligible_intake <- rbeta(n = 1, shape1 = 350, shape2 = 150)
  logodds_eligible_intake   <- log(p_eligible_intake/(1-p_eligible_intake))
  logodds_eligible_intake_1 <- logodds_eligible_intake + rnorm(n = 1, mean = 0.3610133, sd = 0.3610133*0.10)
  p_eligible_intake_1 <- exp(logodds_eligible_intake_1) / (1 + exp(logodds_eligible_intake_1))
  
  p_eligible_consult <- rbeta(n = 1, shape1 = 450, shape2 = 50)
  logodds_eligible_consult   <- log(p_eligible_consult/(1-p_eligible_consult))
  logodds_eligible_consult_1 <- logodds_eligible_consult + rnorm(n = 1, mean = 0.6466272, sd = 0.6466272*0.10)
  p_eligible_consult_1 <- exp(logodds_eligible_consult_1) / (1 + exp(logodds_eligible_consult_1))
  
  # Return the parameters
  pars <- list(
    
    d_death_shape = d_death_coefs['shape'],
    d_death_rate  = exp(d_death_coefs['rate']),
    
    d_intake_shape = exp(d_intake_coefs['shape']),
    d_intake_scale = exp(d_intake_coefs['scale']),
    t_intake = 1/365.25,
    c_intake = 132,
    p_eligible_intake = p_eligible_intake,
    
    d_testing_shape = exp(d_testing_coefs['shape']),
    d_testing_scale = exp(d_testing_coefs['scale']),
    t_testing = 1/365.25,
    c_testing = 623,
    
    d_consult_shape = exp(d_consult_coefs['shape']),
    d_consult_scale = exp(d_consult_coefs['scale']),
    t_consult = 1/365.25,
    c_consult = 184,
    p_eligible_consult = p_eligible_consult,
    
    d_surgery_shape = exp(d_surgery_coefs['shape']),
    d_surgery_scale = exp(d_surgery_coefs['scale']),
    t_surgery = 1/52,
    c_surgery = 8132,
    
    d_recovery_shape = exp(d_recovery_coefs['shape']),
    d_recovery_scale = exp(d_recovery_coefs['scale']),
    t_recovery = 3/12,
    c_recovery = 4576,
    
    d_followup_shape = exp(d_followup_coefs['shape']),
    d_followup_scale = exp(d_followup_coefs['scale']),
    t_followup = 1/365.25,
    c_followup = 132,
    n_followup_rounds = 5,
    
    u_prior_surgery = (1 - rlnorm(n = 1, meanlog = -0.9212659, sdlog = 0.09975135)),
    u_after_surgery = (1 - rlnorm(n = 1, meanlog = -2.30756, sdlog = 0.09975135)),
    
    p_eligible_intake_1  = p_eligible_intake_1,
    p_eligible_consult_1 = p_eligible_consult_1,
    c_surgery_1  = 12716,
    c_recovery_1 = 6329
    
  )
  
  return(pars)
  
}

# The fn_run_DES() function runs the simulation for a certain number of simulated individuals
# specified through the "n_individuals" argument and a set of parameter values specified 
# through the "pars" argument. A random number seed can be set as well through the "seed"
# argument, which is helpful for reproducibility and to reduce unwarranted variation between
# runs with different parameter values. 
#
# If the fn_run_DES() function is used with the deterministic parameter values and random
# number seed, it is basically the same as the analyses performed in the 2_health_economics.R
# script. The benefit, however, is that this function can be used to easily run the model for
# different sets of parameter values, such as in a probabilistic analysis as is demonstrated.
#
# The code within the fn_run_DES() is basically the same as in the 2_health_economics.R
# script, but with some minor changes:
# - comments have been removed to reduce the number of lines of code
# - only the attributes are extracted from the simulation and summarized
# - the mean health economic outcomes are returned rather than patient-level results
# - before returning the results, all objects are removed from the local function environment,
#   except for the output, so that the objects don't accumulate and utilize all the memory if
#   hundreds of runs are performed in the probabilistic analysis.
fn_run_DES <- function(pars, n_individuals, seed = 123) {
  
  with(data = pars, expr = {
    
    fn_eligible_intake <- function() {
      
      out <- if(runif(1) < p_eligible_intake) {0} else {1}
      
      return(out)
      
    }
    
    fn_eligible_consult <- function() {

      out <- if(runif(1) < p_eligible_consult) {0} else {1}
      
      return(out)
      
    }
    
    fn_discount <- function(amount, at = NULL, start = NULL, duration = NULL, rate = 0.03, timefactor = 1) {

      if(is.null(at) & !is.null(start) & !is.null(duration)) {
        
        undisc <- amount * duration;
        disc   <- (1/timefactor) * (amount / (1 + rate) ^ (start*timefactor) / log(1 + rate) * (1 - 1 / (1 + rate) ^ (duration*timefactor)));
        
      } else if(!is.null(at) & is.null(start) & is.null(duration)) {
        
        undisc <- amount
        disc   <- amount / ((1 + rate) ^ (at*timefactor))
        
      } else {
        
        stop('Time argument(s) not specified correctly for fn_discount() function')
        
      }
      
      out    <- c(Undiscounted = undisc, Discounted = disc)
      
      return(out)
      
    }
    
    fn_calculate_impact <- function(CurrentTime, Attrs) {
      
      TimeOfReferral  <- Attrs[1]
      TimeOfSurgery   <- Attrs[2]
      
      TimeOfDeath   <- CurrentTime
      TimeToDeath   <- TimeOfDeath - TimeOfReferral
      TimeToSurgery <- TimeOfSurgery - TimeOfReferral
      
      if(is.na(TimeOfSurgery) == TRUE) {
        CombinedQALYs <- fn_discount(amount = u_prior_surgery, start = TimeOfReferral, duration = TimeToDeath)
      } else {
        QALYs_before  <- fn_discount(amount = u_prior_surgery, start = TimeOfReferral, duration = TimeToSurgery)
        QALYs_after   <- fn_discount(amount = u_after_surgery, start = TimeToSurgery, duration = TimeToDeath - TimeToSurgery)
        
        CombinedQALYs <- QALYs_before + QALYs_after
        
      }
      
      out <- c(TimeOfDeath, TimeToDeath, TimeToSurgery, CombinedQALYs)
      
      return(out)
      
    }
    
    fn_summarise <- fn_summarize <- function(sim_out, keys = NULL) {
      
      if(is.null(keys)) keys <- unique(sim_out$key);
      
      df <- as.data.table(sim_out)[key %in% keys];
      setorder(df, name, time);
      df <- df[, .(value = value[.N]), by = list(name, key)];
      df <- dcast(df, name~key, value.var = "value");
      setcolorder(df, c("name", keys));
      
      return(df)
      
    }
    
    trj_end <- trajectory() %>% 
      set_attribute(
        keys   = c("TimeOfDeath", "TimeToDeath", "TimeToSurgery", "QALYs", "dQALYs"), 
        values = function() fn_calculate_impact(
          CurrentTime = now(.env = sim), 
          Attrs       = get_attribute(.env = sim, keys = c("TimeOfReferral", "TimeOfSurgery"))
        )
      )
    
    trj_main <- trajectory() %>% 

      set_attribute(keys = "TimeOfReferral", values = function() now(.env = sim)) %>% 
      renege_in(t = function() now(.env = sim) + rgompertz(1, d_death_shape, d_death_rate), out = trj_end) %>% 
      
      timeout(task = function() rweibull(1, d_intake_shape, d_intake_scale)) %>% 
      
      seize(resource = "Intake") %>% 
      set_attribute(
        keys   = c('Costs', 'dCosts'), 
        values = function() fn_discount(c_intake, at = now(.env = sim)), 
        mod    = '+'
      ) %>% 
      timeout(task = t_intake) %>% 
      release(resource = "Intake") %>% 
      
      branch(option = function() fn_eligible_intake(), continue = c(F),
             trajectory() %>% 
               set_attribute(keys = "Rejected", values = 1) %>% 
               wait()
      ) %>% 
      
      timeout(task = function() rweibull(1, d_testing_shape, d_testing_scale)) %>% 
      
      seize(resource = "Testing") %>% 
      set_attribute(
        keys   = c('Costs', 'dCosts'), 
        values = function() fn_discount(c_testing, at = now(.env = sim)), 
        mod    = '+'
      ) %>% 
      timeout(task = t_testing) %>% 
      release(resource = "Testing") %>% 
      
      timeout(task = function() rweibull(1, d_consult_shape, d_consult_scale)) %>% 
      
      seize(resource = "Consult") %>% 
      set_attribute(
        keys   = c('Costs', 'dCosts'), 
        values = function() fn_discount(c_consult, at = now(.env = sim)), 
        mod    = '+'
      ) %>% 
      timeout(task = t_consult) %>% 
      release(resource = "Consult") %>%
      
      branch(option = function() fn_eligible_consult(), continue = c(F),
             trajectory() %>% 
               set_attribute(keys = "Rejected", values = 1) %>% 
               wait()
      ) %>% 
      
      timeout(task = function() rweibull(1, d_surgery_shape, d_surgery_scale)) %>% 

      seize(resource = "Surgery") %>% 
      set_attribute(keys = "TimeOfSurgery", values = function() now(.env = sim)) %>% 
      set_attribute(
        keys   = c('Costs', 'dCosts'), 
        values = function() fn_discount(c_surgery, at = now(.env = sim)), 
        mod    = '+'
      ) %>% 
      timeout(task = t_surgery) %>% 
      release(resource = "Surgery") %>%
      
      timeout(task = function() rweibull(1, d_recovery_shape, d_recovery_scale)) %>% 
      
      seize(resource = "Recovery") %>% 
      set_attribute(
        keys   = c('Costs', 'dCosts'), 
        values = function() fn_discount(c_recovery, at = now(.env = sim)), 
        mod    = '+'
      ) %>% 
      timeout(task = t_recovery) %>% 
      release(resource = "Recovery") %>%

      timeout(task = function() rweibull(1, d_followup_shape, d_followup_scale)) %>% 
      
      seize(resource = "FollowUp") %>% 
      set_attribute(keys = "FollowUpCount", values = 1, mod = "+") %>% 
      set_attribute(
        keys   = c('Costs', 'dCosts'), 
        values = function() fn_discount(c_followup, at = now(.env = sim)), 
        mod    = '+'
      ) %>% 
      timeout(task = t_followup) %>% 
      release(resource = "FollowUp") %>%
      rollback(amount = 6, times = n_followup_rounds - 1) %>% 
      wait()
    
    sim <- simmer() %>% 
      add_resource(name = "Intake", capacity = Inf) %>% 
      add_resource(name = "Testing", capacity = Inf) %>% 
      add_resource(name = "Consult", capacity = Inf) %>% 
      add_resource(name = "Surgery", capacity = Inf) %>% 
      add_resource(name = "Recovery", capacity = Inf) %>% 
      add_resource(name = "FollowUp", capacity = Inf) %>% 
      add_generator(name_prefix = "Ind", trajectory = trj_main, mon = 2, 
                    distribution = at(rep(x = 0, times = n_individuals)))
  
    set.seed(seed); sim %>% reset() %>% run();

    df_0 <- fn_summarise(get_mon_attributes(sim))

    p_eligible_intake  <- p_eligible_intake_1
    p_eligible_consult <- p_eligible_consult_1
    c_surgery  <- c_surgery_1
    c_recovery <- c_recovery_1
    
    set.seed(seed); sim %>% reset() %>% run();

    df_1 <- fn_summarise(get_mon_attributes(sim))
    
    dQALYs_0 <- mean(df_0$dQALYs, na.rm = TRUE)
    dQALYs_1 <- mean(df_1$dQALYs, na.rm = TRUE)
    dCosts_0 <- mean(df_0$dCosts, na.rm = TRUE)
    dCosts_1 <- mean(df_1$dCosts, na.rm = TRUE)
    
    out <- c(
      unlist(pars),
      dQALYs_0 = dQALYs_0,
      dQALYs_1 = dQALYs_1,
      dCosts_0 = dCosts_0,
      dCosts_1 = dCosts_0,
      incQALYs = dQALYs_1 - dQALYs_0,
      incCosts = dCosts_1 - dCosts_0,
      ICER     = (dCosts_1 - dCosts_0) / (dQALYs_1 - dQALYs_0),
      NHB_20k  = (dQALYs_1 - dQALYs_0) - (dCosts_1 - dCosts_0) / 20000
    )
    
    rm(list = ls()[ls() != 'out']); gc();
    
    return(out)
    
  })
  
}




## 3. DETERMINISTIC ----

# This section runs the simulation model for the basecase parameters and with the same random
# number seed as used in the 2_health_economics.R script to verify that the fn_runDES() has
# been implemented correctly.

# Loading packages required to run the simulation
library(simmer)           # v4.4.2    
library(flexsurv)         # v2.0      
library(data.table)       # v1.14.0   

# Obtaining the deterministic basecase parameters
pars_deterministic <- fn_pars_deterministic()

# Running the simulation and timing how long it takes to perform the run
# - takes less than 10 sec on a MacBook Pro 2019 2.6 GHz 6-Core Intel Core i7
system.time({
  out_deterministic  <- fn_run_DES(pars = pars_deterministic, n_individuals = 10^4)
})

# Checking whether the results match those obtained through the 2_health_economics.R script:
# - incQALYs = 0.4798688
# - inCosts = 5773.767
# - ICER =  12031.97
out_deterministic[c('incQALYs', 'incCosts', 'ICER')]




## 4. PROBABILISTIC ----

# This section illustrates how the functions can be used to perform a probabilistic analysis.
# The general flow to perform a probabilistic analysis is as follows:
# - Perform "n_runs" probabilistic analysis runs, in each of which:
#   - Sample parameter values (using the fn_pars_probabilistic() function)
#   - Run the simulation with the sampled parameter values (using the fn_run_DES() function)
#   - Record the outcomes from the simulation with the sampled parameters
#   - Continue to the next run
#
# This routine is implemented in a loop below. Rather than using a for-loop, the sapply()
# function is used, which conveniently tries to summarize the outcomes into a matrix, which
# works well for a probabilistic analysis.
#
# Note that the random number seed is set to the run number (i_run) before sampling the
# parameter values to ensure different values are sampled in each run. On the other hand, the
# same random number seed is used to run the simulation model in each run. This is done 
# deliberately to reduce unwarranted variation between runs, as the objective is to 
# quantify the impact that the canges in parameter values have on the outcomes, i.e. the
# parameter uncertainty, not the stochastic uncertainty (which should be simulated out).

# Illustration with a low number of individuals for a low number of runs.
n_runs <- 10

out_probabilistic <- t(sapply(1:n_runs, function(i_run) {
  
  # Set the seed for sampling the parameter values and sample the values using the 
  # fn_pars_probabilistic() function
  set.seed(i_run)
  pars <- fn_pars_probabilistic()
  
  # Run the simulation for the sampled parameter values and return the results
  fn_run_DES(pars = pars, n_individuals = 100, seed = 123)
  
}))


# Given that likely a few hundred of runs are required to obtain stable outcomes from the
# probabilistic analysis, it will take a long time to run e.g. 500 runs of 50,000 individuals
# in sequence. Therefore, the code below illustrates how the analysis can be run in parallel
# using the multiple CPU cores that computers now a days have available. This can be done
# quite easily in R using the base R parallel package.
#
# IMPORTANT! The code below runs the probabilistic analysis in parallel to decrease the run
# time. This means that if your PC has 10 cores and all these are used, 10 instances of the 
# model run at the same time, which may require substantial memory. Therefore, it is 
# important to start with a low number of individuals (and potentially cores) to see how much
# your PC can handle:
# - you can check the utilization of your PC memory using the Task Manager (Windows) or the
#   Activity Monitor (MacOS).
# - results of a probabilistic analysis of 500 runs of 50,000 individuals per run, are
#   saved in the out_probabilistic 500runs 50000individuals.RDS, which can be loaded so that
#   you don't have to run the analysis yourself to analyse the results.
#
# Note that there are different ways that computing cluster can work. For the purpose of
# running a probabilistic analysis the main difference is the memory management. In a "fork" 
# cluster, each thread is run separately but they use the same environment. In a "socket" 
# cluster, each thread has its own environment with all objects required to execute the code,
# which will have to be exported from the main global environment. Therefore, a fork cluster 
# is easier and quicker to set up and more efficient in terms of memory usage, and generally
# preferable. However, fork clusters do not work on Windows machines. Hence, below we show 
# the implementation for both fork and socket clusters.

# Load the package to run the loop in parallel
library(parallel)

n_runs <- 500

# Specify the fork cluster (for Mac users) through the number of cores that are to be used
# - you can use detectCores() to find out how many are available
#detectCores()
cl_fork <- makeForkCluster(nnodes = 6)
  
# Run the analysis using the defined cluster. Not that it is very straightward to make a loop
# specified using an sapply() function (or any other function from the "apply" family) run
# in parallel, as they have a parallel equivalent like parSapply() which only requires the
# cluster to be specified
out_probabilistic <- t(parSapply(cl_fork, 1:n_runs, function(i_run) {
  
  # Set the seed for sampling the parameter values and sample the values using the 
  # fn_pars_probabilistic() function
  set.seed(i_run)
  pars <- fn_pars_probabilistic()
  
  # Run the simulation for the sampled parameter values and return the results
  fn_run_DES(pars = pars, n_individuals = 5*10^4, seed = 123)
  
}))

# It is always smart to immediately save results that took a while to obtain
#saveRDS(object = out_probabilistic, file = 'out_probabilistic.RDS')

# Stop the cluster to free up computer resources
stopCluster(cl_fork)





# Specify the socket cluster (for Windows and Mac users) through the number of cores that are
# to be used. Because this is a socket cluster, we will need to export all the object from
# the global environment that should be available to the threads to execute the code
# - Note that this also includes functions of specific packages, as each thread basically is
#   a new clean instance of R. To prevent exporting all the functions from, for example the
#   simmer package, one by one, we can load the required packages before running the 
#   simulation as demonstrated below. That means we only need to export the objects and 
#   functions that are not included in those packages.
cl_sock <- makePSOCKcluster(names = 6)
clusterExport(cl_sock, c('fn_pars_probabilistic', 'fn_run_DES'))

out_probabilistic <- t(parSapply(cl_sock, 1:n_runs, function(i_run) {
  
  # Load all packages required to run the code below
  library(simmer)
  library(flexsurv)
  library(data.table)
  
  # Set the seed for sampling the parameter values and sample the values using the 
  # fn_pars_probabilistic() function
  set.seed(i_run)
  pars <- fn_pars_probabilistic()
  
  # Run the simulation for the sampled parameter values and return the results
  fn_run_DES(pars = pars, n_individuals = 5*10^4, seed = 123)
  
}))

#saveRDS(object = out_probabilistic, file = 'out_probabilistic.RDS')

stopCluster(cl_sock)


# We can now visualize the result, for example in a incremental cost-effectiveness plane
library(ggplot2)    # v3.3.3

# Turn the matrix with the results of the probabilistic analysis into a data.frame for more
# straightforward analysis and plotting
df_out <- as.data.frame(out_probabilistic)

# Generate the incremental cost-effectiveness plot with the mean outcomes, as well as a 95%
# confidence ellipse
plot_ICE <- ggplot() + 
  geom_point(data = df_out, mapping = aes(x = incQALYs, y = incCosts)) + 
  geom_point(mapping = aes(x = mean(df_out$incQALYs), y = mean(df_out$incCosts)), shape = 16, size = 3, color = 'red') +
  stat_ellipse(data = df_out, mapping = aes(x = incQALYs, y = incCosts), size = 1, color = 'red')

last_plot()



