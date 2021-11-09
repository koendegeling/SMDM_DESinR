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
# This 1_basic_structure.R script implements the overall model structure that reflect the
# pathway that the individuals follow. The subsequent 2_health_economics.R script illustrates
# how this script can be extended to include health economic outcomes for two different 
# strategies. Finally, the 3_probabilistic_analysis.R script demonstrates how the simulation
# model can be analyzed through a probabilistic analysis.
#
# This script contains the following sections:
# 1. INITIALIZATION   sets up the environment for the analysis
# 2. PARAMETERS       defines the parameters used to define the trajectory
# 3. FUNCTIONS        defines functions that are used to implement and run the simulation
# 4. TRAJECTORY       defines the trajectory (i.e., model structure)
# 5. SIMULATION       defines and runs the simulations, including analysis of the outcomes
#
# This script was written using R Studio version 1.4.1106 and R version 4.0.3. The versions
# of the packages used are specified in the INITIALIZATION section.
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
#rm(list = ls()); cat("\14");

# Uncomment to install the used packages
#install.packages(pkgs = c("simmer", "simmer.plot", "flexsurv", "data.table"))

# Loading the packages
library(simmer)           # v4.4.2    implementing and running discrete event simulations
library(simmer.plot)      # v0.1.16   plotting discrete event simulation objects and outputs
library(flexsurv)         # v2.0      Gompertz distribution functions
library(data.table)       # v1.14.0   efficient wrangling of simulation output




## 2. PARAMETERS ----

# This section defines the parameters that are used to define the discrete event simulation.
# These are the parameters that define and are used in the trajectory. The prefix of the 
# parameter name indicates what type of parameter it is:
#   d_    distribution parameter
#   n_    count
#   p_    probability
#   t_    time duration

# Parameters of the Gompertz distribution for background mortality
d_death_shape <- 0.1
d_death_rate  <- 0.005

# Parameters for the intake process: parameters of the Weibull distribution that is used to
# model the time to the intake (d_), the duration of the intake being 1 day (t_), and the 
# probability that an individual is eligible to continue in the process based on the intake
# (p_).
d_intake_shape <- 2.0
d_intake_scale <- 0.3
t_intake <- 1/365.25
p_eligible_intake <- 0.7

# Parameters for the testing process: parameters of the Weibull distribution that is used to
# model the from the intake to the testing procedures (d_), and the duration of the testing
# procedures being 1 day (t_).
d_testing_shape <- 2.5
d_testing_scale <- 0.2
t_testing <- 1/365.25

# Parameters for the final consult before surgery: parameters of the Weibull distribution
# that is used to model the time from the testing procedures to the consult (d_), the
# duration of the consult being 1 day (t_), and the probability that an individual is 
# eligible to continue to the surgery (p_).
d_consult_shape <- 3.0
d_consult_scale <- 0.1
t_consult <- 1/365.25
p_eligible_consult <- 0.9

# Parameters for the surgery: parameters of the Weibull distribution that is used to model
# the time from the final consult to the surgery (d_), and the duration of the surgery being
# 1 week including the hospital stay (t_).
d_surgery_shape <- 1.5
d_surgery_scale <- 0.5
t_surgery <- 1/52

# Parameters for the recovery process: parameters of the Weibull distribution that is used to
# model the time from the surgery to the start of the recovery process (d_), and the duration
# of the recovery process being 2 weeks (t_).
d_recovery_shape <- 2.0
d_recovery_scale <- 2/52
t_recovery <- 3/12

# Parameters for the follow-up process: parameters of the Weibull distribution that is used
# to model the time from the recovery process to the first follow-up check, as well as to
# model the time between follow-up checks (d_), and the duration of a follow-up check being 1
# day (t_).
d_followup_shape  <- 3.0
d_followup_scale  <- 1.0
t_followup <- 1/365.25
n_followup_rounds <- 5




## 3. FUNCTIONS ----

# This section defines several supporting functions that are used to implement and run the
# discrete event simulation. Each function is identified by the prefix 'fn_'. Further 
# information on the functions is provided within them.

fn_eligible_intake <- function() {
  
  # Function to determine whether the individual is eligible to continue for testing at the
  # time of the intake, with the following output value that is used to determine the sub-
  # trajectory in the corresponding branch:
  #   0) continue to testing (0 = skip the branch)
  #   1) not eligible for testing (1 = enter the first sub-trajectory in the branch)
  
  out <- if(runif(1) < p_eligible_intake) {0} else {1}
  
  return(out)
  
}

fn_eligible_consult <- function() {
  
  # Function to determine whether the individual is eligible to continue for surgery at the
  # time of the final consult, with the following output value that is used to determine the
  # sub-trajectory in the corresponding branch:
  #   0) continue to surgery (0 = skip the branch)
  #   1) not eligible for surgery (1 = enter the first sub-trajectory in the branch)
  
  out <- if(runif(1) < p_eligible_consult) {0} else {1}
  
  return(out)
  
}

fn_discount <- function(amount, at = NULL, start = NULL, duration = NULL, rate = 0.03, timefactor = 1) {
  
  # This function returns both the undiscounted and discounted net present values based on an
  # amount (amount), discount rate (rate), and factor (timefactor) to adjust for a potential
  # mismatch between the time unit between the discount rate and time period (e.g., the  
  # default value of 1 applies if the simulation is in years). If the "at" argument is  
  # specified, the amount is discounted at that point in time, whereas the amount is 
  # continuously discounted over the corresponding period of time if the "start" and 
  # "duration" arguments are defined.
  
  # If the "start" and "duration" arguments are specified, use the formula for continuous
  # discounting to discount the the amount over a period of time.
  if(is.null(at) & !is.null(start) & !is.null(duration)) {
    
    undisc <- amount * duration;
    disc   <- (1/timefactor) * (amount / (1 + rate) ^ (start*timefactor) / log(1 + rate) * (1 - 1 / (1 + rate) ^ (duration*timefactor)));
    
    # If only the "at" argument is defined, discount the amount at that specific point in time
  } else if(!is.null(at) & is.null(start) & is.null(duration)) {
    
    undisc <- amount
    disc   <- amount / ((1 + rate) ^ (at*timefactor))
    
  } else {
    
    stop('Time argument(s) not specified correctly for fn_discount() function')
    
  }
  
  out <- c(Undiscounted = undisc, Discounted = disc)
  
  return(out)
  
}

fn_calculate_impact <- function(CurrentTime, Attrs) {
  
  # This function records the time of death, time to death from referral (i.e., overall 
  # survival), and the time between referral and surgery.
  #
  # The function has the following input arguments that are to be provided:
  # - CurrentTime         the simulation time at the time of the function call
  # - Attrs               vector with the following attribute values:
  #   [1] TimeOfReferral  the time at which the individual was referred
  #   [2] TimeOfSurgery   the time at which surgery was performed (can be NA)
  # 
  # The function needs to return the following output in line with the set_attribute call in 
  # the trajectory:
  # - c("TimeOfDeath", "TimeToDeath", "TimeToSurgery")
  
  # Extracting the attribute values
  TimeOfReferral  <- Attrs[1]
  TimeOfSurgery   <- Attrs[2]
  
  # Determining the time-to-events
  TimeOfDeath   <- CurrentTime
  TimeToDeath   <- TimeOfDeath - TimeOfReferral
  TimeToSurgery <- TimeOfSurgery - TimeOfReferral
  
  # Required output: c("TimeOfDeath", "TimeToDeath", "TimeToSurgery")
  out <- c(TimeOfDeath, TimeToDeath, TimeToSurgery)
  
  return(out)
  
}

fn_summarise <- fn_summarize <- function(sim_out, keys = NULL) {
  
  # This function summarizes monitored attribute values extracted from the simmer()
  # environment into a data.frame using the get_mon_attributes() function, which is provided
  # to function through the 'sim_out' argument, to their last recorded value per individual.
  # If no specific attributes are defined through the 'keys' argument, all attributes/keys
  # are summarized. Functions from the data.table package are used due to the potentially 
  # large size of 'sim_out'.
  
  if(is.null(keys)) keys <- unique(sim_out$key);
  
  df <- as.data.table(sim_out)[key %in% keys];
  setorder(df, name, time);
  df <- df[, .(value = value[.N]), by = list(name, key)];
  df <- dcast(df, name~key, value.var = "value");
  setcolorder(df, c("name", keys));
  
  return(df)
  
}




## 4. TRAJECTORY ----

# This section implements the trajectory/model structure through which individuals will be
# simulated. In summary, the pathway for the individuals is as follows:
# - Referral    The model starts with the referral of an individual.
# - Intake      The eligibility of individuals is initially assessed during an intake, where
#               the vast majority of individuals are found to be eligible.
# - Testing     Those considered eligible during the intake receive further tests to 
#               determine their eligibility for surgery.
# - Consult     The final consult determines the final eligibility for surgery, where
#               a small proportion are found not to be eligible based on the test results
# - Surgery     During the surgery, the knee replacement is placed.
# - Recovery    After surgery, individuals go through a period of recovery.
# - Follow-up   After completing the recovery, individuals enter a follow-up process in which
#               they are checked at a regular interval.
#
# Whilst in this process, individuals are at risk of death from background mortality.
# Remember that individuals have a relatively low quality of life, which improves if they
# receive their knee replacement.
#
# Two trajectories are defined to implement the above pathway: trj_end and trj_main. The 
# trj_end trajectory is administrative in that its sole function is to record the final 
# outcomes for an individual once they die. The trj_main trajectory represents the actual
# pathway as describe above. Because the trj_end trajectory is referred to from within the
# trj_main trajectory, it has to be defined first.

# The trj_end determines the outcomes by calling the fn_calculate_impact() function, which 
# requires the time of the simulation at that time and some individual attributes as inputs.
trj_end <- trajectory() %>% 
  set_attribute(
    keys   = c("TimeOfDeath", "TimeToDeath", "TimeToSurgery"), 
    values = function() fn_calculate_impact(
      CurrentTime = now(.env = sim), 
      Attrs       = get_attribute(.env = sim, keys = c("TimeOfReferral", "TimeOfSurgery"))
    )
  )

# The main trajectory implements the general pathway
trj_main <- trajectory() %>% 
  
  # Referral:
  # - recording the time of referral 
  # - setting the moment at which the individual should be transferred to the trj_end
  #   trajectory, regardless of where in the pathway they are at that moment, based on the
  #   background mortality, using the renege_in() function.
  set_attribute(keys = "TimeOfReferral", values = function() now(.env = sim)) %>% 
  renege_in(t = function() now(.env = sim) + rgompertz(1, d_death_shape, d_death_rate), out = trj_end) %>% 
  
  # Time to next event: Intake
  timeout(task = function() rweibull(1, d_intake_shape, d_intake_scale)) %>% 
  
  # Intake
  seize(resource = "Intake") %>% 
  timeout(task = t_intake) %>% 
  release(resource = "Intake") %>% 
  
  # Determine the outcome of the intake, and hence what happens next, using the 
  # fn_eligible_intake() function:
  #   0) continue to testing (i.e., skip the branch)
  #   1) not eligible (i.e., wait until the individual is transferred to trj_end)
  branch(option = function() fn_eligible_intake(), continue = c(F),
         
         # 1) not eligible
         trajectory() %>% 
           set_attribute(keys = "Rejected", values = 1) %>% 
           wait()
         
  ) %>% 
  
  # Time to next event: Testing
  timeout(task = function() rweibull(1, d_testing_shape, d_testing_scale)) %>% 
  
  # Testing
  seize(resource = "Testing") %>% 
  timeout(task = t_testing) %>% 
  release(resource = "Testing") %>% 
  
  # Time to next event: Consult
  timeout(task = function() rweibull(1, d_consult_shape, d_consult_scale)) %>% 
  
  # Final consult
  seize(resource = "Consult") %>% 
  timeout(task = t_consult) %>% 
  release(resource = "Consult") %>%
  
  # Determine the outcome of the final consult, and hence what happens next, using the 
  # fn_eligible_consult() function:
  #   0) continue to testing (i.e., skip the branch)
  #   1) not eligible (i.e., wait until the individual is transferred to trj_end)
  branch(option = function() fn_eligible_consult(), continue = c(F),
         
         # 1) not eligible
         trajectory() %>% 
           set_attribute(keys = "Rejected", values = 1) %>% 
           wait()
         
  ) %>% 
  
  # Time to next event: Surgery
  timeout(task = function() rweibull(1, d_surgery_shape, d_surgery_scale)) %>% 
  
  # Surgery:
  # - also record the time of surgery to calculate the time between referral and surgery, as
  #   well as to track whether the individual has received surgery
  seize(resource = "Surgery") %>% 
  set_attribute(keys = "TimeOfSurgery", values = function() now(.env = sim)) %>% 
  timeout(task = t_surgery) %>% 
  release(resource = "Surgery") %>%
  
  # Time to next event: Recovery
  timeout(task = function() rweibull(1, d_recovery_shape, d_recovery_scale)) %>% 
  
  # Recovery
  seize(resource = "Recovery") %>% 
  timeout(task = t_recovery) %>% 
  release(resource = "Recovery") %>%
  
  # Follow up:
  # - track the number of follow-up visits
  # - after the follow-up visit, wait until the next follow-up visit and roll back
  # - this continues until the individual is transferred to the trj_end trajectory
  
  # Time to next event: Follow up
  timeout(task = function() rweibull(1, d_followup_shape, d_followup_scale)) %>% 
  
  seize(resource = "FollowUp") %>% 
  set_attribute(keys = "FollowUpCount", values = 1, mod = "+") %>% 
  timeout(task = t_followup) %>% 
  release(resource = "FollowUp") %>%
  rollback(amount = 5, times = n_followup_rounds - 1) %>% 
  
  # Wait until moved to traj_end
  wait()

# This was the end of the trj_main trajectory


# Uncomment to plot the trajectory 
#plot(trj_main)




## 5. SIMULATION ----

# This sections defines and runs the simulation. Note that the name 'sim' for the simulation
# environment has been hard-coded into the trajectory, so this name has to be used. Also note
# that we are not considering capacity constraints for now, so the capacity for each resource
# is set to be infinite (i.e., Inf). 
#
# Given that we are considering a health economic evaluation without resource constraints,
# all individuals can simply enter the simulation at time zero, for which the at() function
# is used to define the "distribution" argument of the add_generator() function.

n_individuals <- 10^4

sim <- simmer() %>% 
  add_resource(name = "Intake", capacity = Inf) %>% 
  add_resource(name = "Testing", capacity = Inf) %>% 
  add_resource(name = "Consult", capacity = Inf) %>% 
  add_resource(name = "Surgery", capacity = Inf) %>% 
  add_resource(name = "Recovery", capacity = Inf) %>% 
  add_resource(name = "FollowUp", capacity = Inf) %>% 
  add_generator(name_prefix = "Ind", trajectory = trj_main, mon = 2, 
                distribution = at(rep(x = 0, times = n_individuals)))

# Running the simulation until all events have occurred
set.seed(123); sim %>% reset() %>% run();

# Extracting the simulation output
df_attributes <- get_mon_attributes(sim)
df_arrivals   <- get_mon_arrivals(sim)
df_resources  <- get_mon_resources(sim)  

# Summarizing the recorded attributes using the custom function
df_out <- fn_summarise(df_attributes)


# Inspecting the monitored arrival times
head(df_arrivals)

# Plotting the resource use
head(df_resources)
plot(df_resources, names = c('Intake', 'Consult', 'Surgery'), steps = TRUE)

# Assessing monitored attributes
head(df_out)
df_out[ , .(
  'Probabilty of Surgery' = mean(!is.na(TimeOfSurgery)),
  'Mean Time to Surgery'  = mean(TimeToSurgery, na.rm =  TRUE),
  'Mean Overall Survival' = mean(TimeToDeath)
)]



