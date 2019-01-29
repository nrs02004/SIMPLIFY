## This function generates a function which runs a single, simple randomized two-arm trial...
## to evaluate the effect of withdrawal of a medication on FEV

## n_con: Number of patients on control arm (receiving the medication)
## n_wd: Number of patients on the withdrawal arm
## fev_base: The population mean fev for patients on n_con (this number doesn't actually matter because we are considering a difference)
## fev_delta: The population mean difference between withdrawal and control arm
## fev_sd: The standard deviation of the fev measures
## dropout_avg_t: The time to dropout (Simulation uses an exponential distribution with this average
## level: The level of the test (used for constructing the CI)

gen_one_trial_fev_closure <- function(n_con,
                              n_wd,
                              fev_base,
                              fev_delta = 0,
                              fev_sd,
                              level = 0.05){
  ## This function runs a single trial with the parameters from above
  ## It outputs a confidence interval for the mean difference of FEV
  run_trial <- function(){
    out_control_base <- rnorm(n_con, mean = fev_base, sd = fev_sd)
    out_wd_base <- rnorm(n_wd, mean = fev_base, sd = fev_sd)
  
    out_control_followup <- rnorm(n_wd, mean = fev_base, sd = fev_sd)
    out_wd_followup <- rnorm(n_wd, mean = fev_base + fev_delta, sd = fev_sd)
  
    out_control_del <- out_control_followup - out_control_base
    out_wd_del <- out_wd_followup - out_wd_base
  
    fit <- t.test(out_wd_del,out_control_del, conf.level = 1-level)
    interval.est <- fit$conf.int[1:2]
    return(interval.est)
  }
  return(run_trial)
}


## This function generates a function which single, simple randomized two-arm trial...
## to evaluate the effect of withdrawal of a medication on Pulmonary Exacerbations

## n_con: Number of patients on control arm (receiving the medication)
## n_wd: Number of patients on the withdrawal arm
## prop_PEX_in_6_m: The probability on the control arm of experiencing a pulmonary exacerbation within 6 months (sims use exponential dist)
## prop_PEX_in_6_m_delta: The change in probability from control to withdrawal of PEX within 6 months
## dropout_in_6_m: The probability of dropout within 6 months (sims use exponential dist)
## study_end: study length in years
## level: The level of the test (used for constructing the CI)
gen_one_trial_pex_closure <- function(n_con,
                                      n_wd,
                                      prop_PEX_in_6_m,
                                      prop_PEX_in_6_m_delta,
                                      dropout_in_6_m,
                                      study_end,
                                      level = 0.05){
  
  ## Converting the proportion PEX within 6 months into an exponential rate
  pex_rate_con <- -log(1-prop_PEX_in_6_m)/0.5
  pex_rate_wd <- -log(1-(prop_PEX_in_6_m + prop_PEX_in_6_m_delta))/0.5
  dropout_rate <- -log(1-dropout_in_6_m)/0.5 
  
  ## This function runs a single trial with the parameters from above
  ## It outputs a confidence interval for the log-hazard ratio of the tx effect (on PEX)
  run_trial <- function(){
    pex_time_control <- rexp(n_con, rate = pex_rate_con)
    pex_time_wd <- rexp(n_wd, rate = pex_rate_wd)
    
    dropout_time_control <- pmin(rexp(n_con, rate = dropout_rate), study_end)
    dropout_time_wd <- pmin(rexp(n_wd, rate = dropout_rate), study_end)
    
    tx_group <- c(rep("control", n_con), rep("wd", n_wd))
    outcome <- Surv(time = c(pex_time_control, pex_time_wd),
                    event = c(dropout_time_control > pex_time_control, 
                               dropout_time_wd > pex_time_wd))
    fit <- coxph(outcome~tx_group)
    interval.est <- summary(fit, conf.int = 1-level)$conf.int[c(3,4)]
    return(interval.est)
  }
  return(run_trial)
}


## A function for evaluating if the margin falls below the confidence interval (to evaluate trial success)
eval.trial <- function(interval, margin){
  return(margin < interval[1])
}

## A function to calculate our power for 
eval.design <- function(run_trial, margin, n_trial){
  trials.tests <- replicate(n_trial, eval.trial(run_trial(), margin))
  return(mean(trials.tests))
}
  