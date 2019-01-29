## This script generates power plots for the simple non-adaptive design with various sample sizes and trial lengths/margins



library(tidyr)
library(ggplot2)
library(dplyr)
library(plyr)
library(survival)

source("sim_functions.R")

## FEV SIMULATIONS

fev_base = 80
fev_delta = 0
fev_sd = 7.5
level = 0.05


set.seed(1)

ns <- c(50, 100, 150, 200, 250, 300, 350)
margins_fev <- c(-3.5, -3)
pow_fev <- array(dim = c(length(ns), length(margins_fev)))

for(i in 1:length(ns)){
  for(j in 1:length(margins_fev)){
  run_trial <- gen_one_trial_fev_closure(ns[i],
                                         ns[i],
                                         fev_base,
                                         fev_delta,
                                         fev_sd,
                                         level)

  pow_fev[i,j] <- eval.design(run_trial, margin = margins_fev[j], n_trial=1e4)
  print(paste(i,j))
  }
}

## Reformatting the data to be plotted
plot_dat <- data.frame(pow_fev)
colnames(plot_dat) <- margins_fev
plot_dat$ns <- ns
plot_dat_tall <- plot_dat %>% gather(key = "margin", value = "power", -ns)

p_fev <- ggplot(plot_dat_tall, aes(x = ns, y = power, color = margin)) + 
  geom_line() + 
  xlab("# patients per arm") +
  labs(title = "Power for FEV non-inferiority evaluation")

ggsave("../../plots/FEV_power.pdf",p_fev, width = 6, height = 4)


## PEX SIMULATIONS

prop_PEX_in_6_ms_delta = 0

set.seed(1)
margins_pex <- c(.70, .60)
prop_PEX_in_6_ms <- c(0.35,0.425,0.5)
followups <- c(0.25,0.5) ## followup times in years
dropout_in_6_m <- 0.01
pow_pex <- array(dim = c(length(ns), length(margins_pex), length(prop_PEX_in_6_ms), length(followups)))

for(i in 1:length(ns)){
  for(j in 1:length(margins_pex)){
    for(k in 1:length(prop_PEX_in_6_ms)){
      for(l in 1:length(followups)){
        run_trial <- gen_one_trial_pex_closure(ns[i],
                                               ns[i],
                                               prop_PEX_in_6_ms[k],
                                               prop_PEX_in_6_ms_delta,
                                               dropout_in_6_m,
                                               study_end = followups[l],
                                               level)
    
        pow_pex[i,j,k,l] <- eval.design(run_trial, margin = margins_pex[j], n_trial=1e3)
        print(paste(i,j,k,l))
      }
    }
  }
}

## Formatting data to be plotted
pex_dat_tall <- data.frame(power = NA, n = NA, margin = NA, PEX_prob = NA, followup = NA)
row_ind <- 1
for(i in 1:length(ns)){
  for(j in 1:length(margins_pex)){
    for(k in 1:length(prop_PEX_in_6_ms)){
      for(l in 1:length(followups)){
        out <- c(pow_pex[i,j,k,l],ns[i], margins_pex[j], prop_PEX_in_6_ms[k], followups[l])
        pex_dat_tall[row_ind,] <- out
        row_ind <- row_ind + 1
      }
    }
  }
}

p_pex <- pex_dat_tall %>% mutate(PEX_prob = mapvalues(as.factor(PEX_prob), c("0.35","0.425","0.5"),
                                           c("PEx rate = 0.35","PEx rate = 0.425","PEx rate = 0.5")),
                        followup = as.factor(ifelse(followup == 0.25,"3 months followup", "6 months followup")),
                        margin = as.factor(margin)) %>%
ggplot(aes(x = n, y = power, color = margin)) + 
  geom_line() +
  facet_grid(PEX_prob~followup) +
  xlab("# patients per arm") +
  labs(title = "Power for PEx non-inferiority evaluation")

ggsave("../../plots/PEx_power.pdf",p_pex, width = 6, height = 4)
