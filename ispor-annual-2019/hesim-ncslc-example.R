library("hesim")
library("iviNSCLC")
library("data.table")
library("ggplot2")
library("scales")
theme_set(theme_minimal())

# Some helper functions
## Kind of a bad function because it depends on the global strategies variable
## defined below
strategy_factor <- function(x, rev = FALSE){ 
  if (rev == FALSE){
    x[, strategy_name := factor(strategy_id, 
                                levels = strategies$strategy_id, 
                                labels = strategies$strategy_name)]
  } else{
    x[, strategy_name := factor(strategy_id, 
                                levels = rev(strategies$strategy_id), 
                                labels = rev(strategies$strategy_name))]
  }
}

## Bad function for the same reason
state_factor <- function(x, rev = FALSE){ 
  if (rev == FALSE){
    x[, state_name := factor(state_id, 
                             levels = 1:nrow(tmat), 
                             labels = colnames(tmat))]
  } else{
    x[, state_name := factor(state_id, 
                             levels = rev(1:nrow(tmat)), 
                             labels = colnames(tmat))]
  }  
}

# 1. The decision problem and model structure(s)
## Treatment sequences
strategies <- data.table(strategy_id = 1:2,
                         strategy_name = c("gefitinib", "erlotinib"))

## Target population 
patients <- data.table(patient_id = 1:1000) # Simulating 1000 identical patients

## Model structure (health states)
states <- data.table(state_id = 1:2,
                     state_name = c("Stable", "Progression")) # Non-death health states

##  Model structure (health state transitions)
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- c(states$state_name, "Dead")
transitions <- create_trans_dt(tmat)

### hesim_data
hesim_dat <- hesim_data(patients = patients,
                        strategies = strategies,
                        states = states,
                        transitions = transitions)

# 2. The parameter estimates
## Multi-state model
## We will use the Weibull NMA results. Let's take 1,000 samples from 
## the posterior distribution, which we will use for the PSA
n_samples <- 1000
params_mstate_nma_wei <- iviNSCLC::params_mstate_nma$weibull
params_mstate_nma_wei$n_samples <- n_samples
for (p in c("a0", "a1")){
  params_mstate_nma_wei$coefs[[p]] <- params_mstate_nma_wei$coefs[[p]][1:n_samples, ] 
}

## Utility (based on iviNSCLC::params_utility$state_utility)
utility_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = c(0.7540, 0.6532),
                                       sd = c(0, 0.02223000)),
                            dist = "norm",
                            hesim_data = hesim_dat)


## Costs 
### Inpatient costs (based on iviNSCLC::params_costs_inpt)
hospcost_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                        mean = c(1909, 5805),
                                        sd = c(469, 606)),
                             dist = "norm",
                             hesim_data = hesim_dat)

### Drug costs
discount_lower <- .2
discount_upper <- .3
#### Note: use iviNSCLC::annualized_tx_costs() to compute annualized costs

#### Costs in stable disease based on 1L drug
drugcost_stable <- data.table(strategy_id = strategies$strategy_id,
                              state_id = 1,
                              wac = c(93890.7998, 102895.7950))

#### Costs in progression based on 2L drug (use osimertinib)
drugcost_progression <- data.table(strategy_id = strategies$strategy_id,
                              state_id = 2,
                              wac = c(177955.5223, 177955.5223))

#### Combine costs
##### Note that chemo is only administered for 6 21 day cycles (126 days) so we need
##### to add time intervals
drugcost_tbl <- rbind(drugcost_stable,
                       drugcost_progression)
drugcost_tbl[, max := wac * (1 - discount_lower)]
drugcost_tbl[, min := wac * (1 - discount_upper)]
drugcost_tbl <- stateval_tbl(drugcost_tbl,
                             dist = "unif",
                             hesim_data = hesim_dat)

# 3. Simulation
## 3.1 Constructing the model
### Health state transitions
#### Input data 
transmod_data <- expand(hesim_dat, 
                        by = c("strategies", "patients", "transitions"))
#### As noted above, since we did not fit with flexsurv, we need to do a bit of manual work
#### to set up the covariates so that they are consistent with 
#### iviNSCLC::params_mstate_nma (this is what create_transmod_data() does under the
#### hood in the iviNSCLC package). We are creating transition-strategy interaction
#### terms.
for (p in c("a0", "a1")){
  # gefitinib is the reference arm (i.e., predicts the absolute effect)
  transmod_data[, paste0("gef_s1p1_", p) := ifelse(transition_id == 1, 1, 0)]
  transmod_data[, paste0("gef_s1d_", p) := ifelse( transition_id == 2, 1, 0)]
  transmod_data[, paste0("gef_p1d_", p) := ifelse(transition_id == 3, 1, 0)]
  
  # erlotinib
  transmod_data[, paste0("d_erl_s1p1_", p) := ifelse(strategy_id == 2 & transition_id == 1, 1, 0)]
  transmod_data[, paste0("d_erl_s1d_", p) := ifelse(strategy_id == 2 & transition_id == 2, 1, 0)]
  transmod_data[, paste0("d_erl_p1d_", p) := ifelse(strategy_id == 2 & transition_id == 3, 1, 0)]
}
head(transmod_data[, .(strategy_id, patient_id, transition_id, 
                        gef_s1p1_a0, gef_s1d_a0, gef_p1d_a0,
                        d_erl_s1p1_a0)])

##  We will need to modify the coefficients in params_mstate_nma_wei to reflect
## the variables in transmod_data (this is automated in the iviNSCLC R package).
## Note that if models fit with 
## flexsurv, then it's a bit easier because hesim automatically creates
## the "input matrix" from an input dataset and the fitted multi-state model. 
## It would of course be a nice addition to hesim to allow users to fit NMAs 
## directly from R (rather than JAGS) as it would allow for a closer
## integration with hesim and less manual work.
cols <- grep("^gef|^d_erl",colnames(transmod_data), value = TRUE)
for (p in c("a0", "a1")){
  cols_p <- grep(p, cols, value = TRUE)
  params_mstate_nma_wei$coefs[[p]] <- params_mstate_nma_wei$coefs[[p]][, cols_p]
}
head(params_mstate_nma_wei$coefs$a0[, c("gef_s1p1_a0", "gef_s1d_a0",
                                         "gef_p1d_a0", "d_erl_s1p1_a0")])

#### Model 
transmod <- create_IndivCtstmTrans(object = params_mstate_nma_wei, 
                                   input_data = transmod_data,
                                   trans_mat = tmat, 
                                   clock = "forward") # The NMA is a "clock-forward" multi-state model



### Utility model
utilitymod <- create_StateVals(utility_tbl, n = n_samples) 

### Cost models
hospcostmod <- create_StateVals(hospcost_tbl, n = n_samples)
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
costmods <- list(Hospital = hospcostmod,
                 Drug = drugcostmod)

### Economic model
econmod <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = costmods)

## 3.2 Simulating outcomes
### Disease progression
max_yrs <- 20
max_age <- 100
econmod$sim_disease(max_t = max_yrs * 12, max_age = max_age * 12)
econmod$disprog_[, ':=' (time_start = time_start/12, time_stop = time_stop/12)]
econmod$sim_stateprobs(t = seq(0, 20 , 1/26)) # Biweekly probabilities

#### Disease progression
head(econmod$disprog_)

#### Plot of state probabilities
stprobs <- econmod$stateprobs_[, .(prob_mean = mean(prob),
                                   prob_lower = quantile(prob, .025),
                                  prob_upper = quantile(prob, .975)),
                                by = c("strategy_id", "state_id", "t")]
strategy_factor(stprobs)
state_factor(stprobs)
ggplot(stprobs[t < 15], aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() + facet_wrap(~state_name) + 
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "") +
  theme(legend.position = "bottom")

### QALYs
econmod$sim_qalys(dr = c(0, .03))
mean_qalys <- econmod$qalys_[, .(mean_lys = mean(lys),
                             mean_qalys = mean(qalys)),
                             by = c("strategy_id", "state_id", "dr")]
state_factor(mean_qalys, rev = TRUE)
strategy_factor(mean_qalys)

#### Life-years plot
ggplot(mean_qalys[dr == 0], 
       aes(x = strategy_name, y = mean_lys,
          fill = state_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Mean life-years")

#### (Discounted) QALYs plot
ggplot(mean_qalys[dr == .03], 
       aes(x = strategy_name, y = mean_qalys,
          fill = state_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Mean QALYs")

### Costs
econmod$sim_costs(dr = .03)
mean_costs <- econmod$costs_[, .(mean_costs = mean(costs)),
                             by = c("strategy_id", "category", "dr")]
strategy_factor(mean_costs)
ggplot(mean_costs, 
       aes(x = strategy_name, y = mean_costs, fill = factor(category))) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("") + ylab("Costs") + 
  scale_y_continuous(label = dollar_format())

# 4. Cost-effectiveness analysis
ce_sim <- econmod$summarize()
icea_out <- icea(ce_sim, dr_qalys = .03, dr_costs = .03)
icea_pw_out <- icea_pw(ce_sim, comparator = 1, dr_qalys = .03, dr_costs = .03)

## Cost-effectiveness acceptability frontier
strategy_factor(icea_out$mce)
ggplot(icea_out$mce[best == 1], aes(x = k, y = prob, col = strategy_name)) +
  geom_line() + xlab("Willingness to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, 200000, 25000), label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")
