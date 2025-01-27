source("source/functions.R")
sourceCpp("source/functions.cpp")

# `data' read in as list with 3 primary components, 3 secondary components, and 8 tertiary components (datasets)
# `data$train`: contains info on training data
#               data$train$y contains realized show (0), and no-show(1) outcomes
#               data$train$pred contains predicted probabilities of no-show
#               data$train$lead contains lead times for each patient (number of days before appointment it was scheduled)
#               each of data$train$y, data$train$pred, and data$train$lead are lists of length 8, containing values for 8 different datasets
#               indexed by data$departments
# `data$test`: contains info on test data; structured the same as data$train
# `data$departments`: contains vector of the departments that index secondary list components for data$train and data$test


# Calculate AUC on Validation Set for Each Department
bench_auc(data)

#define parameters
tau <- .3 #decision threshold
slot.prob <- .9 #probability a spot is filled 
ob.prob <- .9 #Probability a high risk patient is overbooked
S <- 20 #Number of appointment Slots
mu <- logN_params(22, 6)$mu #Mean service time per patient; can be a 2xS matrix for individualized distributions; yields log-normal mean
sd <- logN_params(22, 6)$sd #SD service time per patient; can be a 2xS matrix for individualized distributions; yields log-normal SD
block.time <- 30 #Length of each appointment slot/block
otmin.thresh = 15
wait.thresh <- 15
idle.thresh <- 10

#Define which dataset to use
K <- 6

#### Generate PPV, NPV, positive prediction rate, and overall no-show rate from training set
# y: vector of true outcomes (0=show, 1 = no show) from training set
# pred: vector of predicted probabilities from training set
# tau: decision threshold such that `pred` > tau -> no show
metrics <- gen_metrics(y=data$train$y[[K]],pred=data$train$pred[[K]], tau=tau)


#### Generate probability matrix of no show, single show, and double show by appointment slot (this is NOT cumulative; each slot is indepndent in this matrix)
# preds: NULL or Sx3 matrix containing predicted probabilities for a simulated patient day, where each row represents a booking slot 
#        and each column represents probability of no shows, single show or double show for each slot. NA's represent no booking occurred
#        If NULL, then day simulated purely on probability based on object `metrics`, slot.prob, and ob.prob
# complete: if TRUE, then object `preds` assumed to be scheduled day without any further bookings. 
#           if FALSE, then missing slots in object `preds` can continue to be booked with some probability based on metrics, slot.prob, and ob.prob
#           ignored if `preds`= NULL
# S: number of slots in a scheduling  period
# slot.prob: probability an unscheduled slot will be scheduled by the first patient
# ob.prob: probability a slot deemed high-risk will be overbooked
# metrics: object from function `gen_metrics` containing PPV, NPV, PPR, and overall no-show rate
probs <- gen_probs(preds=NULL, complete=F, S=S, slot.prob=slot.prob, ob.prob=ob.prob, metrics)

#### Simulate Distribution
# N: number of replications to simulate distribution
# S: number of appointment slots
# probs: object from gen_probs(); probability matrix for each appointment slot show rates
# mu: log-normal mean service time 
# sd: log-normal sd service time
# block.time: length of service blocks/slots
# print.time: If `TRUE`, computation time is printed 
sim_data <- est_dist(N=10000, S=S, probs=probs, mu=mu, sd=sd, block.time=block.time, print.time=T)


#### Get Estimated Performance Metrics Related to `Costs` from simulated data
# sim_data: is object from function est_dist()
# discrete: if "TRUE", yields costs estimated based only on discrete appointment blocks; ignores distribution of appointment times
# wait.thresh: maximum allowable waiting time after appointment scheduled start time until counted towards costs
# idle.thresh: Minimum allowable observation time needed for appointment slot not to be considered `idle`
# block.time: appointment block sizes in terms of minutes
perf_costs <- est_performance(sim_data, discrete=F, otmin.thresh = 15, wait.thresh = wait.thresh, idle.thresh = idle.thresh, block.time=block.time)

#### optimize a scheduled day
# c1, c2, c3 are costs defined above
# y: training data realized outcomes
# pred: training data predicted outcomes
# S: Number of appointment slots
# slot.prob: probaiblity of a slot being filled
# ob.prob: probability of overbooked slot given first patient high risk
# block.time: length of an appointment slot
# discrete: T or F; whether or not to optimze based on on discrete probabilities of incorporate time 
# nsim: number of simulations for simulated time based data
# grid.start: vector of length 2 giving min and max values to start grid search for tau (defaul is c(0.01, 0.99))
# seed: set random seed if desired

#Define Costs 
c1 <- .2  # cost of each excess waiting patient 
c2 <- .75   # cost of Idle Period
c3 <- 1.5 # cost per 30 min of OT 
costs <- c(c1, c2, c3)


wait.thresh=15
otmin.thresh=15
otprob_upper=.15
excess_wait_upper=5
idle.thresh =10

# get optimized threshold 
tau <- optimize(optimizer="thresh", costs,
                      otprob_upper=otprob_upper, otmin.thresh=otmin.thresh, excess_wait_upper = excess_wait_upper, wait.thresh = wait.thresh, idle.thresh=idle.thresh,
                      y=data$train$y[[K]], pred=data$train$pred[[K]], predicted.probs = NULL,
                      S=S, slot.prob=slot.prob, ob.prob=ob.prob,
                      mu=mu, sd=sd,
                      block.time=block.time, 
                      nsim=20000)
tau=tau[1]
#### Simulates a scheduled day
# optimizer: "thresh" or "costs"... Optimizes day based on pre-defined probability thresholds or costs
# otprob_upper: Upper bound probability for probability of overtime beyond 'otmin_upper" minutes
# otmin_upper: Upper bound for acceptable number of minutes of overtime
# excess_wait_upper: mean number of patients with excess waits deemed acceptable based on wait.thresh
# y: show (0) and no-show (1) data from which to simulate (test data)
# pred: no-show probabilities from which to simulate (test data)
# lead: lead time from which to simulate (test data)
# adaptive: T or F. If T, tau is re-optimized after sequential patient bookings based on frequency defined by 'update.freq'
# update.freq: how often (in terms of # of days) to re-optimize threshold
# begin.ob: The minimum number of slots that must be filled before patients will begin overbooking
# tau: threshold for prediction (only used if 'adaptive'=F)
# slot.prob: probability a slot will be initially filled with first patient
# ob.prob: probability a high-risk patient will be booked against 
# mu: mean service time
# sd: standard deviation of service time
# block.time: length of each service period
# y.train: training data outcomes (0 = show, 1= no-show), used for optimizing threshold
# pred.train: training data predicted probability for no-show
# costs: vector of costs for optimization. See `est_performance` function for more
# wait.thresh = max wait time allowable before excess wait times are counted
# idle.thresh = max idle time allowed in an appointment slot before slot counted as `idle slot`
# nsim= number of replications to simulate distributions of service times by slot

## Note: set adaptive=F and tau>0.999 to simulate no overbooking 

#Optimized once beforehand
results <- NULL
results <- vector(mode="list", length=100)
for(i in 1:10000){
  results[[i]] <- simulate(optimizer="thresh", costs,
                           otprob_upper=otprob_upper, otmin.thresh=otmin.thresh, excess_wait_upper = excess_wait_upper,
                            wait.thresh =waith.thresh, idle.thresh=idle.thresh,
                           y= data$test$y[[K]], pred= data$test$pred[[K]], lead=data$test$lead[[K]], 
                   adaptive=F, update.freq=7, begin.ob=S/2, tau=tau, S=S,
                   slot.prob=slot.prob, ob.prob=ob.prob, mu=mu, sd=sd, block.time=block.time,
                   y.train=data$train$y[[K]], pred.train=data$train$pred[[K]],
                    nsim=100000)
}

analyze(results,otmin.thresh=otmin.thresh, wait.thresh=wait.thresh, idle.thresh=idle.thresh, block.time=block.time, costs=costs)


#Adaptive Optimization
results2 <- NULL
results2 <- vector(mode="list", length=1000)
for(i in 1:1000){
  results2[[i]] <- simulate(optimizer="thresh", costs,
                           otprob_upper=otprob_upper, otmin.thresh = otmin.thresh, excess_wait_upper = excess_wait_upper, wait.thresh=wait.thresh, 
                           y= data$test$y[[K]], pred= data$test$pred[[K]], lead=data$test$lead[[K]], 
                           adaptive=T, update.freq=7, begin.ob=S/2, tau=tau, S=S,
                           slot.prob=slot.prob, ob.prob=ob.prob, mu=mu, sd=sd, block.time=block.time,
                           y.train=data$train$y[[K]], pred.train=data$train$pred[[K]], 
                           idle.thresh=idle.thresh, nsim=10000)
  print(i)
}
res_adapt = analyze(results2, wait.thresh=wait.thresh, idle.thresh=idle.thresh, block.time=block.time, costs=costs)






