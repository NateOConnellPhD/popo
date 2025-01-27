library(ROCR)
library(Rcpp)

#Read in Data
data <- readRDS("data/predictions_051624.RData")
data$train$lead <- lapply(data$train$lead, function(x) ifelse(x<0, 0, x))
data$test$lead <- lapply(data$test$lead, function(x) ifelse(x<0, 0, x))
data$test$y <- lapply(data$test$y, function(x) as.numeric(as.character(x)))
data$train$y <- lapply(data$train$y, function(x) as.numeric(as.character(x)))


##### Converts specified algorithmic Means and SDs to LogNormal Means and SD for distribution sim
# mu: input algorithmic mean 
# sd: input SD
logN_params = function(mu, sd){
   location = log(mu^2/sqrt(mu^2 + sd^2))
   shape = sqrt(log(1+(sd^2/ mu^2)))
   list(mu= location,
   sd =shape)
 }


####### Simulate Service Time Distribution ##########
# 'N'     = number of replications
# 'S'     = number of appointment slots
# 'probs' = vector of length 3 for probabilities of No Show, Single Show, and Double Show or
#           an S x 3 matrix representing No show, single show and double show prob by appointment slot
# 'mu'    = point estiamte of mean service time at the clinic, or 2 X S matrix of service time by each expected patient
#           by appointment slot (S; Columns)
# 'sd'    = point estimate of standard deviation of service time at the clinic, or 2 X S matrix of SD estimates for each patient by slot
# 'block_time' = length of each scheduled visit period
est_dist <- function(N=5000, S=20, probs=NULL, mu, sd, block.time = 30, discrete=F, seed=NULL, print.time = F){
  if(is.null(seed)==F) set.seed(seed)
  
  if(discrete==F){
    #start computation time
    start.time <- Sys.time()
    
    #Define probability matrix to contain probability of show rate by appointment slot
    if(is.matrix(probs)){
      prob_matrix <- probs
    } else{
      prob_matrix <- matrix(rep(probs,S), nrow=S, ncol=3, byrow = T)
    }
    
    #define mu and sd vectors
    if(is.matrix(mu)){
      st_mu <- as.vector(c(mu[1,], mu[2,])) #if supplied as matrix, convert to vector for faster computation
    } else{
      st_mu <- rep(mu, S*2)
    }
    
    if(is.matrix(sd)){
      st_sd <- as.vector(c(sd[1,], sd[2,])) #if supplied as matrix, convert to vector for faster computation
    } else{
      st_sd <- rep(sd, S*2)
    }
    
    estimated_dist <- inner_sim(N, st_mu, st_sd, block.time, prob_matrix)
    service_dist <- estimated_dist$serv
    wait_dist <- estimated_dist$waits
    
    #End computation time
    end.time <- Sys.time()
    
    #Return list object
    return(if(print.time==T){
      list(service_dist =  service_dist,
           wait_dist = wait_dist,
           discrete_probs = rcpp_prob_matrix(probs),
           show_probs = probs,
           simulation_time = print(paste("Distributuion simulated across ", N, " replications in ", round(end.time-start.time,3), " Seconds", sep="")))
    } else {
      list(service_dist = service_dist,
           wait_dist = wait_dist,
           discrete_probs = rcpp_prob_matrix(probs),
           show_probs = probs,
           simulation_time = paste("Distributuion simulated across ", N, " replications in ", round(end.time-start.time,3), " Seconds", sep="")
      )
    })
  } else if(discrete==T){
    return(list(discrete_probs = rcpp_prob_matrix(probs)))
  }
}

# Compute theoretical probability distribution given no show, single show and double show probs
# probs: vector of length 3 containing no show, single show, and double show probabilities
# n: number of time slots in a day
# Returns from Rcpp function `rcpp_prob_matrix`
p_theory <- function(probs, S){
  #initial first row (time slot) based on fixed probabilities of no show, single show or double show
  prob_matrix <- matrix(NA, nrow=S, ncol=3)
  if(is.matrix(probs)){
    prob_matrix <- probs
  } else{
    prob_matrix[1:S,] <- matrix(rep(probs,S), nrow=S, ncol=3, byrow = T)
  }
  
  return(rcpp_prob_matrix(prob_matrix))
}



###### Get AUC Stats for test data based on training data by dataset ##########
# yields speciality, AUC as estimated on validation set, sample size of the training set, and sample size of the test set
bench_auc <- function(output){
  AUCS = n_train = n_test = NULL
  for(i in 1:length(output$departments)){
    pred_stats <- prediction(output$test$pred[[i]], output$test$y[[i]])
    AUCS[i] <- unlist(performance(pred_stats, measure="auc")@y.values)
    n_train[i] = length(unlist(output$train$y[i]))
    n_test[i] = length(unlist(output$test$y[i]))
  }
  
  data.frame(specialty= output$departments, AUC = AUCS, n_train = n_train, n_test = n_test)
}

### Get Lambda, PPR, PPV, and NPV for a given threshold `tau` and dataset `i`
gen_metrics <- function(y, pred, tau=NULL){
  N <- length(y)
  lambda <- sum(y)/N #population no-show rate
  tn <- length(y[pred<tau & y==0])
  fp <- length(y[pred>=tau & y==0])
  tp <- length(y[pred>=tau & y==1])
  fn <- length(y[pred<tau & y==1])
  
  pos_n = (tp+fp)
  neg_n = (tn+fn)
  
  if(pos_n>0 & neg_n>0){
    data.frame(lambda=lambda,
               PPR=(fp+tp)/N, #PPR
               PPV = tp/pos_n, #PPV
               NPV = tn/neg_n, #NPV
               tau=tau
               
    ) 
  } else if(pos_n == 0 & neg_n>0){
    data.frame(lambda=lambda,
               PPR=0, #PPR
               PPV = 0, #PPV
               NPV = tn/neg_n, #NPV
               tau=tau)
  } else if(pos_n>0 & neg_n==0){
    data.frame(lambda=lambda,
               PPR=1, #PPR
               PPV = tp/pos_n, #PPV
               NPV = 0, #NPV
               tau=tau)
  }
}

#### Calculate probability Matrix for appointment Slots 
# `preds`: is a matrix of predicted probabilities from function `generate_day()$prob`
# `day.complete`:  is T or F. Only applicable if `preds` is given. It denotes whether a given day is complete with no further changes
#                or if unfilled slots can still be filled going forward
# `S`: number of slots in a day. Only needed if `preds` is not given
# slot.prob: is probability a slot will be filled if `day.complete=F`
# ob.prob: is probability a booked slot will be overbooked if `day.complete=F`
# metrics: is an object vector from `get_stats` function

gen_probs <- function(preds=NULL, complete=F, S=20, slot.prob=.95, ob.prob=.9, metrics){

  if(is.null(preds)==F){
    S <- nrow(preds)
  }
  
  if(is.null(metrics)){
    return("`metrics` is missing. Must be an object from function `get_stats` needs to be included")
  }
  
  if(is.null(preds)){
    alpha <- rep(metrics$PPR, S)
    slot.prob <- rep(slot.prob, S)
    ob.prob <- rep(ob.prob, S)
    
  } else if(is.null(preds)==F & complete==T){
    
    #calculate probability first patient in slot is predicted 'no-show` (if already filled, probability of 0 or 1)
    alpha <- ifelse(is.na(preds[,1]), 0,
                    ifelse(preds[,1]>=metrics$tau, 1, 0))
    
    #calculate probability of first slot being filled (if already filled, Probability of 1)
    slot.prob <- ifelse(is.na(preds[,1])==T, 0, 1)
    
    #ob.prob
    ob.prob <- ifelse(is.na(preds[,2]), 0, 1)
    
  } else if(is.null(preds)==F & complete==F){
    alpha <- ifelse(is.na(preds[,1]), metrics$PPR,
                    ifelse(preds[,1]>=metrics$tau, 1, 0))
    
    slot.prob <- ifelse(is.na(preds[,1])==T, slot.prob, 1)
    
    ob.prob <- ifelse(is.na(preds[,2]) & is.na(preds[,1]), ob.prob, 
                      ifelse(is.na(preds[,2]) & preds[,1] > metrics$tau, ob.prob, 
                             ifelse(is.na(preds[,2])==F, 1,0)))
  }
  
  probs <- matrix(NA, nrow=S, ncol=3)
  for(i in 1:S){
    no_show <- (1-slot.prob[i]) + slot.prob[i]*((1-alpha[i])*(1-metrics$NPV) + alpha[i]*metrics$PPV*(1-ob.prob[i]) + alpha[i]*metrics$PPV*metrics$lambda*ob.prob[i])
    single_show <- slot.prob[i]*((1-alpha[i])*metrics$NPV + alpha[i]*(1-metrics$PPV)*(1-ob.prob[i]) + 
                               (alpha[i]*(1-metrics$PPV)*metrics$lambda+alpha[i]*metrics$PPV*(1-metrics$lambda))*ob.prob[i])
    double_show <- slot.prob[i]*(alpha[i]*(1-metrics$PPV)*(1-metrics$lambda)*ob.prob[i])
    probs[i,] <- c(no_show, single_show, double_show)
  }
  probs
}


# Compute expected costs associated with current overbooking scheme
# pmat: is object from 'est_dist'
# discrete: if T, then gives estimated costs for discrete blocks; if F, assumes costs based on times
# otmin.thresh: threshold for defining unacceptable length of clinic overtime 
# wait.thresh: Threshold for defining an excess waiting period per patient (default is 10 min; eg. if a patient waits > 10min, they're considered to have an excess wait)
# idle.thresh: minimum amount of service time needed for a slot for it not to be considered idle
# block.time: length of appointment blocks

est_performance <- function(pmat, 
                            otmin.thresh = otmin.thresh,
                            wait.thresh = wait.thresh, 
                            idle.thresh = idle.thresh, 
                            block.time=block.time, 
                            discrete=F ){
  if(discrete==T ){
    n <- nrow(pmat$show_probs)
    list("Expected Idle Slots" = (sum(pmat$discrete_probs[,1])), #Mean number of idle times
         "Expected Excess Waiting Patients" = (sum(apply(pmat$discrete_probs, 1, function(x) x[3:length(x)]%*%seq(from=1, to=length(x)-2, by=1)))), #Excess wait periods across day
         "Expected OT Patients" =as.vector((pmat$discrete_probs[n,][3:length(pmat$discrete_probs[n,])]%*%seq(from=1, to=length(pmat$discrete_probs[n,])-2, by=1))), #Mean number of OT patients after final appointment
         "Expected Total Patient" = (sum(apply(pmat$show_probs, 1, function(x) x[2:length(x)]%*%seq(from=1, to=length(x)-1, by=1))))
    )
  } else if(discrete==F){
   
    perf = rcpp_est_perf(pmat$wait_dist, pmat$service_dist, pmat$show_probs, 
                        otmin.thresh, wait.thresh, idle.thresh, block.time)
    
    list("Expected Idle Slots" = perf$idle_slots,
         "Expected Excess Waiting Patients" = perf$excess_pt_mean,
         "Expected Excess Wait Time Per Excess patient" = perf$excess_pt_min,
         "Expected OT Prob" = perf$ot_prob,
         "Expected OT Minutes" = perf$ot_mins,
         "Expected OT Conditional Minutes" = perf$ot_cond_min,
         "Expected Total Patient" = perf$ttl_pts)
  }
}






#### Optimize predcited schedule ####
# Optimizer: "thresh" or "costs"
# Costs: vector of length 3; representing: cost of each excess waiting patient, cost of Idle Period, and cost per 30 min of OT 
# otprob_upper: Upper bound for acceptable overtime probability for overtime wait > otmin_upper
# otmin_upper: Threshold for defining an unacceptable overtime length
# excess_wait_upper: Upperr bound for acceptable number of patients with excess waits > wait.thresh
# waith.thresh: maximum acceptable waiting time for patients excessively waiting
# y: vector of realize training outcomes
# pred: vector of realized training predicted probabilities
# predicted.probs: matrix of real-time predicted probabilities 
# S: number of appointment slots
# slot.prob: probability an open appointment slot is booked
# ob.prob: probability an already booked slot is overbooked given it is eligible
# mu: log normal mean
# sd: log normal SD
# block.time: length of appointment slots
# idle.thresh: minimum acceptable idle time during an appointment block for that block to be not be considered an idle slot
# discrete: logical; if T then assume discrete time blocks ignoring distributions on patient service times
# nsim: number of iterations to simulate distribution
# grid.start = vector of length 2 defining upper and lower decision threshold bounds to estimate performance over
# model optimizes based on an acceptable over-time probability upper bound, acceptable over_time expected min, and/or acceptable number of excess waits > wait.thresh

optimize <- function(optimizer= "thresh", costs, 
                     otprob_upper=NULL, otmin.thresh, excess_wait_upper = NULL, wait.thresh, idle.thresh,
                     y, pred, predicted.probs = NULL, S, slot.prob, ob.prob, 
                     mu, sd, 
                     block.time,   discrete=F, 
                     nsim=10000, grid.start = c(0.01, .99), seed=NULL){
  if(is.null(seed)==F) set.seed(seed)
  if(optimizer=="thresh"){
    tau <- seq(grid.start[2], grid.start[1], -.01)
    res <- matrix(NA, nrow=length(tau), ncol=8)
    colnames(res) <- c("Tau","Exp Idle Slots", "E(Excess waits)", "E(Waiting Minutes|Excess wait)", "OT Prob", "E(OT MIN)", "E(OT MIN|OT)", "E(Total Patients)")
    for(i in 1:length(tau)){
      metrics <- gen_metrics(y=y, pred=pred, tau=tau[i])
      probs <- gen_probs(preds=predicted.probs, complete=F, S=S, slot.prob=slot.prob, ob.prob=ob.prob, metrics)
      sim_data <- est_dist(N=nsim, S=S, probs=probs, mu=mu, sd=sd, block.time=block.time, discrete=F)
      res[i,1] = tau[i]
      
      res[i,2:8] <- unlist(est_performance(sim_data, 
                                           otmin.thresh=otmin.thresh,
                                           wait.thresh=wait.thresh, 
                                           idle.thresh = idle.thresh, 
                                           block.time=block.time, 
                                           discrete=F))
      if (res[i,3] >excess_wait_upper |res[i,5] > otprob_upper) break
    }
    res = res[i-1,]
    return(res)
  }
  else if(optimizer=="costs"){
    if(discrete==T){
      #Initialize search grid for  tau
      a <- grid.start[1]
      b <- grid.start[2]
      c <- (b-a)/4
      tau <- seq(from=a, to =b, by=c) #taus to search
      #find optimal tau
      total.cost <- NULL
      while(tau[length(tau)/2+1.5]-tau[length(tau)/2-.5]>.02){
        for(i in 1:length(tau)){
          tau <- seq(from=a, to =b, by=c)
          metrics <- gen_metrics(y=y, pred=pred, tau=tau[i])
          probs <- gen_probs(preds=predicted.probs, complete=F, S=S, slot.prob=slot.prob, ob.prob=ob.prob, metrics)
          sim_data <- est_dist(N=nsim, S=S, probs=probs, mu=mu, sd=sd, block.time=block.time, discrete=T)
          perf_costs <- est_performance(sim_data,otmin.thresh=otmin.thresh, wait.thresh=wait.thresh, idle.thresh = idle.thresh, block.time=block.time, discrete=T)
          total.cost[i] <- perf_costs$`Expected Excess Waiting Patients`*costs[1] + perf_costs$`Expected Idle Slots`*costs[2] + perf_costs$`Expected OT Patients`*costs[3]
        }
        a <- min(tau[which.min(total.cost)-1], tau[which.min(total.cost)])
        b <- max(tau[which.min(total.cost)+1],tau[which.min(total.cost)])
        c <- (b-a)/4
      }
    }else if(discrete==F){
      #Initialize search grid for  tau
      a <- grid.start[1]
      b <- grid.start[2]
      c <- (b-a)/4
      tau <- seq(from=a, to =b, by=c)  
      #find optimal tau
      total.cost <- NULL
      while(tau[length(tau)/2+1.5]-tau[length(tau)/2-.5]>.02){
        for(i in 1:length(tau)){
          tau <- seq(from=a, to =b, by=c)
          metrics <- gen_metrics(y=y,pred=pred, tau=tau[i])
          probs <- gen_probs(preds=predicted.probs, complete=F, S=S, slot.prob=slot.prob, ob.prob=ob.prob, metrics)
          sim_data <- est_dist(N=nsim, S=S, probs=probs, mu=mu, sd=sd, block.time=block.time, discrete=F)
          perf_costs <- est_performance(sim_data,otmin.thresh=otmin.thresh, wait.thresh=wait.thresh, idle.thresh = idle.thresh, block.time=block.time, discrete=F)
          total.cost[i] <- perf_costs$`Expected Excess Waiting Patients`*costs[1] + perf_costs$`Expected Idle Slots`*costs[2] + perf_costs$`Expected OT Minutes`*costs[3]/block.time
        }
        a <- min(tau[which.min(total.cost)-1], tau[which.min(total.cost)], na.rm=T)
        b <- max(tau[which.min(total.cost)+1], tau[which.min(total.cost)], na.rm=T)
        c <- (b-a)/4
      }    
    }
    return(list("Optimal Threshold" = tau[which.min(total.cost)], 
                "Estimated Cost" = min(total.cost),
                "Metrics" = unlist(perf_costs)))
  }
}


#### Simulates a scheduled day
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


simulate <- function(optimizer="thresh", costs, 
                     otprob_upper=NULL, otmin.thresh=NULL, excess_wait_upper = NULL, wait.thresh, idle.thresh,
                     y, pred, lead, adaptive=T, update.freq = 7, begin.ob = S/2, tau, S, 
                     slot.prob=NULL, ob.prob=NULL, mu, sd, block.time,
                     y.train, pred.train, 
                     nsim=10000, seed=NULL){
  if(is.null(seed)==F) set.seed(seed)
  
  if(adaptive==T){
    
    #determine if initial spot is filled for each appointment slot based on slot.prob
    spot_filled <- runif(S,0,1)
    spot_filled[spot_filled>=slot.prob] <- NA
    spot_filled[spot_filled<slot.prob] <- 1
    
    #Generate random patients that fill out first slots for each appointment (pulls from supplied data set)
    ids <- round(runif(S,1, length(pred)))
    
    #Generate Sx2 matrix of appointment slots for scheduling period
    a1 <- pred[ids]*spot_filled # estimated probability of no-show
    y1 <- y[ids]*spot_filled # realized outcome for show or no show
    lead1 <- lead[ids]*spot_filled # lead time for scheduling order
    a2 <- y2 <- lead2 <-  rep(NA, S) #empty vectors for second appointment at each time slot
    
    #re-order IDS by lead time
    ids <- ids[order(-lead1)]
    
    #generate frequency pattern for re-optimization
    update.sched <- (seq(from=max(lead1, na.rm=T), to=0, by=-update.freq))
    days <- (rev(sort(unique(lead1))))
    n_days <- length(days)
    
    for(i in 1:n_days){
      i = i+1
      ## Update Threshold (tau) if necessary
      if(i ==1){
        #Initialize threshold estimate with no patients currently scheduled 
        tau <- optimize(optimizer=optimizer, costs, 
                        otprob_upper=otprob_upper, otmin.thresh = otmin.thresh, excess_wait_upper=excess_wait_upper, wait.thresh=wait.thresh, 
                        y=y.train, pred=pred.train, predicted.probs = NULL, S=S, slot.prob=slot.prob, 
                        ob.prob=ob.prob, mu=mu, sd=sd, block.time=block.time, idle.thresh=idle.thresh, 
                        discrete=F, nsim=nsim)[1]
        
        #Define time at which last optimization occurred 
        last_update <- update.sched[1]
        
        #vector to hold current patients in initial appointment by slot
        a1.temp <- rep(NA,S) 
        a1.temp[which(lead1==days[i])] <- a1[which(lead1==days[i])]
        
        
      } else if (length(a1.temp[is.na(a1.temp)==F])<begin.ob){
        
        a1.temp[which(lead1==days[i])] <- a1[which(lead1==days[i])]
        
      } else if(length(a1.temp[is.na(a1.temp)==F])>=begin.ob){
        
        #cycle through optimized updates from the most recent update preceeding the i-1th patient up to (but not including) the one before the i'th patient
        temp_update <- update.sched[update.sched - days[i]<=0]
        temp_update <- update.sched[(length(update.sched)- length(temp_update)): length(update.sched)]
        
        #If last simulated update was before the last scheduled update, then update tau
        update_i <- 1
        update_tau <- ifelse(last_update == temp_update[1], 0, 1)
        
        while(update_tau==1){
          #Updates threshold based on estiamted predicted probability for each slot show-rate
          tau <- optimize(optimizer=optimizer, costs, 
                          otprob_upper=otprob_upper, otmin.thresh=otmin.thresh, excess_wait_upper=excess_wait_upper, wait.thresh=wait.thresh,
                          y=y.train, pred=pred.train, predicted.probs = cbind(a1.temp, a2), S=S, slot.prob=slot.prob, 
                          ob.prob=ob.prob, mu=mu, sd=sd, block.time=block.time,  idle.thresh=idle.thresh, 
                          discrete=F, nsim=nsim, grid.start=c(max(tau-.15,0.01), tau+.15))[1]
          
          #update last updatem time point
          last_update <- temp_update[update_i]
          
          # Allocate overbooking based on uniform probability of it occurring from now until day of appointment 
          time_to_next_update  <- min(temp_update[update_i], days[i]) - ifelse(is.na(temp_update[update_i+1]),-1,temp_update[update_i+1])   #Time till next update
          prob_ob_temp <- time_to_next_update*slot.prob/(min(temp_update[update_i], days[i])+1) #probability of overbooking before next optimization
          
          #Get vectors of ids from which to draw randomly for next update cycle
          max.lead <- min(days[i], temp_update[update_i]) #Calculate time left until appointment
          min.lead <- max.lead-time_to_next_update+1
          potential_obs <- pred[lead <= max.lead & lead >=min.lead] #subsets data to only with patients with lead time <= most recent patient
          potential_y <- y[lead <= max.lead & lead >=min.lead] #subset realized outcomes 
          potential_lead <- lead[lead <= max.lead & lead >=min.lead]
          upper.thresh <- length(potential_obs) #Defines upper threshold for random selection among temp_data
          
          #get vector of patients eligible for over-booking
          a1.update <- a1.temp[is.na(a1.temp)==F & a1.temp>=tau & is.na(a2)==T]
          ob_indicator <- ifelse(runif(length(a1.update))<=prob_ob_temp, 1, NA)
          
          if(length(ob_indicator[is.na(ob_indicator)==F])>0){
            temp_ids <-  round(runif(length(ob_indicator[is.na(ob_indicator)==F]), 1, upper.thresh))
            ob_indicator[is.na(ob_indicator)==F] <- potential_y[temp_ids]
            y2[is.na(a1.temp)==F & a1.temp>=tau & is.na(a2)==T] <- ob_indicator #updates actual outcomes for overbooked appointments
            ob_indicator[is.na(ob_indicator)==F] <- potential_lead[temp_ids]
            lead2[is.na(a1.temp)==F & a1.temp>=tau & is.na(a2)==T] <- ob_indicator # update lead times
            ob_indicator[is.na(ob_indicator)==F] <- potential_obs[temp_ids]
            a2[is.na(a1.temp)==F & a1.temp>=tau & is.na(a2)==T] <- ob_indicator #updates overbooked appointment probabilities
          }
          
          #update time at which re-optimization occurs
          update_i <- update_i+1
          
          #Update keep going counter
          update_tau <- ifelse(length(temp_update)==1 | ifelse(is.na(temp_update[update_i]), -1,temp_update[update_i])  <= ifelse(is.na(days[i+1]), 0, days[i+1]) |  sum(is.na(ob_indicator))==0, 0, 1)
        }
        #Update next slot
        a1.temp[which(lead1==days[i])] <- a1[which(lead1==days[i])]
      }
      
    }
    
    #service time distribution parameters
    if(is.matrix(mu)){
      st_mu <- st_mu
    } else{
      st_mu <- matrix(mu, nrow=2, ncol=S)
    }
    
    if(is.matrix(sd)){
      st_sd <- st_sd
    } else{
      st_sd <- matrix(sd, nrow=2, ncol=S)
    }
    
    times <- (cbind(y1,y2)-1)*-1
    times <- t(matrix(rlnorm(S*2, mean=st_mu, sd=st_sd), nrow=2, byrow=F)*t(times))
    
    #Return List
    list(pred = cbind(a1.temp, a2),
         y = cbind(y1,y2), 
         times = times,
         lead = cbind(lead1, lead2)
    )
    

    
  } else if(adaptive==F){
    #determine if initial spot is filled for each appointment slot 
    spot_filled <- runif(S,0,1)
    spot_filled[spot_filled>=slot.prob] <- NA
    spot_filled[spot_filled<slot.prob] <- 1
    
    #Generate random patients that fill out first slots for each appointment
    ids <- round(runif(S,1, length(pred)))
    
    #Generate Sx2 matrix of appointment slots for period
    a1 <- pred[ids]*spot_filled
    y1 <- y[ids]*spot_filled
    lead1 <- lead[ids]*spot_filled
    a2 <- y2 <- lead2 <-  NULL
    
    for(i in 1:length(a1)){
      
      temp_data <- pred[lead < lead[ids[i]]]
      upper.thresh <- length(temp_data)
      
      if(is.na(a1[i]) | a1[i]<tau | upper.thresh<1){
        a2[i] <- y2[i] <- lead2[i] <-NA
      } else if(a1[i]>=tau){
        ob_indicator <- ifelse(runif(1)>=ob.prob, NA, 1)
        temp_id <- round(runif(1, 1, upper.thresh))
        a2[i] <- ob_indicator*pred[lead < lead[ids[i]]][temp_id]
        y2[i] <- ob_indicator*y[lead < lead[ids[i]]][temp_id]
        lead2[i] <- ob_indicator*lead[lead < lead[ids[i]]][temp_id]
      } 
    }
    
    #service time distribution parameters
    if(is.matrix(mu)){
      st_mu <- st_mu
    } else{
      st_mu <- matrix(mu, nrow=2, ncol=S)
    }
    
    if(is.matrix(sd)){
      st_sd <- st_sd
    } else{
      st_sd <- matrix(sd, nrow=2, ncol=S)
    }
    times <- (cbind(y1,y2)-1)*-1
    times <- t(matrix(rlnorm(S*2, mean=st_mu, sd=st_sd), nrow=2, byrow=F)*t(times))
    
    #Return List
    list(pred = cbind(a1,a2),
         y = cbind(y1,y2), 
         times = times,
         lead= cbind(lead1, lead2)
    )
  }
}




#### Analyze Simulation Results from object of simulation()
# give list of results from simulation() object. Can be from multiple simulations
# wait.thresh and idle.thresh are same as defined in prior functions
wait.thresh=15
otmin.thresh=15
otprob_upper=.10
excess_wait_upper=5
idle.thresh =10
test=results
analyze <- function(test, otmin.thresh, wait.thresh, idle.thresh, block.time, costs){
  if(length(test)==4){
    out <- analyze_sim(test, otmin.thresh=otmin.thresh, wait.thresh=wait.thresh, idle.thresh=idle.thresh, block.time=block.time, costs=costs)
    
    data.frame(idle_slots=out[[1]], "excess_wait"=out[[2]], "mean_excess_wait" = out[[3]], 
               "ot_min"= out[[4]]," total_patients" = out[[5]])
  } else if (length(test)>4){
    out <- lapply(test, function(x) analyze_sim(x, otmin.thresh=otmin.thresh, wait.thresh=wait.thresh, idle.thresh=idle.thresh, block.time=block.time, costs=costs))
    out <- mapply(c, out)
    out <- apply(matrix(unlist(out), nrow=7, byrow=F), 1, mean)
    data.frame(idle_slots=out[1], excess_wait=out[2], mean_excess_wait = out[3], 
               ot_min= out[4], ot_prob =out[5],  total_patients = out[6], costs=out[7])
  }
}


# Used inside function analyze. Analyzes single list output
analyze_sim <- function(test, otmin.thresh, wait.thresh, idle.thresh, block.time, costs){
  S <- nrow(test$times)
  #initialize list to contain distribution
  service_dist <- vector(mode="integer", length=S)
  wait_dist <- vector(mode="integer", length=S*2)
  
  #service Time Matrix (rows for patient, columns for appointment slot)
  service_times <- t(test$times)
  service_times[is.na(service_times)] <- 0
  
  #initialize wait times
  wait_times <- vector(mode="integer", length=S*2)
  
  #Compute Service Times by appointment slot
  wait_times[1:2] <- c(0,(service_times[1,1])*ifelse(service_times[2,1]>0,1,0))
  service_dist[1] <- service_times[1,1] + service_times[2,1]
  
  second_apt_index <- rep(0,S)
  second_apt_index[service_times[2,]>0] <- 1
  for(i in 2:S){
    i2 <- i*2
    carry_over <- max(service_dist[i-1]-block.time,0)
    st_index <- service_times[1,i]
    wait_times[(i2-1)] <- carry_over
    wait_times[i2] <- (st_index+carry_over)*second_apt_index[i]
    service_dist[i] <- st_index + service_times[2,i] + carry_over
  }
  
  idle_slots <- length(service_dist[service_dist<idle.thresh])
  excess_waits <- length(wait_times[wait_times>wait.thresh])
  excess_wait_mean <- sum(wait_times[wait_times>wait.thresh])/length(wait_times[wait_times>wait.thresh])
  excess_wait_mean <- ifelse(is.na(excess_wait_mean)==T, 0, excess_wait_mean)
  ot_min <- max(service_dist[length(service_dist)]-block.time, 0)
  ot <- max(service_dist[length(service_dist)]-block.time-otmin.thresh, 0)
  total_patients <- length(test$y[test$y==0 & is.na(test$y)==F])
  
  cost <- idle_slots*costs[2] + excess_waits*costs[1] + ot_min*costs[3]/block.time
  
  list("Idle Slots"  = idle_slots, 
       "Excess Waiting Patients" = excess_waits,
       "Mean Excess Wait Time" = excess_wait_mean,
       "OT Minutes" = ot_min,
       "OT" = ifelse(ot>0, 1, 0),
       "Total Patients" = total_patients, 
       "Cost" = cost)
}

