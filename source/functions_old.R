### Define inside main loop of est_dist
# cppFunction('List cycle_dist(NumericVector mu, NumericVector sd, int block_time, NumericMatrix probs) {
#             
#             // initialize 
#             int S2 = mu.size();
#             int S = S2/2;
#             int j;
#             double carry_over;
#             double second_apt_index;
#             NumericVector wait_times(S2);
#             NumericVector service_dist_temp(S);
#             NumericVector service_times(S2);
#             NumericVector apt_pts(S2);
#             NumericVector p0 = probs( _ , 0);
#             NumericVector p1 = probs( _ , 1);
#             NumericVector p2 = probs( _ , 2);
#             NumericVector apts(S);
#             NumericVector runis(S);
#             
#             // Random sample from uniform distribution for calculating no show, single show or double show
#             runis = Rcpp::runif(S);
#             for(int i=0; i < S; i++){
#               if(runis[i]<= p0[i]){
#                 apts[i] = 0;
#               }else if(runis[i] > p0[i] && runis[i] <= p1[i] + p0[i]){
#                 apts[i] = 1; 
#               }else if(runis[i] > p1[i]+p0[i]){
#                 apts[i] = 2;
#               }
#             }
#             
#            // Expand apt shows to a longer vector of indicators for whether or not a patient showed or not
#            for(int i = 0; i < S; i++){
#               if(apts[i]>=1){
#                 apt_pts[i] = 1;
#               } else if (apts[i]==0){
#                 apt_pts[i] = 0;
#               }
# 
#               if(apts[i]==2){
#                 apt_pts[S+i] = 1;
#               }else if(apts[i]<2){
#                 apt_pts[S+i] = 0;
#               }
#             }
#             
#             // Calculate service times from normal distribution multipled by indicator apt_pts 
#             for(int i=0; i < S2; i++){
#              service_times[i] = apt_pts[i]*R::rnorm(mu[i], sd[i]);
#             }
#             
#             // Cycle over appointments sequentially and calculate wait time and service times
#             for(int i=0; i < S; i++){
#               j = i+1;
# 
#               if(i==0){
#                 carry_over = 0;
#               } else if(service_dist_temp[i-1]<=block_time){
#                 carry_over = 0;
#               } else{
#                 carry_over = service_dist_temp[i-1]-block_time;
#               }
# 
#               if(service_times[S+i]>0){
#                 second_apt_index = 1;
#               } else{
#                 second_apt_index = 0;
#               }
# 
#               wait_times[i*2] = carry_over;
#               wait_times[i*2+1] = (service_times[i]+carry_over)*second_apt_index;
#               service_dist_temp[i] = service_times[i] + service_times[S+i] + carry_over;
#             }
# 
#             return Rcpp::List::create(Rcpp::Named("waits") = wait_times,
#                                       Rcpp::Named("serv")=service_dist_temp);
# 
#             }')
# 
# ### Inner est_dist simulation loop in Rcpp
# cppFunction('List inner_sim_old(int N, NumericVector mu, NumericVector sd, int block_time, NumericMatrix probs) {
#             
#             // initialize 
#             int S2 = mu.size();
#             int S = S2/2;
#             double carry_over;
#             double second_apt_index;
#             NumericMatrix wait_times(N,S2);
#             NumericMatrix service_dist_temp(N, S);
#             NumericVector service_times(S2);
#             NumericVector apt_pts(S2);
#             NumericVector p0 = probs( _ , 0);
#             NumericVector p1 = probs( _ , 1);
#             NumericVector p2 = probs( _ , 2);
#             NumericVector apts(S);
#             NumericVector runis(S);
#             
#             for(int k=0; k<N; k++){
#                   // Random sample from uniform distribution for calculating no show, single show or double show
#                 runis = Rcpp::runif(S);
#                 for(int i=0; i < S; i++){
#                   if(runis[i]<= p0[i]){
#                     apts[i] = 0;
#                   }else if(runis[i] > p0[i] && runis[i] <= p1[i] + p0[i]){
#                     apts[i] = 1; 
#                   }else if(runis[i] > p1[i]+p0[i]){
#                     apts[i] = 2;
#                   }
#                 }
#                 
#                // Expand apt shows to a longer vector of indicators for whether or not a patient showed or not
#                for(int i = 0; i < S; i++){
#                   if(apts[i]>=1){
#                     apt_pts[i] = 1;
#                   } else if (apts[i]==0){
#                     apt_pts[i] = 0;
#                   }
#     
#                   if(apts[i]==2){
#                     apt_pts[S+i] = 1;
#                   }else if(apts[i]<2){
#                     apt_pts[S+i] = 0;
#                   }
#                 }
#                 
#                 // Calculate service times from normal distribution multipled by indicator apt_pts 
#                 for(int i=0; i < S2; i++){
#                  service_times[i] = apt_pts[i]*R::rnorm(mu[i], sd[i]);
#                 }
#                 
#                 // Cycle over appointments sequentially and calculate wait time and service times
#                 for(int i=0; i < S; i++){
# 
#                   if(i==0){
#                     carry_over = 0;
#                   } else if(service_dist_temp(k,i-1)<=block_time){
#                     carry_over = 0;
#                   } else{
#                     carry_over = service_dist_temp(k,i-1)-block_time;
#                   }
#     
#                   if(service_times[S+i]>0){
#                     second_apt_index = 1;
#                   } else{
#                     second_apt_index = 0;
#                   }
#     
#                   wait_times(k, i*2) = carry_over;
#                   wait_times(k, i*2+1) = (service_times[i]+carry_over)*second_apt_index;
#                   service_dist_temp(k,i) = service_times[i] + service_times[S+i] + carry_over;
#                 }
#               
#             }
#             
# 
#             return Rcpp::List::create(Rcpp::Named("waits") = wait_times,
#                                       Rcpp::Named("serv")= service_dist_temp);
# 
#  }')
# 
# 
# 
# cppFunction('List inner_sim(int N, NumericVector mu, NumericVector sd, int block_time, NumericMatrix probs) {
#             
#             // initialize 
#             int S2 = mu.size();
#             int S = S2/2;
#             double carry_over;
#             double second_apt_index;
#             NumericMatrix wait_times(N,S2);
#             NumericMatrix service_dist_temp(N, S);
#             NumericVector service_times(S2);
#             NumericVector apt_pts(S2);
#             NumericVector p0 = probs( _ , 0);
#             NumericVector p1 = probs( _ , 1);
#             NumericVector p2 = probs( _ , 2);
#             NumericVector runis(S);
#             
#             for(int k=0; k<N; k++){
#             
#                   // Random sample from uniform distribution for calculating no show, single show or double show
#                 runis = Rcpp::runif(S);
#                 
#                 for(int i=0; i < S; i++){
#                   
#                   double r1 = runis[i];
#                   double p0_temp = p0[i];
#                   double p1_temp = p1[i];
#                   if(r1<=  p0_temp){
#                     apt_pts[i] = 0;
#                     apt_pts[S+i] = 0;
#                   }else if(r1 >  p0_temp && r1 <=p1_temp + p0_temp){
#                     apt_pts[i] = 1;
#                     apt_pts[S+i] = 0;
#                   }else if(r1 > p1_temp + p0_temp){
#                     apt_pts[i] = 1;
#                     apt_pts[S+i] = 1;
#                   }
#                   
#                   // Calculate service times from normal distribution multipled by indicator apt_pts
#                   service_times[i] = apt_pts[i]*R::rnorm(mu[i], sd[i]);
#                   service_times[S+i] = apt_pts[S+i]*R::rnorm(mu[S+i], sd[S+i]);
#                   
#                   // calculate wait time and service times
#                   if(i==0){
#                     carry_over = 0;
#                   } else if(service_dist_temp(k,i-1)<=block_time){
#                     carry_over = 0;
#                   } else{
#                     carry_over = service_dist_temp(k,i-1)-block_time;
#                   }
#     
#                   if(service_times[S+i]>0){
#                     second_apt_index = 1;
#                   } else{
#                     second_apt_index = 0;
#                   }
#     
#                   wait_times(k, i*2) = carry_over;
#                   wait_times(k, i*2+1) = (service_times[i]+carry_over)*second_apt_index;
#                   service_dist_temp(k,i) = service_times[i] + service_times[S+i] + carry_over;
# 
#                 }
#  
#             }
#             
# 
#             return Rcpp::List::create(Rcpp::Named("waits") = wait_times,
#                                       Rcpp::Named("serv")= service_dist_temp);
# 
#  }')


#Define Inner for loop in rcpp for est_dist
# cppFunction('List cycle_dist(NumericVector service_times, int block_time) {
#             
#             int S2 = service_times.size();
#             int S = S2/2;
#             int j;
#             double carry_over;
#             double second_apt_index;
#             NumericVector wait_times(S2);
#             NumericVector service_dist_temp(S);
# 
#             for(int i=0; i < S; i++){
#               j = i+1;
#               
#               if(i==0){
#                 carry_over = 0;
#               } else if(service_dist_temp[i-1]<=block_time){
#                 carry_over = 0;
#               } else{
#                 carry_over = service_dist_temp[i-1]-block_time;
#               }
#               
#               if(service_times[S+i]>0){
#                 second_apt_index = 1;
#               } else{
#                 second_apt_index = 0;
#               }
#               
#               wait_times[i*2] = carry_over;
#               wait_times[i*2+1] = (service_times[i]+carry_over)*second_apt_index;
#               service_dist_temp[i] = service_times[i] + service_times[S+i] + carry_over;
#             }
#             
#             return Rcpp::List::create(Rcpp::Named("waits") = wait_times,
#                                       Rcpp::Named("serv")=service_dist_temp);
# 
#             }')

# define random appointment allocation by slot in RCPP for est_dist
# cppFunction('List apts(NumericMatrix probs) {
#             
#             NumericVector p0 = probs( _ , 0);
#             NumericVector p1 = probs( _ , 1);
#             NumericVector p2 = probs( _ , 2);
#             int length_p = p0.size();
#             NumericVector out(length_p);
#             NumericVector runis(length_p);
#             
#             runis = Rcpp::runif(length_p);
# 
#             for(int i=0; i < length_p; i++){
#               if(runis[i]<= p0[i]){
#                 out[i] = 0;
#               }else if(runis[i] > p0[i] && runis[i] <= p1[i] + p0[i]){
#                 out[i] = 1; 
#               }else if(runis[i] > p1[i]+p0[i]){
#                 out[i] = 2;
#               }
#             }
#             return Rcpp::List::create(Rcpp::Named("apts") = out);
#             }')

# 
# est_dist_old <- function(N=5000, S=20, probs=c(.1,.8,.1), mu=24, sd=4, block.time = 30, discrete=F, seed=NULL, print.time = F){
#   if(is.null(seed)==F) set.seed(seed)
#   
#   if(discrete==F){
#     #start computation time
#     start.time <- Sys.time()
#     
#     #Initialize S2 (slightly faster than doing S*2 in multiple places)
#     S2 = S*2
#     
#     #initialize list to contain distribution
#     service_dist <- matrix(NA, nrow=N, ncol=S)
#     wait_dist <- matrix(NA, ncol=S*2, nrow=N)
#     
#     #Define probability matrix to contain probability of show rate by appointment slot
#     if(is.matrix(probs)){
#       prob_matrix <- probs
#     } else{
#       prob_matrix <- matrix(rep(probs,S), nrow=S, ncol=3, byrow = T)
#     }
#     
#     #define mu and sd vectors
#     if(is.matrix(mu)){
#       st_mu <- as.vector(c(mu[1,], mu[2,])) #if supplied as matrix, convert to vector for faster computation
#     } else{
#       st_mu <- rep(mu, S2)
#     }
#     
#     if(is.matrix(sd)){
#       st_sd <- as.vector(c(sd[1,], sd[2,])) #if supplied as matrix, convert to vector for faster computation
#     } else{
#       st_sd <- rep(sd, S2)
#     }
#    
#   
#     #initialize vectors to hold data
#     service_dist_temp <- rep(0, S)
#     #initialize wait times
#     wait_times <- rep(0, S2)
#     
#     for(k in 1:N){
#       #Number of patients to each appointment slot
#       apt_pts = {
#         nprobmat <- nrow(prob_matrix)
#         out <- rep(NA_real_, nprobmat)
#         runis <- runif(nprobmat) 
#         out[runis<prob_matrix[,1]] <- 0L
#         out[runis>=prob_matrix[,1] & runis < (1-prob_matrix[,3])] <- 1L
#         out[is.na(out)] <- 2L
#         out
#       }
#       
#       #service Time Matrix (rows for patient, columns for appointment slot)
#       service_times <- rnorm(S2, mean=st_mu, sd=st_sd)*c(apt_pts>=1,apt_pts>=2)
#       
#       #Compute Service Times by appointment slot
#       wait_times[1:2] <- c(0,(service_times[1])*ifelse(service_times[S+1]>0,1,0))
#       service_dist_temp[1] <- service_times[1] + service_times[S+1]
#       
#       second_apt_index <- rep(0,S)
#       second_apt_index[service_times[(S+1):(S2)]>0] <- 1
#       for(i in 2:S){
#         i2 <- i*2
#         if(service_dist_temp[i-1]<=30){
#           carry_over<-0
#         } else{
#           carry_over <- service_dist_temp[i-1]-block.time
#         }
#         st_index <- service_times[i]
#         wait_times[(i2-1)] <- carry_over
#         wait_times[i2] <- (st_index+carry_over)*second_apt_index[i]
#         service_dist_temp[i] <- st_index + service_times[S+i] + carry_over
#       }
#       service_dist[k,] <- service_dist_temp
#       wait_dist[k,] <- wait_times
#     }
#     
#     #End computation time
#     end.time <- Sys.time()
#     
#     #Return list object
#     return(if(print.time==T){
#       list(service_dist =  lapply(seq_len(ncol(service_dist)), function(x) service_dist[,x]),
#            wait_dist = wait_dist,
#            discrete_probs = p_theory(probs,S),
#            simulation_time = print(paste("Distributuion simulated across ", N, " replications in ", round(end.time-start.time,3), " Seconds", sep="")))
#     } else {
#       list(service_dist = lapply(seq_len(ncol(service_dist)), function(x) service_dist[,x]),
#            wait_dist = wait_dist,
#            discrete_probs = p_theory(probs,S),
#            simulation_time = paste("Distributuion simulated across ", N, " replications in ", round(end.time-start.time,3), " Seconds", sep="")
#       )
#     })
#   } else if(discrete==T){
#     return(list(discrete_probs = p_theory(probs,S)))
#   }
# }

# est_dist_old <- function(N=5000, S=20, probs=c(.1,.8,.1), mu=24, sd=4, block.time = 30, discrete=F, seed=NULL, print.time = F){
#   if(is.null(seed)==F) set.seed(seed)
#   
#   if(discrete==F){
#     #start computation time
#     start.time <- Sys.time()
#     
#     #initialize list to contain distribution
#     service_dist <- matrix(NA, nrow=N, ncol=S)
#     wait_dist <- matrix(NA, ncol=S*2, nrow=N)
#     
#     #Define probability matrix to contain probability of show rate by appointment slot
#     if(is.matrix(probs)){
#       prob_matrix <- probs
#     } else{
#       prob_matrix <- matrix(rep(probs,S), nrow=S, ncol=3, byrow = T)
#     }
#     
#     #define mu and sd vectors
#     if(is.matrix(mu)){
#       st_mu <- as.vector(c(mu[1,], mu[2,])) #if supplied as matrix, convert to vector for faster computation
#     } else{
#       st_mu <- rep(mu, 2*S)
#     }
#     
#     if(is.matrix(sd)){
#       st_sd <- as.vector(c(sd[1,], sd[2,])) #if supplied as matrix, convert to vector for faster computation
#     } else{
#       st_sd <- rep(sd, 2*S)
#     }
#     
#     #initialize vectors to hold data
#     service_dist_temp <- NULL
#     #initialize wait times
#     wait_times <- vector(mode="integer", length=S*2)
#     
#     for(k in 1:N){
#       #Number of patients to each appointment slot
#       apt_pts = {
#         out <- rep(NA_real_, nrow(prob_matrix))
#         runis <- runif(nrow(prob_matrix)) 
#         out[runis<prob_matrix[,1]] <- 0L
#         out[runis>=prob_matrix[,1] & runis < (1-prob_matrix[,3])] <- 1L
#         out[is.na(out)] <- 2L
#         out
#       }
#       
#       #service Time Matrix (rows for patient, columns for appointment slot)
#       service_times <- matrix(rnorm(S*2, mean=st_mu, sd=st_sd)*c(apt_pts>=1,apt_pts>=2), nrow=2, byrow=T)
#       
#       #Compute Service Times by appointment slot
#       wait_times[1:2] <- c(0,(service_times[1,1])*ifelse(service_times[2,1]>0,1,0))
#       service_dist_temp[1] <- service_times[1,1] + service_times[2,1]
#       
#       second_apt_index <- rep(0,S)
#       second_apt_index[service_times[2,]>0] <- 1
#       
#       for(i in 2:S){
#         i2 <- i*2
#         carry_over <- max(service_dist_temp[i-1]-block.time,0)
#         st_index <- service_times[1,i]
#         wait_times[(i2-1)] <- carry_over
#         wait_times[i2] <- (st_index+carry_over)*second_apt_index[i]
#         service_dist_temp[i] <- st_index + service_times[2,i] + carry_over
#       }
#       service_dist[k,] <- service_dist_temp
#       wait_dist[k,] <- wait_times
#     }
#     
#     #End computation time
#     end.time <- Sys.time()
#     
#     #Return list object
#     return(if(print.time==T){
#       list(service_dist =  lapply(seq_len(ncol(service_dist)), function(x) service_dist[,x]),
#            wait_dist = wait_dist,
#            discrete_probs = p_theory(probs,S),
#            simulation_time = print(paste("Distributuion simulated across ", N, " replications in ", round(end.time-start.time,3), " Seconds", sep="")))
#     } else {
#       list(service_dist = lapply(seq_len(ncol(service_dist)), function(x) service_dist[,x]),
#            wait_dist = wait_dist,
#            discrete_probs = p_theory(probs,S),
#            simulation_time = paste("Distributuion simulated across ", N, " replications in ", round(end.time-start.time,3), " Seconds", sep="")
#       )
#     })
#   } else if(discrete==T){
#     return(list(discrete_probs = p_theory(probs,S)))
#   }
# }


# #### Inner est_dist simulation loop in Rcpp
# cppFunction('List rcpp_est_perf(NumericMatrix wait_dist, double wait_thresh, 
#                                 NumericMatrix service_dist, double block_time) {
# 
#             double N = wait_dist( _ , 1).size();
#             double S = wait_dist( 1, _ ).size();
#             double ind; 
#             
#             //excess patients
#             double excess_pt = 0;
#             for(int i=0; i<N; i++){
#               double temp = 0;
#               for(int k= 0; k < S; k++){
#                 if(wait_dist(i,k)>wait_thresh){
#                   ind = 1;
#                 } else{
#                   ind =0;
#                 }
#                 temp += ind;
#               }
#               excess_pt += temp;
#             }
#             double excess_pt_mean = excess_pt/N;
#             
#             
#            //excess minutes
#            double temp=0;
#            double temp_size = 0;
#            for(int i=0; i<N; i++){
#               for(int k=0; k<S; k++){
#                 if(wait_dist(i,k)>wait_thresh){
#                   temp += wait_dist(i,k);
#                   temp_size +=1;
#                 } 
#               }
#            }
#            double excess_min_avg = temp/temp_size;
#            
#            //idle slots
#            double idle_temp=0;
#            for(int i=0; i<N; i++){
#               for(int k=0; k<S/2; k++){
#                 if(service_dist(i,k) == 0 ){
#                   idle_temp += 1;
#                 }
#               }
#            }
#            double idle_slots = idle_temp/N;
#            
#            //OT Minutes
#            NumericVector ot_mins = service_dist(_, S/2-1);
#            double ot_min_sum = 0;
#            for(int i= 0; i<N; i++){
#               if(ot_mins[i] > block_time){
#                 ot_min_sum += ot_mins[i]-block_time; 
#               }
#            }
#            double ot = ot_min_sum/N;
#             
# 
#             return Rcpp::List::create(Rcpp::Named("excess_pt_mean") = excess_pt_mean,
#             Rcpp::Named("excess_pt_min") =excess_min_avg,
#             Rcpp::Named("idle_slots") =idle_slots,
#             Rcpp::Named("ot_mins") =ot); 
#                                       
#             }')
#
#
##### Estimate probability matrix in Rcpp
# cppFunction('NumericMatrix rcpp_prob_matrix(NumericMatrix probs) {
# 
#             double S = probs(_,1).size();
#             
#             //initialize matrix of 0s
#             NumericMatrix m(S,S+3);
#             std::fill( m.begin(), m.end(), 0 ) ;
#             
#             for(int i=0; i<S; i++){
#             if(i>0){
#               m(i,0) = probs(i,0)*(m(i-1,0) + m(i-1, 1));
#               m(i,1) = probs(i,1)*(m(i-1,0)+m(i-1,1)) + probs(i,0)*m(i-1,2);
#               m(i,2) = probs(i,2)*(m(i-1, 0)+m(i-1, 1)) + probs(i,1)*m(i-1,2) + probs(i,0)*m(i-1,3);
#               
#               for(int k=0; k < S-1; k++){
#                 m(i, k+3) = probs(i,2)*m(i-1, k+2) + probs(i,1)*m(i-1, k+3) + probs(i,0)*m(i-1, k+4);
#               }
#             } else if(i==0){
#               m(i,0) = probs(i,0);
#               m(i,1) = probs(i,1);
#               m(i,2) = probs(i,2);
#             }
#            }
#             
# 
#             return m; 
#                                       
#             }')
