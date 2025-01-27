//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>
#include <thread>
using namespace Rcpp;

// [[Rcpp::export]]
void inner_sim_bench(int N, NumericVector mu, NumericVector sd, int block_time, NumericMatrix probs) {
            
 Rcpp::Clock clock;

 clock.tick("entire");

            // initialize 
            int S2 = mu.size();
            int S = S2/2;
            double carry_over;
            double second_apt_index;
            NumericMatrix wait_times(N,S2);
            NumericMatrix service_dist_temp(N, S);
            NumericVector service_times(S2);
            NumericVector apt_pts(S2);
            NumericVector p0 = probs( _ , 0);
            NumericVector p1 = probs( _ , 1);
            NumericVector p2 = probs( _ , 2);
            NumericVector runis(S);
            
            for(int k=0; k<N; k++){
            
                  // Random sample from uniform distribution for calculating no show, single show or double show
                runis = Rcpp::runif(S);
                
                for(int i=0; i < S; i++){
                  
		clock.tick("runi" + std::to_string(k));
                  double r1 = runis[i];
                  if(r1<=  p0[i]){
                    apt_pts[i] = 0;
                    apt_pts[S+i] = 0;
                  }else if(r1 >  p0[i] && r1 <=p1[i] +p0[i]){
                    apt_pts[i] = 1;
                    apt_pts[S+i] = 0;
                  }else if(r1 > p1[i] + p0[i]){
                    apt_pts[i] = 1;
                    apt_pts[S+i] = 1;
                  }
		clock.tock("runi"+ std::to_string(k));

                 clock.tick("st"+ std::to_string(k));
                  // Calculate service times from normal distribution multipled by indicator apt_pts
                  service_times[i] = apt_pts[i]*R::rlnorm(mu[i], sd[i]);
                  service_times[S+i] = apt_pts[S+i]*R::rlnorm(mu[S+i], sd[S+i]);
                  clock.tock("st"+ std::to_string(k));
			
		clock.tick("main"+ std::to_string(k));
                  // calculate wait time and service times
                  if(i==0){
                    carry_over = 0;
                  } else if(service_dist_temp(k,i-1)<=block_time){
                    carry_over = 0;
                  } else{
                    carry_over = service_dist_temp(k,i-1)-block_time;
                  }
    
                  if(service_times[S+i]>0){
                    second_apt_index = 1;
                  } else{
                    second_apt_index = 0;
                  }
    
                  wait_times(k, i*2) = carry_over;
                  wait_times(k, i*2+1) = (service_times[i]+carry_over)*second_apt_index;
                  service_dist_temp(k,i) = service_times[i] + service_times[S+i] + carry_over;
		clock.tock("main"+ std::to_string(k));

                }
 
            }
clock.tock("entire");
            
clock.stop("inner_sim_times");

          

 }
