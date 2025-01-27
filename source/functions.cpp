#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List inner_sim(int N, NumericVector mu, NumericVector sd, int block_time, NumericMatrix probs) {
            
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
                  
                  double r1 = runis[i];
                  double p0_temp = p0[i];
                  double p1_temp = p1[i];
                  if(r1<=  p0_temp){
                    apt_pts[i] = 0;
                    apt_pts[S+i] = 0;
                  }else if(r1 >  p0_temp && r1 <=p1_temp + p0_temp){
                    apt_pts[i] = 1;
                    apt_pts[S+i] = 0;
                  }else if(r1 > p1_temp + p0_temp){
                    apt_pts[i] = 1;
                    apt_pts[S+i] = 1;
                  }
                  
                  // Calculate service times from log normal distribution multipled by indicator apt_pts
                  service_times[i] = apt_pts[i]*R::rlnorm(mu[i], sd[i]);
                  service_times[S+i] = apt_pts[S+i]*R::rlnorm(mu[S+i], sd[S+i]);
                  
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

                }
 
            }
            

            return Rcpp::List::create(Rcpp::Named("waits") = wait_times,
                                      Rcpp::Named("serv")= service_dist_temp);

 }


// [[Rcpp::export]]
List rcpp_est_perf(NumericMatrix wait_dist, 
                   NumericMatrix service_dist, 
                   NumericMatrix probs,
                   double otmin_thresh, 
                   double wait_thresh,
                   double idle_thresh,
                   double block_time) {

            double N = wait_dist( _ , 1).size();
            double S = wait_dist( 1, _ ).size();
            double ind;

            //excess patients
            double excess_pt = 0;
            for(int i=0; i<N; i++){
              double temp = 0;
              for(int k= 0; k < S; k++){
                if(wait_dist(i,k)>wait_thresh){
                  ind = 1;
                } else{
                  ind =0;
                }
                temp += ind;
              }
              excess_pt += temp;
            }
            double excess_pt_mean = excess_pt/N;


           //excess minutes among those exceeding wait threshold
           double temp=0;
           double temp_size = 0;
           for(int i=0; i<N; i++){
              for(int k=0; k<S; k++){
                if(wait_dist(i,k)>wait_thresh){
                  temp += wait_dist(i,k);
                  temp_size +=1;
                }
              }
           }
           double excess_min_avg = temp/temp_size;

           //idle slots
           double idle_temp=0;
           for(int i=0; i<N; i++){
              double temp = 0;
              for(int k=0; k<S/2; k++){
                if(service_dist(i,k) < idle_thresh ){
                  ind = 1 ;
                } else{
                  ind = 0;
                }
                  temp += ind;
              }
               idle_temp += temp;
           }
           double idle_slots = idle_temp/N;

           //OT Expected Minutes and Probability
           NumericVector ot_mins = service_dist(_, S/2-1);
           double ot_min_sum = 0;
           double ot_p = 0;
           for(int i= 0; i<N; i++){
               
              if(ot_mins[i] > (block_time+otmin_thresh)){
                ot_min_sum += ot_mins[i]-block_time;
                ot_p += 1;
                  
                }
                }
           double ot = ot_min_sum/N;
           double ot_cond_min = ot_min_sum/ot_p;
           double ot_prob = ot_p/N;

           
           //total patients expected 
           double ttl_pts = 0;
           double outer_size = probs(_,1).size();
           double inner_size = probs(1,_).size();
           for(int i=0; i< outer_size; i++){
              for(int k=0; k<inner_size; k++){
                ttl_pts+= k*probs(i,k);
              }
           }

            return Rcpp::List::create(Rcpp::Named("excess_pt_mean") = excess_pt_mean,
            Rcpp::Named("excess_pt_min") =excess_min_avg,
            Rcpp::Named("idle_slots") =idle_slots,

	        Rcpp::Named("ot_prob") = ot_prob,
	        Rcpp::Named("ot_mins") =ot,
            Rcpp::Named("ot_cond_mins") =ot_cond_min,
            Rcpp::Named("ttl_pts") = ttl_pts);

            }

// [[Rcpp::export]]
NumericMatrix rcpp_prob_matrix(NumericMatrix probs) {

            double S = probs(_,1).size();
            
            //initialize matrix of 0s
            NumericMatrix m(S,S+3);
            std::fill( m.begin(), m.end(), 0 ) ;
            
            for(int i=0; i<S; i++){
            if(i>0){
              m(i,0) = probs(i,0)*(m(i-1,0) + m(i-1, 1));
              m(i,1) = probs(i,1)*(m(i-1,0)+m(i-1,1)) + probs(i,0)*m(i-1,2);
              m(i,2) = probs(i,2)*(m(i-1, 0)+m(i-1, 1)) + probs(i,1)*m(i-1,2) + probs(i,0)*m(i-1,3);
              
              for(int k=0; k < S-1; k++){
                m(i, k+3) = probs(i,2)*m(i-1, k+2) + probs(i,1)*m(i-1, k+3) + probs(i,0)*m(i-1, k+4);
              }
            } else if(i==0){
              m(i,0) = probs(i,0);
              m(i,1) = probs(i,1);
              m(i,2) = probs(i,2);
            }
           }
            

            return m; 
                                      
            }

