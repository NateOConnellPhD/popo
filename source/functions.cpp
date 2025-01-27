#include <Rcpp.h>
#include <vector>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
List inner_sim(int N, NumericVector mu, NumericVector sd, int block_time, NumericMatrix probs) {
  int S2 = mu.size();
  int S = S2 / 2;
  
  NumericMatrix wait_times(N, S2);
  NumericMatrix service_dist_temp(N, S);

  // Extract probabilities
  NumericVector p0 = probs(_, 0);
  NumericVector p1 = probs(_, 1);
  NumericVector p2 = probs(_, 2);

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> runif(0.0, 1.0);
  std::normal_distribution<> rnorm(0.0, 1.0);

  for (int k = 0; k < N; ++k) {
    std::vector<double> runis(S);
    for (int i = 0; i < S; ++i) {
      runis[i] = runif(rng);
    }

    double carry_over = 0.0;
    for (int i = 0; i < S; ++i) {
      double r1 = runis[i];
      double apt1 = (r1 > p0[i]) ? 1.0 : 0.0;
      double apt2 = (r1 > (p0[i] + p1[i])) ? 1.0 : 0.0;

      double service1 = apt1 * std::exp(mu[i] + sd[i] * rnorm(rng));
      double service2 = apt2 * std::exp(mu[S + i] + sd[S + i] * rnorm(rng));

      carry_over = (i == 0 || service_dist_temp(k, i - 1) <= block_time) ? 0.0 : service_dist_temp(k, i - 1) - block_time;

      wait_times(k, i * 2) = carry_over;
      wait_times(k, i * 2 + 1) = (service1 + carry_over) * apt2;
      service_dist_temp(k, i) = service1 + service2 + carry_over;
    }
  }

  return List::create(Named("waits") = wait_times,
                      Named("serv") = service_dist_temp);
}

// [[Rcpp::export]]
List rcpp_est_perf(NumericMatrix wait_dist,
                   NumericMatrix service_dist,
                   NumericMatrix probs,
                   double otmin_thresh,
                   double wait_thresh,
                   double idle_thresh,
                   double block_time) {
  int N = wait_dist.nrow();
  int S = wait_dist.ncol();

  double excess_pt = 0.0, temp = 0.0, temp_size = 0.0, idle_temp = 0.0;
  NumericVector ot_mins = service_dist(_, S / 2 - 1);

  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < S; ++k) {
      if (wait_dist(i, k) > wait_thresh) {
        excess_pt += 1.0;
        temp += wait_dist(i, k);
        temp_size += 1.0;
      }
      if (k < S / 2 && service_dist(i, k) < idle_thresh) {
        idle_temp += 1.0;
      }
    }
  }

  double excess_pt_mean = excess_pt / N;
  double excess_min_avg = temp / temp_size;
  double idle_slots = idle_temp / N;

  double ot_min_sum = 0.0, ot_p = 0.0;
  for (int i = 0; i < N; ++i) {
    if (ot_mins[i] > (block_time + otmin_thresh)) {
      ot_min_sum += ot_mins[i] - block_time;
      ot_p += 1.0;
    }
  }

  double ot = ot_min_sum / N;
  double ot_cond_min = ot_p > 0 ? ot_min_sum / ot_p : 0.0;
  double ot_prob = ot_p / N;

  double ttl_pts = 0.0;
  int outer_size = probs.nrow();
  int inner_size = probs.ncol();

  for (int i = 0; i < outer_size; ++i) {
    for (int k = 0; k < inner_size; ++k) {
      ttl_pts += k * probs(i, k);
    }
  }

  return List::create(Named("excess_pt_mean") = excess_pt_mean,
                      Named("excess_pt_min") = excess_min_avg,
                      Named("idle_slots") = idle_slots,
                      Named("ot_prob") = ot_prob,
                      Named("ot_mins") = ot,
                      Named("ot_cond_mins") = ot_cond_min,
                      Named("ttl_pts") = ttl_pts);
}

// [[Rcpp::export]]
NumericMatrix rcpp_prob_matrix(NumericMatrix probs) {
  int S = probs.nrow();
  NumericMatrix m(S, S + 3);

  for (int i = 0; i < S; ++i) {
    if (i == 0) {
      m(i, 0) = probs(i, 0);
      m(i, 1) = probs(i, 1);
      m(i, 2) = probs(i, 2);
    } else {
      m(i, 0) = probs(i, 0) * (m(i - 1, 0) + m(i - 1, 1));
      m(i, 1) = probs(i, 1) * (m(i - 1, 0) + m(i - 1, 1)) + probs(i, 0) * m(i - 1, 2);
      m(i, 2) = probs(i, 2) * (m(i - 1, 0) + m(i - 1, 1)) + probs(i, 1) * m(i - 1, 2) + probs(i, 0) * m(i - 1, 3);

      for (int k = 0; k < S - 1; ++k) {
        m(i, k + 3) = probs(i, 2) * m(i - 1, k + 2) + probs(i, 1) * m(i - 1, k + 3) + probs(i, 0) * m(i - 1, k + 4);
      }
    }
  }

  return m;
}

