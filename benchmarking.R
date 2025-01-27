library(Rcpp)
sourceCpp("source/bench_rcpp.cpp")
inner_sim_bench(100, st_mu, st_sd, block.time, prob_matrix)
inner_sim_times

pred
microbenchmark::microbenchmark(
  simple = {mean(pred)},
  better = {sum(pred)/length(pred)},
  times=1000
)
mean(pred)
sum(pred)/length(pred)
