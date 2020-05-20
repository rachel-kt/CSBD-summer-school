library(deSolve)
moments <- data.frame(ode(y = m, times = times, func = model, parms = c))
data_m = data.frame(cbind(data_times, mean_bs, sd_mean, sd_bs, sd_sd))
likl <- function(data, moments){
ll=0
sd_data = moments$X14 - moments$X4^2
   if (sum(c(data_m$sd_mean, sd_data)>0)>1)
      {
        ll =  -1*sum(((data_m$mean_bs - moments$X1)^2)/(data_m$sd_mean^2) + ((data_m$sd_bs - sd_data)^2)/(data_m$sd_sd^2), na.rm = TRUE)
   }
else
  return(ll)   
}
l = likl(data = data_m, moments = moments)

p_sd = 0.05 # proposed sd
current_par = c
# Prior distributionl

prior <- function(param){
  c1 = param[1]
  c2 = param[2]
  l1 = param[3]
  t1 = param[4]
  l2 = param[5]
  t2 = param[6]
  #meanlog <- log(6) - 0.5 * log(1 + 9/(6^2))
  #sdlog <- sqrt(log(1 + 9/(6^2)))
  meanlog <- log(param) - 0.5 * log(1 + p_sd/(param^2))
  sdlog <- sqrt(log(1 + p_sd^2/(param^2)))
  # rlnx2 <- rlnorm(1000, meanlog = meanlog, sdlog = sdlog)
  # hist(rlnx2)
  # mean(rlnx2)
  # sd(rlnx2)
  c1prior = rlnorm(1, mean = meanlog[1], sd = sdlog[1])
  c2prior = rlnorm(1, mean = meanlog[2], sd = sdlog[2])
  l1prior = rlnorm(1, mean = meanlog[3], sd = sdlog[3])
  t1prior = rlnorm(1, mean = meanlog[4], sd = sdlog[4])
  l2prior = rlnorm(1, mean = meanlog[5], sd = sdlog[5])
  t2prior = rlnorm(1, mean = meanlog[6], sd = sdlog[6])
  return(c(c1prior,c2prior,l1prior,t1prior,l2prior,t2prior))
}

mh_iter = 10000
acc_counter = 0
likelihood_rec = matrix(0,mh_iter,1) #l_rec
param_rec = matrix(0,mh_iter, length(c)) #chain
param_rec[1,] = 0.05
temp_par = param_rec
moments <- data.frame(ode(y = m, times = times, func = model, parms = c))
l_old = likl(data = data_m, moments = moments)
likelihood_rec[1] = l_old
for (it in 2:mh_iter) 
{
  if(it%%10==0){print(it)}
  #print(it)
  current_par = param_rec[it-1,]
  proposed_par = prior(current_par)
  plog_back = sum((dnorm(proposed_par, mean = current_par, sd = p_sd, log = T)))
  plog_for = sum((dnorm(current_par, mean = proposed_par, sd = p_sd, log = T)))
  temp_par = proposed_par
  moments <- data.frame(ode(y = m, times = times, func = model, parms = temp_par))
  l_new = likl(data = data_m, moments = moments)
  likelihood_rec[it] = l_new
  acc_ratio = exp(l_new + (plog_back) - (l_old) - plog_for)
  alpha = min(1,acc_ratio)
  if((alpha>=runif(1,0,1)))
  {
    param_rec[it,] = proposed_par
    l_old = l_new
    acc_counter = acc_counter + 1 
  } else
  {
    param_rec[it,] = param_rec[it-1,]
  }
}

iter = seq(1,mh_iter,1)
plot(x = iter, y = param_rec[,1], type = "l")
#abline()
lines(x = iter, y = param_rec[,2], type = "l")
lines(x = iter, y = param_rec[,3], type = "l")
lines(x = iter, y = param_rec[,4], type = "l")
lines(x = iter, y = param_rec[,5], type = "l")
lines(x = iter, y = param_rec[,6], type = "l")
tail(param_rec)
c
tail(likelihood_rec)
head(likelihood_rec)
