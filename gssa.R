############################## initialisation #################################### 
set.seed(100)

c1 = 0.02 # gene activation rate - c1
c2 = 0.009 #gene inactivation rate - c2
l1 = 0.5 # transcription rate - l1
t1 = 0.05 # RNA degradation rate - t1
l2 = 0.1 # RNA translation rate - l2
t2 = 0.005 # protein degradation rate - t2

state_change_vec = c(c1, c2, l1, t1, l2, t2)

#stoich_matrix
R1 = c(1,0,0)
R2 = c(-1,0,0)
R3 = c(0,1,0)
R4 = c(0,-1,0)
R5 = c(0,0,1)
R6 = c(0,0,-1)
stoich_m = rbind(R1,R2,R3,R4,R5,R6)

max_t = 75*60

# N simulations
rm(paths, dist, X_path)
t = 0
X = c(1,0,0)
X_path = c(t,X)
N = 700
for (k in 1:N)
{
  print(paste(k,"th loop running"))
  rm(paths, dist, X_path)
  t = 0
  X = c(1,0,0)
  X_path = c(t,X)
  while(t<=max_t)
  {
   a = c(c1*abs(X[1]-1), c2*X[1], l1*X[1], t1*X[2], l2*X[2], t2*X[3])
    a_sum = cumsum(a)
    r = runif(2,0,1)
    t = t + (1/a_sum[6])*log((1/r[1]))
    {
      j = which(a_sum > (r[2]*a_sum[6]))[1]
      X = X + stoich_m[j,]
     }
    X_path = rbind(X_path, c(t,X))
    X_path
  }
  paths = data.frame(X_path)
  X5 = ceiling(paths$X1)
  paths = cbind(paths, X5)
  dist = data.frame(seq(1:ceiling(max(paths$X1))))
  dist = cbind(dist, matrix(0,ceiling(max(paths$X1)),3))
  colnames(dist) = c("time", "n_gene", "n_rna", "n_protein")
  #sampling time points
  for (i in 1:ceiling(max(paths$X1)))
  {
    j = max(which((paths$X1 - i)<0))
    dist[i,2] = paths[j,2]
    dist[i,3] = paths[j,3]
    dist[i,4] = paths[j,4]
  }
  
  # sampling time points
  if (k==1) 
  {
    sampled_time = dist$time
    gene_lev = sampled_time
    rna_lev = sampled_time
    prot_lev = sampled_time
  }
  else  
  gene_lev = cbind(gene_lev, dist$n_gene)
  rna_lev = cbind(rna_lev, dist$n_rna)
  prot_lev = cbind(prot_lev, dist$n_protein)
  
}

gene_lev = data.frame(gene_lev)
rna_lev = data.frame(rna_lev)
prot_lev = data.frame(prot_lev)

par(mfrow = c(1,2))
path_time = data.frame(paths$X1, paths$X2)
active_run = rle(paths$X2)
active_run = data.frame(active_run$lengths, active_run$values)
colnames(active_run) = c("lengths", "values")
active_run$t = 0
ind2 = 1+active_run[1,1]
ind1 = 1
active_run[1,3] = path_time[ind2,1] - path_time[ind1,1]
for (j in 2:nrow(active_run)) 
{
  ind1 = ind2
  ind2 = ind2 + active_run[j,1]
  active_run[j,3] = path_time[ind2,1] - path_time[ind1,1]
}

time_active = active_run[active_run$values == 1,]
breaks_at = seq(0,max(time_active$t, na.rm = TRUE)+50, by = 50)
time_h = hist(time_active$t, freq = FALSE, main = "Active times", xlab = "times", col = "green", breaks = breaks_at)

par(new=T)
plot(x=time_h$mids, y= c1*exp(-c1*time_h$mids), type = "l", col = "red",xlab = "", ylab = "")

#barplot(time_h$density)

time_inactive = active_run[active_run$values == 0,]
breaks_at = seq(0,max(time_inactive$t, na.rm = TRUE)+50, by = 50)
time_h = hist(time_inactive$t, freq = FALSE, main = "Inactive times", xlab = "times", breaks = breaks_at, col = "red")
par(new=T)
plot(x=time_h$mids, y= c2*exp(-c2*time_h$mids), type = "l", xlab = "", ylab = "", col = "red")
#########################################

par(mfrow = c(1,2))
hist(rna_lev$V11, freq = TRUE, main = "", xlab = "rna copy number", col = "green", breaks = seq(0,max(rna_lev$V11)+1,1))
hist(prot_lev$V11, freq = TRUE, main = "",  xlab = "protein copy number", col = "blue", breaks = seq(0,max(prot_lev$V11)+50,50))

par(mfrow = c(2,2))
# plotting a simulation
tmax = ceiling(max(dist$time/60))
max_y = max(dist[,2:4]) + 10
#plot(1,xlim=c(1,tmax),ylim=c(0,max_y),xlab='time', ylab='number of each species')
plot(x = dist$time/60, y = dist$n_gene, xlim=c(0,tmax),ylim=c(0,max_y), type = "l", xlab = "time", ylab = "number of each species", main = "single simulation", col = "red")
lines(x = dist$time/60, y = dist$n_rna, col = "green")
lines(x = dist$time/60, y = dist$n_protein, col = "blue")
legend(x=0,y=max_y,c("gene", "rna", "protein"),cex=.5,col=c("red", "green","blue"),pch=c("_","_","_"), y.intersp = 0.2, bty = "n")

# plotting protein levels
tmax = ceiling(max(prot_lev$prot_lev/60))
max_y = max(prot_lev[,2:ncol(prot_lev)]) + 10
plot(x = prot_lev$prot_lev/60, 
     y = prot_lev$V2, 
     type = "l", 
     xlim=c(0,tmax),
     ylim=c(0,max_y),
     xlab = "time",
     ylab='protein copy number',
     col = "azure3",
     main = paste(N,"simulations"),
     cex = 0.2)
for (i in 3:ncol(prot_lev)) 
{
  #colum = paste("prot_lev$V",i, sep = "")
  lines(x = prot_lev$prot_lev/60, y = prot_lev[,i], type = "l", col = "azure3")
}

# protein and rna correlation
max_y = max(prot_lev[,2:ncol(prot_lev)])
tmax = max(rna_lev[,2:ncol(rna_lev)])
#plot(1,xlim=c(1,tmax),ylim=c(0,max_y),xlab='time', ylab='protein copy number')
plot(x = rna_lev$V2,
         y = prot_lev$V2,
         type = "p",
         xlim = c(0,tmax),
         ylim = c(0,max_y),
         xlab = "rna copy number",
         ylab ='protein copy number',
         col = "azure3",
         main = "rna vs protein",
         cex = 0.2)
#for (i in 3:ncol(rna_lev))
#{
# lines(x = rna_lev[,i], y = prot_lev[,i], type = "p", col = "azure3")
#}

# protein and rna correlation
# mean_prot = apply(prot_lev[,2:ncol(prot_lev)], 1, mean)
# mean_rna = apply(rna_lev[,2:ncol(rna_lev)], 1, mean)
# max_y = max(mean_prot) + 5
# tmax = max(mean_rna) 
# plot(x = mean_rna, 
#      y = mean_prot, 
#      type = "p", 
#      xlim = c(1,tmax),
#      ylim = c(0,max_y),
#      xlab = "rna copy number",
#      ylab ='protein copy number',
#      col = "azure3",
#      main = "rna vs protein",
#      cex = 0.2)

cor_rna_prot = diag(as.matrix(cor(x = prot_lev[,2:ncol(prot_lev)], y = rna_lev[,2:ncol(rna_lev)])))
plot(cor_rna_prot, xlab = "simulations", ylab = "r", main = "correlation", col = "azure3")

##############################################################################
# ODE Solutions
library(deSolve)
## Chaos in the atmosphere
## def two_state_promoter_moments(m, t, c):
model <- function(t, m, c) {
  with(as.list(c(m, c)), {
    # dX <- a * X + Y * Z
    # dY <- b * (Y - Z)
    # dZ <- -X * Y + c * Y - Z
    dm = seq(0,0, length.out = 14)
    
    c1 = c[1]
    c2 = c[2]
    c3 = c[3]
    c4 = c[4]
    c5 = c[5]
    c6 = c[6]
    
    x_1 = m[1]
    x_2 = m[2]
    x_3 = m[3]
    x_4 = m[4]
    x_11 = m[5]
    x_12 = m[6]
    x_13 = m[7]
    x_14 = m[8]
    x_22 = m[9]
    x_23 = m[10]
    x_24 = m[11]
    x_33 = m[12]
    x_34 = m[13]
    x_44 = m[14]
    
    dm[1] = x_2*c2 - x_1*c1
    dm[2] = x_1*c1 - x_2*c2
    dm[3] = x_2*c3 - x_3*c4
    dm[4] = x_3*c5 - x_4*c6
    dm[5] = x_1*c1 + x_2*c2 - 2*x_11*c1 + 2*x_12*c2
    dm[6] = x_11*c1 - x_2*c2 - x_1*c1 - x_12*c1 - x_12*c2 + x_22*c2
    dm[7] = x_12*c3 - x_13*c1 - x_13*c4 + x_23*c2
    dm[8] = x_13*c5 - x_14*c1 - x_14*c6 + x_24*c2
    dm[9] = x_1*c1 + x_2*c2 + 2*x_12*c1 - 2*x_22*c2
    dm[10] = x_13*c1 + x_22*c3 - x_23*c2 - x_23*c4
    dm[11] = x_14*c1 - x_24*c2 + x_23*c5 - x_24*c6
    dm[12] = x_2*c3 + x_3*c4 + 2*x_23*c3 - 2*x_33*c4
    dm[13] = x_24*c3 + x_33*c5 - x_34*c4 - x_34*c6
    dm[14] = x_3*c5 + x_4*c6 + 2*x_34*c5 - 2*x_44*c6
    list(dm)
  })
}

#inputs
c = c(c1, c2, l1, t1, l2, t2) #proposed state
m = seq(0,0, length.out = 14) 
m[1] = 1
m[5] = 1
#times <- seq(0, 100, by = 0.01)
times <- sampled_time
out <- ode(y = m, times = times, func = model, parms = c)
out = data.frame(out)
moments_ode = cbind(out$X4, out$X14) # MEAN AND VARIANCE OF PROTEIN OPY NUMBER
max_y = max(prot_lev[,2:ncol(prot_lev)]) + 10
tmax = ceiling(max(prot_lev$prot_lev/60))
#plot(1,xlim=c(1,tmax),ylim=c(0,max_y),xlab='t', ylab='protein copy number')
par(mfrow = c(1,1))
plot(x = prot_lev$prot_lev/60, 
     y = prot_lev$V2, 
     type = "l", 
     xlim=c(0,tmax),
     ylim=c(0,max_y),
     xlab = "t",
     ylab='protein copy number',
     col = "azure3")

for (i in 3:ncol(prot_lev)) 
{
  #colum = paste("prot_lev$V",i, sep = "")
  lines(x = prot_lev$prot_lev/60, y = prot_lev[,i], type = "l", col = "azure3")
}

lines(x = sampled_time/60, y = moments_ode[,1], type = "l", col = "blue")
var = moments_ode[,2]-moments_ode[,1]^2
#sd_mom1 = moments_ode[,2]+moments_ode[,1]^2
lines(x = sampled_time/60, y = moments_ode[,1]-sqrt(var), type = "l", col = "red")
lines(x = sampled_time/60, y = moments_ode[,1]+sqrt(var), type = "l", col = "red")

########################################################################################
#data = read.csv("/media/rachel/Part2/CSBD/students_folder/data/YFPDCS2.csv", header = FALSE)

data = prot_lev
times_exp = data$prot_lev
mean_bootstrap = data$prot_lev
sd_bootstrap = data$prot_lev
sec_ord_boot = data$prot_lev
N_boot = 1000
for (i in 1:N_boot) 
{
  index_bs = floor(runif(100, 2, ncol(data)))
  {
    samples = as.matrix(data[,index_bs])  
    mean_t = apply(samples, 1, mean)
    sd_t = apply(samples, 1, sd)
    sec_order = apply(samples^2, 1, mean)
  }
  mean_bootstrap = cbind(mean_bootstrap, mean_t)
  sd_bootstrap = cbind(sd_bootstrap, sd_t)
  sec_ord_boot = cbind(sec_ord_boot, sec_order)
}

mean_bootstrap = data.frame(mean_bootstrap) # mean
mean_bs = apply(mean_bootstrap[,2:ncol(mean_bootstrap)], 1, mean) # mean of mean bootstraps
sd_mean = apply(mean_bootstrap[,2:ncol(mean_bootstrap)], 1, sd) # sd of mean

sd_bootstrap = data.frame(sd_bootstrap) # sd
sd_bs = apply(sd_bootstrap[,2:ncol(sd_bootstrap)], 1, mean) # mean of sd bootstraps
sd_sd = apply(sd_bootstrap[,2:ncol(sd_bootstrap)], 1, sd) # sd of sd

sec_ord_boot = data.frame(sec_ord_boot)
sec_mean_bs = apply(sec_ord_boot[,2:ncol(sec_ord_boot)], 1, mean) # mean of second order means
sec_sd_bs = apply(sec_ord_boot[,2:ncol(sec_ord_boot)], 1, sd) # mean of second order means

mm = max(mean_bootstrap[,2:ncol(mean_bootstrap)])
ms = max(sd_bootstrap[,2:ncol(sd_bootstrap)])

pdf("bootstraps_4.pdf") 
par(mfrow = c(2,1))
########################### plot of bootstrap values for mean #####################
plot(x = mean_bootstrap$mean_bootstrap/60, 
     y = mean_bootstrap$mean_t, 
     type = "l", 
     xlim = c(0,max(mean_bootstrap$mean_bootstrap/60)),
     ylim = c(0,mm),
     xlab = "t",
     ylab = "mean")
for (i in 3:ncol(mean_bootstrap)) 
{
  colum = paste("mean_bootstrap$mean_t.",i, sep = "")
  lines(x = mean_bootstrap$mean_bootstrap/60, y = mean_bootstrap[,i], type = "l", col = "azure3")
}
lines(x = mean_bootstrap$mean_bootstrap/60, y = mean_bs, type = "l", col = "red")
lines(x = mean_bootstrap$mean_bootstrap/60, y = mean_bs - sd_mean, type = "l", col = "blue")
lines(x = mean_bootstrap$mean_bootstrap/60, y = mean_bs + sd_mean, type = "l", col = "blue")

########################### plot of bootstrap values for standard deviation ##################### 
plot(x = sd_bootstrap$sd_bootstrap/60, 
     y = sd_bootstrap$sd_t, 
     type = "l", 
     xlim = c(0,max(sd_bootstrap$sd_bootstrap/60)),
     ylim = c(0,ms),
     xlab = "t",
     ylab = "sd")
for (i in 3:ncol(sd_bootstrap)) 
{
  colum = paste("sd_bootstrap$sd_t.",i, sep = "")
  lines(x = sd_bootstrap$sd_bootstrap/60, y = sd_bootstrap[,i], type = "l", col = "azure3")
}
lines(x = sd_bootstrap$sd_bootstrap/60, y = sd_bs, type = "l", col = "red")
lines(x = sd_bootstrap$sd_bootstrap/60, y = sd_bs - sd_sd, type = "l", col = "blue")
lines(x = sd_bootstrap$sd_bootstrap/60, y = sd_bs + sd_sd, type = "l", col = "blue")
dev.off()

data_times = mean_bootstrap$mean_bootstrap
data_m = data.frame(cbind(data_times, mean_bs, sd_mean, sd_bs, sd_sd))

# var = moments_ode[,2]-moments_ode[,1]^2
# moments_ode_mh = data.frame(out$time, out$X4, var)
# colnames(moments_ode_mh) = c("time", "mean", "var")
 
library(deSolve)
#moments <- data.frame(ode(y = m, times = times, func = model, parms = c))
#data_m = data.frame(cbind(data_times, mean_bs, sd_mean, sd_bs, sd_sd))
data_m = data.frame(cbind(data_times, mean_bs, sd_mean, sd_bs, sd_sd, sec_mean_bs,sec_sd_bs))
likl <- function(data, moments){
  ll=0
  sd_data = moments$X14 - moments$X4^2
  if (sum(c(data_m$sd_mean, sd_data)>0)>1)
  {
    #ll =  -1*sum((((data_m$mean_bs - moments$X4)^2)/(data_m$sd_mean^2) + ((data_m$sd_bs - sd_data)^2)/(data_m$sd_sd^2)), na.rm = TRUE)
    ll =  -1*sum(((((data_m$mean_bs - moments$X4)^2)/data_m$sd_mean^2) + ((data_m$sec_mean_bs - moments$X14)^2)/data_m$sec_sd_bs^2), na.rm = TRUE)
    
  }
  else{}
    return(ll)   
}
#l = likl(data = data_m, moments = moments)

p_sd = 0.01 # proposed sd
current_par = c
# Prior distribution

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

mh_iter = 10
acc_counter = 0
likelihood_rec = matrix(0,mh_iter,1) #l_rec
param_rec = matrix(0,mh_iter, length(c)) #chain
param_rec[1,] = 0.05
#temp_par = param_rec
moments <- data.frame(ode(y = m, times = times, func = model, parms = param_rec[1,]))
l_old = likl(data = data_m, moments = moments)
likelihood_rec[1] = l_old
for (it in 2:mh_iter) 
{
  #if(it%%10==0){print(c(it, l_new, plog_back, l_old, plog_for, acc_ratio))}
  #print(it)
  current_par = param_rec[it-1,]
  print(current_par)
  proposed_par = prior(current_par)
  print(proposed_par)
  mlog = log(current_par) - 0.5 * log(1 + p_sd/(current_par^2))
  print(mlog)
  slog = sqrt(log(1 + p_sd^2/(current_par^2)))
  print(slog)
  plog_back = sum(dlnorm(proposed_par, meanlog = mlog, sdlog = slog, log = T))
  mlog = log(proposed_par) - 0.5 * log(1 + p_sd/(proposed_par^2))
  slog = sqrt(log(1 + p_sd^2/(proposed_par^2)))
  plog_for = sum(dlnorm(current_par, meanlog = mlog, sdlog = slog, log = T))
  temp_par = proposed_par
  moments <- data.frame(ode(y = m, times = times, func = model, parms = temp_par))
  l_new = likl(data = data_m, moments = moments)
  likelihood_rec[it] = l_new
  acc_ratio = exp(l_new + plog_back - l_old - plog_for)
  if(it%%10==0){print(c(it, l_new, plog_back, l_old, plog_for, acc_ratio))}
  alpha = min(1,acc_ratio)
  if((alpha>=runif(1,0,1)))
  {
    param_rec[it,] = proposed_par
    l_old = l_new
    print(c(it, l_new, plog_back, l_old, plog_for, acc_ratio))
    acc_counter = acc_counter + 1 
  } else
  {
    param_rec[it,] = param_rec[it-1,]
  }
}

iter = seq(1,mh_iter,1)
plot(x = iter, y = param_rec[,1], type = "l")
lines(x = iter, y = param_rec[,2], type = "l")
lines(x = iter, y = param_rec[,3], type = "l")
lines(x = iter, y = param_rec[,4], type = "l")
lines(x = iter, y = param_rec[,5], type = "l")
lines(x = iter, y = param_rec[,6], type = "l")
tail(param_rec)
c
tail(likelihood_rec)
head(likelihood_rec)
plot(x = iter, y = likelihood_rec)
