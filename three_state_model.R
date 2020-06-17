############################## initialisation #################################### 
set.seed(14)
# gene_off_1 activation rate - c1
c1 = 0.02
# gene_off_2 activation rate - c2
c2 = 0.03
# gene_on activation rate - c3
c3 = 0.08
# transcription rate - l1
l1 = 0.5
# RNA degradation rate - t1
t1 = 0.05
# RNA translation rate - l2
l2 = 0.6
# protein degradation rate - t2
t2 = 0.005

state_change_vec = c(c1, c2, c3, l1, t1, l2, t2)
#stoich_metric
R1 = c(-1,1,0,0,0)
R2 = c(0,-1,1,0,0)
R3 = c(1,0,-1,0,0)
R4 = c(0,0,0,1,0)
R5 = c(0,0,0,-1,0)
R6 = c(0,0,0,0,1)
R7 = c(0,0,0,0,-1)

stoich_m = rbind(R1,R2,R3,R4,R5,R6,R7)

max_t = 1000*60

# N simulations
t = 0
X = c(0,0,1,0,0)
X_path = c(t,X)
#for (k in 2:100)
{
  t = 0
  X = c(0,0,1,0,0)
  X_path = c(t,X)
  while(t<=max_t)
  {
    a = c(c1*X[1], c2*X[2],c3*X[3], l1*X[3], t1*X[4], l2*X[4], t2*X[5])
    a_sum = cumsum(a)
    r = runif(2,0,1)
    t = t + (1/a_sum[7])*log((1/r[1]))
    {
      j = which(a_sum > (r[2]*a_sum[7]))[1]
      X = X + stoich_m[j,]
    }
    X_path = rbind(X_path, c(t,X))
    X_path
  }
  paths = data.frame(X_path)
  #paths = paths[1:20000,]
  X5 = ceiling(paths$X1)
  paths = cbind(paths, X5)
  dist = data.frame(seq(1:ceiling(max(paths$X1))))
  dist = cbind(dist, matrix(0,ceiling(max(paths$X1)),5))
  colnames(dist) = c("time", "gene_off_1", "gene_off_2", "gene_on", "n_rna", "n_protein")
  #sampling time points
  for (i in 1:ceiling(max(paths$X1)))
  {
    j = max(which((paths$X1 - i)<0))
    dist[i,2] = paths[j,2]
    dist[i,3] = paths[j,3]
    dist[i,4] = paths[j,4]
    dist[i,5] = paths[j,5]
    dist[i,6] = paths[j,6]
  }
  tmax = ceiling(max(dist$time))
  max_y = max(dist[,2:6]) + 10
  plot(1,xlim=c(1,tmax),ylim=c(0,max_y),xlab='t', ylab='n')
  plot(x = dist$time, y = dist$gene_on, xlim=c(1,tmax),ylim=c(0,max_y), type = "l", xlab = "t", ylab = "n", col = "red")
  lines(x = dist$time, y = dist$n_rna, col = "green")
  lines(x = dist$time, y = dist$n_protein, col = "blue")
  legend(x=0,y=max_y,c("gene", "rna", "protein"),cex=.8,col=c("red", "green","blue"),pch=c("_","_","_"))
  
  # sampling time points
  if (k==1) 
  {
    sampled_time = dist$time
    gene_lev = sampled_time
    rna_lev = sampled_time
    prot_lev = sampled_time
  }
  gene_lev = cbind(gene_lev, dist$n_gene)
  rna_lev = cbind(rna_lev, dist$n_rna)
  prot_lev = cbind(prot_lev, dist$n_protein)
  rm(paths, dist, X_path)
}
# gene_lev = data.frame(gene_lev)
# rna_lev = data.frame(rna_lev)
# prot_lev = data.frame(prot_lev)
# 
# # plotting protein levels
# tmax = ceiling(max(prot_lev$prot_lev))
# max_y = max(prot_lev[,2:ncol(prot_lev)]) + 10
# 
# plot(1,xlim=c(1,tmax),ylim=c(0,max_y),xlab='t', ylab='protein copy number')
# plot(x = prot_lev$prot_lev, 
#      y = prot_lev$V2, 
#      type = "l", 
#      xlim=c(1,tmax),
#      ylim=c(0,max_y),
#      xlab = "t",
#      ylab='protein copy number',
#      col = "azure3")
# for (i in 3:ncol(prot_lev)) 
# {
#   #colum = paste("prot_lev$V",i, sep = "")
#   lines(x = prot_lev$prot_lev, y = prot_lev[,i], type = "l", col = "azure3")
# }

#active_run = rle(paths$X2)
active_run = rle(dist$gene_on)
active_run = data.frame(active_run$lengths, active_run$values)
active_run = active_run[active_run[,2]==1,]
colnames(active_run) = c("lengths", "values")

active_run = rle(paths$X4)
active_run = data.frame(active_run$lengths, active_run$values)
active_run = active_run[active_run[,2]==1,]
colnames(active_run) = c("lengths", "values")
#ac_ina_run = rle(paths$X2)
#ac_ina_run = data.frame(ac_ina_run$lengths, ac_ina_run$values)
#colnames(ac_ina_run) = c("lengths", "values")
for (j in nrow(ac_ina_run))
{
  #real_time = paths[ac_ina_run[j]+]
}

#on.time = (1, sum(paths[1:(active_run[1,1]+1)]

range(active_run$lengths)
breaks = seq(min(range(active_run$lengths)),max(range(active_run$lengths)), by=1) 
active_run.cut = cut(active_run$lengths, breaks, right=FALSE)
active_run.freq = table(active_run.cut)
barplot_val = data.frame(ac_tab)
x_bar_height = barplot_val$Freq
x_bar = seq(1,nrow(x_barplot),1)
barplot(table(active_run.cut), ylim = c(0,max(active_run.freq)+10), col = "red", cex.names = 0.5, xlab = "active time", ylab = "frequency")
barplot(height = x_bar_height, ylim = c(0,max(active_run.freq)+10), col = "red", cex.names = 0.5, xlab = "active time", ylab = "frequency")
hist(x = x_bar_height, col = "red", cex.names = 0.5, xlab = "active time", ylab = "frequency", freq = FALSE)

y_bar = (c1*c2/(c2-c1))*(exp(-c1*x_bar)-exp(-c2*x_bar))
lines(x = x_bar, y = y_bar, col = "green")
rm(paths, dist, X_path)


#legend(x=0,y=40,c("gene act", "rna", "protein"),cex=.8,col=c("red", "green","blue"),pch=c("_","_","_"))
par(new=TRUE)

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
moments_ode = cbind(out$X4, out$X14)
plot(1,xlim=c(1,tmax),ylim=c(0,max_y),xlab='t', ylab='protein copy number')
plot(x = prot_lev$prot_lev, 
     y = prot_lev$V2, 
     type = "l", 
     xlim=c(1,tmax),
     ylim=c(0,max_y),
     xlab = "t",
     ylab='protein copy number',
     col = "azure3")
for (i in 3:ncol(mean_bootstrap)) 
{
  #colum = paste("prot_lev$V",i, sep = "")
  lines(x = prot_lev$prot_lev, y = prot_lev[,i], type = "l", col = "azure3")
}

lines(x = sampled_time, y = moments_ode[,1], type = "l", col = "blue")
sd_mom = moments_ode[,2]-moments_ode[,1]*moments_ode[,1]
sd_mom1 = moments_ode[,2]+moments_ode[,1]*moments_ode[,1]
lines(x = sampled_time, y = sqrt(sd_mom), type = "l", col = "red")
lines(x = sampled_time, y = sqrt(sd_mom1), type = "l", col = "red")
########################################################################################
#data = read.csv("/media/rachel/Part2/CSBD/students_folder/data/YFPDCS2.csv", header = FALSE)
data = prot_lev
times_exp = data$prot_lev
mean_bootstrap = data$prot_lev
sd_bootstrap = data$prot_lev
for (i in 1:1000) 
{
  index_bs = floor(runif(100, 1, ncol(data)))
  {
    samples = as.matrix(data[,index_bs])  
    mean_t = apply(samples, 1, mean)
    sd_t = apply(samples, 1, sd)
  }
  mean_bootstrap = cbind(mean_bootstrap, mean_t)
  sd_bootstrap = cbind(sd_bootstrap, sd_t)
}
mean_bootstrap = data.frame(mean_bootstrap) # mean
mean_bs = apply(mean_bootstrap[,2:ncol(mean_bootstrap)], 1, mean) # mean of mean bootstraps
sd_mean = apply(mean_bootstrap[,2:ncol(mean_bootstrap)], 1, sd) # sd of mean

sd_bootstrap = data.frame(sd_bootstrap) # sd
sd_bs = apply(sd_bootstrap[,2:ncol(sd_bootstrap)], 1, mean) # mean of sd bootstraps
sd_sd = apply(sd_bootstrap[,2:ncol(sd_bootstrap)], 1, sd) # sd of sd

mm = max(mean_bootstrap[,2:ncol(mean_bootstrap)])
ms = max(sd_bootstrap[,2:ncol(sd_bootstrap)])

########################### plot of bootstrap values for mean #####################
plot(x = mean_bootstrap$mean_bootstrap, 
     y = mean_bootstrap$mean_t, 
     type = "l", 
     xlim = c(0,max(mean_bootstrap$mean_bootstrap)),
     ylim = c(0,mm),
     xlab = "t",
     ylab = "mean")
for (i in 3:ncol(mean_bootstrap)) 
{
  colum = paste("mean_bootstrap$mean_t.",i, sep = "")
  lines(x = mean_bootstrap$mean_bootstrap, y = mean_bootstrap[,i], type = "l", col = "azure3")
}
lines(x = mean_bootstrap$mean_bootstrap, y = mean_bs, type = "l", col = "red")
lines(x = mean_bootstrap$mean_bootstrap, y = mean_bs - sd_mean, type = "l", col = "blue")
lines(x = mean_bootstrap$mean_bootstrap, y = mean_bs + sd_mean, type = "l", col = "blue")

########################### plot of bootstrap values for standard deviation ##################### 
plot(x = sd_bootstrap$sd_bootstrap, 
     y = sd_bootstrap$sd_t, 
     type = "l", 
     xlim = c(0,max(sd_bootstrap$sd_bootstrap)),
     ylim = c(0,ms),
     xlab = "t",
     ylab = "sd")
for (i in 3:ncol(sd_bootstrap)) 
{
  colum = paste("sd_bootstrap$sd_t.",i, sep = "")
  lines(x = sd_bootstrap$sd_bootstrap, y = sd_bootstrap[,i], type = "l", col = "azure3")
}
lines(x = sd_bootstrap$sd_bootstrap, y = sd_bs, type = "l", col = "red")
lines(x = sd_bootstrap$sd_bootstrap, y = sd_bs - sd_sd, type = "l", col = "blue")
lines(x = sd_bootstrap$sd_bootstrap, y = sd_bs + sd_sd, type = "l", col = "blue")



uncertainties = times_exp
uncertainties = cbind(times_exp, mean_bs, sd_mean, sd_bs, sd_sd)
