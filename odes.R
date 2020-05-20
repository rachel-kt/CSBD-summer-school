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
# gene activation rate - c1
c1 = 0.02
#gene inactivation rate - c2
c2 = 0.5
# transcription rate - l1
l1 = 0.5
# RNA degradation rate - t1
t1 = 0.05
# RNA translation rate - l2
l2 = 0.1
# protein degradation rate - t2
t2 = 0.05
c = c(c1, c2, l1, t1, l2, t2)
m = seq(0,0, length.out = 14) 
m[1] = 1
m[5] = 1
#times <- seq(0, 100, by = 0.01)
times <- times_exp
out <- ode(y = m, times = times, func = model, parms = c)
plot(out)
## add a 3D figure if package scatterplot3D is available
if (require(scatterplot3d))
  scatterplot3d(out[,-1], type = "l")
moments_ode = data.frame(out[,1], out[,5], sqrt(out[,15]-out[,5]^2))
colnames(moments_ode) = c("time", "mean_protein", "sd")
