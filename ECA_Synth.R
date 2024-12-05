###################################################
####       ECA Computation & 3D Plotting      #####
### Giovanni Motta, Haoming Zhao, Josh Friedman ###
####              June 2023                    ####
###################################################
#### Load the Data ####
rm(list=ls())

# Load Libraries
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggnewscale)

# For loading Synthetic Dataset
load("Synth_Data_ECA.rda")
Face_PL <- Face_Synth
times=read.csv("Synth_Times.csv")
pre=read.csv("Synth_Pre.csv")

# For defining parameters
N=dim(Face_PL)[1]
P=dim(Face_PL)[2]
T0=dim(Face_PL)[3]
faces <- Face_PL 
faces <- faces + 1

#### Synchronizing Uneven Time Series ####
### Define Periods
t2=min(times$B - times$A)
t2max=max(times$B - times$A)
t3=min(times$C - times$B)
t3max=max(times$C - times$B)

### Defining Synch Parameters
T <- t2+t3

C2=array(0,c(N,P,t2))
C3=array(0,c(N,P,t3))

h_c=c(0.3, 0.14)
u=c(1:(T-1))/(T-1)
u2=seq(1/t2,1,length.out=t2)
u3=seq(1/t3,1,length.out=t3)

#The number of factors
r=3

# Rescaling / Synchronizing based on Period Times 

for (n in (1:N)){	
  t2a = times[n,2]
  t2b = times[n,3]
  t2ba = t2b-t2a
  t=1:t2ba/t2ba
  #t=1:tmax[n]
  k=matrix(0,t2ba,t2)
  for (p in (1:P)){	
    y=faces[n,p,(t2a+1):t2b]
    for (j in (1:t2)){	
      #m=dnorm((t-u[j])/h)/h
      m=as.numeric(abs((t-u2[j])/h_c[2])<=1)*(rep(1,(t2b-t2a))-((t-u2[j])/h_c[2])^2)
      k[,j]=m/sum(m)}
    C2[n,p,]=y%*%k
  }}

for (n in (1:N)){	
  t3c = times[n,4]
  t3b = times[n,3]
  t3cb = t3c-t3b
  t=1:t3cb/t3cb
  #t=1:tmax[n]
  k=matrix(0,t3cb,t3)
  for (p in (1:P)){	
    y=faces[n,p,(t3b+1):t3c]
    for (j in (1:t3)){	
      #m=dnorm((t-u[j])/h)/h
      m=as.numeric(abs((t-u3[j])/h_c[2])<=1)*(rep(1,(t3c-t3b))- ((t-u3[j])/h_c[2])^2)
      k[,j]=m/sum(m)}
    C3[n,p,]=y%*%k
  }}

library(abind)
C <- abind(C2, C3, along = 3)


## Evolutionary Correspondence Analysis

DN=array(0,c(N,N,T))
DP=array(0,c(P,P,T))
Pi =array(0,c(P,N,T))
Psi=array(0,c(P,N,T))
F=Psi
V=array(0,c(r+1,r+1,T))
hat.phi = array(0,c(P,r+1,T))
hat.lam = array(0,c(N,r+1,T))

#### Check for and replace 0's in dataset with .00001
for (t in 1:T){
  F[,,t]=C[,,t]/sum(C[,,t])
  
  for (p in 1:P){
    DP[p,p,t]=sum(F[p,,t])}
  
  for (n in 1:N){
    DN[n,n,t]=sum(F[,n,t])}}

for (t in 1:T){
  Pi[,,t]=solve(DP[,,t])%*%F[,,t]}

DN2=DN
for (t in 1:T){
  for (n in 1:N){
    if (DN[n,n,t]<0.0001){DN2[n,n,t]=0.0001}
  }}

DN=DN2
for (t in 1:T){
  Psi[,,t]=F[,,t]%*%solve(DN[,,t])
}

######################################
# Extracting loadings and factors 
######################################

V2=array(0,c(P,P,T))

for (t in 1:T){
  Y=solve(sqrt(DP[,,t]))%*%F[,,t]%*%solve(sqrt(DN[,,t]))
  X=Y-sqrt(DP[,,t])%*%matrix(1,P,1)%*%matrix(1,1,N)%*%sqrt(DN[,,t])
  X=Y
  S=t(X)%*%X
  W=eigen(S)$vectors[,1:(r+1)]
  V[,,t]=diag(eigen(S)$values[1:(r+1)])
  V2[,,t]=diag(eigen(S)$values[1:P])
  hat.phi[,,t]=solve(DP[,,t])%*%F[,,t]%*%solve(sqrt(DN[,,t]))%*%W
  hat.lam[,,t]=solve(sqrt(DN[,,t]))%*%W%*%sqrt(V[,,t])
}

tilde.phi=hat.phi
for (j in 1:(r+1)){
  for (t in 1:(T-1)){
    a=matrix(tilde.phi[,j,t+1]-tilde.phi[,j,t],P,1)
    b=matrix(tilde.phi[,j,t+1]+tilde.phi[,j,t],P,1)
    if (norm(a,"F")> norm(b,"F")){tilde.phi[,j,t+1] = - tilde.phi[,j,t+1]}
  }}

######################################
# Trace ratio
######################################


trace.ratio=matrix(0,P-1,T)
for (j in 2:P){
  for (t in 1:T){
    trace.ratio[j-1,t]=sum(V2[2:j,2:j,t])/sum(V2[2:P,2:P,t])
  }}		

x11()
plot(1:T/T,trace.ratio[1,],type="l",ylim=c(0.1,1),xlab="",ylab="",lwd=2)
for (j in 1:(P-2)){
  lines(1:T/T,trace.ratio[j+1,],col=j+1)
}



# Smoothing and plotting the loadings ----

adj.lam=hat.lam

for (j in 1:(r+1)){
  for (t in 1:(T-1)){
    a=matrix(adj.lam[,j,t+1]-adj.lam[,j,t],N,1)
    b=matrix(adj.lam[,j,t+1]+adj.lam[,j,t],N,1)
    if (norm(a,"F")> norm(b,"F")){adj.lam[,j,t+1] = - adj.lam[,j,t+1]}
  }}

tilde.lam=array(0,c(T,N,r+1))
for (i in 1:N){
  tilde.lam[,i,]=t(adj.lam[i,,])}


m = 59
w=(m + 1 - 0:m)/sum(1:(m+1))
w=matrix(w,1,m+1)

for (t in (m+1):T){
  for (j in 1:(r+1)){
    for (i in 1:N){
      y=matrix(adj.lam[i,j,t:(t-m)],m+1,1)
      tilde.lam[t,i,j]=w%*%y}}}

# get condition information
pre <- pre %>%
  mutate(Condition = factor(Condition),
         ID = factor(ID))

# Naming Variables - Dyadic & Ind Datasets (D/-41/82)
AUN <- data.frame(
  AU = as.factor(colnames(Face_PL)),
  FAU = as.factor(c("Inner Brow Raiser","R Outer Brow Raiser","Brow Lowerer",
                    "Upper Lid Raiser","Cheek Raiser","Lid Tightener","Nose Wrinkler",
                    "Upper Lip Raiser","Lip Corner Puller","Dimpler","Lip Corner Depressor",
                    "Chin Raiser","Lip Stretcher","Lip Tightener","Lips Part","Jaw Drop",
                    "Lip Suck","Blink","Gaze X","Gaze Y","Pitch","Roll","Yaw")), 
  Face_Part = as.factor(c("Brow","Brow","Brow","Eye","Mouth","Eye","Nose","Nose","Mouth",
                          "Mouth","Mouth","Jaw","Mouth","Mouth","Mouth","Jaw","Mouth","Eye",
                          "Eye","Eye","Head","Head","Head")))

# Select Row/Col profiles of interest 
rows_N <- rownames(Face_PL)[c(11,27,28,35)]     # max/min Var/SD(F1) ***

rows_P <- c(3,5,18,19)

######################################
# Contribution ratio(N)
######################################

lamf1 <- tilde.lam[,,2]
cumratio = matrix(0,N,T)

for (i in 1:T) {
  for (j in 1:N) {
    numerator <- sum((tilde.lam[i, 1:j, 2])^2)  
    denominator <- sum((tilde.lam[i, , 2])^2)  
    cumratio[j, i] <- numerator / denominator    
  }
}

x11()
par(mar=c(5, 4, 6, 2), xpd=TRUE)
plot((1:T), cumratio[1,], col = 2, type="l", ylim=c(0,1), 
     ylab = "Contribution Ratio (N)",
     xlab = "Synchronized Time (t)")
polygon(c(1:T, rev(1:T)), c(rep(0,times = 847), rev(cumratio[1,])), 
        col = 2, border = NA)
for (i in 2:N){
  lines(1:T,cumratio[i,],col=ifelse(pre$Condition[i] == "Fam", 2, 3), 
        type="l")
  polygon(c(1:T, rev(1:T)), c(cumratio[i - 1,], rev(cumratio[i,])), 
          col = ifelse(pre$Condition[i] == "Fam", 2, 3), border = NA)
}

Contr_N <- cumratio[c(which(pre$ID %in% rows_N), which(pre$ID %in% rows_N)-1),]
T_text <- 0
Y_text <- 0
for(i in 1:length(rows_N)){
  T_text[i] <- 
    which(Contr_N[i,]-Contr_N[i+length(rows_N),]==
            max(Contr_N[i,]-Contr_N[i+length(rows_N),]))
  
  Y_text[i] <- mean(c(Contr_N[i,T_text[i]], Contr_N[i+length(rows_N),T_text[i]]))
}

text(T_text, Y_text, labels = rows_N, font = 2, col = 1)

legend(x = "top", legend = levels(pre$Condition), fill = c(2:3), inset = -0.1,
       ncol = 2, title = "Familiarity Condition")

######################################
# Contribution ratio(P)
######################################

phif1 <- tilde.phi[,2,]
cumratiop = matrix(0,P,T)
for (i in 1:T) {
  for (j in 1:P) {
    numerator <- sum((tilde.phi[1:j, 2, i])^2)  
    denominator <- sum((tilde.phi[, 2, i])^2)  
    cumratiop[j, i] <- numerator / denominator    
  }
}

x11()

par(mar=c(5, 4, 7, 2), xpd=TRUE)
plot(1:T, cumratiop[1,], col = 1, type="l", ylim=c(0,1), 
     ylab = "Contribution Ratio (P)",
     xlab = "Synchronized Time (t)")
polygon(c(1:T, rev(1:T)), c(rep(0,times = 847), rev(cumratiop[1,])), 
        col = 1, border = NA)
for (i in 2:P){
  lines(1:T,cumratiop[i,],col=i, type="l")
  polygon(c(1:T, rev(1:T)), c(cumratiop[i - 1,], rev(cumratiop[i,])), 
                 col = i, border = NA)
}

Contr_P <- cumratiop[c(rows_P, rows_P-1),]
T_text <- 0
Y_text <- 0
for(i in 1:length(rows_P)){
  T_text[i] <- 
    which(Contr_P[i,]-Contr_P[i+length(rows_P),]==
            max(Contr_P[i,]-Contr_P[i+length(rows_P),]))
  
  Y_text[i] <- mean(c(Contr_P[i,T_text[i]], Contr_P[i+length(rows_P),T_text[i]]))
}

text(T_text, Y_text, labels = AUN$AU[rows_P], font = 2, col = 1)

legend(x = "top", legend = AUN$FAU, fill = c(1:23), inset = -0.17,
       ncol = 6, title = "Facial Action Synchrony (in order from bottom to top)")

######################################
# Combined Contribution Plot 
######################################

x11()

cumration = array(0,c(N,T,r-1))
cumratiop = array(0,c(P,T,r-1))

for (k in 1:(r-1)){
  for (i in 1:T) {
    for (j in 1:N) {
      numerator <- sum((tilde.lam[i, 1:j, k+1])^2)  
      denominator <- sum((tilde.lam[i, , k+1])^2)  
      cumration[j, i, k] <- numerator / denominator 
    }
    for (j in 1:P){
      numerator <- sum((tilde.phi[1:j, k+1, i])^2)  
      denominator <- sum((tilde.phi[, k+1, i])^2)  
      cumratiop[j, i, k] <- numerator / denominator
    }
  }
}

library(berryFunctions)

par(mar=c(0, 0, 0, 0), mfcol = c(2,2), oma = c(7,5,1,0), xpd=TRUE)
for(i in 1:(r-1)){
  # Plot N Contribution plots for each factor 
  par(mar = c(0,0,1,0))
  plot(1:T, cumration[1,,i], col = 24, type="l", ylim=c(0,1),
       xlab = "", ylab = "", axes = FALSE, main = ifelse(i==1, "Factor 1", "Factor 2"))
  polygon(c(1:T, rev(1:T)), c(rep(0,times = T), rev(cumration[1,,i])), 
          col = 24, border = NA)
  for (j in 2:N){
    lines(1:T,cumration[j,,i],col=ifelse(pre$Condition[j] == "Unf", 24, 25), 
          type="l")
    polygon(c(1:T, rev(1:T)), c(cumration[j - 1,,i], rev(cumration[j,,i])), 
            col = ifelse(pre$Condition[j] == "Unf", 24, 25), border = NA)
  }
  Contr_N <- cumration[c(which(pre$ID %in% rows_N), which(pre$ID %in% rows_N)-1),,i]
  T_text <- 0;  Y_text <- 0
  for(j in 1:length(rows_N)){
    T_text[j] <- 
      which(Contr_N[j,]-Contr_N[j+length(rows_N),]==
              max(Contr_N[j,]-Contr_N[j+length(rows_N),]))
    Y_text[j] <- mean(c(Contr_N[j,T_text[j]], Contr_N[j+length(rows_N),T_text[j]]))
  }
  textField(T_text, Y_text, labels = rows_N, font = 2, cex = 1.05,
            fill = "white", col = "black", border = "black", 
            field = "rounded", em = 1.2, ex = 1.2)
  if(i == 1) {axis(2)}
  # Plot P Contribution plots for each factor
  par(mar = c(0,0,0,0))
  plot(1:T, cumratiop[1,,i], col = 1, type="l", ylim=c(0,1), 
       xlab = "", ylab = "", axes = FALSE)
  polygon(c(1:T, rev(1:T)), c(rep(0,times = 847), rev(cumratiop[1,,i])), 
          col = 1, border = NA)
  for (k in 2:P){
    lines(1:T,cumratiop[k,,i],col=k, type="l")
    polygon(c(1:T, rev(1:T)), c(cumratiop[k - 1,,i], rev(cumratiop[k,,i])), 
            col = k, border = NA)
  }
  Contr_P <- cumratiop[c(rows_P, rows_P-1),,i]
  T_text <- 0;  Y_text <- 0
  for(k in 1:length(rows_P)){
    T_text[k] <- 
      which(Contr_P[k,]-Contr_P[k+length(rows_P),]==
              max(Contr_P[k,]-Contr_P[k+length(rows_P),]))
    Y_text[k] <- mean(c(Contr_P[k,T_text[k]], Contr_P[k+length(rows_P),T_text[k]]))
  }
  textField(T_text, Y_text, labels = AUN$FAU[rows_P], font = 2, cex = 1.05,
            fill = "white", col = "black", border = "black",
            field = "rounded", ex = 1.2)
  if(i == 1) {axis(1); axis(2)}
  if(i == 2) {axis(1)}
  if(i == 2) {
    mtext("Time (t)",side=1,line=2.5,outer=TRUE,cex=1.3)
    mtext("Contribution Ratio",side=2,line=3,outer=TRUE,cex=1.3)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = "bottom", fill = c(NA, 24:25, 1:23),
       legend = c("Dyad | FAU", "Unfamiliar","Familiar", as.character(AUN$FAU)), 
       ncol = 9, xpd = TRUE, bty = 'n')

######################################
# Average Profile(DP)
######################################


cumpfl = matrix(0,P,T)
for (i in 1:T) {
  for (j in 1:P) {
    numerator <- sum((DP[1:j, 1:j, i]))
    cumpfl[j, i] <- numerator 
  }
}

x11()
plot((1:T)/T, cumpfl[1,], col = 1, type="l", ylim=c(0,1), ylab = "Average Profile")
for (i in 2:P){
  lines(1:T/T,cumpfl[i,],col=i, type="l")
  polygon(c(1:T/T, rev(1:T/T)), c(cumpfl[i - 1,], rev(cumpfl[i,])), col = i, border = NA)
}

x11()
plot((1:T)/T, cumpfl[1,], col = 1, type="l", ylim=c(0,1), ylab = "Average Profile")
for (i in 2:P){
  lines(1:T/T,cumpfl[i,],col=i, type="l")
}

######################################
# Average Profile(DN)
######################################


cumpfln = matrix(0,N,T)
for (i in 1:T) {
  for (j in 1:N) {
    numerator <- sum((DN[1:j, 1:j, i]))
    cumpfln[j, i] <- numerator 
  }
}

x11()
plot((1:T)/T, cumpfln[1,], col = 1, type="l", ylim=c(0,1), ylab = "Average Profile")
for (i in 2:N){
  lines(1:T/T,cumpfln[i,],col=i, type="l")
  polygon(c(1:T/T, rev(1:T/T)), c(cumpfln[i - 1,], rev(cumpfln[i,])), col = i, border = NA)
}

x11()
plot((1:T)/T, cumpfln[1,], col = 1, type="l", ylim=c(0,1), ylab = "Average Profile")
for (i in 2:N){
  lines(1:T/T,cumpfln[i,],col=i, type="l")
}

#####################################################################
#### Extracting and Smoothing Time-Varying means of the factors ####
###################################################################
h2 = 0.08
TT=floor(T/2)
v=seq(0,1,length.out=TT)

mt0=matrix(0,P,TT)
mt1=matrix(0,P,TT)
mt2=matrix(0,P,TT)
mt3=matrix(0,P,TT)
mt4=matrix(0,P,TT)
mt5=matrix(0,P,TT)
j=2
k=3
l=4# the 1st, 2nd and 3rd factors

##Calculate the first factor
X0j=t(tilde.phi[,j,1:(T-1)])
X1j=t(tilde.phi[,j,2:(T)])
Z0j=cbind(rep(1,T-1),X0j)


for (s in 1:TT){
  D=diag(u- s/TT)
  Ws=diag(dnorm((v[s]- u)/h2))/h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm=solve(crossprod(Z0j,crossprod(D*Ws*D,Z0j)))
  Wsm=Ws-crossprod(Ws*D,crossprod(crossprod(D1sm,t(Z0j)),crossprod(Z0j,D*Ws)))
  mt0[,s]= crossprod(X0j,crossprod(Wsm,matrix(1,T-1,1)))/sum(Wsm)
  mt1[,s]= crossprod(X1j,crossprod(Wsm,matrix(1,T-1,1)))/sum(Wsm)
}

rm(Ws,D,D1sm)

##Calculate the second factor

X0k=t(tilde.phi[,k,1:(T-1)])
X1k=t(tilde.phi[,k,2:(T)])
Z0k=cbind(rep(1,T-1),X0k)	

for (s in 1:TT){
  D=diag(u- s/TT)
  Ws=diag(dnorm((v[s]- u)/h2))/h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm=solve(crossprod(Z0k,crossprod(D*Ws*D,Z0k)))
  Wsm=Ws-crossprod(Ws*D,crossprod(crossprod(D1sm,t(Z0k)),crossprod(Z0k,D*Ws)))
  mt2[,s]= crossprod(X0k,crossprod(Wsm,matrix(1,T-1,1)))/sum(Wsm)
  mt3[,s]= crossprod(X1k,crossprod(Wsm,matrix(1,T-1,1)))/sum(Wsm)
}
rm(Ws,D,D1sm)

##Calculate the third factor

X0l=t(tilde.phi[,l,1:(T-1)])
X1l=t(tilde.phi[,l,2:(T)])
Z0l=cbind(rep(1,T-1),X0l)	

for (s in 1:TT){
  D=diag(u- s/TT)
  Ws=diag(dnorm((v[s]- u)/h2))/h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm=solve(crossprod(Z0l,crossprod(D*Ws*D,Z0l)))
  Wsm=Ws-crossprod(Ws*D,crossprod(crossprod(D1sm,t(Z0l)),crossprod(Z0l,D*Ws)))
  mt4[,s]= crossprod(X0l,crossprod(Wsm,matrix(1,T-1,1)))/sum(Wsm)
  mt5[,s]= crossprod(X1l,crossprod(Wsm,matrix(1,T-1,1)))/sum(Wsm)
}

rm(Ws,D,D1sm)

rows_N <- rownames(Face_PL)[c(11,27,28,35)]   

#Select factors in tilde.lam 3rd dim
Loadings <- array(data = tilde.lam[,,2:(r+1)], dim = c(T,N,r)) %>%
  matrix(nrow = dim(.)[1] * dim(.)[2], ncol = dim(.)[3]) %>%
  as.data.frame(.) %>%
  mutate(ID = as.factor(rep(pre$ID, each = T)), 
         Condition = as.factor(rep(pre$Condition, each = T)),
         T = rep(c(1:T), times = N),
         Period = as.factor(ifelse(T <= t2, "Sim", "Com"))) %>%  # for Period 2-3 ONLY 
  set_colnames(c("Factor1", "Factor2", "Factor3", 
                 "ID","Condition", "T_LD", "Period")) %>%
  filter(ID %in% rows_N)

##Select columns interested in from 1~23 here for Loadings and TV-muF here
rows_P <- c(3,5,18,19)
P_selected <- length(rows_P)

# final dataframe for Factor Means
library(reshape2)
mtr <- merge(melt(mt0[rows_P,]), melt(mt2[rows_P,]), by = c("Var1", "Var2")) %>% 
  merge(melt(mt4[rows_P,]), by = c("Var1","Var2")) %>%  # for 3rd Factor
  arrange(Var1, Var2) %>% 
  rename(FAU = Var1, T_FCT = Var2, 
         Factor1 = value.x, Factor2 = value.y, Factor3 = value) %>% 
  mutate(FAU = as.factor(AUN$FAU[rows_P[FAU]]), 
         Period = as.factor(ifelse(T_FCT <= t2/2, "Sim","Com")))

######################################
# Compute Barycenters 
######################################

bc_f <- mtr %>% 
  group_by(FAU, Period) %>% 
  filter(T_FCT == round(median(T_FCT))) %>%           # for median
  bind_rows(mtr %>% group_by(FAU) %>%                 # for shift
              filter(T_FCT == round(t2/2)+1) %>%         # for shift
              mutate(Period = c("Shift"))) %>%        # for shift
  arrange(FAU, T_FCT) %>%
  rename(x = Factor1, y = Factor2, z = Factor3) %>%   # for median
  ungroup()                                           


bc_l <- Loadings %>% 
  group_by(ID, Period) %>% 
  filter(T_LD == round(median(T_LD))) %>%             # for median
  bind_rows(Loadings %>% group_by(ID) %>%                                  # for shift
          filter(T_LD == round(t2)+1) %>%                                  # for shift
          mutate(Period = c("Shift"))) %>%                                 # for shift
  arrange(ID, T_LD) %>%                                                    # for shift
  select(-c(Condition, T_LD)) %>%                     # for median
  rename(x = Factor1, y = Factor2, z = Factor3) %>%   # for median
  merge(pre, by = "ID")

######################################
# Supplementary Profiles(4.27)
######################################
Pi_a <- matrix(0, P, N)
V_a <- matrix(0,3,3)
for (i in 1:P){
  for (j in 1:N){
    Pi_a[i,j] <- median(Pi[i,j,])
  }
}

#Select factor in loops range

V_a[1,1] <- mean(V[2,2,])
V_a[2,2] <- mean(V[3,3,])
V_a[3,3] <- mean(V[4,4,])


V_a[1,1] <- (V_a[1,1])^-0.5
V_a[2,2] <- (V_a[2,2])^-0.5
V_a[3,3] <- (V_a[3,3])^-0.5

lambda_a <- pre %>% select(-ID) %>% 
  mutate(Condition = ifelse(Condition == "Unf", -1, 1))
lambda_a2 <- as.matrix(lambda_a)

#SELECT FACTORS in V_a
Suppro1 <- Pi_a %*% lambda_a2[,c(1,2)] %*% V_a[-3,-3] # "-" here is which factor NOT wanted
Suppro2 <- Pi_a %*% lambda_a2[,c(1,3)] %*% V_a[-3,-3] # "-" here is which factor NOT wanted
Suppro3 <- Pi_a %*% lambda_a2[,c(1,4)] %*% V_a[-3,-3] # "-" here is which factor NOT wanted

#Suppro1 <-Pi_a %*% lambda_a2 %*% V_a[1,1]
supp_P <- Suppro1 %>% data.frame() %>% rename(Condition = X1, TSR = X2) %>%
  mutate(FAU = AUN$FAU) %>% 
  merge(Suppro2 %>% data.frame() %>% rename(Condition = X1, Ars = X2) %>% 
          mutate(FAU = AUN$FAU), by = c("FAU", "Condition")) %>% 
  merge(Suppro3 %>% data.frame() %>% rename(Condition = X1, Val = X2) %>% 
          mutate(FAU = AUN$FAU), by = c("FAU", "Condition")) %>% 
  mutate(across(TSR:Val, ~.x/30), Condition = Condition /2) %>%
  filter(FAU %in% AUN$FAU[rows_P])

######################################
# Plot Loadings(2.13)
######################################
library(ggrepel)
x11()
T_new <- Loadings$T_LD[Loadings$T_LD %% 50 == 1]
#SELECT FACTORS in AESTHETIC PARAMETERS
ggplot(data = Loadings, aes(x = Factor1, y = Factor2, color = Condition)) +
  geom_vline(xintercept = 0, linetype = "longdash") + 
  geom_hline(yintercept = 0, linetype = "longdash") +
  scale_color_manual(values=c("#000000", "#0072B2", 
                              "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#FF8888"), 
                     labels = c("Familiar","Unfamiliar", as.character(AUN$FAU[rows_P])),
                     breaks = c(levels(pre$Condition), as.character(AUN$FAU[rows_P])),
                     name = "Condition & \nFacial Action") +
  scale_shape_manual(values = c(15,17,19,2,6,9),
                     labels = c("Simple", "Shift", "Complex", "Arousal","Valence","TSR-S"), 
                     breaks = c("Sim", "Shift", "Com", "Ars","Val","TSR"), 
                     name = "Period")+ # & \nSupplementary Points") +
  scale_size(guide = "none") + 
  geom_point(data = Loadings %>% filter(T_LD == t2 + 1), 
             aes(Factor1, Factor2, color = "Fam", shape = "Shift"), size = 7) + 
  geom_point(data = Loadings %>% filter(T_LD %in% T_new), aes(shape = Period, size = 3)) +
  geom_path(aes(group = ID, color = Condition), size = 2,
            lineend = "round", linejoin = "round") +
  geom_path(data = bc_f, aes(x=x, y=y, group = FAU, color = FAU), 
            size = 2, alpha = .5, show.legend = FALSE, arrow = arrow(),
            lineend = "round", linejoin = "round", linetype = "solid") +
  geom_point(data = bc_f, aes(x=x, y=y, color = FAU, size = 3, shape = Period), 
             show.legend = TRUE) + 
  geom_text_repel(data = Loadings %>% filter(T_LD %in% T_new),
                   aes(label = T_LD, color = Condition), size = 5,
                   box.padding   = 0.25, point.padding = 0.05,
                   segment.color = 'grey50', max.overlaps = 25, 
                   bg.color = "white", bg.r = 0.15, show.legend = FALSE) +
  #geom_label_repel(data = Loadings %>% filter(T_LD == 1),
  #                aes(label = paste("Dyad", ID, sep = " "), color = "black"), 
  #                fontface = "bold", nudge_x = 0, nudge_y = 0.35, 
  #                arrow = arrow(length = unit(0.015, "npc")), size = 4,
  #                box.padding   = 0.25, point.padding = 0.5, segment.linetype = 6,
  #                segment.color = 'grey50', max.overlaps = 25) +
  geom_label_repel(data = bc_f %>% filter(Period == "Sim"),
                  aes(x = x, y = y, label = FAU, color = FAU), fontface = "bold", 
                  box.padding   = 0.25, point.padding = 0.75, inherit.aes = FALSE,
                  segment.color = 'grey50', max.overlaps = 25, 
                  nudge_x = -.05, nudge_y = -.05, show.legend = FALSE) +
  geom_point(data = supp_P, aes(x=Condition, y=TSR, col= FAU, shape = "TSR"), 
             size = 3, stroke = 2, show.legend = FALSE) +
  geom_point(data = supp_P, aes(x=Condition, y=Val, col= FAU, shape = "Val"), 
             size = 3, stroke = 2, show.legend = FALSE) +
  geom_point(data = supp_P, aes(x=Condition, y=Ars, col= FAU, shape = "Ars"), 
             size = 3, stroke = 2, show.legend = FALSE) +
  xlab("Factor 1: Familiarity") +
  ylab("Factor 2: Affect-Attention Tradeoff") + 
  theme_classic2() +
  theme(axis.title = element_text(size = 24, face = "bold"), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20))
  #ggtitle("Time Varying Smoothed Loadings") + 
  #theme(plot.title = element_text(hjust = 0.5))


######################################
# Plot Factors(2.17)
######################################
x11()
T_new2 <- mtr$T_FCT[mtr$T_FCT %% 25 == 1]
## Plot Factors
#SELECT FACTORS in AESTHETIC PARAMETERS
ggplot(data = mtr, aes(x = Factor1, y = Factor2, color = FAU)) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_color_manual(values=c("#000000", "#0072B2", 
                              "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#FF8888"), 
                     labels = c("Familiar","Unfamiliar", as.character(AUN$FAU[rows_P])),
                     breaks = c(levels(pre$Condition), as.character(AUN$FAU[rows_P])),
                     name = "Condition & \nFacial Action") +
  scale_shape_manual(values = c(15,17,19), 
                     labels = c("Simple", "Shift", "Complex"), 
                     breaks = c("Sim", "Shift", "Com"), name = "Period") +
  scale_size(guide = "none") + 
  geom_point(data = mtr %>% filter(T_FCT == t2/2 + 1), 
             aes(Factor1, Factor2, color = "Fam", shape = "Shift"), size = 7) + 
  geom_point(data = mtr %>% filter(T_FCT %in% T_new2), aes(shape = Period), size = 3) +
  geom_path(aes(group = FAU), size = 2, lineend = "round", linejoin = "round") +
  geom_path(data = bc_l, aes(x = x, y = y, group = ID, color = Condition), 
            size = 2, alpha = .5, arrow = arrow(), show.legend = FALSE,
            lineend = "round", linejoin = "round", linetype = "solid") +
  geom_point(data = bc_l, aes(x = x, y = y, shape = Period, color = Condition), 
             size = 5, alpha = .5) +
  geom_text_repel(data = mtr %>% filter(T_FCT %in% T_new2),
                  aes(label = T_FCT*2-1), size = 5, box.padding = 0.25, 
                  segment.color = 'grey50', max.overlaps = 25, point.padding = 0.05,
                  bg.color = "white", bg.r = 0.15, show.legend = FALSE) +
  geom_label_repel(data = bc_l %>% filter(Period == "Sim"),
                  aes(x=x, y=y, label = paste("Dyad", ID, sep = " "), color = "black"), 
                  fontface = "bold", show.legend = FALSE, 
                  arrow = arrow(length = unit(0.015, "npc")), size = 4,
                  box.padding   = 0.25, point.padding = 0.5, segment.linetype = 6,
                  segment.color = 'grey50', max.overlaps = 25) +
  geom_label_repel(data = mtr %>% filter(T_FCT == 1),
                   aes(label = FAU), fontface = "bold", nudge_x = .05, nudge_y = .05,
                   box.padding   = 0.25, point.padding = 0.75,
                   segment.color = 'grey50', max.overlaps = 25, 
                   show.legend = FALSE) +
  xlab("Factor 1: Familiarity") +
  ylab("Factor 2: Affect-Attention Tradeoff") + 
  #xlim(-.42, .75) + ylim(-.75, .42) + 
  theme_classic2() +
  theme(axis.title = element_text(size = 24, face = "bold"), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20))
  #ggtitle("Time Varying Means of Factors") + 
  #theme(plot.title = element_text(hjust = 0.5)) 


######################################
# TV-VAR - Matrix A(t)(2.19)
######################################

# First Factor VAR
At1=array(0,c(P,P,TT))
for (s in 1:TT){
  D = diag(u - s/TT)
  Ws = diag(dnorm((v[s] - u) / h2)) / h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm = solve(crossprod(Z0j, crossprod(D*Ws*D,Z0j)))
  Wsm = Ws - crossprod(Ws*D, crossprod(crossprod(D1sm, t(Z0j)), crossprod(Z0j, D*Ws)))
  XX0 = X0j - crossprod(matrix(1,1,T-1), matrix(mt0[,s],1,P))
  XX1 = X1j - crossprod(matrix(1,1,T-1), matrix(mt1[,s],1,P))
  At1[,,s] = crossprod(crossprod(XX0, crossprod(Wsm, XX1)), 
                      solve(crossprod(XX0, crossprod(Wsm,XX0))))
}

# Second Factor VAR
At2=array(0,c(P,P,TT))
for (s in 1:TT){
  D = diag(u - s/TT)
  Ws = diag(dnorm((v[s] - u) / h2)) / h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm = solve(crossprod(Z0k, crossprod(D*Ws*D,Z0k)))
  Wsm = Ws - crossprod(Ws*D, crossprod(crossprod(D1sm, t(Z0k)), crossprod(Z0k, D*Ws)))
  XX0 = X0k - crossprod(matrix(1,1,T-1), matrix(mt2[,s],1,P))
  XX1 = X1k - crossprod(matrix(1,1,T-1), matrix(mt3[,s],1,P))
  At2[,,s] = crossprod(crossprod(XX0, crossprod(Wsm, XX1)), 
                      solve(crossprod(XX0, crossprod(Wsm,XX0))))
}


######################################
# Time Invariant VAR(6.12)
######################################

M=dim(hat.phi)[2]
A.hat=array(0,c(P,P,M-1))		# qxq estimated VAR(1)-matrix
S.hat=array(0,c(P,P,M-1))		# qxq estimated error cov
Gf0.hat=S.hat						    # qxq estimated lag-0 cov 
GfA.hat=array(0,c(P^2,P^2,M-1))


for (m in 2:M){
  y=scale(matrix(t(hat.phi[,m,2:(T)]),T-1,P),center=TRUE,scale=FALSE)
  x=scale(matrix(t(hat.phi[,m,1:(T-1)]),T-1,P),center=TRUE,scale=FALSE)	
  A.hat[,,m-1]=t(solve(t(x)%*%x)%*%t(x)%*%y)
  eps.hat=y-x%*%t(A.hat[,,m-1])
  S.hat[,,m-1]=t(eps.hat)%*%eps.hat/(T-P-2)
  Gf0.hat[,,m-1]=t(x)%*%x/(T)
  gfA= kronecker(S.hat[,,m-1],solve(t(x)%*%x))
  GfA.hat[,,m-1]=matrix(gfA,P^2,P^2)
}

# A_t.1
A.hat1 <- A.hat[,,1]
asvar1 <- matrix(0, nrow = P, ncol = P)
for(i in 1:P){
  for (j in 1:P){
    asvar1[i,j] = sqrt(GfA.hat[P*(i-1)+j,P*(i-1)+j,1])}
}

# A_t.2
A.hat2 = A.hat[,,2]
asvar2 <- matrix(0, nrow = P, ncol = P)
for(i in 1:P){
  for (j in 1:P){
    asvar2[i,j] = sqrt(GfA.hat[P*(i-1)+j,P*(i-1)+j,2])}
}

######################################
# Data Cleaning for Plotting 
######################################
# Create variable names string
P_names <- AUN$FAU[rows_P]

# Change all '1's' in this section to '2's' for second factor vis

# Create Confidence Intervals for ggplot
upper <- matrix(A.hat1[rows_P, rows_P] + 1.96 * asvar1[rows_P, rows_P], # 1|2
                nrow = P_selected, ncol = P_selected,
                dimnames = list(c(colnames(Face_PL)[rows_P]), 
                                c(colnames(Face_PL)[rows_P])))

lower <- matrix(A.hat1[rows_P, rows_P] - 1.96 * asvar1[rows_P, rows_P], # 1|2
                nrow = P_selected, ncol = P_selected,
                dimnames = list(c(colnames(Face_PL)[rows_P]), 
                                c(colnames(Face_PL)[rows_P])))

# filter At 
data_list <- At1 %>%                                                    # 1|2
  melt(varnames = c("i","j","T_VAR"), value.name = "Factors_VAR") %>%
  filter(i %in% rows_P & j %in% rows_P) %>% 
  rowwise() %>%
  mutate(iN = colnames(Face_PL)[i], jN = colnames(Face_PL)[j],
         Facial_Action = paste(iN, jN, sep = ","),
         lty = ifelse(Factors_VAR > upper[which(rownames(upper) == iN),
                                          which(rownames(upper) == jN)] | 
                        Factors_VAR < lower[which(rownames(lower) == iN),
                                            which(rownames(lower) == jN)] | 
                        upper[which(rownames(upper) == iN),
                              which(rownames(upper) == jN)] * 
                        lower[which(rownames(lower) == iN),
                              which(rownames(lower) == jN)] > 0,
                      "solid", "dashed")) %>% # add Section label here for lty?  
  group_by(i) %>% 
  group_split()

######################################
# Plot TV-VAR(2.26)
######################################
# Create plots
library(gridExtra)
# Create Task Windows
ann_col_pal <- c("#0072B2", "#D55E00")
ann_T <- data.frame(
  xmin = c(0,t2), xmax = c(t2,T), ymin = c(-Inf, -Inf), ymax = c(Inf, Inf),
  Task = factor(c("Simple", "Complex"))
)

col_pal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")


plot_list <- list()

for (i in seq_along(data_list)){
  
  TV_plot <- ggplot(data_list[[i]], aes(x = 2*T_VAR, y = Factors_VAR, colour = jN)) +
    scale_alpha(guide = 'none') +  
    scale_x_continuous(expand = c(0,0)) + 
    scale_color_manual(values = col_pal, labels = P_names, 
                       breaks = colnames(Face_PL)[rows_P], name = "Caused by:") +
    scale_fill_manual(name = "Task Period:", values = ann_col_pal,
                      breaks = c("Simple", "Complex"), labels = c("Simple","Complex")) + 
    geom_rect(data = ann_T, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = Task),
              alpha = 0.2, inherit.aes = FALSE) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_path(aes(group = Facial_Action, color = jN, alpha = ifelse(lty=="solid",1,.2)), 
              size = 4, lineend = "round", linejoin = "round") +
    labs(title = P_names[i]) + 
    ylim(-5.5,20) + 
    theme_classic2() + 
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size=18), 
          plot.title = element_text(color = col_pal[i], face = "bold", size = 24)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 8)),
           fill = guide_legend(override.aes = list(size = 8)))  
    
  plot_list[[i]] <- TV_plot
}

# plot figure
x11()
VAR <- ggarrange(plotlist = plot_list, common.legend = TRUE, legend = "right",
                 ncol = 2, nrow = 2) 

annotate_figure(VAR, 
                left = text_grob("Granger Causality", 
                                 face = "bold", size = 28, rot = 90), 
                bottom = text_grob("Time (s)", 
                                   face = "bold", size = 28), 
                #top = text_grob("How Nonverbal Synchrony Granger Causes Nonverbal Synchrony Over Time", 
                #                face = "bold", size = 15),
                fig.lab.face = "bold")


######################################
# Plot TIV-VAR(6.13)
######################################
# Create titles
titles_list1 <- list()

for(i in 1:P_selected){
  j = ifelse(i == 1, 2, 1) 
  k = ifelse(i == 1, 3, ifelse(i == 2, 3, 2))
  l = ifelse(i == 1, 4, ifelse(i == 2, 4, ifelse(i == 3, 4, 3)))
  titles_list1[[i]] <- list(bquote(phi[.(P_names[i])]*italic((t))~caused~by~
                                     phi[.(P_names[min(c(i,j,k,l))])]*italic((t-1))),
                            bquote(phi[.(P_names[i])]*italic((t))~caused~by~
                                     phi[.(P_names[min(c(i,j,k,l))+1])]*italic((t-1))),
                            bquote(phi[.(P_names[i])]*italic((t))~caused~by~
                                     phi[.(P_names[min(c(i,j,k,l))+2])]*italic((t-1))),
                            bquote(phi[.(P_names[i])]*italic((t))~caused~by~
                                     phi[.(P_names[min(c(i,j,k,l))+3])]*italic((t-1))))
}

# Create titles 
titles_list2 <- list()


for(i in 1:P_selected){
  j = ifelse(i == 1, 2, 1) 
  k = ifelse(i == 1, 3, ifelse(i == 2, 3, 2))
  l = ifelse(i == 1, 4, ifelse(i == 2, 4, ifelse(i == 3, 4, 3)))
  # for more than 4 rows_P, add another index here, and below
  titles_list2[[i]] <- 
    bquote(How~phi[.(P_names[i])]*italic((t))~is~caused~by~
             phi[.(P_names[min(c(i,j,k,l))])]*italic((t-1))*','~ 
             phi[.(P_names[min(c(i,j,k,l))+1])]*italic((t-1))*','~ 
             phi[.(P_names[min(c(i,j,k,l))+2])]*italic((t-1))*','~and~ 
             phi[.(P_names[min(c(i,j,k,l))+3])]*italic((t-1)))
}

# Create time dimension for plotting
t_seq = c(1:TT)

x11()
par(mfrow=c(4,4))
for (i in 1:4){
  for (j in 1:4){
    index1 = as.numeric(rows_P[i])
    index2 = as.numeric(rows_P[j])
    plot(t_seq*2, At[index1,index2,], col='red', type="l", 
         ylim = c(-10,10),
         ylab = 'VAR coeffcient', xlab = 'T', 
         main = titles_list1[[i]][[j]])
    lines(t_seq*2, rep(A.hat[index1,index2,1],TT), col = 'black')
    lines(t_seq*2, rep(A.hat[index1,index2,1] + 1.96 * asvar[index1,index2],TT), 
          lty='dashed',col = 'green')
    lines(t_seq*2, rep(A.hat[index1,index2,1] - 1.96 * asvar[index1,index2],TT), 
          lty='dashed',col = 'blue')
  }
}

#ylim=c(min(At[index1,index2,])-2,max(At[index1,index2,])+2)

x11()
par(mfrow=c(4,4))
for (i in 1:4){
  for (j in 1:4){
    index1 = as.numeric(rows_P[i])
    index2 = as.numeric(rows_P[j])
    plot(t_seq*2, rep(A.hat[index1,index2,1],TT), col = 'red',type="l", 
         ylim=c(-10,10), 
         ylab = 'VAR coeffcient', xlab = 'T', 
         main = titles_list1[[i]][[j]])
    lines(t_seq*2, rep(A.hat[index1,index2,1] + 1.96 * asvar[index1,index2],TT), 
          lty='dashed',col = 'green')
    lines(t_seq*2, rep(A.hat[index1,index2,1] - 1.96 * asvar[index1,index2],TT), 
          lty='dashed',col = 'blue')
  }
}
