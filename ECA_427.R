#### Load the Data ####
rm(list=ls())
setwd("D:/Google Download/")
#load("3309_X2_ST_Imp.rda")
load("3309_D_X2_Imp.rda")

#Face_PL <- Face_PL[,c(3,5,6,8,18,19:23),]

N=dim(Face_PL)[1]
P=dim(Face_PL)[2]
T0=dim(Face_PL)[3]
faces <- Face_PL 
faces <- faces + 1

#### Synchronizing Uneven Time Series ####
### Define Periods
times=read.csv("N41_Times_Imp.csv")
#times=read.csv("N82_Times.csv")
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
V4=array(0,c(4+1,4+1,T))
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
  V4[,,t]=diag(eigen(S)$values[1:5])
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
plot((1:T)/T, cumratio[1,], col = 1, type="l", ylim=c(0,1), ylab = "Contribution Ratio")
for (i in 2:N){
  lines(1:T/T,cumratio[i,],col=i, type="l")
  polygon(c(1:T/T, rev(1:T/T)), c(cumratio[i - 1,], rev(cumratio[i,])), col = i, border = NA)
}

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
plot((1:T)/T, cumratiop[1,], col = 1, type="l", ylim=c(0,1), ylab = "Contribution Ratio(P)")
for (i in 2:P){
  lines(1:T/T,cumratiop[i,],col=i, type="l")
  polygon(c(1:T/T, rev(1:T/T)), c(cumratiop[i - 1,], rev(cumratiop[i,])), col = i, border = NA)
}

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

# Extracting and Smoothing Time-Varying means of the factors -----

h2 = 0.08
TT=floor(T/2)
v=seq(0,1,length.out=TT)

mt0=matrix(0,P,TT)
mt1=matrix(0,P,TT)
mt2=matrix(0,P,TT)
mt3=matrix(0,P,TT)
mt4=matrix(0,P,TT)
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
}

rm(Ws,D,D1sm)

# Data Cleaning ----
# Naming Variables - Dyadic & Ind Datasets (D/-41/82)
AUN <- data.frame(
  AU = as.factor(colnames(Face_PL)),
  FAU = as.factor(c("Inner Brow Raiser","R Outer Brow Raiser","Brow Lowerer",
                    "Upper Lid Raiser","Cheek Raiser","Lid Tightener","Nose Wrinkler",
                    "Upper Lip Raiser","Lip Corner Puller","Dimpler","Lip Corner_Depressor",
                    "Chin Raiser","Lip Stretcher","Lip Tightener","Lips Part","Jaw Drop",
                    "Lip Suck","Blink","Gaze X","Gaze Y","Pitch","Roll","Yaw")), 
  Face_Part = as.factor(c("Brow","Brow","Brow","Eye","Mouth","Eye","Nose","Nose","Mouth",
                          "Mouth","Mouth","Jaw","Mouth","Mouth","Mouth","Jaw","Mouth","Eye",
                          "Eye","Eye","Head","Head","Head")))
# Naming Variables - Cross Datasets (ST)
#AUN <- data.frame(
#  AU = as.factor(colnames(Face_PL)),
#  FAU = as.factor(c("S_GazeX","T_GazeX","S_GazeY","T_GazeY","S_Yaw","T_Yaw",
#                    "S_Pitch","T_Pitch","S_Roll","T_Roll","S_Brow Lowerer",
#                    "T_Brow Lowerer","S_Cheek Raiser","T_Cheek Raiser","S_Lid Tightener",
#                    "T_Lid Tightener","S_Upper Lip Raiser","T_Upper Lip Raiser",
#                    "S_Blink","T_Blink")))

#AUN <- data.frame(AU=as.factor(colnames(Face_PL)),
#                  FAU = as.factor(c("Brow Lowerer","Cheek Raiser","Lid Tightener",
#                                    "Upper Lip Raiser","Blink","Gaze X","Gaze Y",
#                                    "Pitch","Roll","Yaw")))

##Select rows interested in from 1~41 here for Loadings and TV-muF plots
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggnewscale)

#rows_N <- rownames(Face_PL)[c(7,8,9,16)]       # N = 41; 1 fem/mal teacher for fam/unf
#rows_N <- rownames(Face_PL)[c(1,18,22,37)]     # max/min DTW 
#rows_N <- rownames(Face_PL)[c(11,20,28,32)]    # max/min TSR_u
#rows_N <- rownames(Face_PL)[c(28,29,33,35)]    # max/min/median tilde.lam[,,2:3]
                                                 # 28 = med(F1); 29 = max(F2); 33 = med(F2)
                                                 # 35 = min(F2)&max(F1); 36 = min(F1)
rows_N <- rownames(Face_PL)[c(11,27,28,35)]     # max/min Var/SD(F1) ***
#rows_N <- rownames(Face_PL)[c(5,15,16,40)]     # max/min Euclidean Dist across time

#rows_N <- rownames(Face_PL)[c(11,28,35,36)]    # max/min (F1) & med(F1)
#rows_N <- str_extract(rownames(Face_PL)[c(7,8,9,16)], "\\d\\d") # for N = 41_S_Only
#rows_N <- rownames(Face_PL)[c(13:16)]          # N = 82

# For N41_D_Red
#rows_N <- rownames(Face_PL)[c(2,27,29,41)] # min/max F1 loadings values
#rows_N <- rownames(Face_PL)[c(11,23,27,35)] # min/max F1 loadings variance

# get condition information
pre <- read.csv("N41_supp.csv", header = T) %>% 
  mutate(Condition = factor(Condition),
         ID = factor(ID))
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
#rows_P <- c(1,2,5,6)
#rows_P <- c(1:10)
#rows_P <- c(1,2,19,20) # for ST: 1:2 - GazeX; 11-12 - AU04; 13-14 - AU06; 19-20 - AU45
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
  #summarize(x = mean(Factor1, na.rm = T),            # for mean
  #        y = mean(Factor2, na.rm = T),             # for mean
  #      z = mean(Factor3, na.rm = T)) %>%          # for mean
  rename(x = Factor1, y = Factor2, z = Factor3) %>%   # for median
  ungroup()                                           


bc_l <- Loadings %>% 
  group_by(ID, Period) %>% 
  filter(T_LD == round(median(T_LD))) %>%             # for median
  bind_rows(Loadings %>% group_by(ID) %>%                                  # for shift
          filter(T_LD == round(t2)+1) %>%                                  # for shift
          mutate(Period = c("Shift"))) %>%                                 # for shift
  arrange(ID, T_LD) %>%                                                    # for shift
  #summarize(x = mean(Factor1, na.rm = T),            # for mean
  #         y = mean(Factor2, na.rm = T),             # for mean
  #        z = mean(Factor3, na.rm = T)) %>%          # for mean
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

V_a[1,1] <- median(V[2,2,])
V_a[2,2] <- median(V[3,3,])
V_a[3,3] <- median(V[4,4,])


V_a[1,1] <- (V_a[1,1])^-0.5
V_a[2,2] <- (V_a[2,2])^-0.5
V_a[3,3] <- (V_a[3,3])^-0.5
#lambda_a <- read.csv("N82_supp.csv")[,2:3]
lambda_a <- read.csv("N41_supp.csv")[,2:3]

lambda_a2 <- lambda_a %>% 
  mutate(Condition = ifelse(Condition == "Unfamiliar", -1, 1))

lambda_a2 <- as.matrix(lambda_a2/30)

#SELECT FACTORS in V_a
Suppro <- Pi_a %*% lambda_a2 %*% V_a[-3,-3] # "-" here is which factor NOT wanted
#Suppro1 <-Pi_a %*% lambda_a2 %*% V_a[1,1]
supp_P <- Suppro[rows_P,] %>% data.frame() %>% 
  rename(Condition = X1, TSR = X2) %>% 
  mutate(FAU = colnames(Face_PL)[rows_P])

######################################
# Time-Varying Supplementary Profiles
######################################

TVSup = array(0, c(P,2,T))

for(i in 1:T){
  TVSup[,,i] <- Pi[,,i] %*% lambda_a2 %*% V[2:3,2:3,i]
}
######################################
# Supplementary Profiles(3 times)
######################################

lambda_ac1 <- as.matrix(read.csv("N41_supp_full.csv")[,1:2]/30)
lambda_ac2 <- as.matrix(cbind(read.csv("N41_supp_full.csv")[,1], read.csv("N41_supp_full.csv")[,3])/30)
lambda_ac3 <- as.matrix(cbind(read.csv("N41_supp_full.csv")[,1], read.csv("N41_supp_full.csv")[,4])/30)


#SELECT FACTORS in V_a
Suppro_ac1 <- Pi_a %*% lambda_ac1 %*% V_a[-3,-3] # "-" here is which factor NOT wanted
Suppro_ac2 <- Pi_a %*% lambda_ac2 %*% V_a[-3,-3] # "-" here is which factor NOT wanted
Suppro_ac3 <- Pi_a %*% lambda_ac3 %*% V_a[-3,-3] # "-" here is which factor NOT wanted
#Suppro1 <-Pi_a %*% lambda_a2 %*% V_a[1,1]

supp_ac1 <- Suppro_ac1[rows_P,] %>% data.frame() %>% 
  rename(Condition = X1, TSR_u = X2) %>% 
  mutate(FAU = colnames(Face_PL)[rows_P])

supp_ac2 <- Suppro_ac2[rows_P,] %>% data.frame() %>% 
  rename(Condition = X1, Arousal = X2) %>% 
  mutate(FAU = colnames(Face_PL)[rows_P])

supp_ac3 <- Suppro[rows_P,] %>% data.frame() %>% 
  rename(Condition = X1, Valence = X2) %>% 
  mutate(FAU = colnames(Face_PL)[rows_P])

######################################
# Plot Loadings(2.13)
######################################
x11()
T_new <- Loadings$T_LD[Loadings$T_LD %% 20 == 1]
#SELECT FACTORS in AESTHETIC PARAMETERS
ggplot(data = Loadings, aes(x = Factor1, y = Factor2, color = Condition)) +
  geom_vline(xintercept = 0, linetype = "longdash") + 
  geom_hline(yintercept = 0, linetype = "longdash") +
  scale_color_manual(values=c("#000000", "#F0E442", 
                              "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#FF8888"), 
                     labels = c(levels(pre$Condition), as.character(AUN$FAU[rows_P])),
                     breaks = c(levels(pre$Condition), as.character(AUN$FAU[rows_P])),
                     name = "Condition & \n Facial Muscle") +
  scale_shape_manual(values = c(15,17,19),
                     labels = c("Simple", "Shift", "Complex"), 
                     breaks = c("Sim", "Shift", "Com"), name = "Period") +
  scale_size(guide = "none") + 
  geom_point(data = Loadings %>% filter(T_LD == t2 + 1), 
             aes(Factor1, Factor2, color = "Familiar", shape = "Shift"), size = 5) + 
  geom_point(aes(shape = Period, size = 1)) +
  geom_path(aes(group = ID, color = Condition), size = 1,
            lineend = "round", linejoin = "round") +
  geom_text(data = Loadings %>% filter(T_LD %in% T_new), 
            aes(label = T_LD, color = Condition),
            nudge_x = 0.01, nudge_y = 0.01, size = 3.5) +
  geom_path(data = bc_f, aes(x=x, y=y, group = FAU, color = FAU), 
            size = 1, alpha = .5, arrow = arrow(), show.legend = FALSE,
            lineend = "round", linejoin = "round", linetype = "solid") +
  geom_point(data = bc_f, aes(x=x, y=y, color = FAU, size = 2, shape = Period), 
             show.legend = FALSE) + 
  geom_text(data = Loadings %>% filter(T_LD == 1), 
            aes(label = paste("Dyad", ID, sep = " ")), color = 'black',
            nudge_x = 0.025, nudge_y = -0.02, size = 4, fontface = "bold") +
  geom_text(data = bc_f %>% filter(Period == "Sim"), 
            aes(x=x, y=y, label = FAU, color = FAU),
            size = 4, nudge_x = -0.05, nudge_y = -0.025, fontface = "bold") + 
  geom_point(data = supp_P, aes(x=Condition, y=TSR,
                                col= as.character(AUN$FAU[rows_P])), size = 6) +
  xlab("Factor 1") +
  ylab("Factor 2") + 
  xlim(-0.7,1.15) + ylim(-1.15,0.7) + 
  theme_classic2() +
  ggtitle("Time Varying Smoothed Loadings") + 
  theme(plot.title = element_text(hjust = 0.5))


######################################
# Plot Factors(2.17)
######################################
x11()
T_new2 <- mtr$T_FCT[mtr$T_FCT %% 20 == 1]
## Plot Factors
#SELECT FACTORS in AESTHETIC PARAMETERS
ggplot(data = mtr, aes(x = Factor1, y = Factor2, color = FAU)) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_color_manual(values=c("#000000", "#F0E442", 
                              "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#FF8888"), 
                     labels = c(levels(pre$Condition), as.character(AUN$FAU[rows_P])),
                     breaks = c(levels(pre$Condition), as.character(AUN$FAU[rows_P])),
                     name = "Condition & \n Facial Muscle") +
  scale_shape_manual(values = c(15,17,19), 
                     labels = c("Simple", "Shift", "Complex"), 
                     breaks = c("Sim", "Shift", "Com"), name = "Period") +
  scale_size(guide = "none") + 
  geom_point(data = mtr %>% filter(T_FCT == t2/2 + 1), 
             aes(Factor1, Factor2, color = "Familiar", shape = "Shift"), size = 5) + 
  geom_point(aes(shape = Period), size = 1) +
  geom_path(aes(group = FAU), size = .5, lineend = "round", linejoin = "round") +
  geom_text(data = mtr %>% filter(T_FCT %in% T_new2), aes(label = T_FCT), 
            nudge_x = 0.005, nudge_y = 0.005, size = 3.5, show.legend = FALSE) +
  geom_path(data = bc_l, aes(x = x, y = y, group = ID, color = Condition), 
            size = 1, alpha = .5, arrow = arrow(), show.legend = FALSE,
            lineend = "round", linejoin = "round", linetype = "solid") +
  geom_point(data = bc_l, aes(x = x, y = y, shape = Period, color = Condition), 
             size = 5, alpha = .5) +
  geom_text(data = mtr %>% filter(T_FCT == 1), aes(label = FAU), fontface = "bold", 
            nudge_x = 0.025, nudge_y = 0.025, size = 4, show.legend = FALSE) +
  geom_point(data = supp_P, aes(x=Condition, y=TSR,
                                col= as.character(AUN$FAU[rows_P])), size = 6) +
  geom_text(data = bc_l %>% filter(Period == "Sim"), 
            aes(x = x, y = y, label = paste("Dyad",ID, sep = " ")),
                   size = 4, nudge_x = 0.025, nudge_y = 0.025, color = 'black',
                   fontface = "bold", show.legend = FALSE) +
  xlab("Factor 1") +
  ylab("Factor 2") + 
  theme_classic2() +
  ggtitle("Time Varying Means of Factors") + 
  theme(plot.title = element_text(hjust = 0.5)) 

######################################
# TV-VAR - Matrix A(t)(2.19)
######################################

# time-varying VAR matrix

At=array(0,c(P,P,TT))
for (s in 1:TT){
  D = diag(u - s/TT)
  Ws = diag(dnorm((v[s] - u) / h2)) / h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm = solve(crossprod(Z0j, crossprod(D*Ws*D,Z0j)))
  Wsm = Ws - crossprod(Ws*D, crossprod(crossprod(D1sm, t(Z0j)), crossprod(Z0j, D*Ws)))
  XX0 = X0j - crossprod(matrix(1,1,T-1), matrix(mt0[,s],1,P))
  XX1 = X1j - crossprod(matrix(1,1,T-1), matrix(mt1[,s],1,P))
  At[,,s] = crossprod(crossprod(XX0, crossprod(Wsm, XX1)), 
                      solve(crossprod(XX0, crossprod(Wsm,XX0))))
}

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

asvar <- matrix(0, nrow = P, ncol = P)
for(i in 1:P){
  for (j in 1:P){
    asvar[i,j] = sqrt(GfA.hat[P*(i-1)+j,P*(i-1)+j,1])}
}

######################################
# At2
######################################
A.hat2 = A.hat[,,2]
asvar2 <- matrix(0, nrow = P, ncol = P)
for(i in 1:P){
  for (j in 1:P){
    asvar2[i,j] = sqrt(GfA.hat[P*(i-1)+j,P*(i-1)+j,2])}
}

######################################
# Data Cleaning(2.26)
######################################
# Choose variables
#rows_P <- c(3,5,18,19)
#rows_P <- c(1,2,5,6)
#rows_P <- c(1,2,11:14,19,20) # for ST: 1:2 - GazeX; 11-12 - AU04; 13-14 - AU06; 19-20 - AU45
#P_selected <- length(rows_P)

# Create variable names string
P_names <- colnames(Face_PL)[c(rows_P)]
#P_names <- as.character(AUN$FAU[c(rows_P)])
# Create Confidence Intervals for ggplot
upper <- matrix(A.hat[rows_P, rows_P, 1] + 1.96 * asvar[rows_P, rows_P],
                nrow = P_selected, ncol = P_selected,
                dimnames = list(c(colnames(Face_PL)[rows_P]), 
                                c(colnames(Face_PL)[rows_P])))

lower <- matrix(A.hat[rows_P, rows_P, 1] - 1.96 * asvar[rows_P, rows_P],
                nrow = P_selected, ncol = P_selected,
                dimnames = list(c(colnames(Face_PL)[rows_P]), 
                                c(colnames(Face_PL)[rows_P])))

# filter At 
data_list <- At2 %>% 
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

# Create plots
library(gridExtra)

col_pal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
# For 8 colors (P_selected = 8)
#col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


plot_list <- list()

for (i in seq_along(data_list)){
  
  TV_plot <- ggplot(data_list[[i]], aes(x = 2*T_VAR, y = Factors_VAR, colour = jN)) +
    geom_rect(aes(xmin = 0, xmax = t2, ymin = -Inf, ymax = Inf, 
                  fill = "#0072B2"), alpha = .002, col = NA, show.legend = FALSE) + 
    geom_rect(aes(xmin = t2, xmax = T, ymin = -Inf, ymax = Inf, 
                  fill = "#D55E00"), alpha = .002, col = NA, show.legend = FALSE) + 
    #geom_vline(xintercept = t2, linetype = "longdash") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_path(aes(group = Facial_Action, color = jN, alpha = ifelse(lty=="solid",1,.25)), 
              size = 2, lineend = "round", linejoin = "round") +
    scale_alpha(guide = 'none') +  
    scale_x_continuous(expand = c(0,0)) + 
    scale_color_manual(values = col_pal, labels = P_names, 
                       breaks = colnames(Face_PL)[rows_P], name = "caused by") +
    labs(title = P_names[i]) + 
    #ggtitle(titles_list2[[i]]) +
    ylim(-2,12) + 
    theme_classic2() + 
    theme(legend.title = element_text(size = 14), legend.text=element_text(size=12), 
          plot.title = element_text(color = col_pal[i], face = "bold")) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())  
    
  plot_list[[i]] <- TV_plot
}

# plot figure
x11()
VAR <- ggarrange(plotlist = plot_list, common.legend = TRUE, legend = "top",
                 ncol = 2, nrow = 2) # nrow = 4 for ST

annotate_figure(VAR, 
                left = text_grob("Time-Varying VAR Coefficients", 
                                 face = "bold", size = 12, rot = 90), 
                bottom = text_grob("Time (s)", 
                                   face = "bold", size = 12), 
                top = text_grob("How Nonverbal Synchrony Granger Causes Nonverbal Synchrony Over Time", 
                                face = "bold", size = 15),
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

###
### Extra Plotting Vis ####
library(FactoMineR)
ca1 <- melt(Face_PL) %>% 
  spread(Var2, value) %>% 
  group_by(Var1) %>% 
  summarize(across(AU01:Yaw, ~mean(.x))) %>% #filter(Var1 != 39 & Var1 != 42) %>% 
  select(-Var1) %>%  
  CA()
x11()
plot.CA(ca1)

