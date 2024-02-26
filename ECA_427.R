rm(list=ls())
load("D:/NVB/3309_X2.rda")
N=dim(Face_PL)[1]
P=dim(Face_PL)[2]
faces <- Face_PL

times=read.csv("D:/NVB/N82_Times.csv")
t1=min(times$A)
t1max = max(times$A)
t2=min(times$B - times$A)
t2max=max(times$B - times$A)
t3=min(times$C - times$B)
t3max=max(times$C - times$B)
t4=min(times$D - times$C)
t4max=max(times$D - times$C)

T <- t1+t2+t3+t4-1
#make it 860 by letting t3=493

C1=array(0,c(N,P,t1))
C2=array(0,c(N,P,t2))
C3=array(0,c(N,P,t3-1))
C4=array(0,c(N,P,t4))

h_c=c(0.3, 0.14)


#h=2
u=c(1:T)/T
u1=seq(1/t1,1,length.out=t1)
u2=seq(1/t2,1,length.out=t2)
u3=seq(1/(t3-1),1,length.out=(t3-1))
u4=seq(1/t4,1,length.out=t4)

r=2
#u=1:T

for (n in (1:N)){	
  t=1:t1max/t1max
  #t=1:tmax[n]
  k=matrix(0,t1max,t1)
  for (p in (1:P)){	
    y=faces[n,p,1:t1max]
    for (j in (1:t1)){	
      #m=dnorm((t-u[j])/h)/h
      m=as.numeric(abs((t-u1[j])/h_c[1])<=1)*(rep(1,t1max)- ((t-u1)/h_c[1])^2)
      k[,j]=m/sum(m)}
    C1[n,p,]=y%*%k
  }}

for (n in (1:N)){	
  t=1:t2max/t2max
  #t=1:tmax[n]
  k=matrix(0,t2max,t2)
  for (p in (1:P)){	
    y=faces[n,p,1:t2max]
    for (j in (1:t2)){	
      #m=dnorm((t-u[j])/h)/h
      m=as.numeric(abs((t-u2[j])/h_c[2])<=1)*(rep(1,t2max)- ((t-u2[j])/h_c[2])^2)
      k[,j]=m/sum(m)}
    C2[n,p,]=y%*%k
  }}

for (n in (1:N)){	
  t=1:t3max/t3max
  #t=1:tmax[n]
  k=matrix(0,t3max,t3-1)
  for (p in (1:P)){	
    y=faces[n,p,1:t3max]
    for (j in (1:t3-1)){	
      #m=dnorm((t-u[j])/h)/h
      m=as.numeric(abs((t-u3[j])/h_c[2])<=1)*(rep(1,t3max)- ((t-u3[j])/h_c[2])^2)
      k[,j]=m/sum(m)}
    C3[n,p,]=y%*%k
  }}

for (n in (1:N)){	
  t=1:t4max/t4max
  #t=1:tmax[n]
  k=matrix(0,t4max,t4)
  for (p in (1:P)){	
    y=faces[n,p,1:t4max]
    for (j in (1:t4)){	
      #m=dnorm((t-u[j])/h)/h
      m=as.numeric(abs((t-u4[j])/h_c[1])<=1)*(rep(1,t4max)- ((t-u4[j])/h_c[1])^2)
      k[,j]=m/sum(m)}
    C4[n,p,]=y%*%k
  }}

library(abind)
C <- abind(C1, C2, C3, C4, along = 3)
any(C) == 0
dim(C)

DN=array(0,c(N,N,T))
DP=array(0,c(P,P,T))
Pi =array(0,c(P,N,T))
Psi=array(0,c(P,N,T))
F=Psi
V=array(0,c(r+1,r+1,T))
hat.phi = array(0,c(P,r+1,T))
hat.lam = array(0,c(N,r+1,T))

for (t in 1:T){
  F[,,t]=C[,,t]/sum(C[,,t])
  
  for (p in 1:P){
    DP[p,p,t]=sum(F[p,,t])}
  
  for (n in 1:N){
    DN[n,n,t]=sum(F[,n,t])}}

for (t in 1:T){
  Pi[,,t]=solve(DP[,,t])%*%F[,,t]}

for (t in 1:T){
  Psi[,,t]=F[,,t]%*%solve(DN[,,t])
}

for (t in 1:T){
  for (n in 1:N){
    if (DN[n,n,t]<0.0001){DN[n,n,t]=0.0001}
  }}
any(DN)==0

######################################%%%%%%%%%%%%%%%%%
# Extracting loadings and factors
######################################%%%%%%%%%%%%%%%%%

for (t in 1:T){
  Y=solve(sqrt(DP[,,t]))%*%F[,,t]%*%solve(sqrt(DN[,,t]))
  X=Y-sqrt(DP[,,t])%*%matrix(1,P,1)%*%matrix(1,1,N)%*%sqrt(DN[,,t])
  X=Y
  S=t(X)%*%X
  W=eigen(S)$vectors[,1:(r+1)]
  V[,,t]=diag(eigen(S)$values[1:(r+1)])
  hat.phi[,,t]=solve(DP[,,t])%*%F[,,t]%*%solve(sqrt(DN[,,t]))%*%W
  hat.lam[,,t]=solve(sqrt(DN[,,t]))%*%W%*%sqrt(V[,,t])
}

#---- Smoothing and plotting the loadings ----

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


m=59
w=(m + 1 - 0:m)/sum(1:(m+1))
w=matrix(w,1,m+1)

for (t in (m+1):T){
  for (j in 1:(r+1)){
    for (i in 1:N){
      y=matrix(adj.lam[i,j,t:(t-m)],m+1,1)
      tilde.lam[t,i,j]=w%*%y}}}

library(tidyverse)
library(magrittr)
rows_N <- c(17, 18, 31, 32)
Loadings <- array(data = tilde.lam[,,2:(r+1)], dim = c(T,N,r)) %>%
  matrix(nrow = dim(.)[1] * dim(.)[2], ncol = dim(.)[3]) %>%
  as.data.frame(.) %>%
  mutate(ID = as.factor(floor((1:nrow(.) - 1) / T) + 1), 
         T = rep(c(1:T), times = N)) %>%
  set_colnames(c("Factor1", "Factor2", "ID", "T_LD")) %>%
  filter(ID %in% rows_N)

h2 = 0.08
TT=T/2
v=seq(0,1,length.out=TT)

mt0=matrix(0,P,TT)
mt1=matrix(0,P,TT)
mt2=matrix(0,P,TT)
j=2
k=3# the first and second factor
X0=t(cbind(hat.phi[,j,1],hat.phi[,j,-T]))
X1=t(hat.phi[,j,])	
Z0=cbind(rep(1,T),X0)	

for (s in 1:TT){
  D=diag(u- s/TT)
  Ws=diag(dnorm((v[s]- u)/h2))/h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm=solve(crossprod(Z0,crossprod(D*Ws*D,Z0)))
  Wsm=Ws-crossprod(Ws*D,crossprod(crossprod(D1sm,t(Z0)),crossprod(Z0,D*Ws)))
  mt0[,s]= crossprod(X0,crossprod(Wsm,matrix(1,T,1)))/sum(Wsm)
  mt1[,s]= crossprod(X1,crossprod(Wsm,matrix(1,T,1)))/sum(Wsm)
}

X0=t(cbind(hat.phi[,k,1],hat.phi[,k,-T]))
X1=t(hat.phi[,k,])	
Z0=cbind(rep(1,T),X0)	

for (s in 1:TT){
  D=diag(u- s/TT)
  Ws=diag(dnorm((v[s]- u)/h2))/h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm=solve(crossprod(Z0,crossprod(D*Ws*D,Z0)))
  Wsm=Ws-crossprod(Ws*D,crossprod(crossprod(D1sm,t(Z0)),crossprod(Z0,D*Ws)))
  mt2[,s]= crossprod(X0,crossprod(Wsm,matrix(1,T,1)))/sum(Wsm)
}
rm(Ws,D,D1sm)

rows_P <- c(3, 5, 18, 19)
P_selected <- length(rows_P)
mtr1 <- matrix(0, P_selected, TT)
mtr2 <- matrix(0, P_selected, TT)

for (i in seq_along(rows_P)) {
  mtr1[i,] <- mt0[rows_P[i],]
  mtr2[i,] <- mt2[rows_P[i],]
}

## Convert mtr1 and mtr2 to columns which are easier for ggplot to read
mtr1 <- array(data = t(mtr1), c(TT*P_selected,1)) %>% 
  as.data.frame(.)
mtr2 <- array(data = t(mtr2), c(TT*P_selected,1)) %>% 
  as.data.frame(.)

mtr <- data.frame(cbind(mtr1, mtr2,
                        T_FCT = rep(c(1:TT), times = P_selected), 
                        P_FCT = rep(c('AU03', 'AU05', 'AU18', 'AU19'), each = TT))) %>%
  set_colnames(c("Factor1_F", "Factor2_F", "T_FCT", "Facial_Action")) %>%
  mutate(Facial_Action = factor(Facial_Action))

#Barycenters

bc_f <- mtr %>%
     group_by(Facial_Action) %>%
     summarize(bc_fx = mean(as.numeric(Factor1_F), na.rm = TRUE),
               bc_fy = mean(as.numeric(Factor2_F), na.rm = TRUE))

bc_l <- Loadings %>%
  group_by(ID) %>%
  summarize(bc_lx1 = mean(as.numeric(Factor1[T_LD == median(t1)]), na.rm = TRUE),
            bc_ly1 = mean(as.numeric(Factor2[T_LD == median(t1)]), na.rm = TRUE),
            bc_lx2 = mean(as.numeric(Factor1[T_LD == t1 + round(median(t2))]), na.rm = TRUE),
            bc_ly2 = mean(as.numeric(Factor2[T_LD == t1 + round(median(t2))]), na.rm = TRUE),
            bc_lx3 = mean(as.numeric(Factor1[T_LD == t1 + t2+round(median(t3-1))]), na.rm = TRUE),
            bc_ly3 = mean(as.numeric(Factor2[T_LD == t1 + t2+round(median(t3-1))]), na.rm = TRUE),
            bc_lx4 = mean(as.numeric(Factor1[T_LD == t1 + t2+t3-1+round(median(t4))]), na.rm = TRUE),
            bc_ly4 = mean(as.numeric(Factor2[T_LD == t1 + t2+t3-1+round(median(t4))]), na.rm = TRUE),)

######################################
# Supplementary Profiles(4.27)
######################################

Pi_a <- matrix(0, P, N)
V_a <- matrix(0,2,2)
for (i in 1:P){
  for (j in 1:N){
    Pi_a[i,j] <- mean(Pi[i,j,])
  }
}

for (i in 2:3){
  for (j in 2:3){
    V_a[i-1,j-1] <- mean(V[i,j,])
  }
}

V_a[1,1] <- (V_a[1,1])^-0.5
V_a[2,2] <- (V_a[2,2])^-0.5
lambda_a <- read.csv("D:/NVB/N82_supp.csv")[,2:3]
lambda_a <- lambda_a %>% 
  mutate(Condition = case_when(
    Condition == "Unf" ~ -1,
    Condition == "Fam" ~ 1,
    TRUE ~ NA_integer_  
  ))
lambda_a <- as.matrix(lambda_a)

Suppro <- Pi_a %*% lambda_a %*% V_a
######################################
# Plot Loadings(2.13)
######################################

library(ggpubr)
T_new <- Loadings$T_LD[Loadings$T_LD %% 20 == 1]
#I delete label = T_LD in ggplot() so I could add barycenters by feeding vector
ggplot(data = Loadings, aes(x = Factor1, y = Factor2, 
                                color = ID)) +
  #Set the color changed by T, but every N has DIFFERENT color, 
  #and the value of two colors will never be the same.
  geom_point() +
#set color manually so close color for same status.
  scale_color_manual(values=c("darkblue", "darkred", "blue", "red", 'orange', 'green', 'pink', 'purple'))+
  geom_point(data = bc_f, aes(x = bc_fx, y = bc_fy, color = Facial_Action, fill = Facial_Action), size = 7.5, shape = 21L) +
  geom_path(data = Loadings, aes(group = ID, alpha = T_LD), size = 1,
            lineend = "round", linejoin = "round") +
  geom_text(data = Loadings %>% filter(T_LD %in% T_new), 
            aes(label = T_LD),nudge_x = 0.01, nudge_y = 0.01, size = 3.5)+
  xlab("Factor 1") +
  ylab("Factor 2") + 
  theme_classic2()+
  ggtitle("Time Varying Smoothed Loadings") + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed")


######################################
# Plot Factors(2.17)
######################################

T_new2 <- mtr$T_FCT[mtr$T_FCT %% 20 == 1]
#I delete label = T_FCT in ggplot() so I could add barycenters by feeding vector
ggplot(data = mtr, aes(x = Factor1_F, y = Factor2_F, 
                        colour = Facial_Action)) +
  #Set the color changed by T, but every N has DIFFERENT color, 
  #and the value of two colors will never be the same.
  geom_point() +
  scale_color_manual(values=c('orange', 'green', 'purple', "darkblue", "darkred", "blue", "red", "darkgreen"))+
  geom_point(data = bc_l, aes(x = bc_lx, y = bc_ly, color = ID, fill = ID), size = 7.5, shape = 21L)+
  geom_path(data = mtr, aes(group = Facial_Action, alpha = T_FCT),size = 0.2,
            lineend = "round", linejoin = "round") +
  geom_text(data = mtr %>% filter(T_FCT %in% T_new2),
            aes(label = T_FCT),nudge_x = 0.002, nudge_y = 0.002, size = 3.5)+
  xlab("Factor 1") +
  ylab("Factor 2") + 
  theme_classic2() +
  ggtitle("Time Varying Means of Factors") + 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed")

######################################
# Matrix A(t)(2.19)
######################################

# time-varying VAR matrix

At=array(0,c(P,P,TT))
for (s in 1:TT){
  D=diag(u- s/TT)
  Ws=diag(dnorm((v[s]- u)/h2))/h2
  #Ws=diag(as.numeric(abs((v[s]- u)/h2)<=1))*diag(.75*(1- ((v[s]-u)/h2)^2))/h2
  D1sm=solve(crossprod(Z0,crossprod(D*Ws*D,Z0)))
  Wsm=Ws-crossprod(Ws*D,crossprod(crossprod(D1sm,t(Z0)),crossprod(Z0,D*Ws)))
  rm(Ws,D,D1sm)
  XX0=X0 - crossprod(matrix(1,1,T),matrix(mt0[,s],1,P))
  XX1=X1-crossprod(matrix(1,1,T),matrix(mt1[,s],1,P))
  At[,,s]=crossprod(crossprod(XX0,crossprod(Wsm,XX1)),solve(crossprod(XX0,crossprod(Wsm,XX0))))
}

######################################
# Data Cleaning(2.26)
######################################
At1 <- At[,,1]

for(i in 2:TT){
  At1 <- rbind(At1, At[,,i])
}

At1 <-  as.data.frame(At1) %>%
  mutate(T_VAR = floor((1:nrow(At1) - 1) / P) + 1, 
         P_VAR = rep(c(1:P), times = TT)) %>%
  filter(P_VAR == 3 | P_VAR == 6 | P_VAR == 18 | P_VAR == 19)

At23 <- cbind(At1[,3], At1[,6], At1[,18:19], At1[,24:25])
colnames(At23) <- c("P3", "P6", "P18", "P19")

At_full <- array(0, c(P_selected*TT*P_selected, 3))

for(i in 1:P_selected) {
  At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 1] <- At23[, i]
}

for (i in 1:P_selected) {
  if (i == 1) {
    At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 2] <- rep(c('AU04,AU04','AU04,AU07','AU04,AU19',
                                                                     'AU04,AU20'), times = TT)
  } else if (i == 2) {
    At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 2] <- rep(c('AU07,AU04','AU07,AU07','AU07,AU19',
                                                                     'AU07,AU20'), times = TT)
  } else if (i == 3) {
    At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 2] <- rep(c('AU19,AU04','AU19,AU07','AU19,AU19',
                                                                     'AU19,AU20'), times = TT)
  } else if (i == 4) {
    At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 2] <- rep(c('AU20,AU04','AU20,AU07','AU20,AU19',
                                                                     'AU20,AU20'), times = TT)
  } 
}

for (i in 1:P_selected){
  At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 3] <- rep(1:TT, each = P_selected)
  if (i > 1) {
    At_full[((i-1)*TT*P_selected + 1):(i*TT*P_selected), 3] <- At_full[((i-2)*TT*P_selected + 1):((i-1)*TT*P_selected), 3]
  }
}


colnames(At_full) <- c('Factors_VAR', 'Facial_Action', 'T_VAR')
At_full <- as.data.frame(At_full)

At3 <- At_full[1:(P_selected*TT),]
At6 <- At_full[(P_selected*TT + 1):(2*P_selected*TT),]
At18 <- At_full[(2*P_selected*TT + 1):(3*P_selected*TT),]
At19 <- At_full[(3*P_selected*TT + 1):(4*P_selected*TT),]

type_adjust <- function(df){
  df[,1] <- as.numeric(df[,1])
  df[,3] <- as.numeric(df[,3])
  df[,2] <- factor(df[,2])
  return(df)
}

At3 <- type_adjust(At3) 
At6 <- type_adjust(At6)
At18 <- type_adjust(At18)
At19 <- type_adjust(At19)

######################################
# Plot(2.26)
######################################


data_list <- list(At3, At6, At18, At19)
titles_list <- list(
  expression(paste("How ", phi[3](t), " is influenced by ", 
                   phi[3](t-1), ",", phi[6](t-1), ", ", 
                   phi[18](t-1), ", ", phi[19](t-1), ", and ", 
                   phi[20](t-1))),
  expression(paste("How ", phi[6](t), " is influenced by ", 
                   phi[3](t-1), ", ", phi[6](t-1), ", ", 
                   phi[18](t-1), ", ", phi[19](t-1), ", and ", 
                   phi[20](t-1))),
  expression(paste("How ", phi[18](t), " is influenced by ", 
                   phi[3](t-1), ", ", phi[6](t-1), ", ", 
                   phi[18](t-1), ", ", phi[19](t-1), ", and ", 
                   phi[20](t-1))),
  expression(paste("How ", phi[19](t), " is influenced by ", 
                   phi[3](t-1), ", ", phi[6](t-1), ", ", 
                   phi[18](t-1), ", ", phi[19](t-1), ", and ", 
                   phi[20](t-1)))
)

plot_list <- list()

for (i in seq_along(data_list)) {
  
  p <- ggplot(data_list[[i]], ylim = c(-150,150),aes(x = T_VAR, y = Factors_VAR, 
                                                     label = T_VAR, colour = Facial_Action, )) +
    geom_point() +
    geom_path(aes(group = Facial_Action),size = 1,
              lineend = "round", linejoin = "round") +
    xlab("Time") +
    ylab("Time-varying Coefficients") + 
    labs(color = "i,j stands for how \ni is influenced by j") +
    ggtitle(titles_list[[i]]) +
    theme(plot.title = element_text(size = 10.5)) +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylim(-150,150)
  
  plot_list[[i]] <- p
}

library(gridExtra)

grid.arrange(grobs = plot_list, ncol = 2)
