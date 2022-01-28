library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(doParallel)
library(deSolve)
library(numDeriv)

### PLOT OPTIONS ####

path_figure="Figures/"
path_data="Data/"

theme<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour='grey'),
        panel.grid.major.y = element_line(colour='grey'),
        text = element_text(size=20, family="Times"),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5))

label_gamma<-expression("Asymmetry of interaction strength "*gamma)
label_contribution<-"Relative contribution to the lead eigenvector"
label_correlation<-expression("Correlation of C"["1"]*" and C"["2"])

colour_rooney_original<-scale_colour_manual(name="Species",
                                            values=c("chocolate1","chocolate1","chartreuse4","chartreuse4","red"),
                                            guide = guide_legend(reverse = TRUE),
                                            labels = c(expression("R"["1"]),
                                                       expression("R"["2"]),
                                                       expression("C"["1"]),
                                                       expression("C"["2"]),
                                                       "P"))

line_rooney_original<-scale_linetype_manual(name="Species",
                                            values=c("solid","22","solid","22","solid"),
                                            guide = guide_legend(reverse = TRUE),
                                            labels = c(expression("R"["1"]),
                                                       expression("R"["2"]),
                                                       expression("C"["1"]),
                                                       expression("C"["2"]),
                                                       "P"))

colour_rooney_dispersal<-scale_colour_manual(name="Species",
                                             values=c("chocolate1","chocolate1","chartreuse4","chartreuse4","red","red"),
                                             guide = guide_legend(reverse = TRUE),
                                             labels = c(expression("R"["1"]),
                                                        expression("R"["2"]),
                                                        expression("C"["1"]),
                                                        expression("C"["2"]),
                                                        expression("P"["1"]),
                                                        expression("P"["2"])))

line_rooney_dispersal<-scale_linetype_manual(name="Species",
                                             values=c("solid","22","solid","22","solid","22"),
                                             guide = guide_legend(reverse = TRUE),
                                             labels = c(expression("R"["1"]),
                                                        expression("R"["2"]),
                                                        expression("C"["1"]),
                                                        expression("C"["2"]),
                                                        expression("P"["1"]),
                                                        expression("P"["2"])))

perturbation_C_line<-scale_linetype_manual(values=c("solid","22"),
                                           labels=c(expression("C"["1"]),
                                                    expression("C"["2"])),
                                           name='perturbed\nspecies')

perturbation_P_line<-scale_linetype_manual(values=c("solid","22"),
                                           labels=c(expression("P"["1"]),
                                                    expression("P"["2"])),
                                           name='perturbed\nspecies')

x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

### FUNCTIONS ####

# ODE system of Rooney's model
ODE_rooney_original<-function(B, params){
  #B=(R1,R2,C1,C2,P)
  with(params,{
    dB<-rep(0,5)
    # dR1/dt
    dB[1] = r*B[1]*(p*(Nt-sum(B))+Ns)/(p*(Nt-sum(B))+Ns+b) - gamma*a_CR*B[1]*B[3]/(B[1]+R0_C) - l_R*B[1]
    # dR2/dt
    dB[2] = r*B[2]*((1-p)*(Nt-sum(B))+Ns)/((1-p)*(Nt-sum(B))+Ns+b) - a_CR*B[2]*B[4]/(B[2]+R0_C) - l_R*B[2]
    # dC1/dt
    dB[3] = gamma*e*a_CR*B[1]*B[3]/(B[1]+R0_C) - gamma*a_PC*B[3]^2*B[5]/(B[3]^2+B[4]^2+R0_P*(B[3]+B[4])) - m_C*B[3]
    # dC2/dt
    dB[4] = e*a_CR*B[2]*B[4]/(B[2]+R0_C) - a_PC*B[4]^2*B[5]/(B[3]^2+B[4]^2+R0_P*(B[3]+B[4])) - m_C*B[4]
    # dP/dt
    dB[5] = e*a_PC*B[5]*(gamma*B[3]^2+B[4]^2)/(B[3]^2+B[4]^2+R0_P*(B[3]+B[4])) - m_P*B[5]
    return(dB)
  })
}

# ODE of Rooney's model with dispersal for the numerical simulation of time series
ODE_rooney_dispersal<-function(B, params){
  #B=(Nc,R1,R2,C1,C2,P1,P2)
  with(params,{
    dB<-rep(0,6)
    # dR1/dt
    dB[1] = r*B[1]*(p*(Nt-sum(B))+Ns)/(p*(Nt-sum(B))+Ns+b) - gamma*a_CR*B[1]*B[3]/(B[1]+R0_C) - l_R*B[1]
    # dR2/dt
    dB[2] = r*B[2]*((1-p)*(Nt-sum(B))+Ns)/((1-p)*(Nt-sum(B))+Ns+b) - a_CR*B[2]*B[4]/(B[2]+R0_C) - l_R*B[2]
    # dC1/dt
    dB[3] = gamma*e*a_CR*B[1]*B[3]/(B[1]+R0_C) - gamma*a_PC*B[3]*B[5]/(B[3]+R0_P) - m_C*B[3]
    # dC2/dt
    dB[4] = e*a_CR*B[2]*B[4]/(B[2]+R0_C) - a_PC*B[4]*B[6]/(B[4]+R0_P) - m_C*B[4]
    # dP1/dt
    dB[5] = gamma*e*a_PC*B[3]*B[5]/(B[3]+R0_P) - m_P*B[5] + d*(B[3]/(B[3]+B[4])*B[6] - B[4]/(B[3]+B[4])*B[5])
    # dP2/dt
    dB[6] = e*a_PC*B[4]*B[6]/(B[4]+R0_P) - m_P*B[6] + d*(B[4]/(B[3]+B[4])*B[5] - B[3]/(B[3]+B[4])*B[6])
    return(dB)
  })
}

# function for the ODE solver (returns a list)
ODE_function_rooney_original<-function(t, B, params){
  dB<-ODE_rooney_original(B, params)
  return(list(dB))
}

# function for the ODE solver (returns a list)
ODE_function_rooney_dispersal<-function(t, B, params){
  dB<-ODE_rooney_dispersal(B, params)
  return(list(dB))
}

# Jacobian matrix of the original model of McCann
jacobian_rooney_original<-function(B, params){
  #B=(Nc,R1,R2,C1,C2,P)
  with(params,{
    J<-matrix(0,nrow=6,ncol=6)
    
    # dR1?/dtdNc
    J[2,1] = b*r*B[2]*p/(p*B[1]+Ns+b)^2
    # dR1?/dtdR1
    J[2,2] = r*(p*B[1]+Ns)/(p*B[1]+Ns+b) - gamma*a_CR*B[4]*R0_C/(B[2]+R0_C)^2 - l_R
    # dR1?/dtdR2 = 0
    # dR1?/dtdC1
    J[2,4] = -gamma*a_CR*B[2]/(B[2]+R0_C)
    # dR1?/dtdC2 = 0
    # dR1?/dtdP = 0
    
    # dR2?/dtdNc
    J[3,1] = b*r*B[3]*(1-p)/((1-p)*B[1]+Ns+b)^2
    # dR2?/dtdR1 = 0
    # dR2?/dtdR2
    J[3,3] = r*((1-p)*B[1]+Ns)/((1-p)*B[1]+Ns+b) - a_CR*B[5]*R0_C/(B[3]+R0_C)^2 - l_R
    # dR2?/dtdC1 = 0
    # dR2?/dtdC2
    J[3,5] = -a_CR*B[3]/(B[3]+R0_C)
    # dR2?/dtdP = 0
    
    # dC1?/dtdNc = 0
    # dC1?/dtdR1
    J[4,2] = gamma*e*a_CR*B[4]*R0_C/(B[2]+R0_C)^2
    # dC1?/dtdR2 = 0
    # dC1?/dtdC1
    J[4,4] = gamma*e*a_CR*B[2]/(B[2]+R0_C) - gamma*a_PC*B[4]*B[6]*(2*B[5]^2+R0_P*(B[4]+2*B[5]))/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))^2 - m_C
    # dC1?/dtdC2
    J[4,5] = gamma*a_PC*B[4]^2*B[6]*(2*B[5]+R0_P)/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))^2
    # dC1?/dtdP
    J[4,6] = -gamma*a_PC*B[4]^2/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))
    
    # dC2?/dtdNc = 0
    # dC2?/dtdR1 = 0
    # dC2?/dtdR2
    J[5,3] = e*a_CR*B[5]*R0_C/(B[3]+R0_C)^2
    # dC2?/dtdC1
    J[5,4] = a_PC*B[5]^2*B[6]*(2*B[4]+R0_P)/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))^2
    # dC2?/dtdC2
    J[5,5] = e*a_CR*B[3]/(B[3]+R0_C) - a_PC*B[5]*B[6]*(2*B[4]^2+R0_P*(B[5]+2*B[4]))/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))^2 - m_C
    # dC2?/dtdP
    J[5,6] = -a_PC*B[5]^2/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))
    
    # dP?/dtdNc = 0
    # dP?/dtdR1 = 0
    # dP?/dtdR2 = 0
    # dP?/dtdC1
    J[6,4] = (2*B[4]*B[5]^2*(gamma-1) + gamma*R0_P*B[4]^2 + 2*gamma*R0_P*B[4]*B[5] - R0_P*B[5]^2)*e*a_PC*B[6]/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))^2
    # dP?/dtdC2
    J[6,5] = (2*B[5]*B[4]^2*(1-gamma) + R0_P*B[5]^2 + 2*R0_P*B[4]*B[5] - gamma*R0_P*B[4]^2)*e*a_PC*B[6]/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5]))^2
    # dP?/dtdP
    J[6,6] = e*a_PC*(gamma*B[4]^2+B[5]^2)/(B[4]^2+B[5]^2+R0_P*(B[4]+B[5])) - m_P
    
    # dNc?/dtd...
    for(i in 1:6){
      J[1,i] = -sum(J[2:6,i]) # because Nc=Nt-(R1+R2+C1+C2+P)
    }
    return(J)
  })
}

# Jacobian matrix of the modified model of McCann
jacobian_rooney_dispersal<-function(B,params){
  #B=(Nc,R1,R2,C1,C2,P1,P2)
  with(params,{
    J<-matrix(0,nrow=7,ncol=7)
    
    # same as jacobian_McCann_original
    # dR1?/dtdNc
    J[2,1] = b*r*B[2]*p/(p*B[1]+Ns+b)^2
    # dR1?/dtdR1
    J[2,2] = r*(p*B[1]+Ns)/(p*B[1]+Ns+b) - gamma*a_CR*B[4]*R0_C/(B[2]+R0_C)^2 - l_R
    # dR1?/dtdR2 = 0
    # dR1?/dtdC1
    J[2,4] = -gamma*a_CR*B[2]/(B[2]+R0_C)
    # dR1?/dtdC2 = 0
    # dR1?/dtdP1 = 0
    # dR1?/dtdP2 = 0
    
    # same as jacobian_McCann_original
    # dR2?/dtdNc
    J[3,1] = b*r*B[3]*(1-p)/((1-p)*B[1]+Ns+b)^2
    # dR2?/dtdR1 = 0
    # dR2?/dtdR2
    J[3,3] = r*((1-p)*B[1]+Ns)/((1-p)*B[1]+Ns+b) - a_CR*B[5]*R0_C/(B[3]+R0_C)^2 - l_R
    # dR2?/dtdC1 = 0
    # dR2?/dtdC2
    J[3,5] = -a_CR*B[3]/(B[3]+R0_C)
    # dR2?/dtdP1 = 0
    # dR2?/dtdP2 = 0
    
    # dC1?/dtdNc = 0
    # dC1?/dtdR1
    J[4,2] = gamma*e*a_CR*B[4]*R0_C/(B[2]+R0_C)^2
    # dC1?/dtdR2 = 0
    # dC1?/dtdC1
    J[4,4] = gamma*e*a_CR*B[2]/(B[2]+R0_C) - gamma*a_PC*B[6]*R0_P/(B[4]+R0_P)^2 - m_C
    # dC1?/dtdC2 = 0
    # dC1?/dtdP1
    J[4,6] = -gamma*a_PC*B[4]/(B[4]+R0_P)
    # dC1?/dtdP2 = 0
    
    # dC2?/dtdNc = 0
    # dC2?/dtdR1 = 0
    # dC2?/dtdR2
    J[5,3] = e*a_CR*B[5]*R0_C/(B[3]+R0_C)^2
    # dC2?/dtdC1 = 0
    # dC2?/dtdC2
    J[5,5] = e*a_CR*B[3]/(B[3]+R0_C) - a_PC*B[7]*R0_P/(B[5]+R0_P)^2 - m_C
    # dC2?/dtdP1 = 0
    # dC2?/dtdP2
    J[5,7] = -a_PC*B[5]/(B[5]+R0_P)
    
    # dP1?/dtdNc = 0
    # dP1?/dtdR1 = 0
    # dP1?/dtdR2 = 0
    # dP1?/dtdC1
    J[6,4] = gamma*e*a_PC*B[6]*R0_P/(B[4]+R0_P)^2 + d*B[5]*(B[6]+B[7])/(B[4]+B[5])^2
    # dP1?/dtdC2
    J[6,5] = -d*B[4]*(B[6]+B[7])/(B[4]+B[5])^2
    # dP1?/dtdP1
    J[6,6] = gamma*e*a_PC*B[4]/(B[4]+R0_P) - m_P - d*B[5]/(B[4]+B[5])
    # dP1?/dtdP2
    J[6,7] = d*B[4]/(B[4]+B[5])
    
    # dP2?/dtdNc = 0
    # dP2?/dtdR1 = 0
    # dP2?/dtdR2 = 0
    # dP2?/dtdC1
    J[7,4] = -d*B[5]*(B[6]+B[7])/(B[4]+B[5])^2
    # dP2?/dtdC2
    J[7,5] = e*a_PC*B[7]*R0_P/(B[5]+R0_P)^2 + d*B[4]*(B[6]+B[7])/(B[4]+B[5])^2
    # dP2?/dtdP1
    J[7,6] = d*B[5]/(B[4]+B[5])
    # dP2?/dtdP2
    J[7,7] = e*a_PC*B[5]/(B[5]+R0_P) - m_P - d*B[4]/(B[4]+B[5])
    
    # dNc?/dtd...
    for(i in 1:7){
      J[1,i] = -sum(J[2:7,i]) # because Nc=Nt-(R1+R2+C1+C2+P1+P2)
    }
    return(J)
  })
}

# equilibrium biomasses
equilibrium<-function(params_data, time, B0, ODE_function, i){
  params<-as.list(c(Nt=params_data$Nt[i],
                    Ns=params_data$Ns[i],
                    l_R=params_data$l_R[i],
                    m_C=params_data$m_C[i],
                    m_P=params_data$m_P[i],
                    r=params_data$r[i],
                    b=params_data$b[i],
                    a_CR=params_data$a_CR[i],
                    a_PC=params_data$a_PC[i],
                    R0_C=params_data$R0_C[i],
                    R0_P=params_data$R0_P[i],
                    e=params_data$e[i],
                    p=params_data$p[i],
                    gamma=params_data$gamma[i]))
  TS<-as.data.frame(ode(B0, time, ODE_function, params, method="rk4")) # simulation of time series
  B<-as.numeric(TS[dim(TS)[1],2:dim(TS)[2]]) # keep the last row
  return(B)
}

# time series
time_series<-function(params_data, time, B0, ODE_function, i){
  params<-as.list(c(Nt=params_data$Nt[i],
                    Ns=params_data$Ns[i],
                    l_R=params_data$l_R[i],
                    m_C=params_data$m_C[i],
                    m_P=params_data$m_P[i],
                    r=params_data$r[i],
                    b=params_data$b[i],
                    a_CR=params_data$a_CR[i],
                    a_PC=params_data$a_PC[i],
                    R0_C=params_data$R0_C[i],
                    R0_P=params_data$R0_P[i],
                    e=params_data$e[i],
                    d=params_data$d[i],
                    p=params_data$p[i],
                    gamma=params_data$gamma[i]))
  TS<-as.data.frame(ode(B0, time, ODE_function, params, method="rk4")) # simulation of time series
  if(length(time)>1000){
    TS<-TS[seq(1,dim(TS)[1],floor(dim(TS)[1]/1000)),] # keeps only 1000 data points for the plot
  }
  return(TS)
}

# return an aggregated dataframe ready to plot
time_series_for_plot<-function(TS){
  TS<-melt(TS,
           id.vars = "time",
           variable.name = "species",
           value.name = "biomass")
  return(TS)
}

# T matrix
get_T_matrix<-function(pert,B,z){
  T_matrix<-diag(B*pert)^z
  return(T_matrix)
}

# Lyapunov equation
lyapunov<-function(J,T_matrix,VE){
  TVT<-T_matrix%*%VE%*%t(T_matrix)
  TVT<-matrix(array(TVT),ncol=1)
  kron<-kronecker(J,diag(rep(1,dim(J)[2]))) + kronecker(diag(rep(1,dim(J)[2])),J)
  return(-solve(kron)%*%TVT)
}

# compute de lead eigen value and the lead eigen vector
analysis<-function(params_data, B, ODE_function, i){
  params<-as.list(c(Nt=params_data$Nt[i],
                    Ns=params_data$Ns[i],
                    l_R=params_data$l_R[i],
                    m_C=params_data$m_C[i],
                    m_P=params_data$m_P[i],
                    r=params_data$r[i],
                    b=params_data$b[i],
                    a_CR=params_data$a_CR[i],
                    a_PC=params_data$a_PC[i],
                    R0_C=params_data$R0_C[i],
                    R0_P=params_data$R0_P[i],
                    e=params_data$e[i],
                    d=params_data$d[i],
                    z=params_data$z[i],
                    sigma=params_data$sigma[i],
                    p=params_data$p[i],
                    gamma=params_data$gamma[i]))
  #J<-jacobian_function(B,params)
  J<-jacobian(fun=ODE_function, x=B, method="simple", method.args=list(), params=params) # numerical approximation
  T_matrix<-get_T_matrix(params_data$pert[[i]],B,params_data$z[i])
  V<-lyapunov(J,T_matrix,params_data$VE[[i]])
  C<-cov2cor(matrix(as.numeric(V),length(B),length(B)))
  CV<-c(sum(diag(matrix(V,length(B)))),sum(V))/sum(B)
  eigen<-eigen(J) # eigen values and vectors
  resilience<--(max(Re(eigen$values)))
  lead<-which(Re(eigen$values)==-resilience) # index of the lead eigen value
  if(length(lead)>1){
    lead<-lead[1] # if we have to identical lead eigenvalues
  }
  E<-abs(Re(eigen(J)$vectors[,lead])) # real part of the lead eigenvector
  E<-E/sum(E) # relative contribution of each species to the lead eigenvector
  return(list(B=B,
              V=as.numeric(V),
              C=as.numeric(C),
              CV=CV,
              resilience=resilience,
              E=E))
}

# Create the output dataframe
create_data_rooney_original<-function(params_data,results){
  # biomass
  data_B<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+5))
  names(data_B)[1:dim(params_data)[2]]=names(params_data)
  data_B[,1:dim(params_data)[2]]=params_data
  B_names<-c("R1","R2","C1","C2","P")
  names(data_B)[dim(params_data)[2]+(1:5)]=B_names
  # variance
  data_V<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+5^2))
  names(data_V)[1:dim(params_data)[2]]=names(params_data)
  data_V[,1:dim(params_data)[2]]=params_data
  V_names<-expand.grid(B_names,B_names)
  V_names<-paste("V",V_names$Var1,V_names$Var2,sep = "_")
  names(data_V)[dim(params_data)[2]+(1:length(V_names))]=V_names
  # correlation
  data_C<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+5^2))
  names(data_C)[1:dim(params_data)[2]]=names(params_data)
  data_C[,1:dim(params_data)[2]]=params_data
  C_names<-expand.grid(B_names,B_names)
  C_names<-paste("C",C_names$Var1,C_names$Var2,sep = "_")
  names(data_C)[dim(params_data)[2]+(1:length(C_names))]=C_names
  # aggregated CV
  data_CV<-params_data
  data_CV$CV_pop=0 # average biomass CV
  data_CV$CV_tot=0 # total biomass CV
  # asymptotic resilience
  data_resilience<-params_data
  data_resilience$resilience=0 # real part of the lead eigne value
  # contribution to the lead eigen vector
  data_E<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+5))
  names(data_E)[1:dim(params_data)[2]]=names(params_data)
  data_E[,1:dim(params_data)[2]]=params_data
  E_names<-B_names
  names(data_E)[dim(params_data)[2]+(1:5)]=E_names
  # fill in tha data frames
  for (i in 1:dim(params_data)[1]){
    data_B[i,which(names(data_B)%in%B_names)]<-unlist(results[[i]]$B)
    data_V[i,which(names(data_V)%in%V_names)]<-unlist(results[[i]]$V)
    data_C[i,which(names(data_C)%in%C_names)]<-unlist(results[[i]]$C)
    data_CV[i,which(names(data_CV)%in%c("CV_pop","CV_tot"))]=unlist(results[[i]]$CV)
    data_resilience[i,which(names(data_resilience)%in%c("resilience"))]=unlist(results[[i]]$resilience)
    data_E[i,which(names(data_E)%in%E_names)]<-unlist(results[[i]]$E)
  }
  return(list(data_B=data_B,
              data_V=data_V,
              data_C=data_C,
              data_CV=data_CV,
              data_resilience=data_resilience,
              data_E=data_E,
              B_names=B_names,
              V_names=V_names,
              C_names=C_names,
              E_names=E_names))
}

# Create the output dataframe
create_data_rooney_dispersal<-function(params_data,results){
  # biomass
  data_B<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+6))
  names(data_B)[1:dim(params_data)[2]]=names(params_data)
  data_B[,1:dim(params_data)[2]]=params_data
  B_names<-c("R1","R2","C1","C2","P1","P2")
  names(data_B)[dim(params_data)[2]+(1:6)]=B_names
  # variance
  data_V<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+6^2))
  names(data_V)[1:dim(params_data)[2]]=names(params_data)
  data_V[,1:dim(params_data)[2]]=params_data
  V_names<-expand.grid(B_names,B_names)
  V_names<-paste("V",V_names$Var1,V_names$Var2,sep = "_")
  names(data_V)[dim(params_data)[2]+(1:length(V_names))]=V_names
  # correlation
  data_C<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+6^2))
  names(data_C)[1:dim(params_data)[2]]=names(params_data)
  data_C[,1:dim(params_data)[2]]=params_data
  C_names<-expand.grid(B_names,B_names)
  C_names<-paste("C",C_names$Var1,C_names$Var2,sep = "_")
  names(data_C)[dim(params_data)[2]+(1:length(C_names))]=C_names
  # aggregated CV
  data_CV<-params_data
  data_CV$CV_pop=0 # average biomass CV
  data_CV$CV_tot=0 # total biomass CV
  # asymptotic resilience
  data_resilience<-params_data
  data_resilience$resilience=0 # real part of the lead eigne value
  # contribution to the lead eigen vector
  data_E<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+6))
  names(data_E)[1:dim(params_data)[2]]=names(params_data)
  data_E[,1:dim(params_data)[2]]=params_data
  E_names<-B_names
  names(data_E)[dim(params_data)[2]+(1:6)]=E_names
  # fill in tha data frames
  for (i in 1:dim(params_data)[1]){
    data_B[i,which(names(data_B)%in%B_names)]<-unlist(results[[i]]$B)
    data_V[i,which(names(data_V)%in%V_names)]<-unlist(results[[i]]$V)
    data_C[i,which(names(data_C)%in%C_names)]<-unlist(results[[i]]$C)
    data_CV[i,which(names(data_CV)%in%c("CV_pop","CV_tot"))]=unlist(results[[i]]$CV)
    data_resilience[i,which(names(data_resilience)%in%c("resilience"))]=unlist(results[[i]]$resilience)
    data_E[i,which(names(data_E)%in%E_names)]<-unlist(results[[i]]$E)
  }
  return(list(data_B=data_B,
              data_V=data_V,
              data_C=data_C,
              data_CV=data_CV,
              data_resilience=data_resilience,
              data_E=data_E,
              B_names=B_names,
              V_names=V_names,
              C_names=C_names,
              E_names=E_names))
}

### PARAMETERS #### ----

# parameters from Rooney et al 2006
params_data_original<-expand.grid(Nt=2.8, # Total available nutrient in common pool
                                  Ns=0.6, # External input of nutrient (value not given in the appendix)
                                  l_R=0.4, # Loss rate from resource pool
                                  m_C=0.2, # Natural mortality of consumer
                                  m_P=0.05, # Natural mortality of the predator
                                  r=1, # Intrinsic growth rate of resource
                                  b=1, # Half saturation rate for resource
                                  a_CR=1.6, # Attack rate of consumer on resource
                                  a_PC=1.8, # Attack rate of predator on consumer
                                  R0_C=1, # Half saturation of consumer
                                  R0_P=1, # Attack rate of predator on consumer
                                  e=0.8, # Biomass conversion efficiency
                                  d=10e3, # dispersal rate of top predators
                                  z=0.5, # demographic perturbations
                                  sigma=1e-3, # standard deviation
                                  p=seq(0,1,0.1), # Proportion of energy the predator derives from channel 1
                                  gamma=seq(0.5,2.5,0.02)) # Fast to Slow ratio (asymmetry coefficient)
B0<-c(0.8,0.2,0.1,0.1,0.2) # initial biomasses

################# ----
### ANALYSIS #### ----
################# ----
# ORIGINAL MODEL #### ----
### Time series (preliminary testing) #### ----
time<-seq(0,2000,0.01) # time vector for ODE simulation
index<-which(params_data_original$p==0.5 & params_data_original$gamma==1)
TS<-time_series(params_data_original, time, B0, ODE_function_rooney_original, index)
TS<-time_series_for_plot(TS)

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_original+
  theme+
  xlab("Time")+
  ylab("Biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.5"))

index<-which(params_data_original$p==0.5 & params_data_original$gamma==2.5)
TS<-time_series(params_data_original, time, B0, ODE_function_rooney_original, index)
TS<-time_series_for_plot(TS)

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_original+
  theme+
  xlab("Time")+
  ylab("Biomass")+
  ggtitle(expression(gamma*"=2.5 - "*italic(p)*"=0.5"))

index<-which(as.character(params_data_original$p)=="0.6" & params_data_original$gamma==1) # 0.6 is not precisely defined and params_data_original$p==0.6 returns FALSE
TS<-time_series(params_data_original, time, B0, index)
TS<-time_series_for_plot(TS)

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_original+
  theme+
  xlab("Time")+
  ylab("Biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.6"))

### Biomass calculation (TO DO ONCE) #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
time<-seq(0,2000,0.01) # time vector for ODE simulation

# simulation for biomass
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% equilibrium(params_data, time, B0, ODE_function_rooney_original, i)
stopCluster(cl)

# extract the biomass
data_B<-as.data.frame(matrix(0,dim(params_data)[1],5))
for(i in 1:dim(params_data)[1]){
  data_B[i,]<-results[[i]]
}
names(data_B)<-c("R1","R2","C1","C2","P")
data_B<-cbind(params_data,data_B)

write.table(data_B,paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

### Perturbation of P #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,0,0,1)) # perturbation of predators
VE<-diag(rep(params_data_original$sigma[1]^2,5)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(4:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_original, i)
stopCluster(cl)
list_data<-create_data_rooney_original(params_data,results)

##### biomass #### ----
data_B<-list_data$data_B
data_B<-melt(data_B[,which(names(data_B)%in%c("gamma",list_data$B_names))],
             id.vars = "gamma",
             variable.name = "species",
             value.name = "biomass")

p1<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=species),size=1.5)+
      theme+
      colour_rooney_original+
      line_rooney_original+
      xlab(label_gamma)+
      ylab("Biomass")

##### asymptotic resilience (lead eigenvalue) #### ----
resilience<-list_data$data_resilience

p2<-ggplot(data=resilience)+
      geom_line(aes(gamma,log10(resilience)),size=1.5)+
      theme+
      xlab(label_gamma)+
      ylab("Asymptotic resilience")

# contribution of each species to the lead eigenvector
E<-list_data$data_E
E<-E[which(names(E)%in%c("gamma",names(E)[1:5+dim(params_data)[2]]))]
E<-melt(E,
        id.vars = "gamma",
        variable.name = "species",
        value.name = "contribution")

p3<-ggplot(data=E)+
      geom_line(aes(gamma,contribution,colour=species,linetype=species),size=1.5)+
      theme+
      colour_rooney_original+
      line_rooney_original+
      xlab(label_gamma)+
      ylab(label_contribution)

##### Figure #### ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p2, 0, 0, 1, 1)+
  draw_plot(p3, 1, 0, 1, 1)+
  draw_plot(p1, 2, 0, 1, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,2), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_resilience_rooney_original.pdf",sep=""), graph, width = 20, height = 7, device = cairo_pdf)

##### correlation #### ----
data_C<-list_data$data_C

p1<-ggplot(data=data_C)+
      geom_line(aes(gamma,C_C1_C2),colour="chartreuse4",size=1.5)+
      theme+
      ylim(0.7,1.01)+
      xlab(label_gamma)+
      ylab(lebel_correlation)+
      ggtitle("perturbation of P")

### Perturbation of C1 #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,1,0,0)) # perturbation of C1
VE<-diag(rep(params_data_original$sigma[1]^2,5)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(4:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_original, i)
stopCluster(cl)
list_data<-create_data_rooney_original(params_data,results)

##### correlation #### ----
data_C1<-list_data$data_C
data_C1$model<-"pert_C1"

### Perturbation of C2 #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,0,1,0)) # perturbation of C2
VE<-diag(rep(params_data_original$sigma[1]^2,5)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(4:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_original, i)
stopCluster(cl)
list_data<-create_data_rooney_original(params_data,results)

##### correlation #### ----
data_C2<-list_data$data_C
data_C2$model<-"pert_C2"

data_C<-rbind(data_C1,data_C2)

p2<-ggplot(data=data_C)+
  geom_line(aes(gamma,C_C1_C2,linetype=model),colour="chartreuse4",size=1.5)+
  theme+
  perturbation_line
legend<-get_legend(p2)

p2<-ggplot(data=data_C)+
      geom_line(aes(gamma,C_C1_C2,linetype=model),colour="chartreuse4",size=1.5)+
      theme+theme(legend.position="none")+
      perturbation_C_line+
      ylim(0.7,1.01)+
      xlab(label_gamma)+
      ylab(label_correlation)+
      ggtitle(expression("perturbation of C"["1"]*" or C"["2"]))

##### Final figure correlation #### ----
graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.02, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"figure_correlation_rooney_original.pdf",sep=""), graph, width = 12, height = 5, device = cairo_pdf)

### Time series #### ----
time<-seq(0,50,0.01) # time vector for ODE simulation
params_data<-params_data_original[which(params_data_original$p==0.5),]
index<-which(params_data$gamma==2)

biomass<-read.table(paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(4:0)]
B0<-as.numeric(biomass[index,]) # equilibrium biomass

# perturbation of predators
B0_pert<-B0*c(1,1,1,1,1.5)
TS<-time_series(params_data, time, B0_pert, ODE_function_rooney_original, index)
for(i in 1:dim(TS)[1]){
  TS[i,c(2:6)]<-TS[i,c(2:6)]/B0
}
TS<-time_series_for_plot(TS)

p1<-ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_original+
  theme+
  xlab("Time")+
  ylab("Scaled biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.5 - perturbation of P"))

# perturbation of C1
B0_pert<-B0*c(1,1,1.5,1,1)
TS<-time_series(params_data, time, B0_pert, ODE_function_rooney_original, index)
for(i in 1:dim(TS)[1]){
  TS[i,c(2:6)]<-TS[i,c(2:6)]/B0
}
TS<-time_series_for_plot(TS)

p2<-ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_original+
  theme+
  xlab("Time")+
  ylab("Scaled biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.5 - perturbation of C"["1"]))

# perturbation of C2
B0_pert<-B0*c(1,1,1,1.5,1)
TS<-time_series(params_data, time, B0_pert, ODE_function_rooney_original, index)
for(i in 1:dim(TS)[1]){
  TS[i,c(2:6)]<-TS[i,c(2:6)]/B0
}
TS<-time_series_for_plot(TS)

p3<-ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_original+
  theme+
  xlab("Time")+
  ylab("Scaled biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.5 - perturbation of C"["2"]))

graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,2), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_TS_rooney_original.pdf",sep=""), graph, width = 20, height = 7, device = cairo_pdf)

# DISPERSAL MODEL #### ----
### Time series (preliminary testing) #### ----
# sets initial biomass after the biomass calculated in the original model
biomass<-read.table(paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(4:0)]
biomass$P1<-biomass$P*biomass$C1/(biomass$C1+biomass$C2)
biomass$P2<-biomass$P*biomass$C2/(biomass$C1+biomass$C2)
biomass$P<-NULL
# sets simulations
time<-seq(0,50,0.0001) # time vector for ODE simulation
params_data<-params_data_original[which(params_data_original$p==0.5),]
index<-which(params_data$p==0.5 & params_data$gamma==2)
B0<-as.numeric(biomass[index,])
TS<-time_series(params_data, time, B0, ODE_function_rooney_dispersal, index)
TS<-time_series_for_plot(TS)

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney_dispersal+
  theme+
  xlab("Time")+
  ylab("Biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.5"))

index<-which(params_data_original$p==0.5 & params_data_original$gamma==2.5)
TS<-time_series(params_data_original, time, B0, ODE_function_rooney_original, index)

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney+
  theme+
  xlab("Time")+
  ylab("Biomass")+
  ggtitle(expression(gamma*"=2.5 - "*italic(p)*"=0.5"))

index<-which(as.character(params_data_original$p)=="0.6" & params_data_original$gamma==1) # 0.6 is not precisely defined and params_data_original$p==0.6 returns FALSE
TS<-time_series(params_data_original, time, B0, index)

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species),size=1.5)+
  colour_rooney+
  theme+
  xlab("Time")+
  ylab("Biomass")+
  ggtitle(expression(gamma*"=1 - "*italic(p)*"=0.6"))

### Biomass calculation (TO DO ONCE) #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
time<-seq(0,2000,0.01) # time vector for ODE simulation

# sets initial biomass after the biomass calculated in the original model
biomass<-read.table(paste(path_data,"biomass_rooney_original.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(4:0)]
biomass$P1<-biomass$P*biomass$C1/(biomass$C1+biomass$C2)
biomass$P2<-biomass$P*biomass$C2/(biomass$C1+biomass$C2)
biomass$P<-NULL

# simulation for biomass
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% equilibrium(params_data, time, B0, ODE_function_rooney_dispersal, i)
stopCluster(cl)

# extract the biomass
data_B<-as.data.frame(matrix(0,dim(params_data)[1],6))
for(i in 1:dim(params_data)[1]){
  data_B[i,]<-results[[i]]
}
names(data_B)<-c("R1","R2","C1","C2","P1","P2")
data_B<-cbind(params_data,data_B)

write.table(data_B,paste(path_data,"biomass_rooney_dispersal.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

### Perturbation of P1 #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,0,0,1,0)) # perturbation of predators in patch #1
VE<-diag(rep(params_data$sigma[1]^2,6)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_dispersal.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(5:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_dispersal, i)
stopCluster(cl)
list_data<-create_data_rooney_dispersal(params_data,results)

##### biomass #### ----
data_B<-list_data$data_B
P_tot<-data_B[,which(names(data_B)%in%c("gamma","P1","P2"))]
data_B<-melt(data_B[,which(names(data_B)%in%c("gamma",list_data$B_names))],
             id.vars = "gamma",
             variable.name = "species",
             value.name = "biomass")

p1<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=species),size=1.5)+
      geom_line(data=P_tot,aes(gamma,P1+P2),colour="darkred",linetype="dashed",size=1.5)+
      annotate(geom='text',label='"P"["tot"]*" = P"["tot"]*" + P"["tot"]',x=1.7,y=0.6,size=7,parse=TRUE,family="Times",colour="darkred")+
      theme+
      colour_rooney_dispersal+
      line_rooney_dispersal+
      xlab(label_gamma)+
      ylab("Biomass")

##### asymptotic resilience (lead eigenvalue) #### ----
resilience<-list_data$data_resilience

p2<-ggplot(data=resilience)+
      geom_line(aes(gamma,log10(resilience)),size=1.5)+
      theme+
      xlab(label_gamma)+
      ylab("Asymptotic resilience")

# contribution of each species to the lead eigenvector
E<-list_data$data_E
E<-E[which(names(E)%in%c("gamma",names(E)[1:6+dim(params_data)[2]]))]
E<-melt(E,
        id.vars = "gamma",
        variable.name = "species",
        value.name = "contribution")

p3<-ggplot(data=E)+
      geom_line(aes(gamma,contribution,colour=species,linetype=species),size=1.5)+
      theme+
      colour_rooney_dispersal+
      line_rooney_dispersal+
      xlab(label_gamma)+
      ylab(label_contribution)

##### Figure #### ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p2, 0, 0, 1, 1)+
  draw_plot(p3, 1, 0, 1, 1)+
  draw_plot(p1, 2, 0, 1, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,2), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_resilience_rooney_dispersal.pdf",sep=""), graph, width = 20, height = 7, device = cairo_pdf)

##### correlation #### ----
data_C1<-list_data$data_C
data_C1$model<-"pert_P1"

### Perturbation of P2 #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,0,0,0,1)) # perturbation of P2
VE<-diag(rep(params_data_original$sigma[1]^2,6)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_dispersal.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(5:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_dispersal, i)
stopCluster(cl)
list_data<-create_data_rooney_dispersal(params_data,results)

##### correlation #### ----
data_C2<-list_data$data_C
data_C2$model<-"pert_P2"

data_C<-rbind(data_C1,data_C2)

p1<-ggplot(data=data_C)+
      geom_line(aes(gamma,C_C1_C2,linetype=model),colour="chartreuse4",size=1.5)+
      theme+
      perturbation_P_line+
      ylim(0.7,1.01)+
      xlab(label_gamma)+
      ylab(label_correlation)+
      ggtitle(expression("perturbation of P"["1"]*" or P"["2"]))

### Perturbation of C1 #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,1,0,0,0)) # perturbation of C1
VE<-diag(rep(params_data_original$sigma[1]^2,6)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_dispersal.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(5:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_dispersal, i)
stopCluster(cl)
list_data<-create_data_rooney_dispersal(params_data,results)

##### correlation #### ----
data_C1<-list_data$data_C
data_C1$model<-"pert_C1"

### Perturbation of C2 #### ----
params_data<-params_data_original[which(params_data_original$p==0.5),]
params_data$pert<-list(c(0,0,0,1,0,0)) # perturbation of C2
VE<-diag(rep(params_data_original$sigma[1]^2,6)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_rooney_dispersal.txt",sep=""),sep=";",header=TRUE)
biomass<-biomass[,dim(biomass)[2]-(5:0)]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% analysis(params_data,as.numeric(biomass[i,]), ODE_rooney_dispersal, i)
stopCluster(cl)
list_data<-create_data_rooney_dispersal(params_data,results)

##### correlation #### ----
data_C2<-list_data$data_C
data_C2$model<-"pert_C2"

data_C<-rbind(data_C1,data_C2)

p2<-ggplot(data=data_C)+
      geom_line(aes(gamma,C_C1_C2,linetype=model),colour="chartreuse4",size=1.5)+
      theme+
      perturbation_C_line+
      ylim(0.7,1.01)+
      xlab(label_gamma)+
      ylab(label_correlation)+
      ggtitle(expression("perturbation of C"["1"]*" or C"["2"]))

##### Final figure correlation #### ----
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"figure_correlation_rooney_dispersal.pdf",sep=""), graph, width = 14, height = 5, device = cairo_pdf)
