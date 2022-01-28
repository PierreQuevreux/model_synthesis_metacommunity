library(ggplot2)
library(reshape2)
library(scales)
#library(pracma)
library(doParallel)
library(deSolve)
library(nleqslv) # root finding (function nleqslv)
library(cowplot)
library(tidyr)
library(viridis)
library(magick)
library(numDeriv)

### PLOT OPTIONS ####

path_figure="Figures/"
path_data="Data/"

theme<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour='grey'),
        panel.grid.major.y = element_line(colour='grey'),
        text = element_text(family="serif",size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5))

# Dispersal model

corr_colour_TL_2<-scale_colour_manual(values=c("dodgerblue3","chocolate1"),
                                      labels=c("1","2"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='trophic\nlevel')
corr_colour_TL_2_TS<-scale_colour_manual(values=c("dodgerblue1","chocolate1","dodgerblue4","chocolate4"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='trophic\nlevel')
patch_line<-scale_linetype_manual(values=c("solid","22"),
                                  name='patch')
gamma_line<-scale_linetype_manual(values=c("solid","22"),
                                  name=expression(atop("asymmetry","coefficient "*gamma)))
corr_line_TL_2_TS<-scale_linetype_manual(values=c("solid","solid","22","22"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='Trophic\nlevel')

perturbation_prey_line<-scale_linetype_manual(values=c("solid","22","dotted"),
                                              labels=c("patch #1","patch #2"),
                                              name='perturbation\nof prey in')

pert_colour_TL_2<-scale_colour_manual(values=c("dodgerblue3","chocolate1"),
                                      labels=c("1","2"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='perturbed\nspecies')
fill_resilience<-scale_fill_viridis(name="Asymptotic\nresilience",
                                    trans = "log10")

x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

# ea and ma factors
ma<-expression(italic(ma))
ma001<-expression(italic(ma)*"=0.01")
ma01<-expression(italic(ma)*"=0.1")
ma1<-expression(italic(ma)*"=1")
ma10<-expression(italic(ma)*"=10")
ma100<-expression(italic(ma)*"=100")
ea<-expression(italic("\u03B5")*italic(a))
ea001<-expression(italic("\u03B5")*italic(a)*"=0.01")
ea01<-expression(italic("\u03B5")*italic(a)*"=0.1")
ea1<-expression(italic("\u03B5")*italic(a)*"=1")
ea10<-expression(italic("\u03B5")*italic(a)*"=10")
ea100<-expression(italic("\u03B5")*italic(a)*"=100")

# labels
label_dispersal<-expression("Scaled dispersal rate "*italic(d["i"]))
label_correlation<-"Correlation between the two patches"
label_correlation_intra<-"Correlation between predator and prey"
label_CV<-"Coefficient of variation (CV)"
label_gamma<-expression("Asymmetry of interaction strenght "*gamma)
label_omega<-expression("Production asymmetry "*omega)
label_contribution<-"Relative contribution to the lead eigenvector"
label_contribution_2lines<-expression(atop("Relative contribution","to the lead eigenvector"))
label_resilience<-"Asymptotic resilience"

### FUNCTIONS ####
# Jacobian single patch J(k)
jacobian_single<-function(B,params,nSpecies,CommunityID){
  with(params,{
    asym<-params$asym[[CommunityID]]
    J<-matrix(0,nrow=nSpecies,ncol=nSpecies)
    J[1,1]=D*(g*asym[1]/D-2*B[1])
    if(nSpecies>1){
      J[1,1]=J[1,1]-D*m*a*asym[2]*B[2] # basal species
      J[nSpecies,nSpecies]=D*m^(nSpecies-1)*(-r/D-2*B[nSpecies]+e*a*asym[nSpecies]*B[nSpecies-1]) # top species
      for(i in 1:(nSpecies-1)){
        J[i,i+1]=-D*m^(i-1)*m*a*asym[i+1]*B[i] # effect of predators
      }
      for(i in 2:nSpecies){
        J[i,i-1]=D*m^(i-1)*e*a*asym[i]*B[i] # effect of prey
      }
    }
    if(nSpecies>2){
      for(i in 2:(nSpecies-1)){
        J[i,i]=D*m^(i-1)*(-r/D-2*B[i]+e*a*asym[i]*B[i-1]-m*a*asym[i+1]*B[i+1]) # intermediate species
      }
    }
    return(J)
  })
}

# Jacobian of the intrapatch dynamics J'
jacobian_intra_patch<-function(B,params,nSpecies,nCommunity){
  Z<-matrix(0,nrow=nSpecies,ncol=nSpecies)
  M<-NULL
  L<-NULL
  J<-NULL
  for(i in 1:nCommunity){
    L<-NULL
    J<-jacobian_single(B[nSpecies*(i-1)+c(1:nSpecies)],params,nSpecies,i)
    for(j in 1:nCommunity){
      if(j==i){
        L<-cbind(L,J)
      }
      else{
        L<-cbind(L,Z)
      }
    }
    M<-rbind(M,L)
  }
  return(M)
}

# Jacobian of the dispersal dynamics P' (linear dispersal components)
jacobian_disp<-function(B,params,nSpecies,nCommunity){
  with(params,{
    dim=nSpecies*nCommunity
    P<-matrix(0,nrow=dim,ncol=dim)
    # effects of species i on itself
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=-disp[i]
          }
          else{
            P[row,col]=disp[i]/(nCommunity-1)
          }
        }
      }
    }
    # Finalising the matrix
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        P[(j-1)*nSpecies+i,]=D*m^(i-1)*P[(j-1)*nSpecies+i,]
      }
    }
    return(P)
  })
}

# Biomasses at equilibrium in the symmetric case
equilibrium_symmetric<-function(params,nSpecies,nCommunity){ # compute the biomasses at equilibrium
  with(params,{
    A<-diag(rep(-1,nSpecies*nCommunity))
    for(j in 1:nCommunity){
      for(i in 2:nSpecies){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i-1]=e*a
      }
      for(i in 1:(nSpecies-1)){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i+1]=-m*a
      }
    }
    B<-matrix(r/D,
              nrow = nSpecies*nCommunity,
              ncol = 1,
              byrow = TRUE)
    for(j in 1:nCommunity){
      B[(j-1)*nSpecies+1,1]=-g/D
    }
    C<-matrix(0,
              nrow = nSpecies*nCommunity,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    return(as.numeric(C))
  })
}

# Biomasses at equilibrium in an isolated community (in the case of asymmetry)
equilibrium_isolated<-function(params,nSpecies,communityID){ # compute the biomasses at equilibrium
  with(params,{
    A<-diag(rep(-1,nSpecies))
    for(j in 1:nCommunity){
      for(i in 2:nSpecies){
        A[i,i-1]=e*a*asym[[communityID]][i]
      }
      for(i in 1:(nSpecies-1)){
        A[i,i+1]=-m*a*asym[[communityID]][i+1]
      }
    }
    B<-matrix(r/D,
              nrow = nSpecies,
              ncol = 1,
              byrow = TRUE)
    for(j in 1:nCommunity){
      B[1,1]=-g*asym[[communityID]][1]/D
    }
    C<-matrix(0,
              nrow = nSpecies,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    return(as.numeric(C))
  })
}

# biomass
equilibrium_ode<-function(params_data,B0,nSpecies,nCommunity,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  time<-seq(0,params_data$t_max[i],params_data$t_step[i])
  TS<-as.data.frame(ode(B0, time, ODE_function, params, method="rk4"))
  B<-as.numeric(TS[dim(TS)[1],2:dim(TS)[2]])
  V<-0#lyapunov(J+P,T,VE,nSpecies,nCommunity)
  C<-0#cov2cor(matrix(as.numeric(V),nSpecies*nCommunity,nSpecies*nCommunity))
  M1<-0#dispersal_importance(B,params,nSpecies,nCommunity)
  CV<-0#c(sum(diag(matrix(V,nSpecies*nCommunity))),sum(V))/sum(B)
  Press<-0
  resilience<-0
  E<-0
  return(list(B=B,
              V=as.numeric(V),
              C=as.numeric(C),
              M1=M1,
              CV=CV,
              Press=as.numeric(Press),
              resilience=resilience,
              E=E))
}

# T matrix
T_matrix<-function(pert,B,z,nSpecies,nCommunity){
  T<-matrix(0,nSpecies*nCommunity,nSpecies*nCommunity)
  coord=NULL
  for(i in 1:length(pert)){
    coord=(pert[[i]][2]-1)*nSpecies+pert[[i]][1] # pert->list of vector containing the trophic level and the patch of the perturbed species (e.g. (1,2) species 1 in patch 2)
    T[coord,coord]=1
  }
  T<-T*diag(B)^z
  return(T)
}

# Lyapunov equation
lyapunov<-function(J,T,VE,nSpecies,nCommunity){
  TVT<-T%*%VE%*%t(T)
  TVT<-matrix(array(TVT),ncol=1)
  kron<-kronecker(J,diag(rep(1,nSpecies*nCommunity))) + kronecker(diag(rep(1,nSpecies*nCommunity)),J)
  return(-solve(kron)%*%TVT)
}

# Parallelised analytic calculation 
analytic<-function(params_data,B,nSpecies,nCommunity,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  J<-jacobian_intra_patch(B,params,nSpecies,nCommunity)
  P<-jacobian_disp(B,params,nSpecies,nCommunity)
  VE<-params_data$VE[[i]]
  T<-T_matrix(params_data$pert[[i]],B,params_data$z[i],nSpecies,nCommunity)
  V<-lyapunov(J+P,T,VE,nSpecies,nCommunity)
  C<-cov2cor(matrix(as.numeric(V),nSpecies*nCommunity,nSpecies*nCommunity))
  M1<-0#dispersal_importance(B,params,nSpecies,nCommunity)
  CV<-c(sum(diag(matrix(V,nSpecies*nCommunity))),sum(V))/sum(B)
  #Press<--solve(diag(B))%*%solve(J+P)%*%diag(B)
  Press<-0
  eigen<-eigen(J+P) # eigen values and eigen vectors
  resilience<-max(abs(Re(eigen$values))) # lead eigen value
  lead<-which(abs(Re(eigen$values))==resilience) # index of the lead eigen value
  E<-abs(Re(eigen$vectors[,lead])) # lead eigen vector
  E<-E/sum(E) # relative contribution of each element
  return(list(B=B,
              V=as.numeric(V),
              C=as.numeric(C),
              M1=M1,
              CV=CV,
              Press=as.numeric(Press),
              resilience=resilience,
              E=E))
}

# Create the output dataframe
create_data<-function(params_data,results,nSpecies,nCommunity){
  # biomass
  data_B<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_B)[1:dim(params_data)[2]]=names(params_data)
  data_B[,1:dim(params_data)[2]]=params_data
  B_names<-expand.grid(c("B"),c(1:nSpecies),c(1:nCommunity))
  B_names$Var1<-paste(B_names$Var1,B_names$Var2,sep = "_")
  B_names<-paste(B_names$Var1,B_names$Var3,sep = "")
  names(data_B)[(dim(params_data)[2]+1):dim(data_B)[2]]=B_names
  # variance
  data_V<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_V)[1:dim(params_data)[2]]=names(params_data)
  data_V[,1:dim(params_data)[2]]=params_data
  V_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  V_names<-paste(V_names$Var1,V_names$Var2,sep = "")
  V_names<-expand.grid(c("V"),V_names,V_names)
  V_names<-paste(V_names$Var1,V_names$Var2,V_names$Var3,sep = "_")
  names(data_V)[(dim(params_data)[2]+1):dim(data_V)[2]]=V_names
  # correlation
  data_C<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_C)[1:dim(params_data)[2]]=names(params_data)
  data_C[,1:dim(params_data)[2]]=params_data
  C_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  C_names<-paste(C_names$Var1,C_names$Var2,sep = "")
  C_names<-expand.grid(c("C"),C_names,C_names)
  C_names<-paste(C_names$Var1,C_names$Var2,C_names$Var3,sep = "_")
  names(data_C)[(dim(params_data)[2]+1):dim(data_C)[2]]=C_names
  # relative importance of dispersal M1
  data_M1<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_M1)[1:dim(params_data)[2]]=names(params_data)
  data_M1[,1:dim(params_data)[2]]=params_data
  M1_names<-expand.grid(c("M"),c(1:nSpecies),c(1:nCommunity))
  M1_names$Var1<-paste(M1_names$Var1,M1_names$Var2,sep = "_")
  M1_names<-paste(M1_names$Var1,M1_names$Var3,sep = "")
  names(data_M1)[(dim(params_data)[2]+1):dim(data_M1)[2]]=M1_names
  # aggregated CV
  data_CV<-params_data
  data_CV$CV_pop=0 # average biomass CV
  data_CV$CV_tot=0 # total biomass CV
  # response to a press perturbation
  data_P<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_P)[1:dim(params_data)[2]]=names(params_data)
  data_P[,1:dim(params_data)[2]]=params_data
  P_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  P_names<-paste(P_names$Var1,P_names$Var2,sep = "")
  P_names<-expand.grid(c("P"),P_names,P_names)
  P_names<-paste(P_names$Var1,P_names$Var2,P_names$Var3,sep = "_")
  names(data_P)[(dim(params_data)[2]+1):dim(data_P)[2]]=P_names
  # asymptotic resilience
  data_resilience<-params_data
  data_resilience$resilience=0 # real part of the lead eigne value
  # contribution to the lead eigen vector
  data_E<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_E)[1:dim(params_data)[2]]=names(params_data)
  data_E[,1:dim(params_data)[2]]=params_data
  E_names<-expand.grid(c("E"),c(1:nSpecies),c(1:nCommunity))
  E_names$Var1<-paste(E_names$Var1,E_names$Var2,sep = "_")
  E_names<-paste(E_names$Var1,E_names$Var3,sep = "")
  names(data_E)[(dim(params_data)[2]+1):dim(data_E)[2]]=E_names
  # fill in tha data frames
  for (i in 1:dim(params_data)[1]){
    data_B[i,which(names(data_B)%in%B_names)]<-unlist(results[[i]]$B)
    data_V[i,which(names(data_V)%in%V_names)]<-unlist(results[[i]]$V)
    data_C[i,which(names(data_C)%in%C_names)]<-unlist(results[[i]]$C)
    data_M1[i,which(names(data_M1)%in%M1_names)]<-unlist(results[[i]]$M1)
    data_CV[i,which(names(data_CV)%in%c("CV_pop","CV_tot"))]=unlist(results[[i]]$CV)
    data_P[i,which(names(data_P)%in%P_names)]<-unlist(results[[i]]$Press)
    data_resilience[i,which(names(data_resilience)%in%c("resilience"))]=unlist(results[[i]]$resilience)
    data_E[i,which(names(data_E)%in%E_names)]<-unlist(results[[i]]$E)
  }
  return(list(data_B=data_B,
              data_V=data_V,
              data_C=data_C,
              data_M1=data_M1,
              data_CV=data_CV,
              data_P=data_P,
              data_resilience=data_resilience,
              data_E=data_E,
              B_names=B_names,
              V_names=V_names,
              C_names=C_names,
              M1_names=M1_names,
              P_names=P_names,
              E_names=E_names))
}

# ODE of the system for the numerical simulation of time series
ODE_function_TS<-function(t, B, params){
  with(params,{
    dB<-rep(0,nSpecies*nCommunity)
    # intra-patch dynamics
    ID=0
    for(j in 1:nCommunity){
      ID=(j-1)*nSpecies
      dB[ID+1]=D*B[ID+1]*(g*asym[[j]][1]/D - B[ID+1])
      if(nSpecies>1){
        for(i in 2:nSpecies){
          dB[ID+i]=D*m^(i-1)*B[ID+i]*(-r/D - B[ID+i] + e*a*asym[[j]][i]*B[ID+i-1])
        }
        for(i in 1:(nSpecies-1)){
          dB[ID+i]=dB[ID+i]+D*m^(i-1)*B[ID+i]*(- m*a*asym[[j]][i+1]*B[ID+i+1])
        }
      }
    }
    # dispersal dynamics
    d_matrix<-matrix(0,nSpecies,nCommunity) # matrix containing the dispersal terms
    B_matrix<-matrix(B,nSpecies,nCommunity) # matrix containg the biomasses
    # final equations
    for(i in 1:nSpecies){
      d_matrix[i,]=D*m^(i-1)*disp[i]*B_matrix[i,]
    }
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        dB[(j-1)*nSpecies+i]=dB[(j-1)*nSpecies+i]-(1+1/(nCommunity-1))*d_matrix[i,j]+sum(d_matrix[i,])/(nCommunity-1)
      }
    }
    return(list(dB))
  })
}

# Make a table to plot a correlation matrix
table_for_matrix<-function(table,nparams){
  table<-melt(table,
              id.vars = names(table)[1:nparams],
              variable.name = "species",
              value.name = "correlation")
  table<-table %>% separate(species,c(NA,"species_1","species_2"),sep="_")
  table$species_1<-as.factor(table$species_1)
  table$species_2<-as.factor(table$species_2)
  table<-table %>% separate(species_1,c("species_1","community_1"),sep=1)
  table<-table %>% separate(species_2,c("species_2","community_2"),sep=1)
  return(table)
}

# Make a table ready to use for ggplot
table_for_plot<-function(table,nparams,value.name){
  table<-melt(table,
              id.vars = names(table)[1:nparams],
              variable.name = "species",
              value.name = value.name)
  table<-table %>% separate(species,into=c(NA,NA,"species","community"),sep=c(1,2,3))
  return(table)
}

# set asymmetry for two patches
set_asymmetry<-function(nSpecies,gamma,omega){
  asym_1<-c(omega,rep(gamma,nSpecies-1)) # patch with modified growth rate and attack rate
  asym_2<-rep(1,nSpecies) # reference patch
  return(list(list(asym_1,asym_2)))
}

# time series
time_series<-function(params_data,B0,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  time<-seq(0,params_data$t_max[i],params_data$t_step[i])
  TS<-as.data.frame(rk(B0, time, ODE_function_TS, params, method="rk45dp7"))
  return(TS)
}

# return an aggregated dataframe ready to plot
time_series_for_plot<-function(TS){
  if(length(time)>1000){
    TS<-TS[seq(0,length(time),floor(length(time)/1000)),] # keeps only 1000 points
  }
  TS<-melt(TS,
           id.vars = "time",
           variable.name = "species",
           value.name = "biomass")
  TS$community=1
  TS$community[which(TS$species%in%c(3,4))]=2
  TS$species[which(TS$species%in%c(1,3))]=1
  TS$species[which(TS$species%in%c(2,4))]=2
  TS$species<-as.factor(TS$species)
  TS$community<-as.factor(TS$community)
  return(TS)
}

# time series after a pulse perturbation
time_series_perturbation<-function(params_data,B0,nSpecies,nCommunity,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  time<-seq(0,params_data$t_max[i],params_data$t_step[i])
  # perturbation of predators in patch #1
  B<-B0
  coord=NULL
  pert<-params_data$pert[[i]] # selects the list of perturbed species
  for(j in 1:length(pert)){
    coord=(pert[[j]][2]-1)*nSpecies+pert[[j]][1] # pert->list of vector containing the trophic level and the patch of the perturbed species (e.g. (1,2) species 1 in patch 2)
    B[coord]=B[coord]*params_data$pert_factor[i] # generate the pulse perturbation
  }
  TS<-as.data.frame(rk(B, time, ODE_function_TS, params, method="rk45dp7"))
  for(j in 1:(nSpecies*nCommunity)){
    TS[,j+1]<-TS[,j+1]/B0[j] # rescale the biomasses by their value at equilibrium
  }
  if(length(time)>1000){
    TS<-TS[seq(0,length(time),floor(length(time)/1000)),] # keeps only 1000 points
  }
  TS<-melt(TS,
           id.vars = "time",
           variable.name = "species",
           value.name = "biomass")
  #levels(TS$species)<-c(1:4)
  TS$community=1
  TS$community[which(TS$species%in%c(3,4))]=2
  TS$species[which(TS$species%in%c(1,3))]=1
  TS$species[which(TS$species%in%c(2,4))]=2
  TS$species<-as.factor(TS$species)
  TS$community<-as.factor(TS$community)
  TS$ea=params_data$ea[i]
  TS$ma=params_data$ma[i]
  TS$gamma=params_data$gamma[i]
  TS$omega=params_data$omega[i]
  TS$d=params_data$d[i]
  TS$model=params_data$model[i]
  return(TS)
}

# get the params list
get_params<-function(params_data,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i],
                    gamma=params_data$gamma[i],
                    omega=params_data$omega[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  return(params)
}

### PARAMETERS ####
d_min=-5
d_max=5
d_step=0.1
d_interval=10^(seq(d_min,d_max,d_step))

g=1
r=0
D=1
e=0.65
m=c(0.0065,0.065,0.65,6.5,65)
a=c(1/6.5,1/0.65,1/0.065)
sigma=1e-3
z=0.5

pert=list(c(1,1)) # c(species, patch)
disp=list(c(1,1)) # c(species1, species2,..., nSpecies)
asym=list(c(1.5,1.5),c(1,1)) # asymmetry of attack rate in patch #1 and #2

params_data_original<-expand.grid(simu_ID=0,
                                  g=g,
                                  r=r,
                                  D=D,
                                  e=e,
                                  m=m,
                                  a=a,
                                  sigma=sigma,
                                  z=z)
params_data_original$ma=params_data_original$m*params_data_original$a
params_data_original<-params_data_original[params_data_original$ma>0.05 & params_data_original$ma<=15,]
params_data_original$ea=params_data_original$e*params_data_original$a
params_data_original$ea<-as.factor(params_data_original$ea)
params_data_original$ma<-as.factor(params_data_original$ma)
levels(params_data_original$ma)<-c(ma01,ma1,ma10)
params_data_original$ma = factor(params_data_original$ma,levels(params_data_original$ma)[c(3,2,1)])
levels(params_data_original$ea)<-c(ea01,ea1,ea10)

########################## ----
# MAIN TEXT ############## ----
########################## ----
# preliminary test for the ODE integration #### ----
nSpecies=2
nCommunity=2
gamma=seq(1,5) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e3 # high dispersal

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         model="pert_11")
params_data$omega=params_data$gamma
params_data<-merge(params_data_original,params_data)
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$lambda<-params_data$e*params_data$a^2*params_data$m
params_data$pert_factor=1
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# initial biomasses computed analytically before ODE resolution (for d>>1, B_21=B_22)
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity) # B_11 B_21 B_12 B_22
B0[,1]=params_data$g/params_data$D*(2*params_data$omega+params_data$lambda*params_data$omega*(params_data$gamma^2+1)-params_data$lambda*params_data$gamma*(1+params_data$omega*params_data$gamma))/(2+params_data$lambda*(params_data$gamma^2+1))
B0[,2]=2*params_data$e*params_data$a*params_data$g/params_data$D*(1+params_data$omega*params_data$gamma)/(2+params_data$lambda*(params_data$gamma^2+1))/2
B0[,3]=params_data$g/params_data$D*(2+params_data$lambda*(params_data$gamma^2+1)-params_data$lambda*(1+params_data$omega*params_data$gamma))/(2+params_data$lambda*(params_data$gamma^2+1))
B0[,4]=B0[,2]
B0[,1][B0[,1]<1e-6]=NA
B0[,3][B0[,3]<1e-6]=NA

# integration time
params_data$t_max=10000
params_data$t_step=1

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve","reshape2")) %dopar% time_series_perturbation(params_data,B0[i,],nSpecies,nCommunity,i)#time_series(params_data,B0[i,],i)
stopCluster(cl)

#TS<-results[[1]]
#TS<-time_series_for_plot(TS)

TS<-NULL
for(i in 1:dim(params_data)[1]){
  databis<-results[[i]]
  TS<-rbind(TS,databis)
}

ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_grid(~gamma, labeller=label_parsed)+
  corr_colour_TL_2+
  patch_line+
  theme+
  xlab("Time")+
  ylab("Biomass")

# biomass calculation - gamma - omega (TO DO ONCE) #### ----
nSpecies=2
nCommunity=2
gamma=seq(1,5,0.2) # asymmetry coefficient
d=1e3 # high dispersal

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         model="pert_11")
params_data$omega=params_data$gamma
params_data<-merge(params_data_original,params_data)
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$lambda<-params_data$e*params_data$a^2*params_data$m
params_data$simu_ID<-seq(1,dim(params_data)[1])
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution (for d>>1, B_21=B_22)
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity) # B_11 B_21 B_12 B_22
B0[,1]=params_data$g/params_data$D*(2*params_data$omega+params_data$lambda*params_data$omega*(params_data$gamma^2+1)-params_data$lambda*params_data$gamma*(1+params_data$omega*params_data$gamma))/(2+params_data$lambda*(params_data$gamma^2+1))
B0[,2]=2*params_data$e*params_data$a*params_data$g/params_data$D*(1+params_data$omega*params_data$gamma)/(2+params_data$lambda*(params_data$gamma^2+1))/2
B0[,3]=params_data$g/params_data$D*(2+params_data$lambda*(params_data$gamma^2+1)-params_data$lambda*(1+params_data$omega*params_data$gamma))/(2+params_data$lambda*(params_data$gamma^2+1))
B0[,4]=B0[,2]
B0[,1][B0[,1]<1e-6]=NA
B0[,3][B0[,3]<1e-6]=NA

# integration time
params_data$t_max=3000
params_data$t_max[params_data$gamma>2]=5000
params_data$t_max[params_data$gamma>3]=10000
params_data$t_max[params_data$gamma>4]=15000
params_data$t_step=1

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve")) %dopar% time_series(params_data,B0[i,],i)
stopCluster(cl)

data_B<-NULL
for(i in 1:dim(params_data)[1]){
  TS<-results[[i]]
  data_B<-rbind(data_B,TS[dim(TS)[1],2:dim(TS)[2]])
}
data_B<-cbind(params_data,data_B)
data_B<-data_B[,-which(names(data_B)%in%c("pert","asym","disp"))]

write.table(data_B,paste(path_data,"biomass_asymmetry.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

# parameters and simulations - gamma - omega #### ----
nSpecies=2
nCommunity=2
gamma=seq(1,5,0.2) # asymmetry coefficient
d=1e3 # high dispersal

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         model="pert_11")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(1,2))), # perturbation of prey in patch 2
                                           disp=list(c(0,1)), # dispersal of predators
                                           gamma=gamma,
                                           d=d,
                                           model="pert_12"))
params_data$omega=params_data$gamma
params_data<-merge(params_data_original,params_data)
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

biomass<-read.table(paste(path_data,"biomass_asymmetry.txt",sep=""),sep=";",header=TRUE)
biomass<-as.matrix(biomass[,dim(biomass)[2]-seq(3,0)])
biomass<-rbind(biomass,biomass)

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analytic(params_data,biomass[i,],nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,results,nSpecies,nCommunity)

## correlation inter-patch #### ----
data_C<-list_data$data_C
databis<-data_C[,which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22'))]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
      perturbation_prey_line+
      corr_colour_TL_2+
      theme+
      ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation)

ggsave(paste(path_figure,"figure_asymmetry.pdf",sep=""),p1, width = 7, height = 5, device = cairo_pdf)

## biomass #### ----
# biomasses without dispersal
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies*nCommunity)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:(nSpecies*nCommunity))

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  biomass[i,c(1,2)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,1)
  biomass[i,c(3,4)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,2)
}
biomass<-biomass[biomass$model=="pert_11",]

# raw biomass
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      y_axis_log10+
      xlab(label_gamma)+
      ylab("Biomass")

# scaled biomass
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B[,dim(data_B)[2]-(3:0)]<-data_B[,dim(data_B)[2]-(3:0)]/biomass[,dim(biomass)[2]-(3:0)] # rescales biomass by the biomass without dispersal
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p2<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      xlab(label_gamma)+
      ylab("Scaled biomass")

## time series ----
nSpecies=2
nCommunity=2
gamma=c(1,3) # asymmetry coefficient
d=1e3

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         model="pert_11")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(1,2))), # perturbation of prey in patch 2
                                           disp=list(c(0,1)), # dispersal of predators
                                           gamma=gamma,
                                           d=d,
                                           model="pert_12"))
params_data<-merge(params_data_original,params_data)
params_data$omega=params_data$gamma
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$lambda<-params_data$e*params_data$a^2*params_data$m
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=list(list(rep(params_data$gamma[i],nSpecies), # patch with modified attack rate
                                rep(1,nSpecies))) # reference patch
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# initial biomasses computed analytically before ODE resolution
B0<-read.table(paste(path_data,"biomass_asymmetry.txt",sep=""),sep=";",header=TRUE)
B0<-B0[B0$gamma==1 | B0$gamma==3,]
B0<-B0[,dim(B0)[2]+c(-3:0)]
names(B0)<-c(1:4)
B0<-as.matrix(B0) # B_11 B_21 B_12 B_2

# integration time and pulse perturbations depending on parameters
params_data$t_max=300
params_data$t_step=1
params_data$pert_factor=1.1

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve","reshape2")) %dopar% time_series_perturbation(params_data,B0[i,],nSpecies,nCommunity,i)
#results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve","reshape2")) %dopar% time_series(params_data,as.numeric(biomass[i,]),nSpecies,nCommunity,i)
stopCluster(cl)

TS<-NULL
for(i in 1:dim(params_data)[1]){
  databis<-results[[i]]
  TS<-rbind(TS,databis)
}
levels(TS$model)<-c(expression("prey perturbed"*" in patch #1"),expression("prey perturbed"*" in patch #2"))
TS$gamma<-as.factor(TS$gamma)
levels(TS$gamma)<-c(expression(gamma*" = 1"),expression(gamma*" = 3"))

p3<-ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p3)

p3<-ggplot(data=TS)+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(gamma~model, labeller=label_parsed)+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      xlab("Time")+
      ylab("Scaled biomass")

# final figure ----
graph<-ggdraw(xlim = c(0, 3.45), ylim = c(0, 1.55)) +
  draw_plot(p1, 0.05, 0.75, 1, 0.75)+
  draw_plot(p2, 0.05, 0, 1, 0.75)+
  draw_plot(p3, 1.2, 0, 2, 1.5)+
  draw_plot(legend, 3.3, 0.5, 0.05, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,0,1.1), c(1.55,0.85,1.55), size = 30)
ggsave(paste(path_figure,"supp_asymmetry.pdf",sep=""), graph, width = 14, height = 6, device = cairo_pdf)
