# R version 3.6.3 (2020-02-29)
library(ggplot2) # ggplot2_3.3.5
library(reshape2) # reshape2_1.4.4
library(scales) # scales_1.1.1
#library(pracma)
library(doParallel) # doParallel_1.0.16 iterators_1.0.13 foreach_1.5.1
library(deSolve) # deSolve_1.30
library(nleqslv) # root finding (function nleqslv) nleqslv_3.3.2
library(cowplot) # cowplot_1.1.1
library(tidyr) # tidyr_1.1.4
library(viridis) # viridis_0.6.2 viridisLite_0.4.0
library(magick) # magick_2.7.3
library(numDeriv) # to evaluate the Jacobian matrix numDeriv_2016.8-1.1

#sessionInfo()
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.7       pillar_1.6.4     compiler_3.6.3   plyr_1.8.6       tools_3.6.3      lifecycle_1.0.1  tibble_3.1.6     gtable_0.3.0    
# [9] pkgconfig_2.0.3  rlang_0.4.12     DBI_1.1.2        gridExtra_2.3    withr_2.4.3      dplyr_1.0.7      stringr_1.4.0    generics_0.1.1  
# [17] vctrs_0.3.8      grid_3.6.3       tidyselect_1.1.1 glue_1.6.0       R6_2.5.1         fansi_0.5.0      purrr_0.3.4      magrittr_2.0.1  
# [25] codetools_0.2-16 ellipsis_0.3.2   assertthat_0.2.1 colorspace_2.0-2 utf8_1.2.2       stringi_1.7.6    munsell_0.5.0    crayon_1.4.2    

### PLOT OPTIONS ####

path_figure="Figures/"
path_data="Data/"

# names of dataframe columns
table_names<-list(B_names=c("B_11","B_21","B_12","B_22"),
                  V_names=c("V_11_11","V_21_11","V_12_11","V_22_11","V_11_21","V_21_21","V_12_21","V_22_21","V_11_12","V_21_12","V_12_12","V_22_12","V_11_22","V_21_22","V_12_22","V_22_22"),
                  C_names=c("C_11_11","C_21_11","C_12_11","C_22_11","C_11_21","C_21_21","C_12_21","C_22_21","C_11_12","C_21_12","C_12_12","C_22_12","C_11_22","C_21_22","C_12_22","C_22_22"),
                  E_names=c("E_11","E_21","E_12","E_22"))

theme<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour='grey'),
        panel.grid.major.y = element_line(colour='grey'),
        text = element_text(size=20,family="Times"),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5))

# Dispersal model

corr_colour_TL_4<-scale_colour_manual(values=c("dodgerblue3","chocolate1","chartreuse4","red"),
                                      labels=c("1","2","3","4"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='trophic\nlevel')
corr_colour_top<-scale_colour_manual(values=c("chartreuse4","red"),
                                     labels=c("3","4"),
                                     guide = guide_legend(reverse = TRUE),
                                     name='trophic\nlevel')
patch_line<-scale_linetype_manual(values=c("solid","22"),
                                  name='patch')
patch_line_pert<-scale_linetype_manual(values=c("solid","22"),
                                       name='species\nperturbed\nin patch')
fill_colour_TL_4<-scale_fill_manual(values=c("dodgerblue3","chocolate1","chartreuse4","red"),
                                    labels=c("1","2","3","4"),
                                    guide = guide_legend(reverse = TRUE),
                                    name='Trophic\nlevel')
corr_colour_grad<-scale_fill_gradient2(low = "red",
                                       mid = "white",
                                       high = "blue",
                                       midpoint = 0,
                                       limits = c(-1,1),
                                       name="Correlation\ncoefficient")

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
label_CV<-"Coefficient of variation (CV)"
label_gamma<-expression("Asymmetry of interaction strength "*gamma)
label_omega<-expression("Asymmetry of biomass production "*omega)
label_contribution<-"Relative contribution to the lead eigenvector"
label_resilience<-"Asymptotic resilience"

### FUNCTIONS ####

# set asymmetry for two patches
set_asymmetry<-function(nSpecies,gamma,omega){
  asym_1<-c(omega,rep(gamma,nSpecies-1)) # patch with modified growth rate and attack rate
  asym_2<-rep(1,nSpecies) # reference patch
  return(list(list(asym_1,asym_2)))
}

# ODE of the system
ODE_model<-function(B, params){
  with(params,{
    dB<-rep(0,8)
    # dB1/dt patch #1
    dB[1] = D*(B[1]*(asym[[1]][1]*g/D - B[1] - asym[[1]][2]*m*a*B[2]) + disp[1]*(B[5] - B[1]))
    # dB2/dt patch #1
    dB[2] = m*D*(B[2]*(-r/D - B[2] + asym[[1]][2]*e*a*B[1] - asym[[1]][3]*m*a*B[3]) + disp[2]*0)
    # dB3/dt patch #1
    dB[3] = m*D*(B[3]*(-r/D - B[3] + asym[[1]][3]*e*a*B[2] - asym[[1]][4]*m*a*B[4]) + disp[3]*0)
    # dB4/dt patch #1
    dB[4] = m*D*(B[4]*(-r/D - B[4] + asym[[1]][4]*e*a*B[3]) + disp[4]*(B[8]*S_prey/(B[7]+S_prey) - B[4]*S_prey/(B[3]+S_prey)))
    # dB1/dt patch #1
    dB[5] = D*(B[5]*(asym[[2]][1]*g/D - B[5] - asym[[2]][2]*m*a*B[6]) + disp[1]*(B[1] - B[5]))
    # dB2/dt patch #1
    dB[6] = m*D*(B[6]*(-r/D - B[6] + asym[[2]][2]*e*a*B[5] - asym[[2]][3]*m*a*B[7]) + disp[2]*0)
    # dB3/dt patch #1
    dB[7] = m*D*(B[7]*(-r/D - B[7] + asym[[2]][3]*e*a*B[6] - asym[[2]][4]*m*a*B[8]) + disp[3]*0)
    # dB4/dt patch #1
    dB[8] = m*D*(B[8]*(-r/D - B[8] + asym[[2]][4]*e*a*B[7]) + disp[4]*(B[4]*S_prey/(B[3]+S_prey) - B[8]*S_prey/(B[7]+S_prey)))
    return(dB)
  })
}

# ODE function for the ODE solver
ODE_function<-function(t, B, params){
  dB<-ODE_model(B, params)
  return(list(dB))
}

# Jacobian matrix
jacobian_model<-function(B,params){
  with(params,{
    J<-matrix(0,nrow=4,ncol=4)
    # dB11?/dtdb11
    J[1,1] = D*(asym[[1]][1]*g/D - 2*B[1] - asym[[1]][2]*m*a*B[2] - disp[1]*B[2]/(B[2]+S_pred))
    # dB11?/dtdb21
    J[1,2] = D*(- asym[[1]][2]*m*a*B[1] - disp[1]*B[1]*S_pred/(B[2]+S_pred)^2)
    # dB11?/dtdb12
    J[1,3] = D*disp[1]*B[4]/(B[4]+S_pred)
    # dB11?/dtdb22
    J[1,4] = D*disp[1]*B[3]*S_pred/(B[4]+S_pred)^2
    
    # dB21?/dtdb11
    J[2,1] = m*D*(asym[[1]][2]*e*a*B[2] + disp[2]*B[2]*S_prey/(B[1]+S_prey)^2)
    # dB21?/dtdb21
    J[2,2] = m*D*(-r/D - 2*B[2] + asym[[1]][2]*e*a*B[1] - disp[2]*S_prey/(B[1]+S_prey))
    # dB21?/dtdb12
    J[2,3] = -m*D*disp[2]*B[4]*S_prey/(B[3]+S_prey)^2
    # dB21?/dtdb22
    J[2,4] = m*D*disp[2]*S_prey/(B[3]+S_prey)
    
    # dB12?/dtdb11
    J[3,1] = D*disp[1]*B[2]/(B[2]+S_pred)
    # dB12?/dtdb21
    J[3,2] = D*disp[1]*B[1]*S_pred/(B[2]+S_pred)^2
    # dB12?/dtdb12
    J[3,3] = D*(asym[[2]][1]*g/D - 2*B[3] - asym[[2]][2]*m*a*B[4] - disp[1]*B[4]/(B[4]+S_pred))
    # dB12?/dtdb22
    J[3,4] = D*(- asym[[2]][2]*m*a*B[3] - disp[1]*B[3]*S_pred/(B[4]+S_pred)^2)
    
    # dB22?/dtdb11
    J[4,1] = -m*D*disp[2]*B[2]*S_prey/(B[1]+S_prey)^2
    # dB22?/dtdb21
    J[4,2] = m*D*disp[2]*S_prey/(B[1]+S_prey)
    # dB22?/dtdb12
    J[4,3] = m*D*(asym[[2]][2]*e*a*B[4] + disp[2]*B[4]*S_prey/(B[3]+S_prey)^2)
    # dB22?/dtdb22
    J[4,4] = m*D*(-r/D - 2*B[4] + asym[[2]][2]*e*a*B[3] - disp[2]*S_prey/(B[3]+S_prey))
    
    return(J)
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

# Biomass
equilibrium_ode<-function(params_data,B0,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i],
                    S_pred=params_data$S_pred[i],
                    S_prey=params_data$S_prey[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  time<-seq(0,params_data$t_max[i],params_data$t_step[i])
  TS<-as.data.frame(ode(B0, time, ODE_function, params, method="rk4"))
  B<-TS[dim(TS)[1],2:dim(TS)[2]]
  names(B)<-c("1_1","2_1","3_1","4_1","1_2","2_2","3_2","4_2")
  B<-cbind(params_data[i,],B)
  return(B)
}

# Parallelised function solving the ODE
time_series<-function(params_data, B0, i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i],
                    S_pred=params_data$S_pred[i],
                    S_prey=params_data$S_prey[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  time<-seq(0,params_data$t_max[i],params_data$t_step[i]) # time vector for ODE simulation
  TS<-as.data.frame(ode(B0, time, ODE_function, params, method="rk4"))
  if(length(time)>1000){
    TS<-TS[seq(1,dim(TS)[1],floor(dim(TS)[1]/1000)),] # keeps only 1000 data points for the plot
  }
  return(TS)
}

# Sets time series data frame ready for plot
time_series_for_plot<-function(params_data, TS, i){
  TS<-melt(TS,
           id.vars = "time",
           variable.name = "species",
           value.name = "biomass")
  TS<-merge(params_data[i,],TS)
  TS$community=1
  TS$community[which(TS$species%in%c(5:8))]=2
  TS$species[which(TS$species%in%c(1,5))]=1
  TS$species[which(TS$species%in%c(2,6))]=2
  TS$species[which(TS$species%in%c(3,7))]=3
  TS$species[which(TS$species%in%c(4,8))]=4
  TS$species<-as.factor(TS$species)
  TS$community<-as.factor(TS$community)
  return(TS)
}

# Parallelised analytic calculation 
linear_analysis<-function(params_data,biomass,nSpecies,nCommunity,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i],
                    S_pred=params_data$S_pred[i],
                    S_prey=params_data$S_prey[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  B<-as.numeric(biomass[i,])
  #J<-jacobian_model(B,params)
  J<-jacobian(fun=ODE_model, x=B, method="simple", method.args=list(), params=params) # numerical approximation
  VE<-params_data$VE[[i]]
  T<-T_matrix(params_data$pert[[i]],B,params_data$z[i],nSpecies,nCommunity)
  V<-lyapunov(J,T,VE,nSpecies,nCommunity)
  C<-cov2cor(matrix(as.numeric(V),nSpecies*nCommunity,nSpecies*nCommunity))
  CV<-c(sum(diag(matrix(V,nSpecies*nCommunity))),sum(V))/sum(B)
  eigen<-eigen(J) # eigen values and eigen vectors
  resilience<--max((Re(eigen$values))) # lead eigen value
  lead<-which(-(Re(eigen$values))==resilience) # index of the lead eigen value
  E<-abs(Re(eigen$vectors[,lead])) # lead eigen vector
  E<-E/sum(E) # relative contribution of each element
  return(list(B=B,
              V=as.numeric(V),
              C=as.numeric(C),
              CV=CV,
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
  # aggregated CV
  data_CV<-params_data
  data_CV$CV_pop=0 # average biomass CV
  data_CV$CV_tot=0 # total biomass CV
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
# biomasses in an isolated food chain #### ----
nSpecies=4
nCommunity=1
gamma=seq(0.5,5,0.5)
omega=c(1)
params_data<-expand.grid(gamma=gamma,
                         omega=omega)
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}

data_B<-params_data
data_B$'1'=0
data_B$'2'=0
data_B$'3'=0
data_B$'4'=0
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  data_B[i,which(names(data_B)%in%c("1","2","3","4"))]<-equilibrium_isolated(params,nSpecies,1)
}
data_B$omega<-as.factor(data_B$omega)
levels(data_B$omega)<-c(expression(omega*"=1"),
                        expression(omega*"=5"),
                        expression(omega*"=10"))

data_B<-melt(data_B[,which(names(data_B)%in%c("gamma","omega","1","2","3","4"))],
             id.vars = c("gamma","omega"),
             variable.name = "species",
             value.name = "biomass")

ggplot(data=data_B)+
  geom_line(aes(gamma,biomass,colour=species),size=1.5)+
  facet_grid(.~omega, labeller=label_parsed)+
  corr_colour_TL_4+
  theme +
  y_axis_log10+
  xlab(label_gamma)+
  ylab("Biomass")

### FULL MODEL #### ----
# preliminary test for the ODE integration #### ----
nSpecies=4
nCommunity=2
gamma=2 # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=500
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data,B0[i,],i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  data_TS<-rbind(data_TS,time_series_for_plot(params_data, results[[i]], i))
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_grid(ma~ea, labeller=label_parsed)+
  corr_colour_TL_4+
  patch_line+
  theme+
  y_axis_log10_short+
  xlab("Time")+
  ylab("Biomass")

# biomass calculation - gamma - omega (TO DO ONCE) #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e3, # high sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=500
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve")) %dopar% equilibrium_ode(params_data,B0[i,],i)
stopCluster(cl)

data_B<-NULL
for(i in 1:dim(params_data)[1]){
  data_B<-rbind(data_B,results[[i]])
}
data_B<-data_B[,-which(names(data_B)%in%c("pert","asym","disp"))]

write.table(data_B,paste(path_data,"biomass_disp14_full_main.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

# analysis - gamma - omega - load biomasses #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(4,1))), # perturbations of top predators in patch 1
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1,
                         model="pert_41")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]

data_B<-read.table(paste(path_data,"biomass_disp14_full_main.txt",sep=""),sep=";",header=TRUE)
biomass<-na.omit(data_B) # remove simulations without coexistance
params_data<-params_data[biomass$simu_ID,] # remove simulations without coexistance
biomass<-biomass[,-((nSpecies*nCommunity-1):0)+dim(biomass)[2]] # selects the biomasses

# adds perturbations of other species
params_data_42<-params_data
params_data_42$pert=list(list(c(4,2))) # perturbations of top predator in patch 2
params_data_42$model="pert_42"
params_data_31<-params_data
params_data_31$pert=list(list(c(3,1))) # perturbations of predators in patch 1
params_data_31$model="pert_31"
params_data_32<-params_data
params_data_32$pert=list(list(c(3,2))) # perturbations of predators in patch 2
params_data_32$model="pert_32"
params_data<-rbind(params_data,params_data_42,params_data_31,params_data_32)
rm(params_data_42,params_data_31,params_data_32)
biomass<-rbind(biomass,biomass,biomass,biomass)

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% linear_analysis(params_data,biomass,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,results,nSpecies,nCommunity)

## biomass #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B[data_B$model=="pert_31",])+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme
legend_1<-get_legend(p1)

p1<-ggplot(data=data_B[data_B$model=="pert_31",])+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      corr_colour_TL_4+
      patch_line+
      theme + theme(legend.position="none")+
      y_axis_log10+
      xlab(label_gamma)+
      ylab("Biomass")

## biomass scaled by biomass without dispersal #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
# reference biomass without dispersal
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
for (i in 1:dim(data_B)[1]){
  data_B[i,4:dim(data_B)[2]]=data_B[i,4:dim(data_B)[2]]/B0[i,]
}
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p2<-ggplot(data=data_B[data_B$model=="pert_31",])+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      geom_hline(yintercept=1,linetype='dashed',size=1)+
      annotate(geom='text',label="biomass higher\nwith dispersal",x=3.5, y=1.25,size=7,family="Times")+
      annotate(geom='text',label="biomass lower\nwith dispersal",x=3.5, y=0.75,size=7,family="Times")+
      corr_colour_TL_4+
      patch_line+
      theme + theme(legend.position="none")+
      #y_axis_log10+
      xlab(label_gamma)+
      ylab("Scaled biomass")

## final graph #### ----
graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1.05, 0, 1, 1)+
  draw_plot(legend_1, 2.05, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"figure_biomass_disp_14.pdf",sep=""), graph, width = 14, height = 5, device = cairo_pdf)

## asymptotic resilience #### ----
data_resilience<-list_data$data_resilience

p2<-ggplot(data=data_resilience[data_resilience$model=="pert_31",])+
      geom_line(aes(gamma,log10(resilience)),size=1.5)+
      theme+
      xlab(label_gamma)+
      ylab(label_resilience)

## eigenvector #### ----
data_E<-list_data$data_E # lead eigen vector
data_E<-data_E[,which(names(data_E)%in%c("model","gamma","omega",list_data$E_names))]
data_E<-table_for_plot(data_E,3,"contribution")

# p1<-ggplot(data=databis[databis$gamma==1,])+
#   geom_line(aes(d,contribution,colour=species,linetype=community),size=1.5)+
#   corr_colour_TL_2+
#   patch_line+
#   theme
# legend_1<-get_legend(p1)

p3<-ggplot(data=data_E[data_E$model=="pert_31",])+
      geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
      corr_colour_TL_4+
      patch_line+
      theme + #theme(legend.position = "none")+
      #x_axis_log10_short+
      #ylim(0,1)+
      xlab(label_gamma)+
      ylab(label_contribution)

## final graph #### ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  #draw_plot(legend_1, 2, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,1,1), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_resilience_disp_14.pdf",sep=""), graph, width = 17, height = 6, device = cairo_pdf)

## correlation #### ----
data_C<-list_data$data_C
# levels(data_C$model)<-c("perturbation of\nspecies 4 in patch #1",
#                         "perturbation of\nspecies 4 in patch #2",
#                         "perturbation of\nspecies 3 in patch #1",
#                         "perturbation of\nspecies 3 in patch #2")
data_C$pert_species="perturbation of species 4"
data_C$pert_species[which(data_C$model%in%c("pert_31","pert_32"))]="perturbation of species 3"
data_C$pert_community="1"
data_C$pert_community[which(data_C$model%in%c("pert_42","pert_32"))]="2"
data_C<-data_C[,c(which(names(data_C)%in%c("gamma","omega","pert_species","pert_community","C_31_32","C_41_42")))]
data_C<-data_C[,c(1,2,5,6,3,4)]
data_C<-table_for_plot(data_C,4,"correlation")

p1<-ggplot(data=data_C)+
  geom_line(aes(gamma,correlation,colour=species,linetype=pert_community),size=1.5)+
  corr_colour_top+
  patch_line_pert+
  theme
legend_1<-get_legend(p1)

p1<-ggplot(data=data_C)+
      geom_line(aes(gamma,correlation,colour=species,linetype=pert_community),size=1.5)+
      facet_wrap(~pert_species)+
      corr_colour_top+
      patch_line_pert+
      theme+theme(legend.position = "none")+
      #ylim(0.2,1)+
      xlab(label_gamma)+
      ylab(label_correlation)

## biomass CV in perturbed patch #### ----
data_V<-list_data$data_V
data_B<-list_data$data_B
# selects the variance and the biomass of species 3 in patch #1 when perturbation occurs in patch #1
data_V$V_3<-data_V$V_31_31
data_V$V_3[data_V$model=="pert_32" | data_V$model=="pert_42"]<-data_V$V_32_32[data_V$model=="pert_32" | data_V$model=="pert_42"]
data_V$B_3<-data_B$B_31
data_V$B_3[data_V$model=="pert_32" | data_V$model=="pert_42"]<-data_B$B_32[data_B$model=="pert_32" | data_B$model=="pert_42"]
# selects the variance and the biomass of species 4 in patch #1 when perturbation occurs in patch #1
data_V$V_4<-data_V$V_41_41
data_V$V_4[data_V$model=="pert_32" | data_V$model=="pert_42"]<-data_V$V_42_42[data_V$model=="pert_32" | data_V$model=="pert_42"]
data_V$B_4<-data_B$B_41
data_V$B_4[data_V$model=="pert_32" | data_V$model=="pert_42"]<-data_B$B_42[data_B$model=="pert_32" | data_B$model=="pert_42"]
# CV
data_V$CV_3<-sqrt(data_V$V_3)/data_V$B_3
data_V$CV_4<-sqrt(data_V$V_4)/data_V$B_4

data_V$pert_species="perturbation of species 4"
data_V$pert_species[which(data_V$model%in%c("pert_31","pert_32"))]="perturbation of species 3"
data_V$pert_community="1"
data_V$pert_community[which(data_V$model%in%c("pert_42","pert_32"))]="2"
data_V<-data_V[,c(which(names(data_V)%in%c("gamma","omega","pert_species","pert_community","CV_3","CV_4")))]
data_V<-data_V[,c(1,2,5,6,3,4)]
data_V<-table_for_plot(data_V,4,"CV")
data_V$species<-data_V$community

p2<-ggplot(data=data_V)+
      geom_line(aes(gamma,CV,colour=species,linetype=pert_community),size=1.5)+
      facet_wrap(~pert_species)+
      corr_colour_top+
      patch_line_pert+
      theme+theme(legend.position = "none")+
      #ylim(0.2,1)+
      xlab(label_gamma)+
      ylab("CV in the perturbed patch")

## final graph #### ----
graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend_1, 2.02, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"figure_correlation_CV_disp_14.pdf",sep=""), graph, width = 14, height = 6, device = cairo_pdf)

# time series of a pulse perturbation #### ----
nSpecies=4
nCommunity=2
gamma=5 # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model=c("pred_31","pert_32"))
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-list_data$data_B
B0<-B0[B0$gamma==gamma & (B0$model=="pert_31" | B0$model=="pert_32"),]
B0_ref<-B0[,-((nSpecies*nCommunity-1):0)+dim(B0)[2]] # reference biomass for the rescaling
# sets perturbations
B0$B_31[B0$model=="pert_31"]=0.9*B0$B_31[B0$model=="pert_31"]
B0$B_32[B0$model=="pert_32"]=0.9*B0$B_32[B0$model=="pert_32"]
B0<-B0[,-((nSpecies*nCommunity-1):0)+dim(B0)[2]] # selects the biomasses
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=50
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data, as.numeric(B0[i,]), i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  TS<-results[[i]]
  for(j in 1:dim(TS)[1]){
    TS[j,c(2:dim(TS)[2])]<-TS[j,c(2:dim(TS)[2])]/as.numeric(B0_ref[i,])
  }
  TS<-time_series_for_plot(params_data, TS, i)
  data_TS<-rbind(data_TS,TS)
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_wrap(~model)+
  corr_colour_TL_4+
  patch_line+
  theme+
  #y_axis_log10_short+
  xlab("Time")+
  ylab("Scaled biomass")

########################## ----
# APPENDIX ############## ----
########################## ----
### DSIPERSAL OF SPECIES 4 #### ----
# preliminary test for the ODE integration #### ----
nSpecies=4
nCommunity=2
gamma=5 # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(0,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=300
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data,B0[i,],i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  data_TS<-rbind(data_TS,time_series_for_plot(params_data, results[[i]], i))
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_grid(ma~ea, labeller=label_parsed)+
  corr_colour_TL_4+
  patch_line+
  theme+
  y_axis_log10_short+
  xlab("Time")+
  ylab("Biomass")

# biomass calculation - gamma - omega (TO DO ONCE) #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(0,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e3, # high sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=300
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve")) %dopar% equilibrium_ode(params_data,B0[i,],i)
stopCluster(cl)

data_B<-NULL
for(i in 1:dim(params_data)[1]){
  data_B<-rbind(data_B,results[[i]])
}
data_B<-data_B[,-which(names(data_B)%in%c("pert","asym","disp"))]

write.table(data_B,paste(path_data,"biomass_disp4_full_main.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

# analysis - gamma - omega - load biomasses #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(4,1))), # perturbations of top predators in patch 1
                         disp=list(c(0,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1,
                         model="pert_41")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]

data_B<-read.table(paste(path_data,"biomass_disp4_full_main.txt",sep=""),sep=";",header=TRUE)
biomass<-na.omit(data_B) # remove simulations without coexistance
params_data<-params_data[biomass$simu_ID,] # remove simulations without coexistance
biomass<-biomass[,-((nSpecies*nCommunity-1):0)+dim(biomass)[2]] # selects the biomasses

# adds perturbations of other species
params_data_42<-params_data
params_data_42$pert=list(list(c(4,2))) # perturbations of top predator in patch 2
params_data_42$model="pert_42"
params_data_31<-params_data
params_data_31$pert=list(list(c(3,1))) # perturbations of predators in patch 1
params_data_31$model="pert_31"
params_data_32<-params_data
params_data_32$pert=list(list(c(3,2))) # perturbations of predators in patch 2
params_data_32$model="pert_32"
params_data<-rbind(params_data,params_data_42,params_data_31,params_data_32)
rm(params_data_42,params_data_31,params_data_32)
biomass<-rbind(biomass,biomass,biomass,biomass)

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% linear_analysis(params_data,biomass,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,results,nSpecies,nCommunity)

## biomass #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B[data_B$model=="pert_31",])+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position="none")+
  y_axis_log10+
  xlab(label_gamma)+
  ylab("Biomass")+
  ggtitle("Dispersal of species 4")

## biomass scaled by biomass without dispersal #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
# reference biomass without dispersal
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
for (i in 1:dim(data_B)[1]){
  data_B[i,4:dim(data_B)[2]]=data_B[i,4:dim(data_B)[2]]/B0[i,]
}
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

ggplot(data=data_B[data_B$model=="pert_31",])+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position="none")+
  #y_axis_log10+
  xlab(label_gamma)+
  ylab("Scaled biomass")+
  ggtitle("Dispersal of species 1 and 4")

## asymptotic resilience #### ----
data_resilience<-list_data$data_resilience

p2<-ggplot(data=data_resilience[data_resilience$model=="pert_31",])+
  geom_line(aes(gamma,log10(resilience)),size=1.5)+
  theme+
  xlab(label_gamma)+
  ylab(label_resilience)+
  ggtitle("Dispersal of species 4")

## eigen vector #### ----
data_E<-list_data$data_E # lead eigen vector
data_E<-data_E[,which(names(data_E)%in%c("model","gamma","omega",list_data$E_names))]
data_E<-table_for_plot(data_E,3,"contribution")

# p1<-ggplot(data=databis[databis$gamma==1,])+
#   geom_line(aes(d,contribution,colour=species,linetype=community),size=1.5)+
#   corr_colour_TL_2+
#   patch_line+
#   theme
# legend_1<-get_legend(p1)

p3<-ggplot(data=data_E[data_E$model=="pert_31",])+
  geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position = "none")+
  #x_axis_log10_short+
  #ylim(0,1)+
  xlab(label_gamma)+
  ylab(label_contribution)+
  ggtitle("Dispersal of species 4")

## final graph #### ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  #draw_plot(legend_1, 2, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,1,1), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_resilience_disp_4.pdf",sep=""), graph, width = 17, height = 6, device = cairo_pdf)

## correlation #### ----
data_C<-list_data$data_C
levels(data_C$model)<-c("species 4 in patch #1",
                        "species 4 in patch #2",
                        "species 3 in patch #1",
                        "species 3 in patch #2")

p1<-ggplot(data=data_C)+
  geom_line(aes(gamma,C_31_32),size=1.5)+
  facet_wrap(~model)+
  theme + #theme(legend.position="none")+
  ylim(0.2,1)+
  xlab(label_gamma)+
  ylab(label_correlation)+
  ggtitle("Dispersal of species 4")

ggsave(paste(path_figure,"figure_correlation_disp_4.pdf",sep=""), p1, width = 10, height = 10, device = cairo_pdf)

# time series of a pulse perturbation #### ----
nSpecies=4
nCommunity=2
gamma=5 # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(0,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model=c("pred_31","pert_32"))
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-list_data$data_B
B0<-B0[B0$gamma==5 & (B0$model=="pert_31" | B0$model=="pert_32"),]
B0_ref<-B0[,-((nSpecies*nCommunity-1):0)+dim(B0)[2]] # reference biomass for the rescaling
# sets perturbations
B0$B_31[B0$model=="pert_31"]=1.2*B0$B_31[B0$model=="pert_31"]
B0$B_32[B0$model=="pert_32"]=1.2*B0$B_32[B0$model=="pert_32"]
B0<-B0[,-((nSpecies*nCommunity-1):0)+dim(B0)[2]] # selects the biomasses
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=50
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data, as.numeric(B0[i,]), i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  TS<-results[[i]]
  for(j in 1:dim(TS)[1]){
    TS[j,c(2:dim(TS)[2])]<-TS[j,c(2:dim(TS)[2])]/as.numeric(B0_ref[i,])
  }
  TS<-time_series_for_plot(params_data, TS, i)
  data_TS<-rbind(data_TS,TS)
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_wrap(~model)+
  corr_colour_TL_4+
  patch_line+
  theme+
  #y_axis_log10_short+
  xlab("Time")+
  ylab("Scaled biomass")

### DSIPERSAL OF SPECIES 1 #### ----
# preliminary test for the ODE integration #### ----
nSpecies=4
nCommunity=2
gamma=5 # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,0)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=2000
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data,B0[i,],i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  data_TS<-rbind(data_TS,time_series_for_plot(params_data, results[[i]], i))
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_grid(ma~ea, labeller=label_parsed)+
  corr_colour_TL_4+
  patch_line+
  theme+
  y_axis_log10_short+
  xlab("Time")+
  ylab("Biomass")

# biomass calculation - gamma - omega (TO DO ONCE) #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,0)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e3, # high sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=2000
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve")) %dopar% equilibrium_ode(params_data,B0[i,],i)
stopCluster(cl)

data_B<-NULL
for(i in 1:dim(params_data)[1]){
  data_B<-rbind(data_B,results[[i]])
}
data_B<-data_B[,-which(names(data_B)%in%c("pert","asym","disp"))]

write.table(data_B,paste(path_data,"biomass_disp1_full_main.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

# analysis - gamma - omega - load biomasses #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(4,1))), # perturbations of top predators in patch 1
                         disp=list(c(1,0,0,0)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1,
                         model="pert_41")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]

data_B<-read.table(paste(path_data,"biomass_disp1_full_main.txt",sep=""),sep=";",header=TRUE)
biomass<-na.omit(data_B) # remove simulations without coexistance
params_data<-params_data[biomass$simu_ID,] # remove simulations without coexistance
biomass<-biomass[,-((nSpecies*nCommunity-1):0)+dim(biomass)[2]] # selects the biomasses

# adds perturbations of other species
params_data_42<-params_data
params_data_42$pert=list(list(c(4,2))) # perturbations of top predator in patch 2
params_data_42$model="pert_42"
params_data_31<-params_data
params_data_31$pert=list(list(c(3,1))) # perturbations of predators in patch 1
params_data_31$model="pert_31"
params_data_32<-params_data
params_data_32$pert=list(list(c(3,2))) # perturbations of predators in patch 2
params_data_32$model="pert_32"
params_data<-rbind(params_data,params_data_42,params_data_31,params_data_32)
rm(params_data_42,params_data_31,params_data_32)
biomass<-rbind(biomass,biomass,biomass,biomass)

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% linear_analysis(params_data,biomass,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,results,nSpecies,nCommunity)

## biomass #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B[data_B$model=="pert_31",])+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position="none")+
  y_axis_log10+
  xlab(label_gamma)+
  ylab("Biomass")+
  ggtitle("Dispersal of species 1")

## asymptotic resilience #### ----
data_resilience<-list_data$data_resilience

p2<-ggplot(data=data_resilience[data_resilience$model=="pert_31",])+
  geom_line(aes(gamma,log10(resilience)),size=1.5)+
  theme+
  xlab(label_gamma)+
  ylab(label_resilience)+
  ggtitle("Dispersal of species 1")

## eigen vector #### ----
data_E<-list_data$data_E # lead eigen vector
data_E<-data_E[,which(names(data_E)%in%c("model","gamma","omega",list_data$E_names))]
data_E<-table_for_plot(data_E,3,"contribution")

# p1<-ggplot(data=databis[databis$gamma==1,])+
#   geom_line(aes(d,contribution,colour=species,linetype=community),size=1.5)+
#   corr_colour_TL_2+
#   patch_line+
#   theme
# legend_1<-get_legend(p1)

p3<-ggplot(data=data_E[data_E$model=="pert_31",])+
  geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position = "none")+
  #x_axis_log10_short+
  #ylim(0,1)+
  xlab(label_gamma)+
  ylab(label_contribution)+
  ggtitle("Dispersal of species 1")

## final graph #### ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  #draw_plot(legend_1, 2, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,1,1), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_resilience_disp_1.pdf",sep=""), graph, width = 17, height = 6, device = cairo_pdf)

## correlation #### ----
data_C<-list_data$data_C
levels(data_C$model)<-c("species 4 in patch #1",
                        "species 4 in patch #2",
                        "species 3 in patch #1",
                        "species 3 in patch #2")

p1<-ggplot(data=data_C)+
  geom_line(aes(gamma,C_31_32),size=1.5)+
  facet_wrap(~model)+
  theme + #theme(legend.position="none")+
  #ylim(0.2,1)+
  xlab(label_gamma)+
  ylab(label_correlation)+
  ggtitle("Dispersal of species 1")

ggsave(paste(path_figure,"figure_correlation_disp_1.pdf",sep=""), p1, width = 10, height = 10, device = cairo_pdf)

### FULL MODEL NO DENSITY-DEPENDENT DISPERSAL #### ----
# preliminary test for the ODE integration #### ----
nSpecies=4
nCommunity=2
gamma=5 # asymmetry coefficient
omega=5 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e3, # low sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=100
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data,B0[i,],i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  data_TS<-rbind(data_TS,time_series_for_plot(params_data, results[[i]], i))
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_grid(ma~ea, labeller=label_parsed)+
  corr_colour_TL_4+
  patch_line+
  theme+
  y_axis_log10_short+
  xlab("Time")+
  ylab("Biomass")

# biomass calculation - gamma - omega (TO DO ONCE) #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e3, # high sensitivity to predator abundance
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=500
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve")) %dopar% equilibrium_ode(params_data,B0[i,],i)
stopCluster(cl)

data_B<-NULL
for(i in 1:dim(params_data)[1]){
  data_B<-rbind(data_B,results[[i]])
}
data_B<-data_B[,-which(names(data_B)%in%c("pert","asym","disp"))]

write.table(data_B,paste(path_data,"biomass_disp14_full_main.txt",sep=""),sep=";",qmethod = "double",row.names = FALSE)

# analysis - gamma - omega - load biomasses #### ----
nSpecies=4
nCommunity=2
gamma=seq(0.5,5,0.5) # asymmetry coefficient
omega=1
d=1e3 # scaled dispersal rate

params_data<-expand.grid(pert=list(list(c(4,1))), # perturbations of top predators in patch 1
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1,
                         model="pert_41")
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
params_data$simu_ID<-seq(1,dim(params_data)[1])
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]

data_B<-read.table(paste(path_data,"biomass_disp14_full_main.txt",sep=""),sep=";",header=TRUE)
biomass<-na.omit(data_B) # remove simulations without coexistance
params_data<-params_data[biomass$simu_ID,] # remove simulations without coexistance
biomass<-biomass[,-((nSpecies*nCommunity-1):0)+dim(biomass)[2]] # selects the biomasses

# adds perturbations of other species
params_data_42<-params_data
params_data_42$pert=list(list(c(4,2))) # perturbations of top predator in patch 2
params_data_42$model="pert_42"
params_data_31<-params_data
params_data_31$pert=list(list(c(3,1))) # perturbations of predators in patch 1
params_data_31$model="pert_31"
params_data_32<-params_data
params_data_32$pert=list(list(c(3,2))) # perturbations of predators in patch 2
params_data_32$model="pert_32"
params_data<-rbind(params_data,params_data_42,params_data_31,params_data_32)
rm(params_data_42,params_data_31,params_data_32)
biomass<-rbind(biomass,biomass,biomass,biomass)

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("numDeriv")) %dopar% linear_analysis(params_data,biomass,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,results,nSpecies,nCommunity)

## biomass #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B[data_B$model=="pert_31",])+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position="none")+
  y_axis_log10+
  xlab(label_gamma)+
  ylab("Biomass")+
  ggtitle("Dispersal of species 1 and 4")

## biomass scaled by biomass without dispersal #### ----
data_B<-list_data$data_B
data_B<-data_B[,c(which(names(data_B)%in%c("gamma","omega","model",list_data$B_names)))]
# reference biomass without dispersal
B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  for(j in 1:nCommunity){
    B0[i,(1:nSpecies)+nSpecies*(j-1)]<-equilibrium_isolated(params,nSpecies,j)
  }
}
for (i in 1:dim(data_B)[1]){
  data_B[i,4:dim(data_B)[2]]=data_B[i,4:dim(data_B)[2]]/B0[i,]
}
data_B<-table_for_plot(data_B,3,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

ggplot(data=data_B[data_B$model=="pert_31",])+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position="none")+
  #y_axis_log10+
  xlab(label_gamma)+
  ylab("Scaled biomass")+
  ggtitle("Dispersal of species 1 and 4")

## asymptotic resilience #### ----
data_resilience<-list_data$data_resilience

p2<-ggplot(data=data_resilience[data_resilience$model=="pert_31",])+
  geom_line(aes(gamma,log10(resilience)),size=1.5)+
  theme+
  xlab(label_gamma)+
  ylab(label_resilience)+
  ggtitle("Dispersal of species 1 and 4")

## eigen vector #### ----
data_E<-list_data$data_E # lead eigen vector
data_E<-data_E[,which(names(data_E)%in%c("model","gamma","omega",list_data$E_names))]
data_E<-table_for_plot(data_E,3,"contribution")

# p1<-ggplot(data=databis[databis$gamma==1,])+
#   geom_line(aes(d,contribution,colour=species,linetype=community),size=1.5)+
#   corr_colour_TL_2+
#   patch_line+
#   theme
# legend_1<-get_legend(p1)

p3<-ggplot(data=data_E[data_E$model=="pert_31",])+
  geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_4+
  patch_line+
  theme + #theme(legend.position = "none")+
  #x_axis_log10_short+
  #ylim(0,1)+
  xlab(label_gamma)+
  ylab(label_contribution)+
  ggtitle("Dispersal of species 1 and 4")

## final graph #### ----
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  #draw_plot(legend_1, 2, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,1,1), c(1,1,1), size = 30)
ggsave(paste(path_figure,"figure_resilience_disp_14.pdf",sep=""), graph, width = 17, height = 6, device = cairo_pdf)

## correlation #### ----
data_C<-list_data$data_C
levels(data_C$model)<-c("species 4 in patch #1",
                        "species 4 in patch #2",
                        "species 3 in patch #1",
                        "species 3 in patch #2")

p1<-ggplot(data=data_C)+
  geom_line(aes(gamma,C_31_32),size=1.5)+
  facet_wrap(~model)+
  theme +
  ylim(0.2,1)+
  xlab(label_gamma)+
  ylab(label_correlation)+
  ggtitle("Dispersal of species 1 and 4")

ggsave(paste(path_figure,"figure_correlation_disp_14.pdf",sep=""), p1, width = 10, height = 10, device = cairo_pdf)

# time series of a pulse perturbation #### ----
nSpecies=4
nCommunity=2
gamma=5 # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e3 #d_interval

params_data<-expand.grid(pert=list(list(c(2,1))), # irrelevant here
                         disp=list(c(1,0,0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         S0_prey=1e-3, # high sensitivity to prey abundance
                         S0_pred=1e-3, # low sensitivity to predator abundance
                         model=c("pred_31","pert_32"))
params_data<-merge(params_data_original,params_data)
# selects ea=10 and ma=10
params_data<-params_data[params_data$ea==paste(ea10) & params_data$ma==paste(ma10),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
# initial biomasses computed analytically before ODE resolution
B0<-list_data$data_B
B0<-B0[B0$gamma==5 & (B0$model=="pert_31" | B0$model=="pert_32"),]
B0_ref<-B0[,-((nSpecies*nCommunity-1):0)+dim(B0)[2]] # reference biomass for the rescaling
# sets perturbations
B0$B_31[B0$model=="pert_31"]=1.2*B0$B_31[B0$model=="pert_31"]
B0$B_32[B0$model=="pert_32"]=1.2*B0$B_32[B0$model=="pert_32"]
B0<-B0[,-((nSpecies*nCommunity-1):0)+dim(B0)[2]] # selects the biomasses
# density-depend dispersal sensitivity S (set according to the biomass in the reference patch #2)
params_data$S_pred=params_data$S0_pred*B0[,6]
params_data$S_prey=params_data$S0_prey*B0[,7]
# integration time
params_data$t_max=100
params_data$t_step=0.001

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("nleqslv","deSolve","reshape2")) %dopar% time_series(params_data, as.numeric(B0[i,]), i)
stopCluster(cl)

data_TS<-NULL
for(i in 1:dim(params_data)[1]){
  TS<-results[[i]]
  for(j in 1:dim(TS)[1]){
    TS[j,c(2:dim(TS)[2])]<-TS[j,c(2:dim(TS)[2])]/as.numeric(B0_ref[i,])
  }
  TS<-time_series_for_plot(params_data, TS, i)
  data_TS<-rbind(data_TS,TS)
}

ggplot(data=data_TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_wrap(~model)+
  corr_colour_TL_4+
  patch_line+
  theme+
  #y_axis_log10_short+
  xlab("Time")+
  ylab("Scaled biomass")
