#####################################################################################
##################### Bootstrapping para age y transit time #########################
#####################################################################################

rm(list=ls()) 


if(!require(readxl)){install.packages("readxl")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(rlist)){install.packages("rlist")}
if(!require(FME)){install.packages("FME")}
if(!require(SoilR)){install.packages("SoilR")}
if(!require(forecast)){install.packages("forecast")}
if(!require(bayestestR)){install.packages("bayestestR")}


#########################################################################################

#--------------------- Radiocarbon curves


hua <- data.frame(
  year = Hua2021$`SHZone1-2`[ , 1],
  delta14C = Hua2021$`SHZone1-2`[ , 2]
)

yrs=seq(1965,2019,by=1/4) # A series of years by quarters

nz2=spline(hua, xout=yrs) #Spline interpolation of the NH_Zone 2 dataset at a quaterly basis
nhz2=ts((nz2$y-1),start=1965, freq=4) #Transformation into a time-series object

m=ets(nhz2) #Fits an exponential smoothing state space model to the time series
f2=forecast(m,h=11*4) #Uses the fitted model to forecast 11 years into the future

bc=data.frame(year=c(hua[-dim(hua)[1],1],
                     seq(tsp(f2$mean)[1],tsp(f2$mean)[2], by=1/tsp(f2$mean)[3])),
              
              Delta14C=c(hua[-dim(hua)[1],2],as.numeric(f2$mean)))


preandpost <- bind.C14curves(prebomb=IntCal13, postbomb=Hua2021$`SHZone1-2`,time.scale="AD")

preandpost <- dplyr::filter(preandpost, Year.AD > -3000 & Year.AD < 1945)

preandpost <- preandpost[,1:2]
colnames(preandpost) <- c("year", "Delta14C")

AtmFc1 <- rbind(preandpost[,1:2], bc)

AtmFc <- BoundFc(map=AtmFc1, format="Delta14C", lag=3) 


rm(bc, f2, hua, m, nhz2, nz2, preandpost)


#########################################################################################

# R system
propagation_R <- function(data, num_iter, input_vector){
  
  pb <- txtProgressBar(min = 0, max = num_iter, style = 3)

  # Definir el número de iteraciones
  n_iter <- num_iter
  
  # Crear un vector para almacenar los valores p de cada iteración
  iter <- numeric(n_iter)
  Iter_number <- numeric(n_iter)
  System_age <- numeric(n_iter)
  POM_age <- numeric(n_iter)
  MAOM_age <- numeric(n_iter)
  Transit_time <- numeric(n_iter)
  
  
  # Iniciar el loop de iteraciones
  for (i in 1:n_iter){
    # Dividir cada unidad experimental en subconjuntos de tamaño 'subsampling' usando purrr
    
    subset_par <- data.frame(k1= sample(data$p1, 1),
                             k2= sample(data$p2, 1),
                             alfa= sample(data$p3, 1))
    
    
    
    #-------------- A matrix and inputs
    
    ks <- subset_par[1:2]
    
    A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
    #de esa matriz.
    
    ## Add transfer coefficients to A matrix:
    alpha_2_1 <- subset_par[3]
    
    A[2,1] <- as.numeric(alpha_2_1*subset_par[1])
    
    u <-matrix((input_vector), ncol=1)
    
    #---------------- Age and transit time ----------------------
    

    Sist_age <- systemAge(A=A, u=u)
    Trans_time <- transitTime(A=A, u=u)
    
    Iter_number[i] <- i
    System_age[i] <- as.numeric(Sist_age$meanSystemAge)
    POM_age[i] <- as.numeric(Sist_age$meanPoolAge[1])
    MAOM_age[i] <- as.numeric(Sist_age$meanPoolAge[2])
    Transit_time[i] <- as.numeric(Trans_time$meanTransitTime)
    
    
    salida <- as.data.frame(cbind(
      Iter_number, System_age, POM_age, MAOM_age, Transit_time
    ))
    
    setTxtProgressBar(pb, i)
    
  }
  
  return(salida)
}

#########################################################################################

# CC system
propagation_CC <- function(data, num_iter, input_vector){
  
  pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
  
  # Definir el número de iteraciones
  n_iter <- num_iter
  
  # Crear un vector para almacenar los valores p de cada iteración
  iter <- numeric(n_iter)
  Iter_number <- numeric(n_iter)
  System_age <- numeric(n_iter)
  POM_age <- numeric(n_iter)
  MAOM_age <- numeric(n_iter)
  Transit_time <- numeric(n_iter)
  
  
  # Iniciar el loop de iteraciones
  for (i in 1:n_iter){
    # Dividir cada unidad experimental en subconjuntos de tamaño 'subsampling' usando purrr
    
    subset_par <- data.frame(k1= sample(data$p1, 1),
                             k2= sample(data$p4, 1),
                             alfa= sample(data$p5, 1))
    
    
    
    #-------------- A matrix and inputs
    
    ks <- subset_par[1:2]
    
    A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
    #de esa matriz.
    
    ## Add transfer coefficients to A matrix:
    alpha_2_1 <- subset_par[3]
    
    A[2,1] <- as.numeric(alpha_2_1*subset_par[1])
    
    u <-matrix((input_vector), ncol=1)
    
    #---------------- Age and transit time ----------------------
    
    
    Sist_age <- systemAge(A=A, u=u)
    Trans_time <- transitTime(A=A, u=u)
    
    Iter_number[i] <- i
    System_age[i] <- as.numeric(Sist_age$meanSystemAge)
    POM_age[i] <- as.numeric(Sist_age$meanPoolAge[1])
    MAOM_age[i] <- as.numeric(Sist_age$meanPoolAge[2])
    Transit_time[i] <- as.numeric(Trans_time$meanTransitTime)
    
    
    salida <- as.data.frame(cbind(
      Iter_number, System_age, POM_age, MAOM_age, Transit_time
    ))
    
    setTxtProgressBar(pb, i)
    
  }
  
  return(salida)
}


######################## R_ Sistem #######################################################

#------------------------ Cargar MCMC result

MCMC_R_results <- list.load('lists_MCMC/MCMC_R.rdata')

pars_df_R <- as.data.frame((MCMC_R_results$pars))


age_transit_dist_R <- propagation_R(pars_df_R, 10000, c(11.45 * 0.45, 0))




export <- list("age_transit_dist_R" = age_transit_dist_R) #Crea el objeto "salida" (elegir cualquier nombre) en el que se listan las

write.xlsx(export, "outputs_age_tt_error_prop/age_transit_dist_R.xlsx", colWidths = c("auto", "auto"),overwrite=TRUE) #Exporta en arhivo .xlsx la lista creada anteriormente.



######################## CC_ Sistem #######################################################

#------------------------ Cargar MCMC result

MCMC_CC_results <- list.load('lists_MCMC/MCMC_CC.rdata')

pars_df_CC <- as.data.frame((MCMC_CC_results$pars))


age_transit_dist_CC <- propagation_CC(pars_df_CC, 10000, c(6.37 * 0.45, 0))




export <- list("age_transit_dist_CC" = age_transit_dist_CC) #Crea el objeto "salida" (elegir cualquier nombre) en el que se listan las

write.xlsx(export, "outputs_age_tt_error_prop/age_transit_dist_CC.xlsx", colWidths = c("auto", "auto"),overwrite=TRUE) #Exporta en arhivo .xlsx la lista creada anteriormente.



##########################################################################################





R_age_and_TT_R <- read_excel("outputs_age_tt_error_prop/age_transit_dist_R.xlsx", sheet = "age_transit_dist_R")

resumen <- function(vector){
  data.frame(Mean_Age = mean(vector),
             Mean_Age_SD = sd(vector),
             Q1 = as.numeric(quantile(vector, probs = c(0.25))),
             Median = median(vector),
             Q3 = as.numeric(quantile(vector, probs = c(0.75)))
  )
}


R_System <- rbind(
  cbind(Variable="Mean Age" ,resumen(R_age_and_TT_R$System_age)),
  cbind(Variable="POM Age" ,resumen(R_age_and_TT_R$POM_age)),
  cbind(Variable="MAOM Age" ,resumen(R_age_and_TT_R$MAOM_age)),
  cbind(Variable="Transit time" ,resumen(R_age_and_TT_R$Transit_time))
)


R_age_and_TT_CC <- read_excel("outputs_age_tt_error_prop/age_transit_dist_CC.xlsx", sheet = "age_transit_dist_CC")

CC_System <- rbind(
  cbind(Variable="Mean Age" ,resumen(R_age_and_TT_CC$System_age)),
  cbind(Variable="POM Age" ,resumen(R_age_and_TT_CC$POM_age)),
  cbind(Variable="MAOM Age" ,resumen(R_age_and_TT_CC$MAOM_age)),
  cbind(Variable="Transit time" ,resumen(R_age_and_TT_CC$Transit_time))
)



str(R_age_and_TT_R)

ci_System_age_R <- ci(R_age_and_TT_R$System_age, method = "ETI") #Credible interval alfa R System
ci_System_age_CC<- ci(R_age_and_TT_CC$System_age, method = "ETI") #Credible interval alfa R System



