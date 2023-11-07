#########################################################################################
############################## Model CC SensVar plot ####################################
#########################################################################################

rm(list=ls()) 

if(!require(readxl)){install.packages("readxl")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(data.table)){install.packages("data.table")}
if(!require(SoilR)){install.packages("SoilR")}
if(!require(forecast)){install.packages("forecast")}
if(!require(FME)){install.packages("FME")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(rlist)){install.packages("rlist")}
if(!require(bayestestR)){install.packages("bayestestR")}
if(!require(cowplot)){install.packages("cowplot")}
if(!require(rlang)){install.packages("rlang")}




#########################################################################################

MCMC_R_results <- list.load('lists_MCMC/MCMC_R.rdata')
MCMC_CC_results <- list.load('lists_MCMC/MCMC_CC.rdata')



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

AtmFc_R <- BoundFc(map=AtmFc1, format="Delta14C", lag=3) 
AtmFc_CC <- BoundFc(map=AtmFc1, format="Delta14C", lag=0.5) 


rm(bc, f2, hua, m, nhz2, nz2, preandpost)



#########################################################

#----------- Model 0 - 1963

times <- seq(0, 1963, by=1) # Vector de tiempos 

F0 <- ConstFc(values=c(-50, -150), format="Delta14C") #

pools <- c(x1="Labile",
           x2="Stable"
)


C0 <- c(x1= 53 * 0.1, 
        x2= 53 * 0.9)



inputVector <- c(11.45 *0.45,   #x1
                 0    #x2
)



pars <- MCMC_R_results$bestpar


ks <- pars[1:2]

A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
#de esa matriz.

## Add transfer coefficients to A matrix:
alpha_2_1 <- pars[3]

A[2,1] <- alpha_2_1*pars[1]



Model14 <- Model_14(t = times, A = A, ivList = C0, inputFluxes = inputVector, 
                    initialValF = F0, inputFc = AtmFc_R, c14DecayRate = -0.0001209681)



#----------------- Obtención de salidas del modelo


a <- data.frame(
  time=times,
  getF14(Model14))

b <- data.frame(
  time=times,
  getC(Model14))

c <- getReleaseFlux14(Model14)

d <- data.frame(time=times,
                R_F14 =  getF14R(Model14)
)

e <- data.frame(
  time = times,
  C_14_mean = getF14C(Model14)
)


a <- cbind(a, d[2], e[2])

colnames(b) <- c("time", pools); colnames(a) <-  c("time", pools, "Efflux", "Bulk soil")



#########################################################################################
#------ Outputs to take as inputs

C_14 <- a %>%
  gather('Labile' , 'Stable', 'Efflux', 'Bulk soil',
         key= "Pool", value = "C_14")

C_14$Pool <- factor(C_14$Pool)


C_stock <- b %>% 
  gather('Labile' , 'Stable',
         key= "Pool", value = "C_stock")

C_stock$Pool <- factor(C_stock$Pool)



#########################################################################################


#########################################################################################

#---------------- Modelo a partir de 1964



# Sacar los datos de punto de partida de la simulación anterior para representar a partir de 1964 el sistema de CC.


C_14_1963 <- C_14 %>% 
  dplyr::filter(time==1963)

C_stock_1963 <- C_stock %>% 
  dplyr::filter(time==1963)

F0 <- ConstFc(values=c(C_14_1963[1,3], C_14_1963[2,3]), format="Delta14C") #

C0 <- c(x1= C_stock_1963[1,3], 
        x2= C_stock_1963[2,3])


inputVector <- c(6.37 *0.45,   #x1
                 0    #x2
)


pars <- MCMC_CC_results$pars



fun_error_prop <- function(pars){
  
  #------------ Fase 1964 - 1990
  
  times <- seq(1964, 1990, by=1) # Vector de tiempos 
  
  # Define compartmental matrix
  ks <- pars[1:2]
  
  A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
  #de esa matriz.
  
  ## Add transfer coefficients to A matrix:
  alpha_2_1 <- pars[3]
  
  A[2,1] <- alpha_2_1*pars[1]
  
  
  Model14 <- Model_14(t = times, A = A, ivList = C0, inputFluxes = inputVector, 
                      initialValF = F0, inputFc = AtmFc_CC, c14DecayRate = -0.0001209681)
  
  C_pools <- getC(Model14) ; colnames(C_pools) <-  c("C_POM", "C_MAOM")
  
  C_14_pools <- getF14(Model14) ;  colnames(C_14_pools) <-  c("C_14_POM", "C_14_MAOM")
  
  C_total <- rowSums(getC(Model14)) # stock de SOC; suma todos los pools en cada año
  
  C_14_bulk <- getF14C(Model14) # delta14C del suelo en bulk. Promedio ponderado por masa del delta14C de todos los pools. 
  
  C_14_incub <- getF14R(Model14) # Firma delta14C en CO2. Promedio ponderado de la firma delta14C todos los pools#(getF14) por lo 
  #que se oxida en cada momento de cada pool (getReleaseFlux)
  
  model_output_1 <- data.frame(time=times, 
                               C_pools,
                               C_ha_mean= C_total,
                               C_14_pools,
                               C_14_bulk = C_14_bulk,
                               C_14_CO2 = C_14_incub)
  
  

  
  #------------ Fase 1991 - 2022
  
  F0 <- ConstFc(values=c(model_output_1[c(model_output_1$time==1990), c("C_14_POM")],
                         model_output_1[c(model_output_1$time==1990), c("C_14_MAOM")]), format="Delta14C") #
  
  
  C0 <- c(x1= model_output_1[c(model_output_1$time==1990), c("C_POM")], 
          x2= model_output_1[c(model_output_1$time==1990), c("C_MAOM")])
  
  
  times <- seq(1990, 2022, by=1) # Vector de tiempos 
  
  ks <- pars[c(1, 4)]
  
  A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
  #de esa matriz.
  
  ## Add transfer coefficients to A matrix:
  alpha_2_1 <- pars[5]
  
  A[2,1] <- alpha_2_1*pars[1]
  
  
  Model14 <- Model_14(t = times, A = A, ivList = C0, inputFluxes = inputVector, 
                      initialValF = F0, inputFc = AtmFc_CC, c14DecayRate = -0.0001209681)
  
  C_pools <- getC(Model14) ; colnames(C_pools) <-  c("C_POM", "C_MAOM")
  
  C_14_pools <- getF14(Model14) ;  colnames(C_14_pools) <-  c("C_14_POM", "C_14_MAOM")
  
  C_total <- rowSums(getC(Model14)) # stock de SOC; suma todos los pools en cada año
  
  C_14_bulk <- getF14C(Model14) # delta14C del suelo en bulk. Promedio ponderado por masa del delta14C de todos los pools. 
  
  C_14_incub <- getF14R(Model14) # Firma delta14C en CO2. Promedio ponderado de la firma delta14C todos los pools#(getF14) por lo 
  #que se oxida en cada momento de cada pool (getReleaseFlux)
  
  model_output_2 <- data.frame(time=times, 
                               C_pools,
                               C_14_pools,
                               C_ha_mean= C_total,
                               C_14_bulk = C_14_bulk,
                               C_14_CO2 = C_14_incub)
  
  model_output_2$time <- as.numeric(model_output_2$time)
  
  model_output_2 <- model_output_2 %>% 
    dplyr::filter(time>1990)
  
  
  
  model_output_final <- rbind(model_output_1, model_output_2)
  
  
  return(data.frame(model_output_final))
  
  

}



num <- 3000

sens_C_POM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_POM")))
sens_C_MAOM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_MAOM")))
sens_C_ha_mean <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_ha_mean")))
sens_C14_bulk <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_14_bulk")))
sens_C14_POM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_14_POM")))
sens_C14_MAOM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_14_MAOM")))
sens_C_14_CO2 <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_14_CO2")))



export <- list("sens_C_POM" = sens_C_POM,
               "sens_C_MAOM" = sens_C_MAOM,
               "sens_C_ha_mean" = sens_C_ha_mean,
               "sens_C14_bulk" = sens_C14_bulk,
               "sens_C14_POM" = sens_C14_POM,
               "sens_C14_MAOM" = sens_C14_MAOM,
               "sens_C_14_CO2" = sens_C_14_CO2
) 


write.xlsx(export, "sens_var_outputs/sesns_var_output_CC.xlsx", colWidths = c("auto", "auto"),overwrite=TRUE) #Exporta en arhivo .xlsx la lista creada anteriormente.


