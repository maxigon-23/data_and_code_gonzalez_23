#########################################################################################
################################### Model R SensVar  ####################################
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

summary(MCMC_R_results)
summary(MCMC_CC_results)


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



#--------------------- MODEL

times <- seq(0, 2022, by=1) # Vector de tiempos 

F0 <- ConstFc(values=c(-50, -150), format="Delta14C") # 14C initial values

pools <- c(x1="Labile",
           x2="Stable"
)

C0 <- c(x1= 53 * 0.1, 
        x2= 53 * 0.9)

inputVector <- c(11.45 * 0.45,   #x1
                 0    #x2
)




fun_error_prop <- function(pars){

  ks <- pars[1:2]

  A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
  #de esa matriz.

  ## Add transfer coefficients to A matrix:
  alpha_2_1 <- pars[3]

  A[2,1] <- alpha_2_1*pars[1]


  Model14 <- Model_14(t = times, A = A, ivList = C0, inputFluxes = inputVector,
                      initialValF = F0, inputFc = AtmFc, c14DecayRate = -0.0001209681)


  C_pools <- getC(Model14) ; colnames(C_pools) <-  c("C_POM", "C_MAOM")

  C14_pools <- getF14(Model14)

  salida <- data.frame(time=times,
                       C_pools,
                       C_ha_mean= rowSums(getC(Model14)),
                       C14_bulk = getF14C(Model14),
                       C14_POM = C14_pools[,1],
                       C14_MAOM = C14_pools[,2],
                       C_14_CO2 = getF14R(Model14) )

  return(salida)

}


pars <- MCMC_R_results$pars



num <- 1000

sens_C_POM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_POM")))
sens_C_MAOM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_MAOM")))
sens_C_ha_mean <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_ha_mean")))
sens_C14_bulk <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C14_bulk")))
sens_C14_POM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C14_POM")))
sens_C14_MAOM <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C14_MAOM")))
sens_C_14_CO2 <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_14_CO2")))



export <- list("sens_C_POM" = sens_C_POM,
               "sens_C_MAOM" = sens_C_MAOM,
               "sens_C_ha_mean" = sens_C_ha_mean,
               "sens_C14_bulk" = sens_C14_bulk,
               "sens_C14_POM" = sens_C14_POM,
               "sens_C14_MAOM" = sens_C14_MAOM,
               "sens_C_14_CO2" = sens_C_14_CO2
) 

write.xlsx(export, "sens_var_outputs/sesns_var_output_R.xlsx", colWidths = c("auto", "auto"),overwrite=TRUE) 



