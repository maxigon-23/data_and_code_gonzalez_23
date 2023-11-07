#####################################################################################
########################## Efficiencies calculations ################################
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

MCMC_R_result <- list.load('lists_MCMC/MCMC_R.rdata')

MCMC_CC_result <- list.load('lists_MCMC/MCMC_CC.rdata')

summary(MCMC_CC_result)
summary(MCMC_R_result)

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



# R system
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
    
    C_release <- getReleaseFlux(Model14) #This function computes carbon release from each pool of the given model as funtion of time
    colnames(C_release) <-  c("C_POM_release", "C_MAOM_release")
    
    C_total <- rowSums(getC(Model14)) # stock de SOC; suma todos los pools en cada año
    
    C_14_bulk <- getF14C(Model14) # delta14C del suelo en bulk. Promedio ponderado por masa del delta14C de todos los pools. 
    
    C_14_incub <- getF14R(Model14) # Firma delta14C en CO2. Promedio ponderado de la firma delta14C todos los pools#(getF14) por lo 
    #que se oxida en cada momento de cada pool (getReleaseFlux)
    
    model_output <- data.frame(time=times, 
                               C_pools,
                               C_ha_mean= C_total,
                               C_14_bulk = C_14_bulk,
                               C_14_incub = C_14_incub,
                               C_release)
    
    model_output <-  model_output %>% 
      filter(time>=1964)
    
    
    model_output$total_release <- model_output$C_POM_release +  model_output$C_MAOM_release
    model_output$stabilization_flux <-  model_output$C_POM * alpha_2_1*pars[1]
    model_output$mineraliz_prop <- model_output$total_release / model_output$C_ha_mean
    model_output$mineraliz_MAOM_prop <- model_output$C_MAOM_release / model_output$C_MAOM
    
    model_output$stabilization_efficiency <- model_output$stabilization_flux / inputVector[1]


    

  
  return(model_output)
    
    
}

#########################################################################################

            

pars <- MCMC_R_result$pars


num <- 1000

sens_stabilization_flux <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("stabilization_flux")))
sens_stabilization_efficiency <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("stabilization_efficiency")))
sens_total_release <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("total_release")))
sens_mineraliz_prop <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("mineraliz_prop")))
sens_C_MAOM_release <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_MAOM_release")))
sens_mineraliz_MAOM_prop <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("mineraliz_MAOM_prop")))

sens_C_POM_release <- summary(sensRange(num=num, func=fun_error_prop, parInput=pars, sensvar=c("C_POM_release")))



export <- list("sens_stabilization_flux" = sens_stabilization_flux,
               "sens_stabilization_efficiency" = sens_stabilization_efficiency,
               "sens_total_release" = sens_total_release,
               "sens_C_POM_release" = sens_C_POM_release,
               "sens_mineraliz_prop" = sens_mineraliz_prop,
               "sens_C_MAOM_release" = sens_C_MAOM_release,
               "sens_mineraliz_MAOM_prop" = sens_mineraliz_MAOM_prop
) #Crea el objeto "salida" (elegir cualquier nombre) en el que se listan las

write.xlsx(export, "sens_var_outputs/sesns_var_output_R_tasas.xlsx", colWidths = c("auto", "auto"),overwrite=TRUE) #Exporta en arhivo .xlsx la lista creada anteriormente.

