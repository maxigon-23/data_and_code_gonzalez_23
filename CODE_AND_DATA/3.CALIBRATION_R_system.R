#########################################################################################
############################ Crop - Pasture Rotation Model ##############################
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



#########################################################################################

#--------------------- Measured data

#-------- SOC anual stock
data_obs_C <- read_excel("DATA.xlsx", sheet = "C_anual")

data_obs_C <- data_obs_C %>% 
  filter(treatment=="R" & depth=="0-20") %>% 
  select(time, C_ha_mean, C_ha_SE)

data_obs_C$time <- as.numeric(data_obs_C$time)

#-------- bulk soil radiocarbon
data_obs_14_bulk <- read_excel("DATA.xlsx", sheet = "C14_bulk")

data_obs_14_bulk <- data_obs_14_bulk %>%
  group_by(treatment, time) %>% 
  summarise(C_14_bulk = weighted.mean(C_14_mean, stock_C),
            C_14_bulk_SE = weighted.mean(C_14_SE, stock_C)) %>% 
  filter(treatment=="R")

data_obs_14_bulk <- data_obs_14_bulk[ , c("time", "C_14_bulk", "C_14_bulk_SE")]
data_obs_14_bulk <- as.data.frame(data_obs_14_bulk)
data_obs_14_bulk$time <- as.numeric(data_obs_14_bulk$time)

#-------- incubation radiocarbon
data_obs_14_incub <- read_excel("DATA.xlsx", sheet = "C_14_incub")

data_obs_14_incub <- data_obs_14_incub %>%
  group_by(treatment, time) %>% 
  summarise(C_14_incub = weighted.mean(C_14_incub_mean, relative),
            C_14_incub_SE = weighted.mean(C_14_incub_SE, relative)) %>%
  filter(treatment=="R")

data_obs_14_incub <- data_obs_14_incub[ , c("time", "C_14_incub", "C_14_incub_SE")]
data_obs_14_incub <- as.data.frame(data_obs_14_incub)
data_obs_14_incub$time <- as.numeric(data_obs_14_incub$time)

#-------- Fraction partition
data_obs_C$C_POM <- data_obs_C$C_ha_mean * 0.135
data_obs_C$C_MAOM <- data_obs_C$C_ha_mean * ( 1 - 0.135)



#########################################################################################

#--------------------- Radiocarbon curves


hua <- data.frame(
  year = Hua2021$`SHZone1-2`[ , 1],
  delta14C = Hua2021$`SHZone1-2`[ , 2]
)

yrs=seq(1965,2019,by=1/4) 

nz2=spline(hua, xout=yrs)
nhz2=ts((nz2$y-1),start=1965, freq=4) 

m=ets(nhz2) 
f2=forecast(m,h=11*4)

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

#--------------------- Model preparation


times <- seq(0, 2021, by=1) 

F0 <- ConstFc(values=c(-50, -150), format="Delta14C") # 14C initial values

pools <- c(x1="Labile",
           x2="Stable"
)

C0 <- c(x1= 53 * 0.1, 
        x2= 53 * 0.9)

inputVector <- c(11.45 * 0.45,   #x1
                 0    #x2
)


#########################################################################################

#--------------------- Model autoparametrization



autoparam <- function(pars){
  # Define compartmental matrix
  ks <- pars[1:2]
  
  A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
  #de esa matriz.
  
  ## Add transfer coefficients to A matrix:
  alpha_2_1 <- pars[3]
  
  A[2,1] <- alpha_2_1*pars[1]
  
  
  Model14 <- Model_14(t = times, A = A, ivList = C0, inputFluxes = inputVector, 
                      initialValF = F0, inputFc = AtmFc, c14DecayRate = -0.0001209681)
  
  C_pools <- getC(Model14) ; colnames(C_pools) <-  c("C_POM", "C_MAOM")
  
  C_total <- rowSums(getC(Model14)) # stock de SOC; suma todos los pools en cada año
  
  C_14_bulk <- getF14C(Model14) # delta14C del suelo en bulk. Promedio ponderado por masa del delta14C de todos los pools. 
  
  C_14_incub <- getF14R(Model14) # Firma delta14C en CO2. Promedio ponderado de la firma delta14C todos los pools#(getF14) por lo 
  #que se oxida en cada momento de cada pool (getReleaseFlux)
  
  model_output <- data.frame(time=times, 
                             C_pools,
                             C_ha_mean= C_total,
                             C_14_bulk = C_14_bulk,
                             C_14_incub = C_14_incub)
  
  model_output <-  model_output %>% 
    filter(time>=1964)
  
  
  return(data.frame(model_output))
  
}





model_cost <- function(pars){
  
  model_output <- autoparam(pars)
  
  # Extrae la predicción del modelo para el tiempo correspondiente
  predicted_C_14_incub <- model_output[model_output$time==2021, "C_14_incub"]
  
  # Extrae la observación correspondiente
  observed_C_14_incub <- data_obs_14_incub["C_14_incub"]  # Reemplaza "value_column" con el nombre de la columna de tus datos
  
  # Calcula el coste como la diferencia absoluta entre la observación y la predicción
  cost1 <- abs(predicted_C_14_incub - observed_C_14_incub)
  
  cost2 <- modCost(model = model_output[, c("time", "C_POM")],
                   obs = as.data.frame(data_obs_C[, c("time", "C_POM")]),
                   cost = cost1)
  
  cost3<- modCost(model = model_output[, c("time", "C_MAOM")],
                   obs = as.data.frame(data_obs_C[, c("time", "C_MAOM")]),
                   cost = cost2)
  
  
  cost4<- modCost(model = model_output[, c("time", "C_ha_mean")],
                  obs = as.data.frame(data_obs_C[, c("time", "C_ha_mean", "C_ha_SE")]),
                  err = "C_ha_SE",
                  cost = cost3)
  
  
  cost5 <- modCost(model = model_output[model_output$time==2008 | model_output$time==2021, c("time", "C_14_bulk")],
                   obs = as.data.frame(data_obs_14_bulk),
                   err= "C_14_bulk_SE",
                   cost = cost4)
  

  return(cost5)
  
}


inipars <- c(
  1, # k1
  0.0018 ,  #k2
  0.0065 # a2_1
)



C_fit <- modFit(f=model_cost, p=inipars, method = "Marq", upper = c(2.5, 0.2, 0.2),
                lower = c(0.2,0.0005,0.0001))

C_fit$par

var0 <- C_fit$var_ms

MCMC1 <- modMCMC(f = model_cost, p = C_fit$par, niter = 25000, var0 = var0, jump = C_fit$par*0.05,
                wvar0=1, updatecov=50,
                upper = c(2.5, 0.05, 0.05),
                lower = c(0.2, 0,0))

list.save(MCMC1, 'lists_MCMC/MCMC_R.rdata')



