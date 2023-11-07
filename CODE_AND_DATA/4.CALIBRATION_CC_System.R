#########################################################################################
########################## Continuous cropping system Model #############################
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
  filter(treatment=="CC" & depth=="0-20") %>% 
  select(time, C_ha_mean, C_ha_SE)

data_obs_C$time <- as.numeric(data_obs_C$time)

#-------- bulk soil radiocarbon
data_obs_14_bulk <- read_excel("DATA.xlsx", sheet = "C14_bulk")

data_obs_14_bulk <- data_obs_14_bulk %>%
  group_by(treatment, time) %>% 
  summarise(C_14_bulk = weighted.mean(C_14_mean, stock_C),
            C_14_bulk_SE = weighted.mean(C_14_SE, stock_C)) %>% 
  filter(treatment=="CC")

data_obs_14_bulk <- data_obs_14_bulk[ , c("time", "C_14_bulk", "C_14_bulk_SE")]
data_obs_14_bulk <- as.data.frame(data_obs_14_bulk)
data_obs_14_bulk$time <- as.numeric(data_obs_14_bulk$time)
data_obs_14_bulk <- data.frame(time=c(2008, 2021),
                               C_14_bulk = c(mean(data_obs_14_bulk$C_14_bulk)),
                               C_14_bulk_SE =c(mean(data_obs_14_bulk$C_14_bulk_SE)))


#-------- incubation radiocarbon
data_obs_14_incub <- read_excel("DATA.xlsx", sheet = "C_14_incub")

data_obs_14_incub <- data_obs_14_incub %>%
  group_by(treatment, time) %>% 
  summarise(C_14_incub = weighted.mean(C_14_incub_mean, relative)) %>% 
  filter(treatment=="CC")

data_obs_14_incub <- data_obs_14_incub[ , c("time", "C_14_incub")]
data_obs_14_incub <- as.data.frame(data_obs_14_incub)
data_obs_14_incub$time <- as.numeric(data_obs_14_incub$time)

#-------- Fraction partition
data_fracciones <- data.frame(time = 1964,
                              C_POM = 52.6098147596899 * 0.135,
                              C_MAOM = 52.6098147596899 * (1 - 0.135))


data_fracciones <- rbind(data_fracciones,
                         data.frame(
                           time = 2021,
                           C_POM = 40.9458888 * 0.116,
                           C_MAOM = 40.9458888 * (1 - 0.116))
)




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

AtmFc_R <- BoundFc(map=AtmFc1, format="Delta14C", lag=3) 
AtmFc_CC <- BoundFc(map=AtmFc1, format="Delta14C", lag=0.5) 


rm(bc, f2, hua, m, nhz2, nz2, preandpost)


#########################################################################################



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


MCMC_R_result <- list.load('lists_MCMC/MCMC_R.rdata')


pars <- MCMC_R_result$bestpar


ks <- pars[1:2]

A <- -1 * diag(ks) 


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

#--------------- Plots


C_14 <- a %>%
  gather('Labile' , 'Stable', 'Efflux', 'Bulk soil',
         key= "Pool", value = "C_14")

C_14$Pool <- factor(C_14$Pool)



C_stock <- b %>% 
  gather('Labile' , 'Stable',
         key= "Pool", value = "C_stock")

C_stock$Pool <- factor(C_stock$Pool)



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



#########################################################################################

#--------------------- Model autoparametrization


autoparam <- function(pars){
  
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
                             C_14_pools,
                             C_ha_mean= C_total,
                             C_14_bulk = C_14_bulk,
                             C_14_incub = C_14_incub)
  

  
  #------------ Fase 1990 - 2022
  
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
                             C_14_incub = C_14_incub)
  
  model_output_2$time <- as.numeric(model_output_2$time)
  
  model_output_2 <- model_output_2 %>% 
    dplyr::filter(time>1990)
  
  
  
  model_output_final <- rbind(model_output_1, model_output_2)
  
  
  return(data.frame(model_output_final))
  
}





model_cost <- function(pars){
  
  model_output <- autoparam(pars)
  
  # Extrae la predicción del modelo para el tiempo correspondiente
  predicted_C_14_incub <- model_output[model_output$time==2021, "C_14_incub"]
  
  # Extrae la observación correspondiente
  observed_C_14_incub <- data_obs_14_incub["C_14_incub"]  # Reemplaza "value_column" con el nombre de la columna de tus datos
  
  # Calcula el coste como la diferencia absoluta entre la observación y la predicción
  cost1 <- abs(predicted_C_14_incub - observed_C_14_incub)
  
  cost2 <- modCost(model = model_output[ , c("time", "C_ha_mean")],
                   obs = as.data.frame(data_obs_C[, c("time", "C_ha_mean", "C_ha_SE")]),
                   err = "C_ha_SE",
                   cost = cost1)
  
  cost3 <- modCost(model = model_output[model_output$time==2008 | model_output$time==2021, c("time", "C_POM")],
                   obs = as.data.frame(data_fracciones[, c("time", "C_POM")]),
                   cost = cost2)
  
  cost4 <- modCost(model = model_output[model_output$time==2008 | model_output$time==2021, c("time", "C_MAOM")],
                  obs = as.data.frame(data_fracciones[, c("time", "C_MAOM")]),
                  cost = cost3)
  
  cost5 <- modCost(model = model_output[model_output$time==2008 | model_output$time==2021, c("time", "C_14_bulk")],
                   err = "C_14_bulk_SE",
                   obs = as.data.frame(data_obs_14_bulk),
                   cost = cost4)
  
  
  return(cost5)
  
}




inipars <- c(
  1.2752462387, # k1 etapa 1
  0.0096, #k2 etapa 1
  0.001, #alfa etapa 1
  0.0024, #k2 etapa 2
  0.006) #alfa etapa 2



C_fit <- modFit(f=model_cost, p=inipars, method = "Marq",  upper = c(2.5, 0.08, 0.08, 0.08, 0.08),
                lower = c(0.2, 0, 0, 0, 0))

C_fit$par

var0 <- C_fit$var_ms

MCMC1 <- modMCMC(f = model_cost, p = C_fit$par, niter = 25000, var0 = var0, jump = C_fit$par*0.05,
                 wvar0=1, updatecov=50,
                 upper = c(2.5, 0.08, 0.08, 0.08, 0.08),
                 lower = c(0.2, 0, 0, 0, 0))


list.save(MCMC1, 'lists_MCMC/MCMC_CC.rdata')



