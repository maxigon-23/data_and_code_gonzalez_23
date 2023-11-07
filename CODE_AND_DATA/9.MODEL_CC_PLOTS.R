#########################################################################################
############################### Model CC Plots ##########################################
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


#------------- Parametrization results

MCMC_R_results <- list.load('lists_MCMC/MCMC_R.rdata')
MCMC_CC_results <- list.load('lists_MCMC/MCMC_CC.rdata')


#--------------------- Measured data

#-------- SOC anual stock
data_obs_C <- read_excel("DATA.xlsx", sheet = "C_anual")

data_obs_C <- data_obs_C %>% 
  filter(treatment=="CC" & depth=="0-20") %>% 
  select(time, C_ha_mean, C_ha_SE)

data_obs_C$time <- as.numeric(data_obs_C$time)


#-------- bulk soil radiocarbon
data_obs_14_bulk <- read_excel("DATA.xlsx", sheet = "C_14_0_20")

data_obs_14_bulk <- data_obs_14_bulk %>%
  dplyr::filter(treatment=="CC") %>% 
  select(time, C_14_bulk_mean, C_14_bulk_SE)

data_obs_14_bulk <- data_obs_14_bulk[ , c("time", "C_14_bulk_mean", "C_14_bulk_SE")]
data_obs_14_bulk <- as.data.frame(data_obs_14_bulk)
data_obs_14_bulk$time <- as.numeric(data_obs_14_bulk$time)



#-------- incubation radiocarbon
data_obs_14_incub <- read_excel("DATA.xlsx", sheet = "C_14_0_20")

data_obs_14_incub <- data_obs_14_incub %>%
  dplyr::filter(treatment=="CC" & time =="2021") %>% 
  select(time, C_14_incub_mean, C_14_incub_SE)

data_obs_14_incub <- data_obs_14_incub[ , c("time", "C_14_incub_mean", "C_14_incub_SE")]
data_obs_14_incub <- as.data.frame(data_obs_14_incub)
data_obs_14_incub$time <- as.numeric(data_obs_14_incub$time)

#-------- Fraction partition

data_obs_C_fractions <- read_excel("DATA.xlsx", sheet = "Fractions_0_20")

data_obs_C_fractions <- data_obs_C_fractions %>% 
  dplyr::filter(treatment == "CC") %>%  
  summarise(time = time,
            POM_mean = C_stock_POM_mean,
            POM_SE = C_stock_POM_SE,
            MAOM_mean = C_stock_MAOM_mean,
            MAOM_SE = C_stock_MAOM_SE)



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


#########################################################################################


#########################################################################################

#----------- Run model




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


inputVector <- c(6.37 * 0.45,   #x1
                 0    #x2
)



#########################################################################################


autoparam_R <- function(pars){
  
  times <- seq(0, 2022, by=1) # Vector de tiempos 
  
  F0 <- ConstFc(values=c(-50, -150), format="Delta14C") # 14C initial values
  
  pools <- c(x1="Labile",
             x2="Stable"
  )
  
  C0 <- c(x1= 53 * 0.1, 
          x2= 53 * 0.9)
  
  inputVector <- c(11.45 *0.45 ,   #x1
                   0    #x2
  )
  # Define compartmental matrix
  ks <- pars[1:2]
  
  A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
  #de esa matriz.
  
  ## Add transfer coefficients to A matrix:
  alpha_2_1 <- pars[3]
  
  A[2,1] <- alpha_2_1*pars[1]
  
  
  Model14 <- Model_14(t = times, A = A, ivList = C0, inputFluxes = inputVector, 
                      initialValF = F0, inputFc = AtmFc_R, c14DecayRate = -0.0001209681)
  
  C_pools <- getC(Model14) ; colnames(C_pools) <-  c("C_POM", "C_MAOM")
  
  C_total <- rowSums(getC(Model14)) # stock de SOC; suma todos los pools en cada año
  
  C_14_pools <- getF14(Model14) ;  colnames(C_14_pools) <-  c("C_14_POM", "C_14_MAOM")
  
  C_14_bulk <- getF14C(Model14) # delta14C del suelo en bulk. Promedio ponderado por masa del delta14C de todos los pools. 
  
  C_14_incub <- getF14R(Model14) # Firma delta14C en CO2. Promedio ponderado de la firma delta14C todos los pools#(getF14) por lo 
  #que se oxida en cada momento de cada pool (getReleaseFlux)
  
  model_output <- data.frame(time=times, 
                             C_pools,
                             C_14_pools,
                             C_ha_mean= C_total,
                             C_14_bulk = C_14_bulk,
                             C_14_incub = C_14_incub)
  
  
  
  return(data.frame(model_output))
  
}

autoparam_CC <- function(pars){
  
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
  
  model_output_2 <- model_output_2 %>% 
    dplyr::filter(time>1990)
  
  
  
  model_output_final <- rbind(model_output_1, model_output_2)
  
  
  model_output_final <- rbind(model_output_1, model_output_2)
  
  
  return(data.frame(model_output_final))
  
}



#########################################################################################



pars_CC <- c(median(MCMC_CC_results$pars[,1]),
             median(MCMC_CC_results$pars[,2]),
             median(MCMC_CC_results$pars[,3]),
             median(MCMC_CC_results$pars[,4]),
             median(MCMC_CC_results$pars[,5]))


salida_CC <- autoparam_CC(pars_CC)



pars_R <- c(median(MCMC_R_results$pars[,1]),
          median(MCMC_R_results$pars[,2]),
          median(MCMC_R_results$pars[,3]))


salida_R <- autoparam_R(pars_R) %>% 
  filter(time<1964)


salida <- rbind(salida_R, salida_CC)

#########################################################################################



sens_C_POM <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C_POM") %>% 
                      dplyr::filter(x <= 1964),
                    read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C_POM")%>% 
                      dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C_POM",
         Mean = salida$C_POM)


sens_C_MAOM <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C_MAOM") %>% 
                       dplyr::filter(x <= 1964),
                     read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C_MAOM")%>% 
                       dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C_MAOM",
         Mean = salida$C_MAOM)


sens_C_ha_mean <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C_ha_mean") %>% 
                          dplyr::filter(x <= 1964),
                        read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C_ha_mean")%>% 
                          dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C_bulk",
         Mean = salida$C_ha_mean)


sens_C14_bulk <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C14_bulk") %>% 
                         dplyr::filter(x <= 1964),
                       read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C14_bulk")%>% 
                         dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C14_bulk",
         Mean = salida$C_14_bulk)


sens_C14_POM <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C14_POM") %>% 
                        dplyr::filter(x <= 1964),
                      read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C14_POM")%>% 
                        dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C14_POM",
         Mean = salida$C_14_POM)


sens_C14_MAOM <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C14_MAOM") %>% 
                         dplyr::filter(x <= 1964),
                       read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C14_MAOM")%>% 
                         dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C14_MAOM",
         Mean = salida$C_14_MAOM)


sens_C_14_CO2 <- rbind(read_excel("sens_var_outputs/sesns_var_output_R.xlsx", sheet = "sens_C_14_CO2") %>% 
                         dplyr::filter(x <= 1964),
                       read_excel("sens_var_outputs/sesns_var_output_CC.xlsx", sheet = "sens_C_14_CO2")%>% 
                         dplyr::filter(x > 1964) ) %>% 
  rename(time = x) %>% 
  mutate(fraction = "C_14_CO2",
         Mean = salida$C_14_incub)


C_stocks_fractions <- rbind(sens_C_ha_mean, sens_C_POM, sens_C_MAOM)
C_14_fractions <- rbind(sens_C14_bulk, sens_C14_POM, sens_C14_MAOM, sens_C_14_CO2)

rm(sens_C_POM, sens_C_MAOM, sens_C_ha_mean, sens_C14_bulk, sens_C14_POM, sens_C14_MAOM, sens_C_14_CO2)


# 
# #--------- Pegar data atmosferica
# 
# str(AtmFc1)
# 
# Atm_C14 <- AtmFc1 %>% 
#   dplyr::filter(year>0 & year<=2021)
# 
# 

# levels(factor(C_14_fractions$fraction))


C_stocks_fractions$fraction <- factor(C_stocks_fractions$fraction, 
                                      levels = c("C_bulk", "C_MAOM", "C_POM"),
                                      labels = c("Bulk",  "MAOM", "POM"))                    

str(C_stocks_fractions)

point_size_1 <- 0.9
point_size_2 <- 1.5


plot_stocks <- ggplot() +
  scale_x_continuous(limits = c(1950, 2022),
                     breaks = seq(from = 1950, to = 2022, by = 10)) +
  geom_ribbon(data=C_stocks_fractions, aes(x=time, ymin=q50 - Sd, ymax=q50 + Sd, fill=fraction), alpha =0.5) +
  geom_ribbon(data=C_stocks_fractions, aes(x=time, ymin=q05, ymax=q95, fill=fraction), alpha =0.35) +
  geom_line(data=C_stocks_fractions, aes(x=time, y=q50, color=fraction), size=0.4) + theme_bw() +
  scale_fill_manual(values = c("tan4","olivedrab4","royalblue1"))+
  scale_color_manual(values = c("tan4","olivedrab4","royalblue1")) +
  geom_point(data= data_obs_C, aes(x= time, y=C_ha_mean), color="tan4", size=point_size_1)+ 
  geom_errorbar(data= data_obs_C, aes(x=time, y= C_ha_mean, ymin=C_ha_mean-C_ha_SE  , ymax=C_ha_mean+C_ha_SE ),
                width=0.8, color="tan4")+
  geom_point(data= data_obs_C_fractions, aes(x= time, y=POM_mean), color="blue4", size=point_size_2)+ 
  geom_errorbar(data= data_obs_C_fractions, aes(x=time, y= POM_mean, ymin=POM_mean-POM_SE  , ymax=POM_mean+POM_SE ),
                width=0.8, color="blue4")+
  geom_point(data= data_obs_C_fractions, aes(x= time, y=MAOM_mean), color="darkgreen", size=point_size_2)+ 
  geom_errorbar(data= data_obs_C_fractions, aes(x=time, y= MAOM_mean, ymin=MAOM_mean-MAOM_SE  , ymax=MAOM_mean+MAOM_SE ),
                width=0.8, color="darkgreen")+   theme(legend.position = "bottom") + 
  ylab(expression(C~stock~(Mg~ha^{-1}))) +
  xlab("Time (years)")+  guides(color=guide_legend(title="Fraction"),
                                fill=guide_legend(title="Fraction"))

levels(C_14_fractions$fraction)
C_14_fractions$fraction <- factor(C_14_fractions$fraction, 
                                  levels = c("C14_bulk", "C14_MAOM", "C14_POM", "C_14_CO2"),
                                  labels = c("Bulk",  "MAOM", "POM", "Efflux"))

colnames(AtmFc1) <- c("time", "Atmosphere")
str(C_14_fractions)

Atm_modificado <- data.frame(time = AtmFc1$time,
                             Mean = NA,
                             Sd = NA,
                             Min = NA,
                             Max = NA,
                             q05 = NA,
                             q25 = NA,
                             q50 = AtmFc1$Atmosphere,
                             q75 = NA,
                             q95 = NA,
                             fraction = "Atmosphere") %>% 
  dplyr::filter(time>1950 & time<=2022)
  
C_14_fractions_plot <- rbind(C_14_fractions, Atm_modificado)

str(C_14_fractions_plot)



plot_C14 <- ggplot() +
  scale_x_continuous(limits = c(1950, 2022),
                     breaks = seq(from = 1950, to = 2022, by = 10)) +
  scale_y_continuous(limits = c(-120, 700),
                     breaks = seq(from = -100, to = 700, by = 100)) +
  geom_ribbon(data=C_14_fractions_plot, aes(x=time, ymin=q50 - Sd, ymax=q50 + Sd, fill=fraction), alpha =0.4) +
  geom_ribbon(data=C_14_fractions_plot, aes(x=time, ymin=q05, ymax=q95, fill=fraction), alpha =0.25) +
  geom_line(data=C_14_fractions_plot, aes(x=time, y=q50, color=fraction), size=0.4) +
  # scale_color_manual(values = c("black","blue", "coral2", "olivedrab4"))+
  # scale_fill_manual(values = c("black","blue", "coral2", "olivedrab4"))+
  theme_bw()+   theme(legend.position = "bottom") + 
  ylab(expression(paste(Delta^{14}, "C (\u2030)"))) +
  xlab("Time (years)")+  guides(color=guide_legend(title="Fraction"),
                                fill=guide_legend(title="Fraction"))+
  # geom_line(data=AtmFc1, aes(x=time, y=Atmosphere, color="Atmosphere"), size=1) +
  scale_fill_manual(values = c("brown","black","blue", "tan2", "olivedrab4"))+
  scale_color_manual(values = c("brown","black","blue", "tan2", "olivedrab4")) +
  geom_point(data= data_obs_14_bulk, aes(x= time, y=C_14_bulk_mean ), color="brown", size=point_size_2)+ 
  geom_errorbar(data= data_obs_14_bulk, aes(x=time, y= C_14_bulk_mean , ymin=C_14_bulk_mean -C_14_bulk_SE  , ymax=C_14_bulk_mean +C_14_bulk_SE ),
                width=1.4, color="brown")+
  geom_point(data= data_obs_14_incub, aes(x= time, y=C_14_incub_mean ), color="tan2", size=point_size_2)+ 
  geom_errorbar(data= data_obs_14_incub, aes(x=time, y= C_14_incub_mean , ymin=C_14_incub_mean -C_14_incub_SE    , ymax=C_14_incub_mean +C_14_incub_SE   ),
                width=1.4, color="tan2")


leyend_size <-  theme(legend.text = element_text(size = 7), 
                      legend.title = element_text(size = 9),
                      legend.key.size = unit(0.35, "cm"),
                      legend.spacing.y = unit(0.17, "cm"))

plot_grid <- ggarrange(plot_stocks + leyend_size, plot_C14 + leyend_size, 
                       ncol = 1, nrow = 2,
                       labels = c("(a)", "(b)"),
                       font.label=list(color="black",size=8))

ggsave(plot=plot_grid, "plots/grid_CC_model_output.png",device = "png" ,units="in", width=4, height=6.64, dpi=650)



