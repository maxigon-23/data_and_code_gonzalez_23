#####################################################################################
########################### Age y transit time plots ################################
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
if(!require(ggpubr)){install.packages("ggpubr")}


#########################################################################################

#---------------- Upload model parametrizations

MCMC_R_results <- list.load('lists_MCMC/MCMC_R.rdata')
MCMC_CC_results <- list.load('lists_MCMC/MCMC_CC.rdata')

summary(MCMC_R_results)
summary(MCMC_CC_results)



#########################################################################################


#--------------- Model R

inputVector <- c(11.45 *0.45,   #x1
                 0    #x2
)

u <-matrix((inputVector), ncol=1)

pars_R <- c(median(MCMC_R_results$pars[,1]),
            median(MCMC_R_results$pars[,2]),
            median(MCMC_R_results$pars[,3]))


ks <- pars_R[1:2]

A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
#de esa matriz.

## Add transfer coefficients to A matrix:
alpha_2_1 <- pars_R[3]

A[2,1] <- alpha_2_1*pars_R[1]

R_Sist_age <- systemAge(A=A, u=u, a = seq(0, 3000, by=0.1), q = c(0.05, 0.5, 0.95))
R_Trans_time <- transitTime(A=A, u=u, a = seq(0, 3000, by=0.1), q = c(0.05, 0.5, 0.95))

#########################################################################################


#########################################################################################

#--------------- Model CC


inputVector <- c(6.37 *0.45,   #x1
                 0    #x2
)

u <-matrix((inputVector), ncol=1)

pars_CC <- c(median(MCMC_CC_results$pars[,1]),
             median(MCMC_CC_results$pars[,4]),
             median(MCMC_CC_results$pars[,5]))



ks <- pars_CC[1:2]

A <- -1 * diag(ks) # Crea automáticamente una matrix poniendo el vector en la diagonal
#de esa matriz.

## Add transfer coefficients to A matrix:
alpha_2_1 <- pars_CC[3]

A[2,1] <- alpha_2_1*pars_CC[1]

CC_Sist_age <- systemAge(A=A, u=u, a = seq(0, 3000, by=0.1), q = c(0.05, 0.5, 0.95))
CC_Trans_time <- transitTime(A=A, u=u, a = seq(0, 3000, by=0.1), q = c(0.05, 0.5, 0.95))


#########################################################################################

#########################################################################################

#-------- Unificar resultados



Ages_TT <- rbind(
  data.frame(
    time=  seq(0, 3000, by=0.1),
    System_Age = R_Sist_age$systemAgeDensity,
    POM_Age = R_Sist_age$poolAgeDensity[,1],
    MAOM_Age = R_Sist_age$poolAgeDensity[,2],
    Transit_time = R_Trans_time$transitTimeDensity,
    System = "R"
  ),
  data.frame(
    time=  seq(0, 3000, by=0.1),
    System_Age = CC_Sist_age$systemAgeDensity,
    POM_Age = CC_Sist_age$poolAgeDensity[,1],
    MAOM_Age = CC_Sist_age$poolAgeDensity[,2],
    Transit_time = CC_Trans_time$transitTimeDensity,
    System = "CC"
  )
)



##


##########################################################################################





R_age_and_TT_R <- read_excel("outputs_age_tt_error_prop/age_transit_dist_R.xlsx", sheet = "age_transit_dist_R")

resumen <- function(vector){
  data.frame(Mean_Age = mean(vector),
             Mean_Age_SD = sd(vector),
             Q1 = as.numeric(quantile(vector, probs = c(0.25))),
             Median = median(vector),
             Q3 = as.numeric(quantile(vector, probs = c(0.75))),
             LI = as.numeric(ci(vector, method = "ETI"))[2],
             LS = as.numeric(ci(vector, method = "ETI"))[3]
  )
}


R_System <- data.frame(
  Sistema = "R",
    rbind(
    cbind(Variable="Mean Age" ,resumen(R_age_and_TT_R$System_age)),
    cbind(Variable="POM Age" ,resumen(R_age_and_TT_R$POM_age)),
    cbind(Variable="MAOM Age" ,resumen(R_age_and_TT_R$MAOM_age)),
    cbind(Variable="Transit time" ,resumen(R_age_and_TT_R$Transit_time))
    )
)

CC_age_and_TT_CC <- read_excel("outputs_age_tt_error_prop/age_transit_dist_CC.xlsx", sheet = "age_transit_dist_CC")

CC_System <-  data.frame(
  Sistema = "CC",
    rbind(
    cbind(Variable="Mean Age" ,resumen(CC_age_and_TT_CC$System_age)),
    cbind(Variable="POM Age" ,resumen(CC_age_and_TT_CC$POM_age)),
    cbind(Variable="MAOM Age" ,resumen(CC_age_and_TT_CC$MAOM_age)),
    cbind(Variable="Transit time" ,resumen(CC_age_and_TT_CC$Transit_time))
  )
)


resumen_Age_TT <- rbind(R_System, CC_System)

str(R_age_and_TT_R)





resumen_mean_age <- resumen_Age_TT %>% 
  filter(Variable=="Mean Age")

resumen_MAOM_age <- resumen_Age_TT %>% 
  filter(Variable=="MAOM Age")


resumen_POM_age_R <- resumen_Age_TT %>% 
  filter(Variable=="POM Age" & Sistema =="R")

resumen_POM_age_CC <- resumen_Age_TT %>% 
  filter(Variable=="POM Age" & Sistema =="CC")

resumen_TT_R <- resumen_Age_TT %>% 
  filter(Variable=="Transit time" & Sistema =="R")

resumen_TT_CC <- resumen_Age_TT %>% 
  filter(Variable=="Transit time" & Sistema =="CC")
  
System_Age <- ggplot()+
  scale_y_continuous(limits = c(0, 0.008),
                       breaks = seq(from = 0, to = 0.008, by = 0.002))+
  theme_bw() +
  geom_vline(data= resumen_mean_age, aes(xintercept = Mean_Age, color = Sistema), linetype = "dashed",  size = 1) +  # Línea para el valor medio
  geom_vline(data= resumen_mean_age, aes(xintercept = LI, color = Sistema), linetype = "dotted",  size = 1, alpha = 0.5)+
  geom_vline(data= resumen_mean_age, aes(xintercept = LS, color = Sistema), linetype = "dotted",  size = 1, alpha = 0.5)+
  geom_rect(data = resumen_mean_age, aes(xmin = LI, xmax = LS, ymin = -Inf, ymax = Inf, fill = Sistema), alpha = 0.1) + # Área sombreada
  theme_bw() +
  guides(fill = FALSE, linetype = FALSE, color = FALSE)+
  geom_line(data= Ages_TT, aes(x= time, y=System_Age, color=System))+ scale_x_continuous(limits = c(0, 700),
                                                                                                   breaks = seq(from = 0, to = 700, by = 100))+
  geom_text(data= resumen_mean_age, aes(x=Mean_Age, y=Inf, color=Sistema, label=paste("Mean = ", base::signif(Mean_Age, digits=4))), size= 2.7, vjust=1.8, hjust=0.45) +
  scale_fill_manual(values = c("coral1","turquoise3")) +
  scale_color_manual(values = c("coral1","turquoise3")) +  ylab("Density distribution") +
  xlab("System age (years)")





POM_Age <- ggplot()+
  geom_line(data= Ages_TT, aes(x= time, y=POM_Age, color=System))+ scale_x_continuous(limits = c(0, 7),
                                                                          breaks = seq(from = 0, to = 7, by = 1))+
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(from = 0, to = 1, by = 0.25))+
  theme_bw()+
  geom_vline(data= resumen_POM_age_R, aes(xintercept = Mean_Age),color = "turquoise3", linetype = "dashed",  size = 1) +  # Línea para el valor medio
  geom_text(data= resumen_POM_age_R, aes(x=Mean_Age, y=Inf, label=paste("Mean = ", base::signif(Mean_Age, digits=3))), color = "turquoise3",  size= 2.7, vjust=1.8, hjust=-0.1) +
  geom_vline(data= resumen_POM_age_CC, aes(xintercept = Mean_Age),color = "coral1", linetype = "dashed",  size = 1) +  # Línea para el valor medio
  geom_text(data= resumen_POM_age_CC, aes(x=Mean_Age, y=Inf, label=paste("Mean = ", base::signif(Mean_Age, digits=3))), color = "coral1",  size= 2.7, vjust=4.5, hjust=-0.1) +
  scale_fill_manual(values = c("coral1","turquoise3")) +
  scale_color_manual(values = c("coral1","turquoise3")) +  ylab("Density distribution") +
  xlab("POM age (years)")






MAOM_Age <- ggplot()+
  scale_y_continuous(limits = c(0, 0.008),
                     breaks = seq(from = 0, to = 0.008, by = 0.002))+
  theme_bw() +
  geom_vline(data= resumen_MAOM_age, aes(xintercept = Mean_Age, color = Sistema), linetype = "dashed",  size = 1) +  # Línea para el valor medio
  geom_vline(data= resumen_MAOM_age, aes(xintercept = LI, color = Sistema), linetype = "dotted",  size = 1, alpha = 0.5)+
  geom_vline(data= resumen_MAOM_age, aes(xintercept = LS, color = Sistema), linetype = "dotted",  size = 1, alpha = 0.5)+
  geom_rect(data = resumen_MAOM_age, aes(xmin = LI, xmax = LS, ymin = -Inf, ymax = Inf, fill = Sistema), alpha = 0.1) + # Área sombreada
  theme_bw() +
  guides(fill = FALSE, linetype = FALSE, color = FALSE)+
  geom_line(data= Ages_TT, aes(x= time, y=MAOM_Age, color=System))+ scale_x_continuous(limits = c(0, 800),
                                                                                         breaks = seq(from = 0, to = 800, by = 100))+
  geom_text(data= resumen_MAOM_age, aes(x=Mean_Age, y=Inf, color=Sistema, label=paste("Mean = ", base::signif(Mean_Age, digits=4))), size= 2.7, vjust=1.8, hjust=0.45) +
  scale_fill_manual(values = c("coral1","turquoise3")) +
  scale_color_manual(values = c("coral1","turquoise3")) +  ylab("Density distribution") +
  xlab("MAOM age (years)")
  



Transit_time <- ggplot(data= Ages_TT)+
  geom_line(aes(x= time, y=Transit_time, color=System))+ scale_x_continuous(limits = c(0, 12),
                                                                       breaks = seq(from = 0, to = 12, by = 2))+
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(from = 0, to = 1, by = 0.25))+
  theme_bw()+
  
  geom_vline(data= resumen_TT_R, aes(xintercept = Mean_Age),color = "turquoise3", linetype = "dashed",  size = 1) +  # Línea para el valor medio
  geom_text(data= resumen_TT_R, aes(x=Mean_Age, y=Inf, label=paste("Mean = ", base::signif(Mean_Age, digits=5))), color = "turquoise3",  size= 2.7, vjust=1.8, hjust=0.45) +
  geom_vline(data= resumen_TT_CC, aes(xintercept = Mean_Age),color = "coral1", linetype = "dashed",  size = 1) +  # Línea para el valor medio
  geom_text(data= resumen_TT_CC, aes(x=Mean_Age, y=Inf, label=paste("Mean = ", base::signif(Mean_Age, digits=3))), color = "coral1",  size= 2.7, vjust=4.5, hjust=0.45) +
  scale_fill_manual(values = c("coral1","turquoise3")) +
  scale_color_manual(values = c("coral1","turquoise3")) +  ylab("Density distribution") +
  xlab("Transit time (years)")





############


arrange <- ggarrange(System_Age,
                     POM_Age,
                     MAOM_Age,
                     Transit_time,
                     ncol = 2, nrow = 2, common.legend = TRUE,
                     legend = "bottom",
                     labels = c("(a)", "(b)", "(c)", "(d)"),
                     font.label=list(color="black",size=8) )
                     
                     




ggsave(plot=arrange, "plots/age_and_TT.png", device = "png", units="cm", width=18, height=15.5, dpi=400)


R_Sist_age$meanSystemAge
CC_Sist_age$meanSystemAge

#---------------- Mean and Median transit time
R_Trans_time$quantiles
R_Trans_time$meanTransitTime

CC_Trans_time$quantiles
CC_Trans_time$meanTransitTime


#--------------- Calculos con edades
Ages_TT_R <- Ages_TT %>% 
  dplyr::filter(System == "R")

Syst_age_total_R <- sum(Ages_TT_R$System_Age)
POM_age_total_R <- sum(Ages_TT_R$POM_Age)
MAOM_age_total_R <- sum(Ages_TT_R$MAOM_Age)


Ages_TT_R %>% 
  filter(time<=59) %>% 
  summarise(prop= sum(System_Age)/Syst_age_total_R)

Ages_TT_R %>% 
  filter(time<=59) %>% 
  summarise(prop= sum(POM_Age)/POM_age_total_R)

Ages_TT_R %>% 
  filter(time>=100) %>% 
  summarise(prop= sum(MAOM_Age)/MAOM_age_total_R)

Ages_TT_R %>% 
  filter(time>=500) %>% 
  summarise(prop= sum(MAOM_Age)/MAOM_age_total_R)

Ages_TT_R %>% 
  filter(time>=1000) %>% 
  summarise(prop= sum(MAOM_Age)/MAOM_age_total_R)

