#########################################################################################
################################ Parameter distributions ################################
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

pars_R <-  data.frame(MCMC_R_results$pars)
pars_CC <-  data.frame(MCMC_CC_results$pars)

colnames(pars_R) <- c("k1", "k2", "alfa")
pars_R <- data.frame(pars_R, System = "R")

colnames(pars_CC) <- c("k1", "k2_1", "alfa_1", "k2_2", "alfa_2")
pars_CC <- data.frame(pars_CC, System = "CC")

ci_k1_R <- as.numeric(ci(pars_R$k1, method = "ETI")) #Credible interval alfa R System
ci_k2_R <- as.numeric(ci(pars_R$k2, method = "ETI")) #Credible interval alfa R System
ci_alfa_R <- as.numeric(ci(pars_R$alfa, method = "ETI")) #Credible interval alfa R System

ci_k1_CC <- as.numeric(ci(pars_CC$k1, method = "ETI")) #Credible interval alfa R System
ci_k2_1_CC <- as.numeric(ci(pars_CC$k2_1, method = "ETI")) #Credible interval alfa R System
ci_k2_2_CC <- as.numeric(ci(pars_CC$k2_2, method = "ETI")) #Credible interval alfa R System
ci_alfa_1_CC <- as.numeric(ci(pars_CC$alfa_1, method = "ETI")) #Credible interval alfa R System
ci_alfa_2_CC <- as.numeric(ci(pars_CC$alfa_2, method = "ETI")) #Credible interval alfa R System



#########################################################################################


plot_histogram <- function(data, par,x_axis_title, color){
  
  median_val <- median(data[[par]])
  
  mean_val <- mean(data[[par]])
  
  ggplot(data=data, aes(x = data[[par]])) +
    geom_histogram(aes( y= after_stat(count / sum(count))), 
                   color=color, fill = color, alpha = 0.15,  bins = 60) +
    geom_density(aes(y=after_stat(count / sum(count))), size=1.2,color="black", fill = color, alpha=0.5, n=60, adjust=3) +
    geom_vline(aes(xintercept=mean_val), 
               color="blue", linetype="dotted", size= 0.8) +
    geom_text(aes(x=mean_val, y=Inf, label=paste("Mean = ", base::signif(mean_val, digits=3))), size= 2.7, color="blue", vjust=1.8, hjust=-0.1) +
    geom_vline(aes(xintercept=median_val), 
               color="red", linetype="dashed", size= 0.8) +
    geom_text(aes(x=median_val, y=Inf, label=paste("Median = ", base::signif(median_val, digits=3))), size= 2.7, color="red", vjust=4, hjust=-0.1) +
    theme_bw() +
    labs(y= "Density", x= paste(x_axis_title)) +
    guides(fill="none")
}


#---------- k1 plot

k1_R <- plot_histogram(pars_R, "k1", x_axis_title="k1", "turquoise3"); k1_CC <- plot_histogram(pars_CC, "k1", x_axis_title="k1", "coral") 


scaling_x_k1 <-  scale_x_continuous(breaks=seq(0, 2, by = 0.5), limits = c(0, 1.5))

k1 <- ggarrange(k1_R + scaling_x_k1 , 
          k1_CC + scaling_x_k1 + theme(axis.title.y =element_blank()), 
          ncol = 2, nrow = 1,
          common.legend = TRUE, legend="bottom",
          labels = c("R system", "CC system"),
          label.x = 0.7, # Cambia para ajustar la posición en el eje x
          label.y = 0.29, # Cambia para ajustar la posición en el eje y
          font.label = list(size = 7,  color = "black", face = "bold")) # Ajusta según tus preferencias) 


ci_k1_R <- ci(pars_R$k1, method = "ETI") #Credible interval k1 R System
ci_k1_CC <- ci(pars_CC$k1, method = "ETI") #Credible interval k1 CC System

IC_k1 <- data.frame(System = c("R" , "CC"),
                    LI = c(as.numeric(ci_k1_R[2]), as.numeric(ci_k1_CC[2])),
                    LS = c(as.numeric(ci_k1_R[3]), as.numeric(ci_k1_CC[3])),
                    mean = c(mean(pars_R$k1), mean(pars_CC$k1)) ) # General table for Credible intervals in both systems


IC_k1$System = factor(IC_k1$System, levels=c('R','CC'))


inset_k1 <- ggplot(data= IC_k1, aes(x= System, color=System)) +
  geom_linerange(aes(ymin= LI, ymax= LS),
                 position = position_dodge(width = 3), lwd = 1) +
  geom_pointrange(aes(x = System, y = mean, ymin = LI,
                      ymax = LS), position = position_dodge(width = 1/2),
                  shape = 15, size=0.7, fill = "WHITE", lwd = 1/2)+ theme_bw() + labs( y="", x="") + guides(color="none")+
  scale_color_manual(values=c("R"="turquoise3", "CC"="coral")) 


k1_plus_comparison <- ggdraw() +
  draw_plot(k1, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(inset_k1, x = 0.79, y = 0.28, width = 0.19, height = 0.52)





#---------- k2 plot


k2_R <- plot_histogram(pars_R, "k2", x_axis_title="k2",  "turquoise3"); k2_1_CC <- plot_histogram(pars_CC, "k2_1", x_axis_title="k2 (Period 1)","coral"); 
k2_2_CC <- plot_histogram(pars_CC, "k2_2", x_axis_title="k2 (Period 2)", "coral4")

scaling_x_k2_R <-  scale_x_continuous(breaks=seq(0.0014, 0.0030, by = 0.0004), limits = c(0.0014, 0.0024))
scaling_x_k2_CC <-  scale_x_continuous(breaks=seq(0.0014, 0.0154, by = 0.007), limits = c(0.0010, 0.017))

k2 <- ggarrange(k2_R + scaling_x_k2_R , 
                k2_1_CC + scaling_x_k2_CC + theme(axis.title.y=element_blank()),
                k2_2_CC + scaling_x_k2_CC + theme(axis.title.y=element_blank()),
                 ncol = 3, nrow = 1,
          common.legend = TRUE, legend="bottom",
          labels = c("R system", "CC system P1", "CC system P2"),
          label.x = 0.47, # Cambia para ajustar la posición en el eje x
          label.y = 0.84, # Cambia para ajustar la posición en el eje y
          font.label = list(size = 7, face = "bold", color = "black")) # Ajusta según tus preferencias)


ci_k2_R <- ci(pars_R$k2, method = "ETI") #Credible interval k2 R System
ci_k2_1_CC <- ci(pars_CC$k2_1, method = "ETI") #Credible interval k2 CC System
ci_k2_2_CC <- ci(pars_CC$k2_2, method = "ETI") #Credible interval k2 CC System

IC_k2 <- data.frame(System = c("R" , "CC P 1", "CC P 2"),
                    LI = c(as.numeric(ci_k2_R[2]), as.numeric(ci_k2_1_CC[2]), as.numeric(ci_k2_2_CC[2])),
                    LS = c(as.numeric(ci_k2_R[3]), as.numeric(ci_k2_1_CC[3]), as.numeric(ci_k2_2_CC[3])),
                    mean = c(mean(pars_R$k2), mean(pars_CC$k2_1), mean(pars_CC$k2_2)))  # General table for Credible intervals in both systems


IC_k2$System = factor(IC_k2$System, levels=c('R','CC P 1', 'CC P 2' ))


inset_k2 <- ggplot(data= IC_k2, aes(x= System, color=System)) +
  geom_linerange(aes(ymin= LI, ymax= LS),
                 position = position_dodge(width = 3), lwd = 1) +
  geom_pointrange(aes(x = System, y = mean, ymin = LI,
                      ymax = LS), position = position_dodge(width = 1/2),
                  shape = 15, size=0.7, fill = "WHITE", lwd = 1/2)+ theme_bw() + labs( y="", x="") + guides(color="none")+
  scale_color_manual(values=c("R"="turquoise3", "CC P 1"="coral", "CC P 2" ="coral4" )) +
  scale_x_discrete(labels = c("R" = "R", "CC P 1" = "CC\nP 1", "CC P 2" = "CC\nP 2"))


k2_plus_comparison <- ggdraw() +
  draw_plot(k2, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(inset_k2 + theme(axis.text.x = element_text(size = 5.5),
                             axis.text.y = element_text(size = 5.5)), x = 0.16, y = 0.22, width = 0.16, height = 0.55)




#------- alfa plot



alfa_R <- plot_histogram(pars_R, "alfa",x_axis_title= "alpha", color="turquoise3")+
   labs(y= "Density", x= expression(alpha), "Period 1") ;

alfa_1_CC <- plot_histogram(pars_CC, "alfa_1",x_axis_title="Alfa P1", color="coral") +
  labs(y= "Density", x= expression(alpha~(Period~1))) 

alfa_2_CC <- plot_histogram(pars_CC, "alfa_2",x_axis_title="Alfa P2", color= "coral4")+
  labs(y= "Density", x= expression(alpha~(Period~2))) 

scaling_x_alfa_R <-  scale_x_continuous(breaks=seq(0.012, 0.020, by = 0.002), limits = c(0.012, 0.020))
scaling_x_alfa_CC <-  scale_x_continuous(breaks=seq(0, 0.12, by = 0.03), limits = c(0, 0.12))

alfa <- ggarrange(alfa_R + scaling_x_alfa_R  , 
                  alfa_1_CC + scaling_x_alfa_CC  + theme(axis.title.y=element_blank()),
                  alfa_2_CC + scaling_x_alfa_CC  + theme(axis.title.y=element_blank()),
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend="bottom",
                  labels = c("R system", "CC system P1", "CC system P2"),
                  label.x = 0.55, # Cambia para ajustar la posición en el eje x
                  label.y = 0.29, # Cambia para ajustar la posición en el eje y
                  font.label = list(size = 7, face = "bold", color = "black") # Ajusta según tus preferencias
)


ci_alfa_R <- ci(pars_R$alfa, method = "ETI") #Credible interval alfa R System
ci_alfa_1_CC <- ci(pars_CC$alfa_1, method = "ETI") #Credible interval alfa CC System
ci_alfa_2_CC <- ci(pars_CC$alfa_2, method = "ETI") #Credible interval alfa CC System

IC_alfa <- data.frame(System = c("R" , "CC P 1", "CC P 2"),
                    LI = c(as.numeric(ci_alfa_R[2]), as.numeric(ci_alfa_1_CC[2]), as.numeric(ci_alfa_2_CC[2])),
                    LS = c(as.numeric(ci_alfa_R[3]), as.numeric(ci_alfa_1_CC[3]), as.numeric(ci_alfa_2_CC[3])),
                    mean = c(mean(pars_R$alfa), mean(pars_CC$alfa_1), mean(pars_CC$alfa_2)))  # General table for Credible intervals in both systems


IC_alfa$System = factor(IC_alfa$System, levels=c('R','CC P 1', 'CC P 2' ))


inset_alfa <- ggplot(data= IC_alfa, aes(x= System, color=System)) +
  geom_linerange(aes(ymin= LI, ymax= LS),
                 position = position_dodge(width = 3), lwd = 1) +
  geom_pointrange(aes(x = System, y = mean, ymin = LI,
                      ymax = LS), position = position_dodge(width = 1/2),
                  shape = 15, size=0.7, fill = "WHITE", lwd = 1/2)+ theme_bw() + labs( y="", x="") +
  scale_color_manual(values=c("R"="turquoise3", "CC P 1"="coral", "CC P 2" ="coral4" ),
                     labels = c("R", "CC\nP 1", "CC\nP 2")) + guides(color="none") +
  scale_x_discrete(labels = c("R" = "R", "CC P 1" = "CC\nP 1", "CC P 2" = "CC\nP 2"))



alfa_plus_comparison <- ggdraw() +
  draw_plot(alfa, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(inset_alfa + theme(axis.text.x = element_text(size = 5.5),
                             axis.text.y = element_text(size = 5.5)), x = 0.491, y = 0.284, width = 0.16, height = 0.55)




############


arrange <- ggarrange(k1_plus_comparison,
          k2_plus_comparison,
          alfa_plus_comparison,
          ncol = 1, nrow = 3, common.legend = TRUE

          
)

ggsave(plot=arrange, "plots/parameter_distributions.png", device = "png", units="cm", width=19, height=21, dpi=400)


