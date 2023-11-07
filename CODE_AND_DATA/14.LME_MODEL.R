#-------------------------------- Statistics ----------------------------------------------


rm(list=ls()) 




############################################################################################
#############################   Observed data preparation    ###############################
############################################################################################



if(!require(readxl)){install.packages("readxl")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(plyr)){install.packages("plyr")}
if(!require(data.table)){install.packages("data.table")}
if(!require(lme4)){install.packages("lme4")}
if(!require(lmerTest)){install.packages("lmerTest")}



stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}


#-------- SOC anual stock

data <- read_excel("DATA.xlsx", sheet = "t_tests") %>% 
  filter(treatment=="R" | treatment=="CC" & depth=="0-20") %>% 
  select(year, treatment, plot, depth, C_stock)


f <- c("year", "treatment", "plot", "depth")
data[, f] <- lapply(data[ , f] , factor)





model <- lmer(data= data, formula = C_stock ~ year + treatment + (1| plot), na.action=na.omit)

summary(model)


summary(model)

car::Anova(model)


qqnorm(resid(model))
plot(fitted(model), resid(model))
abline(0,0)

