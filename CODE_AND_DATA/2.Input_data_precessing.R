#########################################################################################
############################### Input data preparation ##################################
#########################################################################################


rm(list=ls()) 

if(!require(readxl)){install.packages("readxl")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(data.table)){install.packages("data.table")}


#########################################################################################

#--------------------- Data

#-------- SOC anual stock
inputs <- read_excel("DATA.xlsx", sheet = "C_inputs")


a_factor <- c("Period",  "Plot", "Treatment")
inputs[, a_factor] <- lapply(inputs[ , a_factor] , factor)

inputs <- inputs[, c(1:9, 11:20)]

resumen <-  inputs %>% 
  filter(Treatment=="2" | Treatment=="5") %>% 
  group_by(Plot, Period) %>% 
  summarise(
    Treatment = Treatment,
    rango_años = paste(min(Year),"-",max(Year)),
    Total_input_sum = sum(Total_input),
    Dif_years = max(Year) - min(Year),
    Anual_input = Total_input_sum / Dif_years
  ) %>% 
  filter(row_number() == 1) %>%
  ungroup()


resumen_agrupado_dry_matter <- resumen %>% 
  group_by(Treatment) %>% 
  summarise(
    Anual_input_mean = mean(Anual_input),
    Anual_input_sd = sd(Anual_input)
  )

resumen_agrupado_C <- resumen %>% 
  group_by(Treatment) %>% 
  summarise(
    Anual_input_mean = mean(Anual_input) * 0.45,
    Anual_input_sd = sd(Anual_input) * 0.45
  )



