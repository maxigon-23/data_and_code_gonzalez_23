##########################################################################################
################################### T tests #############################################
#########################################################################################



rm(list=ls()) 

if(!require(readxl)){install.packages("readxl")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(data.table)){install.packages("data.table")}
if(!require(purrr)){install.packages("purrr")}
if(!require(plyr)){install.packages("plyr")}


#---------- Funcion de calculo de error estandar
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}



#----------- Aplicacion de prueba t por profundidad
t_test_depth <- function(tabla, var_obs) {
  # Obtener los niveles únicos de 'depth'
  unique_depth <- unique(tabla$depth)
  
  # Crear una lista para almacenar los resultados de la prueba t
  resultados <- list()
  
  # Iterar sobre los niveles de 'depth'
  for (d in unique_depth) {
    # Filtrar los datos por nivel de 'depth'
    subset_data <- tabla %>% filter(depth == d)
    
    # Extraer las variables necesarias
    var_observada <- subset_data[[var_obs]]
    tratamiento <- subset_data$treatment
    
    # Realizar la prueba t
    t_test <- t.test(var_observada ~ tratamiento)
    
    # Almacenar los resultados en la lista
    resultados[[paste(d)]] <- t_test
  }
  
  
  #nombres_depth <- names(resultados)
  p_valores <- sapply(resultados, function(resultado) resultado$p.value)
  
  means <- as.data.frame(sapply(resultados, function(resultado) resultado$estimate))
  means <- t(means)
  colnames(means) <- c("mean_CC", "mean_R")
  
  x <- as.data.frame(cbind(depth=rownames(means), p_valores, means))
  rownames(x) <- NULL
  
  x <- x %>%
    gather('mean_CC' ,'mean_R',  key= "treatment", value = "mean")
  
  tabla <- x %>%
    group_by(depth) %>%
    mutate(signif = ifelse(p_valores > 0.05, 'a', 'b'))
  
  tabla_gathered <- tabla %>%
    mutate(signif = ifelse(treatment == "mean_CC" & signif == "a", "a",
                           ifelse(treatment == "mean_R" & signif == "a", "a",
                                  ifelse(treatment == "mean_CC" & signif == "b", "a",
                                         ifelse(treatment == "mean_R" & signif == "b", "b", NA))))) %>% 
    dplyr::select(depth, treatment, mean, p_valores, signif) %>% 
    mutate(treatment = revalue(treatment, c("mean_CC" = "CC", "mean_R" = "R")))
  
  return(tabla_gathered)
}



#--------------------- Measured data

data <- read_excel("DATA.xlsx", sheet = "t_tests")


#--------------------- Incubation data


incub_0_20_por_plot <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>%
  dplyr::group_by(plot, treatment) %>% 
  dplyr::summarise(C_14_incub = weighted.mean(C_14_incub, incub_rate),
                   incub_rate = mean(incub_rate)
            )

resumen_incub_por_depth <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>%
  dplyr::group_by(treatment, depth) %>% 
  dplyr::summarise(C_14_incub_mean = mean(C_14_incub),
                   C_14_incub_SE = stderr(C_14_incub),
                   incub_rate_mean = mean(incub_rate),
                   incub_rate_SE = stderr(incub_rate),
  )


resumen_incub_0_20 <- incub_0_20_por_plot %>% 
  group_by(treatment) %>% 
  dplyr::summarise(C_14_incub_mean = mean(C_14_incub ) ,
            C_14_incub_SE = stderr(C_14_incub),
            incub_rate_mean = mean(incub_rate),
            incub_rate_SE = stderr(incub_rate)
            )


incub <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>% 
select(year, plot, depth, treatment, C_14_incub, incub_rate)


t_test_depth(incub, "C_14_incub")
t_test_depth(incub, "incub_rate")

t.test(data = incub_0_20_por_plot, C_14_incub ~ treatment )
t.test(data = incub_0_20_por_plot, incub_rate ~ treatment )

            
#--------------------- Bulk radiocarbon data

bulk_C14__0_20_por_plot <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2008" | year=="2021") %>%
  dplyr::group_by(year, plot, treatment) %>% 
  dplyr::summarise(C_14_bulk = weighted.mean(C_14_bulk, C_stock)
                   )

bulk_C14__0_20_por_plot_2008 <- bulk_C14__0_20_por_plot %>% 
  dplyr::filter(year == "2008")

bulk_C14__0_20_por_plot_2021 <- bulk_C14__0_20_por_plot %>% 
  dplyr::filter(year == "2021")

resumen_bulk_0_20 <- bulk_C14__0_20_por_plot %>% 
  group_by(year, treatment) %>% 
  dplyr::summarise(C_14_bulk_mean = mean(C_14_bulk ) ,
                   C_14_bulk_SE = stderr(C_14_bulk)
  )

resumen_bulk_por_depth <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2008" | year=="2021") %>%
  dplyr::group_by(year, treatment, depth) %>% 
  dplyr::summarise(C_14_bulk_mean = mean(C_14_bulk),
                   C_14_bulk_SE = stderr(C_14_bulk)
  )


C14_bulk_2008 <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2008") %>% 
  select(year, plot, depth, treatment, C_14_bulk)

C14_bulk_2021 <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>% 
  select(year, plot, depth, treatment, C_14_bulk)



t_test_depth(C14_bulk_2008, "C_14_bulk")
t_test_depth(C14_bulk_2021, "C_14_bulk")

t.test(data = bulk_C14__0_20_por_plot_2008, C_14_bulk ~ treatment )
t.test(data = bulk_C14__0_20_por_plot_2021, C_14_bulk ~ treatment )




#--------------------- POM data


POM__0_20_por_plot <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>%
  dplyr::group_by(plot, treatment) %>% 
  dplyr::summarise(C_stock_POM = sum(C_stock_POM)
  )


resumen_POM_0_20 <- POM__0_20_por_plot %>% 
  group_by( treatment) %>% 
  dplyr::summarise(C_stock_POM_mean = mean(C_stock_POM ) ,
                   C_stock_POM_SE = stderr(C_stock_POM)
  )

resumen_POM_por_depth <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>%
  dplyr::group_by(year, treatment, depth) %>% 
  dplyr::summarise(C_stock_POM_mean = mean(C_stock_POM),
                   C_stock_POM_SE = stderr(C_stock_POM)
  )


POM_data <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>% 
  select(year, plot, depth, treatment, C_stock_POM)


t_test_depth(POM_data, "C_stock_POM")

t.test(data = POM__0_20_por_plot, C_stock_POM ~ treatment )



#--------------------- MAOM data


MAOM__0_20_por_plot <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>%
  dplyr::group_by(plot, treatment) %>% 
  dplyr::summarise(C_stock_MAOM = sum(C_stock_MAOM)
  )


resumen_MAOM_0_20 <- MAOM__0_20_por_plot %>% 
  group_by( treatment) %>% 
  dplyr::summarise(C_stock_MAOM_mean = mean(C_stock_MAOM ) ,
                   C_stock_MAOM_SE = stderr(C_stock_MAOM)
  )

resumen_MAOM_por_depth <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>%
  dplyr::group_by(year, treatment, depth) %>% 
  dplyr::summarise(C_stock_MAOM_mean = mean(C_stock_MAOM),
                   C_stock_MAOM_SE = stderr(C_stock_MAOM)
  )


MAOM_data <- data %>% 
  filter(depth=="0-10" | depth=="10-20") %>% 
  filter(year=="2021") %>% 
  select(year, plot, depth, treatment, C_stock_MAOM)


t_test_depth(MAOM_data, "C_stock_MAOM")

t.test(data = MAOM__0_20_por_plot, C_stock_MAOM ~ treatment )



#------------ Data frame info 14C (0-20 cm)



resumen_POM_0_20
resumen_MAOM_0_20

resumen_bulk_0_20
resumen_incub_0_20






