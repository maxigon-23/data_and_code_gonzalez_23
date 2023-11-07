################################# Cimate plots #######################################

rm(list=ls()) 
ls() 

if(!require(readtext)){install.packages("readtext")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(readxl)){install.packages("readxl")}
if(!require(plyr)){install.packages("plyr")}
if(!require(xts)){install.packages("xts")}
if(!require(lubridate)){install.packages("lubridate")}
if(!require(Rmisc)){install.packages("Rmisc")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(cowplot)){install.packages("cowplot")}

######################################################################################

temp <- read_excel("DATA.xlsx", sheet = "CLIMATE")

precip <- read_excel("DATA.xlsx", sheet = "precip_68_2022")

temp$T_mean <- (temp$T_Max + temp$T_Min ) /2

#----------------- Creación de datos promedio -----------------------------

#Genera x ordenando los datos de tabla1 por fecha (ascendente) y a la columna de fecha
#le asigna un formato POSIXct para que luego pueda ser leido por la funcion que genera el objeto xts:
x <- precip[order(as.Date(precip$fecha, format = "%Y/%m/%d")),]

#Creamos el objeto xts a partir de la tabla con los datos:
tabla_xts <- xts(x = x, order.by = x$fecha)

tabla_xts$fecha <- NULL

storage.mode(tabla_xts) <- "numeric"

precip_month <- split(tabla_xts$precipitacion, f="months")
precip_month_sum <- lapply(precip_month, FUN = cumsum)
precip_MES <- do.call(rbind, precip_month_sum)


#Ahora me quedo con el ultimo dia de cada mes (el que tiene la acumulacion del mes):
#creando un vector de index_mes con las posiciones de los ultimos dias de cada mes para despues usarlo de
#index de corte:
index_mes <- endpoints(precip_MES, on = "months")
precip_MES1 <- precip_MES[index_mes]
names(precip_MES1) <- "Precipitacion"


#-------- Evapotranspiracion

evapo_month <- split(tabla_xts$evapotranspiracion, f="months")
evapo_monthh_sum <- lapply(evapo_month, FUN = cumsum)
evapo__MES <- do.call(rbind, evapo_monthh_sum)


index_mes <- endpoints(evapo__MES, on = "months")
evapo__MES1 <- evapo__MES[index_mes]
names(evapo__MES1) <- "Evapotranspiracion"



#------------------------------



#------------------ T_Max
x <- temp[order(as.Date(temp$fecha, format = "%Y/%m/%d")),]

#Creamos el objeto xts a partir de la tabla con los datos:
tabla_xts_temp <- xts(x = x, order.by = x$fecha)

#Borramos la columna fecha para que no este duplicada:
tabla_xts_temp$fecha <- NULL

storage.mode(tabla_xts_temp) <- "numeric"

monthly_split <- split(tabla_xts_temp$T_Max, f = "months")

means <- lapply(monthly_split, FUN = mean)
index <- endpoints(tabla_xts_temp, on = "months")
tabla_prueba <- tabla_xts_temp[index]
index1 <- index(tabla_prueba)

temps_monthly <- as.xts(as.numeric(means), order.by = index1)
names(temps_monthly) <- "TMAX_MEAN"


#------------------ T_Min

monthly_split <- split(tabla_xts_temp$T_Min, f = "months")

means <- lapply(monthly_split, FUN = mean)
index <- endpoints(tabla_xts_temp, on = "months")
tabla_prueba <- tabla_xts_temp[index]
index1 <- index(tabla_prueba)

temps_monthly_MIN <- as.xts(as.numeric(means), order.by = index1)
names(temps_monthly_MIN) <- "TMIN_MEAN"

temps_monthly <- cbind(temps_monthly, temps_monthly_MIN$TMIN_MEAN)




#------------------ T_Mean

monthly_split <- split(tabla_xts_temp$T_mean, f = "months")

means <- lapply(monthly_split, FUN = mean)
index <- endpoints(tabla_xts_temp, on = "months")
tabla_prueba <- tabla_xts_temp[index]
index1 <- index(tabla_prueba)

temps_monthly_MEAN <- as.xts(as.numeric(means), order.by = index1)
names(temps_monthly_MEAN) <- "T_MEAN"

temps_monthly <- cbind(temps_monthly, temps_monthly_MEAN$T_MEAN)



#------------------------------------------------------------------

a <- index(temps_monthly)

b <- month(as.POSIXct(a, format="%Y/%m/%d"))
b <- factor(b)


y <- cbind(b, temps_monthly)

xx <- data.frame(date=index(y), coredata(y))
xx$b <- factor(xx$b)

RESUMEN_TMAX <- ddply(xx,~b,summarise,mean=mean(TMAX_MEAN),sd=sd(TMAX_MEAN))
names(RESUMEN_TMAX) <- c("MES", "TMAX_MEAN", "TMAX_SD")

RESUMEN_TMIN <- ddply(xx,~b,summarise,mean=mean(TMIN_MEAN),sd=sd(TMIN_MEAN))
names(RESUMEN_TMIN) <- c("MES", "TMIN_MEAN", "TMIN_SD")

RESUMEN_TMEAN <- ddply(xx,~b,summarise,mean=mean(T_MEAN),sd=sd(T_MEAN))
names(RESUMEN_TMEAN) <- c("MES", "TMEAN_MEAN", "TMEAN_SD")


precip_MES1$b <- month(as.POSIXct(index(precip_MES1), format="%Y/%m/%d"))
precip_MES1$b <- factor(precip_MES1$b)

precip_MES1_df <- as.data.frame(precip_MES1)

evapo__MES1$b <- month(as.POSIXct(index(evapo__MES1), format="%Y/%m/%d"))
evapo__MES1$b <- factor(evapo__MES1$b)

evapo__MES1_df <- as.data.frame(evapo__MES1)


RESUMEN_PREC <- ddply(precip_MES1_df,~b,summarise,mean=mean(Precipitacion),sd=sd(Precipitacion))
names(RESUMEN_PREC) <- c("MES", "PREC_MEAN", "PREC_SD")

RESUMEN_EVAPO <- ddply(evapo__MES1_df,~b,summarise,mean=mean(Evapotranspiracion, na.rm=TRUE),sd=sd(Evapotranspiracion, na.rm=TRUE))
names(RESUMEN_EVAPO) <- c("MES", "EVAPO_MEAN", "EVAPO_SD")


RESUMEN <- cbind(RESUMEN_TMAX, RESUMEN_TMIN$TMIN_MEAN,
                   RESUMEN_TMIN$TMIN_SD, RESUMEN_TMEAN$TMEAN_MEAN, RESUMEN_TMEAN$TMEAN_SD, 
                 RESUMEN_PREC$PREC_MEAN, RESUMEN_PREC$PREC_SD, RESUMEN_EVAPO$EVAPO_MEAN,
                 RESUMEN_EVAPO$EVAPO_SD)

names(RESUMEN) <- c("Mes", "TMAX_MEAN", "TMAX_SD", "TMIN_MEAN", "TMIN_SD",
                    "TMEAN_MEAN", "TMEAN_SD",
                    "PREC_MEAN", "PREC_SD", "EVAPO_MEAN", "EVAPO_SD")



#---------------- Grafico datos promedio ----------------------------------------------


RESUMEN$Mes <- as.numeric(RESUMEN$Mes)

RESUMEN <- RESUMEN[order(RESUMEN$Mes),]
RESUMEN$Mes <- factor(RESUMEN$Mes)

RESUMEN$Mes <- revalue(RESUMEN$Mes, c("1"="Jan", "2"="Feb", "3"="Mar", "4"="Apr", "5"="May",
                                  "6"="Jun", "7"="Jul", "8"="Aug", "9"="Sep", "10"="Oct",
                                  "11"="Nov", "12"="Dec"))


coeff <- 0.17

p1 <- ggplot(data= RESUMEN) +
  geom_bar(aes(Mes, y= PREC_MEAN), width = 0.7, alpha= 0.6 , stat = "identity")   +
  geom_line(aes(x=Mes, y=TMAX_MEAN/coeff, group=1), stat="identity") +
  geom_line(aes(x=Mes, y=TMIN_MEAN/coeff, group=1), stat="identity") + 
  geom_point(aes(x=Mes, y=TMAX_MEAN/coeff)) +
  geom_point(aes(x=Mes, y=TMIN_MEAN/coeff)) +
  theme_set(theme_bw())


p1

p2 <- p1 +  scale_y_continuous(
        # Features of the first axis
    name = "Rainfall and evapotranspiration (mm)",
        # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Temperature (°C)") ) +
    scale_x_discrete(name = "Month")


p2 <- p2 +    annotate(
    "text", label = "Max.\n temperature",
    x = "Nov", y = 28.5/coeff, size = 2.5, colour = "black"
  ) +  annotate(
    "text", label = "Min.\n Temperature",
    x = "Nov", y = 16/coeff, size = 2.5, colour = "black"
  )


p2

clima_prom <-  data.frame(T.MAX= mean(RESUMEN$TMAX_MEAN), T.MIN= mean(RESUMEN$TMIN_MEAN), PREC= sum(RESUMEN$PREC_MEAN))

clima_prom$T.mean <- (RESUMEN_prom$T.MAX + RESUMEN_prom$T.MIN)/2

clima_prom

#ggsave("plots/RESUMEN_promedio_La_estanzuela.jpg",device = "jpg" ,units="in", width=7, height=4, dpi=250)
  



#--------------------- Con evapotranspiración 

RESUMEN_gathered <- RESUMEN %>%
  gather('PREC_MEAN' , 'EVAPO_MEAN', key= "PREC_EVAPO", value = "Lamina") %>% 
  dplyr::select(Mes,PREC_EVAPO,Lamina)





coeff <- 0.12


RESUMEN_gathered$PREC_EVAPO <- factor(RESUMEN_gathered$PREC_EVAPO, levels = c("PREC_MEAN", "EVAPO_MEAN"))

p1 <- ggplot() +
  geom_bar(
    data = RESUMEN_gathered,
    aes(Mes, y = Lamina, fill = PREC_EVAPO),
    width = 0.7, 
    alpha = 0.6, 
    stat = "identity", 
    position = position_dodge2(width = 0.9)
  )  +
geom_line(data=RESUMEN, aes(x=Mes, y=TMAX_MEAN/coeff, group=1), stat="identity") +
  geom_line(data=RESUMEN, aes(x=Mes, y=TMIN_MEAN/coeff, group=1), stat="identity") + 
  geom_point(data=RESUMEN, aes(x=Mes, y=TMAX_MEAN/coeff)) +
  geom_point(data=RESUMEN, aes(x=Mes, y=TMIN_MEAN/coeff)) +
  theme_set(theme_bw()) +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c(PREC_MEAN="#333BFF", EVAPO_MEAN="#CC6600"),
                    labels = c("Rainfall", "Evapotranspiration"),
                    name=NULL)


p1

p2 <- p1 +  scale_y_continuous(
  # Features of the first axis
  name = "Rainfall and \nevapotrasnpiration (mm)",
  # Add a second axis and specify its features
  sec.axis = sec_axis(~.*coeff, name="Temperature (°C)") ) +
  scale_x_discrete(name = "Month")


p3 <- p2 +    annotate(
  "text", label = "Max.\n temperature",
  x = "Nov", y = 28.5/coeff, size = 2.5, colour = "black"
) +  annotate(
  "text", label = "Min.\n Temperature",
  x = "Nov", y = 16/coeff, size = 2.5, colour = "black"
)


p3

clima_prom <-  data.frame(T.MAX= mean(RESUMEN$TMAX_MEAN), T.MIN= mean(RESUMEN$TMIN_MEAN), PREC= sum(RESUMEN$PREC_MEAN),
                          EVAPO = sum(RESUMEN$EVAPO_MEAN))


clima_prom

#ggsave("plots/RESUMEN_promedio_La_estanzuela.jpg",device = "tiff" ,units="in", width=7, height=4, dpi=250)


#----------------------------------------

#------------------------------------------------------------------------------------------------------

#---------------------- Grafico en el tiempo -----------------------------------------------------------



precip_quarters <- split(precip_MES1$Precipitacion, f="quarters")
precip_quarters_sum <- lapply(precip_quarters, FUN = cumsum)
precip_QUARTERS <- do.call(rbind, precip_quarters_sum)

#index de corte:
index_quarters <- endpoints(precip_QUARTERS, on = "quarters")
precip_QUARTERS <- precip_QUARTERS[index_quarters]
names(precip_QUARTERS) <- "Precipitacion"


dataframe_precipitacion <- as.data.frame(precip_QUARTERS)

dataframe_precipitacion <- data.frame(fecha= index(precip_QUARTERS), dataframe_precipitacion)


ggplot(dataframe_precipitacion, aes(x = fecha, y = Precipitacion)) +
  geom_line() +
  labs(x = "Fecha", y = "Precipitación") +
  theme_minimal()



#----------------------------------------------------------------

temps_monthly <- data.frame(fecha=index(temps_monthly), temps_monthly)

coeff <- 0.2

precipitacion <- ggplot() +
  geom_line(data=dataframe_precipitacion, aes(x = fecha, y = Precipitacion)) +
  labs(x = "Time (years)", y = "Rainfall (mm)") +
  theme_bw() + geom_hline(yintercept = mean(mean(dataframe_precipitacion$Precipitacion, na.rm = TRUE)), linetype = "dashed", color = "red")+
  scale_x_datetime(date_labels = "%Y", date_breaks = "10 years") + scale_y_continuous(position = "right")




n <- -2
p3 <- p2 +    annotate(
  "text", label = "Max.\n temperature",
  x = "Oct", y = 26.9/coeff, size = 2.5, colour = "black"
) +  annotate(
  "text", label = "Min.\n Temperature",
  x = "Oct", y = 15.8/coeff, size = 2.5, colour = "black"
) + theme(
  legend.spacing.x = unit(0.2, "cm"), 
  legend.spacing.y = unit(0.2, "cm")
) + theme(
  legend.margin = margin(n, n, n, n),
  legend.text = element_text(size = 8)
)+
  theme(legend.key.size = unit(0.3, 'cm'))

p3

plot_grid(p3, precipitacion, ncol = 2, align = "hv", labels = c("(a)", "(b)"))

ggsave("plots/grid_temperatura_precipitacion_horizontal.png",device = "png" ,units="in", width=9, height=3.2, dpi=250)


