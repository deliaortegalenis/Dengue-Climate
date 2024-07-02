library(INLA)
library(data.table)
library(tidyverse)
library(sf)
library(sp)
library(spdep)
library(dlnm)
library(tsModel)
library(RColorBrewer)
library(raster)
library(rgdal)
library(maps)
library(tidyr)


#load shape file for Valle
suramerica=shapefile("C:/Users/david.arango/OneDrive - PUJ Cali/Desktop/DELIA/scripts_final_paper/SUDAMERICA_ADM2/sudamerica_adm2.shp")
colombia=suramerica[which(suramerica$ADM0=="COLOMBIA"),]
valle=colombia[colombia$ADM1=="Valle del Cauca",]
valle=valle[valle$AREA>0.003,]


# Creating adjacency matrix
library(spdep)
nb.map <- poly2nb(valle) 
nb2INLA("valle.adj", nb.map)
g <- inla.read.graph(filename = "valle.adj")

#load data 
library(readxl)
data <- read_excel("C:/Users/david.arango/OneDrive - PUJ Cali/Desktop/DELIA/scripts_final_paper/base_final_3_vr3")
data <- data[order(data$Año),]
data$tmax2 <- 1.65+(1.16*data$tmax2)
data$tmin2 <- (-13.43)+(1.81*data$tmin2)
data$tprom2 <- (-6.33)+(1.44*data$tprom2)

#creating time variable
library(dplyr)
data <- data %>%
  group_by(cod_mun) %>%
  mutate(time = row_number())

# maximun lag for climate variables
nlag = 6

# Minimum temperature 
lag_tmin <- tsModel::Lag(data$tmin2, group = data$cod_mun, k = 0:nlag)
# Maximum temperature 
lag_tmax <- tsModel::Lag(data$tmax2, group = data$cod_mun, k = 0:nlag)
# Precipitation
lag_prec <- tsModel::Lag(data$prec2, group = data$cod_mun, k = 0:nlag)
# Humidity
lag_hmd <- tsModel::Lag(data$hum2, group = data$cod_mun, k = 0:nlag)
# Maximum temperature Range
lag_tmaxR <- tsModel::Lag(data$tmaxR, group = data$cod_mun, k = 0:nlag)
# Mean temperature 
lag_tprom <- tsModel::Lag(data$tprom2, group = data$cod_mun, k = 0:nlag)
# Mean temperature Range
lag_tpromR <- tsModel::Lag(data$tpromR, group = data$cod_mun, k = 0:nlag)

# nino12
lag_nino12 <- tsModel::Lag(data$nino12, group = data$cod_mun, k = 0:nlag)

# Remove year 2000 from lagged climate variables
lag_tmin <- lag_tmin[data$Año > 2000,]
lag_tmax <- lag_tmax[data$Año > 2000,]
lag_prec <- lag_prec[data$Año > 2000,]
lag_hmd  <- lag_hmd[data$Año > 2000,]
lag_tmaxR <- lag_tmaxR[data$Año > 2000,]
lag_tprom <- lag_tprom[data$Año > 2000,]
lag_tpromR <- lag_tpromR[data$Año > 2000,]
lag_nino12 <- lag_nino12[data$Año > 2000,]


# remove year 2000 from dengue dataframe
data <- data[data$Año > 2000,]

# define dimensions
data$time <- data$time - 12
# total number of months
ntime <- length(unique(data$time))
# total number of years
nyear <- length(unique(data$Año))
# total number of municipio
nmuni <- length(unique(data$cod_mun))

# cross-basis matrix (combining nonlinear exposure and lag functions)
lagknot = equalknots(0:nlag, 2)

# Tmin
basis_tmin <- crossbasis(lag_tmin,
                         argvar = list(fun = "ns", knots = equalknots(data$tmin2, 2)),
                         arglag = list(fun = "ns", knots = nlag/2))

# Tmax
basis_tmax <- crossbasis(lag_tmax,
                         argvar = list(fun = "ns", knots = equalknots(data$tmax2, 2)),
                         arglag = list(fun = "ns", knots = nlag/2))

# Prec
basis_prec <- crossbasis(lag_prec,
                         argvar = list(fun = "ns", knots = equalknots(data$prec2, 2)),
                         arglag = list(fun = "ns", knots = nlag/2))

# Hmd
basis_hmd <- crossbasis(lag_hmd,
                        argvar = list(fun = "ns", knots = equalknots(data$hum2, 2)),
                        arglag = list(fun = "ns", knots = nlag/2))


# TmaxR
basis_tmaxR <- crossbasis(lag_tmaxR,
                          argvar = list(fun = "ns", knots = equalknots(data$tmaxR, 2)),
                          arglag = list(fun = "ns", knots = nlag/2))

# Tprom
basis_tprom <- crossbasis(lag_tprom,
                          argvar = list(fun = "ns", knots = equalknots(data$tprom2, 2)),
                          arglag = list(fun = "ns", knots = nlag/2))

# TpromR
basis_tpromR <- crossbasis(lag_tpromR,
                           argvar = list(fun = "ns", knots = equalknots(data$tpromR, 2)),
                           arglag = list(fun = "ns", knots = nlag/2))

# nino12
basis_nino12 <- crossbasis(lag_nino12,
                           argvar = list(fun = "ns", knots = equalknots(data$nino12, 2)),
                           arglag = list(fun = "ns", knots = nlag/2))

# cross-basis matrix (evaluating linear exposure and lag functions)
lagknot = equalknots(0:nlag, 2)

# Tmin
basis_tmin <- crossbasis(lag_tmin,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "ns", knots = nlag/2))

# Tmax
basis_tmax <- crossbasis(lag_tmax,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "ns", knots = nlag/2))

# Prec
basis_prec <- crossbasis(lag_prec,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "ns", knots = nlag/2))

# Hmd
basis_hmd <- crossbasis(lag_hmd,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns", knots = nlag/2))


# TmaxR
basis_tmaxR <- crossbasis(lag_tmaxR,
                          argvar = list(fun = "lin"),
                          arglag = list(fun = "ns", knots = nlag/2))

# Tprom
basis_tprom <- crossbasis(lag_tprom,
                          argvar = list(fun = "lin"),
                          arglag = list(fun = "ns", knots = nlag/2))

# TpromR
basis_tpromR <- crossbasis(lag_tpromR,
                           argvar = list(fun = "lin"),
                           arglag = list(fun = "ns", knots = nlag/2))


colnames(basis_tmin) = paste0("basis_tmin.", colnames(basis_tmin))
colnames(basis_tmax) = paste0("basis_tmax.", colnames(basis_tmax))
colnames(basis_prec) = paste0("basis_prec.", colnames(basis_prec))
colnames(basis_hmd)  = paste0("basis_hmd.", colnames(basis_hmd))
colnames(basis_tmaxR)  = paste0("basis_tmaxR.", colnames(basis_tmaxR))
colnames(basis_tprom)  = paste0("basis_tprom.", colnames(basis_tprom))
colnames(basis_tpromR)  = paste0("basis_tpromR.", colnames(basis_tpromR))
colnames(basis_nino12)  = paste0("basis_nino12.", colnames(basis_nino12))


# create municipio index 
data$muni_index <- rep(1:nmuni, ntime)

# create year index
data$year_index <- data$Año - 2000 

# set data for models
Y  <- data$casos # response variable
N  <- length(Y) # total number of data points
E  <- data$Total_General/10^5 # model offset so that response is equivalent to an incidence rate per 100,000 people
T1 <- data$mes_final # for random effect to account for annual cycle (seasonality)
T2 <- data$year_index # for random effect to account for inter-annual variability
S1 <- data$muni_index # for municipal spatial random effect
Vu <- data$pob_urbana # include level of urbanisation (% pop living in urban areas)
Vw <- data$acueducto # include % of households with water system
Vx <- data$alcan # include % of households with sewer
Vs <- data$estrato_1 # include % of households with low socioeconomic strata 

# create dataframe for model testing
df <- data.frame(Y, E, T1, T2, S1, Vu, Vw, Vx, Vs)

# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# inla model function

# include formula and set defaults for data, family (to allow other prob dist models e.g. Poisson) and config (to allow for sampling)
mymodel <- function(formula, data = df, family = "nbinomial", config = FALSE)
  
{
  model <- inla(formula = formula, data = data, family = family, offset = log(E),
                control.inla = list(strategy = 'adaptive'), 
                control.compute = list(dic = TRUE, waic=TRUE, config = config, 
                                       cpo = TRUE, return.marginals = FALSE),
                control.fixed = list(correlation.matrix = TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                verbose = FALSE)
  model <- inla.rerun(model)
  return(model)
}


# baseline model
baseformula <- Y ~ 1 + f(T1,  replicate=S1, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) +
  f(S1, model = "bym2", replicate = T2, graph = g, 
    scale.model = TRUE, hyper = precision.prior) 

model <- mymodel(baseformula, family = "poisson")
model$dic$dic
model$waic$waic

# define formulas by updating the baseline formula with different combinations of Tmin, Tmax, Prec and Hmd cross-basis functions
formula0.1 <- update.formula(baseformula, ~. + basis_tmin)
formula0.2 <- update.formula(baseformula, ~. + basis_tmax)
formula0.3 <- update.formula(baseformula, ~. + basis_prec)
formula0.4 <- update.formula(baseformula, ~. + basis_hmd)
formula0.5 <- update.formula(baseformula, ~. + basis_tmaxR)
formula0.6 <- update.formula(baseformula, ~. + basis_tprom) 
formula0.7 <- update.formula(baseformula, ~. + basis_tpromR) 
formula0.8 <- update.formula(baseformula, ~. + basis_nino12) 
formula0.9 <- update.formula(baseformula, ~. + basis_tmax + basis_prec)
formula0.10 <- update.formula(baseformula, ~. + basis_tmin + basis_prec)
formula0.11 <- update.formula(baseformula, ~. + basis_tprom + basis_prec)
formula0.12 <- update.formula(baseformula, ~. + basis_tmaxR + basis_prec)
formula0.13 <- update.formula(baseformula, ~. + basis_tpromR + basis_prec)
formula0.14 <- update.formula(baseformula, ~. + basis_tmax + basis_hmd)
formula0.15 <- update.formula(baseformula, ~. + basis_tmin + basis_hmd)
formula0.16 <- update.formula(baseformula, ~. + basis_tprom + basis_hmd)
formula0.17 <- update.formula(baseformula, ~. + basis_tmaxR + basis_hmd)
formula0.18 <- update.formula(baseformula, ~. + basis_tpromR + basis_hmd)
formula0.19 <- update.formula(baseformula, ~. + basis_hmd + basis_prec)
formula0.20 <- update.formula(baseformula, ~. + basis_tmax + basis_nino12)
formula0.21 <- update.formula(baseformula, ~. + basis_tmin + basis_nino12)
formula0.22 <- update.formula(baseformula, ~. + basis_tprom + basis_nino12)
formula0.23 <- update.formula(baseformula, ~. + basis_tmaxR + basis_nino12)
formula0.24 <- update.formula(baseformula, ~. + basis_tpromR + basis_nino12)
formula0.25 <- update.formula(baseformula, ~. + basis_hmd + basis_nino12)
formula0.26 <- update.formula(baseformula, ~. + basis_prec + basis_nino12)
formula0.27  <- update.formula(baseformula, ~. + basis_tmax + basis_hmd + basis_prec)
formula0.28  <- update.formula(baseformula, ~. + basis_tmin + basis_hmd + basis_prec)
formula0.29  <- update.formula(baseformula, ~. + basis_tprom + basis_hmd + basis_prec)
formula0.30  <- update.formula(baseformula, ~. + basis_tmaxR + basis_hmd + basis_prec)
formula0.31  <- update.formula(baseformula, ~. + basis_tpromR + basis_hmd + basis_prec)
formula0.32  <- update.formula(baseformula, ~. + basis_tmax + basis_prec + basis_nino12)
formula0.33  <- update.formula(baseformula, ~. + basis_tmin + basis_prec + basis_nino12)
formula0.34  <- update.formula(baseformula, ~. + basis_tprom + basis_prec + basis_nino12)
formula0.35  <- update.formula(baseformula, ~. + basis_tmaxR + basis_prec + basis_nino12)
formula0.36  <- update.formula(baseformula, ~. + basis_tpromR + basis_prec + basis_nino12)


model0.1 <- mymodel(formula0.1, df)
model0.2 <- mymodel(formula0.2, df)
model0.3 <- mymodel(formula0.3, df)
model0.4 <- mymodel(formula0.4, df)
model0.5 <- mymodel(formula0.5, df)
model0.6 <- mymodel(formula0.6, df)
model0.7 <- mymodel(formula0.7, df)
model0.8 <- mymodel(formula0.8, df)
model0.9 <- mymodel(formula0.9, df)
model0.10 <- mymodel(formula0.10, df)
model0.11 <- mymodel(formula0.11, df)
model0.12 <- mymodel(formula0.12, df)
model0.13 <- mymodel(formula0.13, df)
model0.14 <- mymodel(formula0.14, df)
model0.15 <- mymodel(formula0.15, df)
model0.16 <- mymodel(formula0.16, df)
model0.17 <- mymodel(formula0.17, df)
model0.18 <- mymodel(formula0.18, df)
model0.19 <- mymodel(formula0.19, df)
model0.20 <- mymodel(formula0.20, df)
model0.21 <- mymodel(formula0.21, df)
model0.22 <- mymodel(formula0.22, df)
model0.23 <- mymodel(formula0.23, df)
model0.24 <- mymodel(formula0.24, df)
model0.25 <- mymodel(formula0.25, df)
model0.26 <- mymodel(formula0.26, df)
model0.27 <- mymodel(formula0.27, df)
model0.28 <- mymodel(formula0.28, df)
model0.29 <- mymodel(formula0.29, df)
model0.30 <- mymodel(formula0.30, df)
model0.31 <- mymodel(formula0.31, df)
model0.32 <- mymodel(formula0.32, df)
model0.33 <- mymodel(formula0.33, df)
model0.34 <- mymodel(formula0.34, df)
model0.35 <- mymodel(formula0.35, df)
model0.36 <- mymodel(formula0.36, df)

model0.1$dic$dic
model0.2$dic$dic
model0.3$dic$dic
model0.4$dic$dic
model0.5$dic$dic
model0.6$dic$dic
model0.7$dic$dic
model0.8$dic$dic
model0.9$dic$dic
model0.10$dic$dic
model0.11$dic$dic
model0.12$dic$dic
model0.13$dic$dic
model0.14$dic$dic
model0.15$dic$dic
model0.16$dic$dic
model0.17$dic$dic
model0.18$dic$dic
model0.19$dic$dic
model0.20$dic$dic
model0.21$dic$dic
model0.22$dic$dic
model0.23$dic$dic
model0.24$dic$dic
model0.25$dic$dic
model0.26$dic$dic
model0.27$dic$dic
model0.28$dic$dic
model0.29$dic$dic
model0.30$dic$dic
model0.31$dic$dic
model0.32$dic$dic
model0.33$dic$dic
model0.34$dic$dic
model0.35$dic$dic
model0.36$dic$dic

model$waic$waic
model0.1$waic$waic
model0.2$waic$waic
model0.3$waic$waic
model0.4$waic$waic
model0.5$waic$waic
model0.6$waic$waic
model0.7$waic$waic
model0.8$waic$waic
model0.9$waic$waic
model0.10$waic$waic
model0.11$waic$waic
model0.12$waic$waic
model0.13$waic$waic
model0.14$waic$waic
model0.15$waic$waic
model0.16$waic$waic
model0.17$waic$waic
model0.18$waic$waic
model0.19$waic$waic
model0.20$waic$waic
model0.21$waic$waic
model0.22$waic$waic
model0.23$waic$waic
model0.24$waic$waic
model0.25$waic$waic
model0.26$waic$waic
model0.27$waic$waic
model0.28$waic$waic
model0.29$waic$waic
model0.30$waic$waic
model0.31$waic$waic
model0.32$waic$waic
model0.33$waic$waic
model0.34$waic$waic
model0.35$waic$waic
model0.36$waic$waic


# define formula by updating the baseline formula with different combinations of socioeconomic variables
formula0.37 <- update.formula(baseformula, ~. + Vu)
formula0.38 <- update.formula(baseformula, ~. + Vw)
formula0.39 <- update.formula(baseformula, ~. + Vx)
formula0.40 <- update.formula(baseformula, ~. + Vs)
formula0.41 <- update.formula(baseformula, ~. + Vu + Vw)
formula0.42 <- update.formula(baseformula, ~. + Vu + Vx)
formula0.43 <- update.formula(baseformula, ~. + Vu + Vs)
formula0.44 <- update.formula(baseformula, ~. + Vw + Vs)
formula0.45 <- update.formula(baseformula, ~. + Vx + Vs)
formula0.46 <- update.formula(baseformula, ~. + Vw + Vx)
formula0.47 <- update.formula(baseformula, ~. + Vu + Vx + Vs)
formula0.48 <- update.formula(baseformula, ~. + Vs + Vu + Vw)


model0.37 <- mymodel(formula0.37, df)
model0.38 <- mymodel(formula0.38, df)
model0.39 <- mymodel(formula0.39, df)
model0.40 <- mymodel(formula0.40, df)
model0.41 <- mymodel(formula0.41, df)
model0.42 <- mymodel(formula0.42, df)
model0.43 <- mymodel(formula0.43, df)
model0.44 <- mymodel(formula0.44, df)
model0.45 <- mymodel(formula0.45, df)
model0.46 <- mymodel(formula0.46, df)
model0.47 <- mymodel(formula0.47, df)
model0.48 <- mymodel(formula0.48, df)



model0.37$dic$dic
model0.38$dic$dic
model0.39$dic$dic
model0.40$dic$dic
model0.41$dic$dic
model0.42$dic$dic
model0.43$dic$dic
model0.44$dic$dic
model0.45$dic$dic
model0.46$dic$dic
model0.47$dic$dic
model0.48$dic$dic


model0.37$waic$waic
model0.38$waic$waic
model0.39$waic$waic
model0.40$waic$waic
model0.41$waic$waic
model0.42$waic$waic
model0.43$waic$waic
model0.44$waic$waic
model0.45$waic$waic
model0.46$waic$waic
model0.47$waic$waic
model0.48$waic$waic

# baseformula with climatic variables
baseformula1 <- Y ~ 1 + f(T1, replicate=S1, model = "rw1", cyclic = TRUE, constr = TRUE,
                          scale.model = TRUE,  hyper = precision.prior) +
  f(S1, model = "bym2", replicate = T2, graph = g, 
    scale.model = TRUE, hyper = precision.prior) + basis_tprom + basis_hmd + basis_prec

# best fitting model0.29 + water network
formula1.2 <- update.formula(baseformula1, ~. + Vw)

model1.1 <- model0.29
model1.2 <- mymodel(formula1.2, df)

model1.1$dic$dic
model1.1$waic$waic


