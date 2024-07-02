
##facet_wrap dengue incidence rate
data2=data %>% 
  group_by(Año, mes_final, cod_mun, Municipio) %>%
  # calculate municipal level incidence rate
  summarise(cases = sum(casos),
            pop = sum(Total_General)) %>% 
  mutate(var = cases / pop * 10^5) 


data2=data.frame(data2)
ggplot(data2,aes(x = mes_final, y = Año, fill = var)) + 
  geom_raster() + facet_wrap(~Municipio) + 
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Dengue incidence rate", colours=rev(brewer.pal(9, "RdYlGn")), trans = "log1p", breaks = c(0, 10, 50, 100), labels = c(0, 10, 50, 100) ) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct"))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

#Tprom
ggplot(data,aes(x = mes_final, y = Año, fill = tprom2)) + 
  geom_raster() + facet_wrap(~Municipio) + 
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Mean temperature", colours=rev(brewer.pal(9, "RdYlGn")), trans = "log1p", breaks = c(0, 20, 25, 28), labels = c(0, 20, 25, 28) ) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct"))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

#Tmáx
ggplot(data,aes(x = mes_final, y = Año, fill = tmax2)) + 
  geom_raster() + facet_wrap(~Municipio) + 
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Maximun temperature", colours=rev(brewer.pal(9, "RdYlGn")), trans = "log1p", breaks = c(24, 28, 32, 34), labels = c(24, 28, 32, 34) ) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct"))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))


#Prec
ggplot(data,aes(x = mes_final, y = Año, fill = prec2)) + 
  geom_raster() + facet_wrap(~Municipio) + 
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Precipitation (mm)", colours=rev(brewer.pal(9, "RdYlGn")), trans = "log1p", breaks = c(0, 100, 500, 1000), labels = c(0, 100, 500, 1000) ) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct"))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

###############################################
##Relative risk mean temperature overall lags##

# extract full coef and vcov and create indicators for each term
coef <- model0.29$summary.fixed$mean
vcov <- model0.29$misc$lincomb.derived.covariance.matrix

# create indicators for the terms associated with Tprom crossbasis
indt <- grep("basis_tprom", model0.29$names.fixed)

# extract predictions from the Tprom DLNM centred on overall mean Tprom (22.6 deg C)
predt <- crosspred(basis_tprom, coef = coef[indt], vcov=vcov[indt,indt],
                   model.link = "log", bylag = 0.25, cen = round(mean(data$tprom2), 0))

# define x values (lag, by lag)
lagbylag <- seq(0, nlag, 0.25)

# plot overall response across all lags 
plot(predt,"overall", xlab = expression(paste("Mean Temperature °C")), 
     ylab = "Relative risk", main = "", 
     ylim = c(range(predt$allRRlow,predt$allRRhigh)))
box()

##########################
# Contour plot

model <- model0.29

# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix

# find position of the terms associated with Tprom crossbasis
indt <- grep("basis_tprom", model$names.fixed)

# extract predictions from the Tprom DLNM centred on overall mean Tprom (22,66)
predt <- crosspred(basis_tprom, coef = coef[indt], vcov=vcov[indt,indt],
                   model.link = "log", bylag = 0.25, cen = round(mean(data$tprom2), 0)) 

y <- predt$predvar
x <- seq(0, nlag, 0.25)
z <- t(predt$matRRfit)

pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(z, 20)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))

filled.contour(x,y,z,
               xlab = "Lag (months)", ylab = expression(paste("Mean Temperature (°C)")),
               col = cols,levels = levels,
               plot.axes = { axis(1, at = 0:nlag, c(0:nlag)) 
                 axis(2)},
               key.title = title(main = "RR"))


# Plot precipitation output 
model <- model0.29

# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix

# find position of the terms associated with prec crossbasis
indt1 <- grep("basis_prec", model$names.fixed)

# extract predictions from the Tmin DLNM centred on overall mean prec (152.02 mm)
predt1 <- crosspred(basis_prec, coef = coef[indt1], vcov=vcov[indt1,indt1],
                    model.link = "log", bylag = 0.25, cen = round(mean(data$prec2), 0)) 

y1 <- predt1$predvar
x1 <- seq(0, nlag, 0.25)
z1 <- t(predt1$matRRfit)

pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(z1, 20)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))

filled.contour(x1,y1,z1,
               xlab = "Lag", ylab = expression(paste("Precipitation (mm)")),
               col = cols,levels = levels,
               plot.axes = { axis(1, at = 0:nlag, c(0:nlag)) 
                 axis(2)},
               key.title = title(main = "RR"))


# relative humidity output 
model <- model0.29

# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix

# find position of the terms associated with hum crossbasis
indt2 <- grep("basis_hmd", model$names.fixed)

# extract predictions from the Tmin DLNM centred on mean hmd
predt2 <- crosspred(basis_hmd, coef = coef[indt2], vcov=vcov[indt2,indt2],
                    model.link = "log", bylag = 0.25, cen = round(mean(data$hum2), 0)) 

y2 <- predt2$predvar
x2 <- seq(0, nlag, 0.25)
z2 <- t(predt2$matRRfit)

pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(z1, 20)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))

filled.contour(x2,y2,z2,
               xlab = "Lag", ylab = expression(paste("Relative humidity")),
               col = cols,levels = levels,
               plot.axes = { axis(1, at = 0:nlag, c(0:nlag)) 
                 axis(2)},
               key.title = title(main = "RR"))
