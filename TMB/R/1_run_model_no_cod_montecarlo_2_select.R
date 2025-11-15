rm(list = ls())

# -- load packages --
library(TMB)
library(tidyverse)
theme_set(theme_bw() +
            theme(
              strip.background = element_rect(fill = "transparent", color = "transparent"),
              strip.text = element_text(size = 11)
            ))
source("simulate_consumption.R")

# -- load data --
load("data.RData")

# compile & load the updated C++ (with fixed kappa and no logeps_h)
compile("Consumption_selective.cpp")
dyn.load(dynlib("Consumption_selective"))

# prep
data$area <- data$paperarea

# ---- parameters MUST match the C++: x, logb_alpha, logb_h, BH, logeps_alpha, logSigma
parameters <- list(
  x = 0,
  logb_alpha = 0,
  logb_h = 0,
  BH = 0,
  logeps_alpha = 0,     # <--- added (still present in C++)
  logSigma = log(20)
)

par.list <- list(
  age3 = parameters,
  age4 = parameters,
  age5 = parameters,
  age6 = parameters
)

# map settings (unchanged): map BH or x to NA if you want them fixed
map.list <- list(
  age3 = list(x = factor(NA), BH = factor(NA)),
  age4 = list(x = factor(NA), BH = factor(NA)),  # , logb_h = factor(NA)),
  age5 = list(x = factor(NA), BH = factor(NA)),  # , logb_h = factor(NA)),
  age6 = list(x = factor(NA), BH = factor(NA))
)

# per-age switch for using logb_h (1) or BH (0)
selective <- c(1, 1, 1, 1)

test  <- dplyr::filter(data, year > 2013)
data2 <- dplyr::filter(data, year <= 2013)

shift <- 0.1

simmall <- parest.all <- results.all <- pred.all <- MCall <- tibble()
area <- data2$area[1]
ages <- 3:6
nsim <- 2000
set.seed(12345)

for (i in seq_along(ages)) {
  age <- ages[i]
  
  maxTemp <- max(data2$meantemp[data2$age == age])
  minTemp <- min(data2$meantemp[data2$age == age])
  
  temp <- data2$meantemp[data2$age == age & data2$area == area]
  
  data.list <- list(
    ny   = length(unique(data2$year)),
    T_0  = minTemp - shift,
    T_l  = maxTemp + shift,
    P    = data2$cod.maturingbiomass[data2$age == age & data2$area == area]/1e3,
    B    = data2$cap.maturingbiomass[data2$age == age & data2$area == area]/1e3,
    Temp = temp,
    f    = data2$EmpC[data2$age == age & data2$area == area],
    
    predP    = test$cod.maturingbiomass[test$age == age & test$area == area]/1e3,
    predB    = test$cap.maturingbiomass[test$age == age & test$area == area]/1e3,
    predTemp = test$meantemp[test$age == age & test$area == area],
    
    logbhTRUE = selective[i]
  )
  
  obj <- MakeADFun(
    data = data.list,
    parameters = par.list[[i]],
    DLL = "Consumption_selective",
    map = map.list[[i]],
    silent = TRUE
  )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  rep <- sdreport(obj)
  ss  <- summary(rep, p.value = TRUE)
  
  # ---- Pull scalar parameter estimates by NAME (safer than slicing)
  # Expected ADREPORT scalars in C++: x, b_alpha, b_h, eps_alpha, sigma
  wanted <- c("x", "b_alpha", "b_h", "eps_alpha", "sigma")
  have <- intersect(wanted, rownames(ss))
  scal <- as_tibble(ss[have, , drop = FALSE], rownames = "parname")
  
  # Build parest tibble in your preferred naming
  # Keep the four you later use downstream: x, balpha, bh, sigma
  # (eps_alpha is estimated but not used in simulateConsumption par vector)
  pick <- c("x", "b_alpha", "b_h", "sigma")
  scal4 <- scal[match(pick, scal$parname), ]
  names(scal4)[2:4] <- c("Estimate", "Std. Error", "z value")  # keep sdreport naming
  
  parest <- tibble(
    Estimate   = scal4$Estimate,
    `Std. Error` = scal4$`Std. Error`,
    `z value`  = scal4$`z value`,
    `P-value`  = scal4$`Pr(>|z|)`,
    par = c("x", "balpha", "bh", "sigma"),
    age = age
  )
  
  parest.all <- bind_rows(parest.all, parest)
  
  # ---- Pull vectors by name (alpha, h, C) for the fitting years
  alpha_hat <- ss[rownames(ss) == "alpha", 1]
  alpha_se  <- ss[rownames(ss) == "alpha", 2]
  h_hat     <- ss[rownames(ss) == "h", 1]
  h_se      <- ss[rownames(ss) == "h", 2]
  C_hat     <- ss[rownames(ss) == "C", 1]
  C_se      <- ss[rownames(ss) == "C", 2]
  
  results <- tibble(
    conv    = opt$convergence,
    B       = data.list$B,
    P       = data.list$P,
    year    = unique(data2$year),
    age     = age,
    area    = area,
    temp    = data.list$Temp,
    minTemp = data.list$T_0,
    maxTemp = data.list$T_l,
    f       = data.list$f,
    
    C       = C_hat,
    seC     = C_se,
    lowC    = C_hat - 2 * C_se,
    highC   = C_hat + 2 * C_se,
    
    alpha   = alpha_hat,
    sealpha = alpha_se,
    lowAlpha = alpha_hat - 2 * alpha_se,
    highAlpha = alpha_hat + 2 * alpha_se,
    
    h       = h_hat,
    seh     = h_se,
    lowH    = h_hat - 2 * h_se,
    highH   = h_hat + 2 * h_se
  )
  
  # ---- Monte Carlo sims (uses b_a, b_h, sigma)
  sims <- sapply(
    X = 1:nsim, FUN = simulateConsumption,
    temp = data.list$predTemp,
    par  = unlist(c(parest$Estimate[parest$par == "balpha"],
                    parest$Estimate[parest$par == "bh"],
                    parest$Estimate[parest$par == "sigma"])),
    B = data.list$predB,
    tmin = minTemp,
    tmax = maxTemp,
    simtemp = TRUE
  )
  
  simsall <- sapply(
    X = 1:nsim, FUN = simulateConsumption,
    temp = data$meantemp[data$age == age & data$area == area],
    par  = unlist(c(parest$Estimate[parest$par == "balpha"],
                    parest$Estimate[parest$par == "bh"],
                    parest$Estimate[parest$par == "sigma"])),
    B = data$cap.maturingbiomass[data$age == age & data$area == area]/1e3,
    tmin = minTemp,
    tmax = maxTemp,
    simtemp = FALSE
  )
  
  simall <- data.frame(
    year = 2001:2017,
    age = age,
    simsall
  )
  names(simall) <- c("year", "age", 1:nsim)
  simall <- pivot_longer(simall, cols = 3:ncol(simall), names_to = "sim", values_to = "consum")
  simmall <- bind_rows(simmall, simall)
  
  MCall <- bind_rows(MCall, tibble(
    year = unique(data$year),
    age = age,
    area = area,
    f   = data$EmpC[data$age == age & data$area == area],
    C   = rowMeans(simsall),
    seC = apply(simsall, 1, sd),
    lowC  = apply(simsall, 1, quantile, prob = .025),
    highC = apply(simsall, 1, quantile, prob = .975)
  ))
  
  pred <- tibble(
    year = unique(test$year),
    age = age,
    area = area,
    f   = test$EmpC[test$age == age & test$area == area],
    C   = rowMeans(sims),
    seC = apply(sims, 1, sd),
    lowC  = apply(sims, 1, quantile, prob = .025),
    highC = apply(sims, 1, quantile, prob = .975)
  )
  
  results.all <- bind_rows(results.all, results)
  pred.all    <- bind_rows(pred.all, pred)
}

cat("---------\nDid all converge?\n", ifelse(all(results.all$conv == 0), "- YES", "- NO"), "\n---------\n")

(parest.all <- parest.all[, c("age", "par", "Estimate", "Std. Error", "z value", "P-value")])

write.table(dplyr::filter(parest.all, par != "x"),
            "../Mahmood_model_parameteres.txt", quote = FALSE, row.names = FALSE)

library(xtable)
names(parest.all)[1] <- "Age"
names(parest.all)[2] <- "Parameter"
names(parest.all)[4] <- "Std.Err"
names(parest.all)[5] <- "Z-value"

parest.all$`p-value` <- ifelse(parest.all$`P-value` < 0.05,
                               paste0("\\color{green} ", round(parest.all$`P-value`, 3)),
                               paste0("\\color{red} ", round(parest.all$`P-value`, 3)))
parest.all <- subset(parest.all, Parameter != "x")
parest.all$stars <- gtools::stars.pval(parest.all$`P-value`)
names(parest.all)[ncol(parest.all)] <- ""
names(parest.all) <- c("Age", "Parameter", "Estimate", "Std.Err", "Z-value", "P-value", "P", "")

parest.all$Parameter[parest.all$Parameter == "balpha"] <- "$b_a$"
parest.all$Parameter[parest.all$Parameter == "bh"]     <- "$b_h$"
parest.all$Parameter[parest.all$Parameter == "sigma"]  <- "$\\sigma$"

parest.all$Estimate <- as.character(signif(as.numeric(parest.all$Estimate), 3))
parest.all$Std.Err  <- as.character(signif(as.numeric(parest.all$Std.Err), 3))
parest.all$`Z-value` <- as.character(signif(as.numeric(parest.all$`Z-value`), 3))
parest.all$`P-value` <- as.character(round(as.numeric(parest.all$`P-value`), 2))
parest.all$`P-value`[parest.all$`P-value` == "0"] <- "0.00"

parest.all$Age[seq(1, nrow(parest.all) - 2, 3)] <- paste0("\\multirow{3}{*}{", parest.all$Age[seq(1, nrow(parest.all) - 2, 3)], "}")
parest.all$Age[-seq(1, nrow(parest.all) - 2, 3)] <- ""

print(xtable(parest.all[, -which(names(parest.all) == "P")], digits = 5,
             label = "tab:parameter_estimates_by_age",
             align = c("r", "c", rep("l", ncol(parest.all) - 2)),
             caption = "Model parameter estimates based on observations for years 2001-2013. Total consumption and mean temperature of the area I, II and III is used."),
      include.rownames = FALSE,
      file = "../Tables/Estimates_by_age_2013.tex",
      hline.after = c(-1, 0, 3, 6, 9, 12),
      sanitize.text.function = function(x) x)

pred.all <- bind_rows(
  pred.all,
  results.all[results.all$year == 2013, c("year", "age", "area", "f", "C", "seC", "lowC", "highC")]
)

results.all$agelab <- paste0("Cod age: ", results.all$age, " yrs")
pred.all$agelab    <- paste0("Cod age: ", pred.all$age, " yrs")

ggplot(results.all, aes(x = year)) +
  geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = .5, fill = "grey") +
  facet_wrap(~agelab, strip.position = "top", scales = "free_y") +
  geom_line(data = pred.all, aes(x = year, y = C), col = "red") +
  geom_ribbon(data = pred.all, aes(x = year, ymin = lowC, ymax = highC), fill = "red", alpha = .2) +
  scale_x_continuous(breaks = seq(2001, 2017, 2), name = "Year") +
  geom_line(aes(y = C), lty = 2) +
  geom_point(data = pred.all, aes(x = year, y = f), col = "black") +
  geom_point(aes(y = f), col = "red") +
  scale_y_continuous(name = "Consumption per cod (kg/month)")

ggsave("../Figures/montecarlo_simulated_future_eps0.pdf", device = "pdf", dpi = "retina", width = 10, height = 6)


# after fitting
nll <- opt$objective
k   <- length(obj$par)  # number of free parameters
n   <- length(data.list$f)

AIC <- 2 * nll + 2 * k
BIC <- 2 * nll + log(n) * k

obs  <- data.list$f
pred <- ss[rownames(ss) == "C", 1]

RMSE <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
MAE  <- mean(abs(obs - pred), na.rm = TRUE)

lowC <- results$lowC
highC <- results$highC
Coverage <- mean(obs >= lowC & obs <= highC, na.rm = TRUE)




MCall$agelab = paste0("Cod age: ", MCall$age, " yrs" )
MCall$color <- ifelse(MCall$year <=2013, "Train", "Test")
MCall[nrow(MCall)+1:4, ] <- MCall[MCall$year == 2013,]
MCall[nrow(MCall)-0:3, "color"] <- "Test"
ggplot(MCall, aes(x = year))+
  geom_ribbon(aes(ymin = lowC, ymax = highC, fill =color), alpha = .3) +
#  theme_bw() + 
  facet_wrap(~agelab) + 
  #geom_line(data = pred.all, aes(x = year, y = C), col = "red") + 
  #geom_ribbon(data = pred.all, aes(x= year, ymin = lowC, ymax = highC), fill = "red", alpha = .2)+
  scale_x_continuous(breaks = seq(2001,2017,2), name = "Year") +
  geom_line(aes(y=C), lty = 2) +
  geom_point(aes(y = f), col = "red")+
  #geom_point(data = pred.all, aes(x = year, y = f), col = "black")+ 
  scale_y_continuous(name = "Consumption per cod (tonnes/month)")+
  theme(legend.title= element_blank(),
        legend.position = "none") +
#  ggsave("../Consumption_report_no_cod/figs/1_consumption_per_year_and_age_2013_montecarlosim_allyears.pdf", width = 10, height = 6)
  ggsave("../Figures/montecarlo_simulated_CI_all_years.pdf", device = "pdf", dpi = "retina", width = 10, height = 6)

# -- Prediction --
test



ggplot(results.all, aes(x = temp, y = C, col = year)) + geom_point() + facet_wrap(~age)
ggplot(results.all, aes(x = f, y = C)) + geom_point()+  geom_abline(intercept = 0,slope = 1)+
  theme_bw()+
  facet_wrap(~age)
names(results.all)[10:11] <- c("Empirical consum", "Model prediction")
ggplot(pivot_longer(results.all, cols = 10:11, names_to = "series", values_to = "value"),  
       aes(x = B, y = value,col = series)) +
  geom_point()+geom_line() + theme_bw() +
  facet_wrap(~age) +
  xlab("Capelin maturing biomass")+
  ylab("Consumption per cod (tonnes/month)") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12)) 
ggsave("../Consumption_report_no_cod/figs/Empirical_vs_model_by_B.pdf", width = 10, height = 6)
ggplot(results.all, aes( x= B, y = temp)) + geom_point() + 
  facet_wrap(~age)
names(results.all)[10:11] <- c("f", "C")

ggplot(results.all, aes(x = year)) +
  #geom_ribbon(aes(ymin = lowAlpha, ymax = highAlpha), alpha = .5, fill = "grey")+
  geom_point(aes( y = alpha))+  geom_line(aes(y=alpha),lty = 2)+
  theme_bw()+ 
  scale_x_continuous(breaks = seq(2001,2017,2), name = "Year") +
  scale_y_continuous(name = "Year") +
  facet_wrap(~age, scales = "free")
ggsave("../Consumption_report_no_cod/figs/alpha_by_years.pdf", width = 10, height = 6)



ggplot(results.all, aes(x = year)) +
  #geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = .5, fill = "grey")+
  geom_point(aes( y = h))+  geom_line(aes(y=h),lty = 2)+
  theme_bw()+ 
  scale_x_continuous(breaks = seq(2001,2017,2), name = "Year") +
  facet_wrap(~age, scales = "free")
ggsave("../Consumption_report_no_cod/figs/h_by_years.pdf", width = 10, height = 6)

ggplot(results.all, aes(x = temp)) +
  #geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = .5, fill = "grey")+
  geom_point(aes(y = h))+  geom_line(aes(y = h),lty = 2)+
  theme_bw()+ 
  facet_wrap(~age, scales = "free_y") +
  ggsave("../Consumption_report_no_cod/figs/h_by_temp.pdf", width = 10, height = 6)
ggplot(results.all, aes(x = temp)) +
  #geom_ribbon(aes(ymin = lowH, ymax = highH), alpha = .5, fill = "grey")+
  geom_point(aes(y = alpha))+  geom_line(aes(y = alpha),lty = 2)+
  theme_bw()+ 
  facet_wrap(~age, scales = "free_y")+
  ggsave("../Consumption_report_no_cod/figs/alpha_by_temp.pdf", width = 10, height = 6)


tmp.data <- data %>% group_by(age) %>% summarize(minTemp = min(temp)-shift,
                                                 maxTemp = max(temp)+shift)

ggplot(rbind(data,test), aes(x = year, y = temp)) + geom_point(aes(col = area))+
  geom_hline(data= tmp.data, aes(yintercept = minTemp), lty = 2, col = "skyblue", lwd = .9) +
  geom_hline(data= tmp.data, aes(yintercept = maxTemp), lty = 2, col = "skyblue", lwd = .9) +
  geom_line(aes(y = meantemp), lty = 2, lwd = .7)+ geom_point(aes(y = meantemp)) +
  facet_wrap(~age) + theme_bw() +
  scale_color_discrete(name = "Area")+
  theme(legend.position = "top")+
  scale_x_continuous(breaks = seq(2001,2017,2), name  = "Years")+
  scale_y_continuous(breaks = 0:8, name = "Temperature (?C)")
ggsave("../Consumption_report_no_cod/figs/temp_by_year_with_limits.pdf", width = 10, height = 6)

filter(parest.all, par == "x")


dat2 <- tibble(
  temp = seq(minTemp-shift, maxTemp+shift, 0.01),
  fun  = (temp-(minTemp-shift)) * sqrt(maxTemp+shift - temp)
)
ggplot(dat2, aes( x= temp)) + geom_line(aes(y = fun), col = "skyblue", lwd = 1) +
  geom_vline(xintercept = minTemp-shift, lty = 2, col = "darkblue")+
  geom_vline(xintercept = maxTemp+shift, lty = 2, col = "darkblue")+
  geom_rug(data = data, aes(x = meantemp), sides = "b") + 
  scale_y_continuous(name = expression((T-T[min])(T[max]-T)^0.5), seq(0,3,0.1))+ 
  scale_x_continuous(name = "T = temperature",breaks = c(seq(4,5,.5), minTemp-shift, maxTemp+shift),
                     minor_breaks = seq(0,10,.125),
                     labels = c(paste0(seq(4,5,.5),"?C"), expression(T[min]), expression(T[max])))
ggsave("../Figures/temperature_function.pdf", device = "pdf", dpi = "retina", width = 9, height = 4)

# ---------------------------------------------
# Predictions 
# ---------------------------------------------
# keep temperature fixed, vary B
Bmean <- mean(data$cap.maturingbiomass[data$age == 3 & data$area == "I"])
meanTemp <- mean(data$temp)
B <- seq(min(data$cap.maturingbiomass)/1e3, max(data$cap.maturingbiomass)/1e3,
         length.out = 200)
Temp <- meanTemp
ages <- 3:6
dat <- data.frame()
for(i in 1:length(ages)){
  a = as.numeric(parest.all[parest.all$age == ages[i] & parest.all$par == "balpha","Estimate"]) *
    (Temp - (minTemp-shift))*sqrt(maxTemp + shift - Temp)
  h = as.numeric(parest.all[parest.all$age == ages[i] & parest.all$par == "bh","Estimate"]) *
    (Temp - (minTemp-shift))*sqrt(maxTemp + shift - Temp) 
  dat2 <- data.frame(
    age = ages[i],
    Temp= Temp,
    B   = B,
    C   = a * B/(1+a*h*B)#unlist(sapply(B, function(B) a * B/(1+a*h*B)))
  )
  dat <- rbind(dat,dat2)
  }
(p1 <- ggplot(dat, aes(x = B, y = C, col = factor(age)) ) + geom_line(lwd  = .9)+
  scale_y_continuous(name = "Model predicted consumption", breaks = seq(0,500,50)) +
  scale_x_continuous(name = "Capelin SSB (in thousand tonnes)", breaks = seq(0,2500,500)) +
  scale_color_discrete(name = "Cod age") + 
    ggtitle(paste0("T = ", round(meanTemp,2), "?C"))
)
B <- Bmean
Temp <- seq(minTemp, maxTemp, length.out = 200)
ages <- 3:6
dat <- data.frame()
for(i in 1:length(ages)){
  a = as.numeric(parest.all[parest.all$age == ages[i] & parest.all$par == "balpha","Estimate"]) *(Temp - (minTemp-shift))*sqrt(maxTemp + shift - Temp)
  h =as.numeric( parest.all[parest.all$age == ages[i] & parest.all$par == "bh","Estimate"]) *
    (Temp - (minTemp-shift))*sqrt(maxTemp + shift - Temp) 
  dat2 <- data.frame(
    age = ages[i],
    Temp= Temp,
    B   = B,
    C   = a * B/(1+a*h*B)#unlist(sapply(1:length(Temp), function(j) a[j] * B/(1+a[j]*h[j]*B)))
  )
  dat <- rbind(dat,dat2)
}
(p2 <- ggplot(dat, aes(x = Temp, y = C, col = factor(age)) ) + geom_line(lwd  = .9)+
  scale_y_continuous(name = "Model predicted consumption", breaks = seq(0,4000,500)) +
  scale_x_continuous(name = "Temperature (?C)", breaks = seq(0,10,.5)) +
  scale_color_discrete(name = "Cod age")+
    ggtitle(paste0("B = ", round(Bmean/1e3,0), " thousand tonnes")))
 # geom_text(x = Inf, y = Inf, label = paste0("B = ", round(Bmean/1e3,0), " thousand tonnes"), 
 #           hjust = 1, vjust = 1.5,
 #           col = "black"))

library(ggpubr)
theme_set(theme_bw())
ggarrange(p1 + theme(legend.position = "none"),
          p2 + theme(axis.title.y = element_blank()), ncol = 2, widths = c(1,1.2))
ggsave("../Figures/Model_predictions.pdf", width = 10, height = 4, device = "pdf", dpi = "retina")
# ---------------------------------------------------------------
# --------------------------------------------------------------
tmp <- tibble(
  year = rep(2001:2017, 4),
  temp = c(data$temp[data$area == "I" & data$age == 3],
           data$temp[data$area == "II" & data$age == 3],
           data$temp[data$area == "III" & data$age == 3],
           data$meantemp[data$area == "I" & data$age == 3]),
  area = rep(c("I","II","III","Mean"), each = 17),
  fun = (temp-minTemp + shift)*sqrt(maxTemp + shift - temp)
)
ggplot(tmp, aes(x = year, y = fun, col = area)) + geom_point() + geom_line(lwd =1) +
  geom_text(aes(label = paste0(round(temp,1), "?C")), vjust = -1, hjust = 1, size = 3,
            show.legend = FALSE) +
  scale_x_continuous(breaks = seq(2001,2017,2), name = "Year") + 
  scale_y_continuous(name = expression((T-T[min])(T[max]-T)^0.5), 
                     breaks = seq(0,3.5,.25))+
  scale_color_discrete(name = "Area")
  ggsave("../Figures/temperature_function_by_year.pdf", device = "pdf", dpi = "retina",
         width = 10, height = 4)
# ------------------------------------------------------------------
  # ---------------------------------------------
  # Predictions 
  # ---------------------------------------------
  # keep temperature fixed, vary B
  Bmean <- mean(data$cap.maturingbiomass[data$age == 3 & data$area == "I"])
  meanTemp <- mean(data$temp)
  B <- seq(min(data$cap.maturingbiomass,0)/1e3, max(data$cap.maturingbiomass)/1e3,
           length.out = 200)
  
  ages <- 3:6
  areas <- unique(data$area)
  maxTemp <- max(data2$meantemp)
  minTemp <- min(data2$meantemp)
  Temp <- seq(minTemp,maxTemp, length.out = 50)
  dat <-data.frame()
  for(i in 1:length(ages)){
      for(t in 1:length(Temp)){
        a = as.numeric(parest.all[parest.all$age == ages[i] & parest.all$par == "balpha","Estimate"]) *
          (Temp[t] - (minTemp-shift))*sqrt(maxTemp + shift - Temp[t])
        h = as.numeric(parest.all[parest.all$age == ages[i] &  parest.all$par == "bh","Estimate"]) *
          (Temp[t] - (minTemp-shift))*sqrt(maxTemp + shift - Temp[t]) 
        dat2 <- tibble(
          i = i,
          age = ages[i],
          agelab = paste0("Cod age: ", ages[i], " yrs"),
          Temp= Temp[t],
          TempC = paste0(round(Temp,1), "?C"),
          minTemp = minTemp,
          maxTemp = maxTemp,
          B   = B,
          C   = a * B/(1+a*h*B)#unlist(sapply(B, function(B) a * B/(1+a*h*B)))
        )
        dat <- rbind(dat,dat2)
      }
    
  }
  ggplot(dat, aes(x = B, y = C, col = Temp, group = factor(Temp)) ) + 
    geom_line(lwd  = .9)+ 
    facet_wrap(~agelab, scales = "fixed") +
    #scale_y_continuous(name = "Model predicted consumption", breaks = seq(0,500,50)) +
    #scale_x_continuous(name = "Capelin SSB (in thousand tonnes)", breaks = seq(0,2500,500)) +
    scale_color_viridis_c(name = "Tempe-\nrature", 
                          breaks = seq(3,6,.5), 
                          labels = paste0(seq(3,6,.5), "?C")
    ) + 
    scale_x_continuous(name = "Capelin Biomass (10^3 tonnes)", breaks = seq(0,3000,500),
                       limits = c(0,2500))+
    scale_y_continuous(name = "Consumption per cod (tonnes/month)", breaks = seq(0,1,.2))+
    guides(color = guide_colorbar(barheight = 25))
  ggsave("../Figures/age_prediction_curves.pdf", device = "pdf", dpi = "retina",
         width = 10, height = 6)
  