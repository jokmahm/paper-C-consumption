# -- load packages --
library(TMB)
library(tidyverse)

# -- load data --
load("data.RData")
sdtemp <- sd(data$temp-data$meantemp)
theme_set(theme_bw()+
            theme(strip.background = element_rect(fill = "transparent", color = "transparent"),
                  strip.text = element_text(size = 11)))

compile("Consumption_selective_eps.cpp")
dyn.load(dynlib("Consumption_selective_eps"))
data$area <- data$paperarea
parameters <- list(x=0,
                   logb_alpha = 0,
                   logb_h = 0, 
                   BH = 0,
                   theta_eps_alpha =0,
                   #logeps_h =0,
                   logSigma = log(5e-2)
)

# -- all areas --

map <- list(x = factor(NA), 
            BH = factor(NA))

test <- filter(data, year >2013)
data <- filter(data, year <=2013)

shift <- 0.1

parest.all <- results.all <- pred.all <-tibble()
maxTemp <- max(data$meantemp)
minTemp <- min(data$meantemp)
area <- data$area[1]
ages <- 3:6
nsim <- 200
set.seed(123412)
for(sim in 1:nsim){
for(i in 1:length(ages)){
  age <- ages[i]
  #maxTemp <- max(data$meantemp[data$age == age])
  #minTemp <- min(data$meantemp[data$age == age])
  # temp <- rnorm(length(unique(data$year)), 
  #               mean = data$meantemp[data$age == age & data$area == area],
  #               sd=sd(data$meantemp[data$age == age & data$area == area]))
  temp <- runif(length(unique(data$year)), 
                min = minTemp,
                max =maxTemp)
  data.list <- list(ny = length(unique(data$year)),
               T_0 = minTemp - shift,
               T_l = maxTemp + shift,
               P = data$cod.maturingbiomass[data$age == age & data$area == area]/1e3,
               B = data$cap.maturingbiomass[data$age == age & data$area == area]/1e3,
               Temp =  temp,
               f = data$EmpC[data$age == age & data$area == area],
               predP = test$cod.maturingbiomass[test$age == age & test$area == area]/1e3,
               predB = test$cap.maturingbiomass[test$age == age & test$area == area]/1e3,
               predTemp = test$meantemp[test$age == age & test$area == area],
               logbhTRUE = 1)
  obj     <- MakeADFun(data.list, parameters,  DLL = "Consumption_selective_eps", silent = TRUE,
                       map = map)#logb_alpha=factor(NA)))#NULL) # list(x=factor(NA)))
  opt     <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  rep <- sdreport(obj)
  ss <- summary(rep, p.value = TRUE)[-(1:(7-length(map))),]
  
  parest <-as_tibble(ss[1:5, ])
  parest$par = c("x", "balpha", "bh", "eps_a", "sigma")
  parest$age <- age
  parest$sim = sim
  parest$conv = opt$convergence
  parest.all <- rbind(parest.all, parest)
  results <- tibble(
    conv = opt$convergence,
    B = data.list$B,
    P = data.list$P, 
    year = unique(data$year),
    age = age,
    area = area,
    temp = data.list$Temp,
    minTemp = data.list$T_0,
    maxTemp = data.list$T_l,
    f = data.list$f,
    C= ss[rownames(ss) == "C",1],
    seC= ss[rownames(ss) == "C",2],
    lowC = pmax(C -2 *seC,0),#*sd(f), 
    highC = C+2*seC,#*sd(f),
    alpha= ss[rownames(ss) == "alpha",1],
    sealpha= ss[rownames(ss) == "alpha",2],
    lowAlpha = pmax(alpha-2 * sealpha,0),
    highAlpha = alpha+2 * sealpha,
    h= ss[rownames(ss) == "h",1],
    seh= ss[rownames(ss) == "h",2],
    lowH = pmax(h-2 * seh,0),
    highH = h+2 * seh,
    sim = sim
  )
  pred <- tibble(
    year = unique(test$year),
    age = age,
    area = area,
    f =  test$EmpC[test$age == age & test$area == area],
    C= ss[rownames(ss) == "predC",1],
    seC= ss[rownames(ss) == "predC",2],
    lowC = pmax(C -2 *seC,0),
    highC = C+2*seC,
    sim = sim,
    conv = opt$convergence
  )
  results.all <- rbind(results.all, results)
  pred.all <- rbind(pred.all, pred)
}
}
cat("---------\nDid all converge?\n", ifelse(all(results.all$conv==0), "- YES","- NO"), 
    "\n---------\n")

(parest.all <- parest.all[,c(8,7,6,5,1:3)])
parest.all <- filter(parest.all, conv ==0)

library(xtable)
names(parest.all)[ncol(parest.all)] <- "P-value"

parest.all$'p-value' <- ifelse(parest.all$'P-value' < 0.05, 
                               paste0("\\color{green} ", round(parest.all$'P-value',3)),
                               paste0("\\color{red} ", round(parest.all$'P-value',3)))
# print(xtable(parest.all[,-(ncol(parest.all)-1)], digits= 5,
#              label = "tab1",
#              caption = "Model parameter estimates based on observations for years 2001-2013. The years 2014-2017 are used for model validation. P-value is red/green if it is larger/smaller than 0.05."),
#       include.rownames = FALSE, file = "../Consumption_report_no_cod/by_age_2013.tex",
#       hline.after = c(-1,0,4,8,12,16),
#       sanitize.text.function = function(x)x)
#parest.all <- filter(parest.all, !(par == "balpha" & Estimate > 1e2))
ggplot(filter(parest.all, par !="x"), aes(x = Estimate)) + facet_wrap(par~age, scales = "free", ncol = 4)+
  stat_density() +theme_bw() +
  ggsave("../Consumption_report_no_cod/figs/kde_montecarlo.pdf", width = 10, height = 6)

pred.all <- filter(pred.all, conv ==0, 
                   !is.nan(seC))
results.all <- filter(results.all, conv ==0, 
                      !is.nan(seC))

pred.all <- rbind(pred.all, 
results.all[results.all$year == 2013, c("year", "age", "area","f", "C", "seC", "lowC", "highC","sim", "conv")]
)

names(pred.all)
results.all$agelab = paste0("Cod age: ", results.all$age, " yrs")

pred.all$agelab = paste0("Cod age: ", pred.all$age, " yrs")
ggplot(results.all, aes(x = year, group =sim))+
  geom_ribbon(aes(ymin = lowC, ymax = highC), alpha = .5, fill = "grey") +
  facet_wrap(~agelab, strip.position = "top", scales = "free_y") + 
  geom_line(data = pred.all, aes(x = year, y = C), col = "red") + 
  geom_ribbon(data = pred.all, aes(x= year, ymin = lowC, ymax = highC), fill = "red", alpha = .2)+
  scale_x_continuous(breaks = seq(2001,2017,2), name = "Year") +
   geom_line(aes(y=C), lty = 2) +
  
  geom_point(data = pred.all, aes(x = year, y = f), col = "black")+ 
  geom_point(aes(y = f), col = "red")+
  scale_y_continuous(name = "Consumption per cod (kg/month)") 
#  ggsave("../Consumption_report_no_cod/figs/1_consumption_per_year_and_age_2013_montecarlo.pdf", width = 10, height = 6)
  ggsave("../Figures/eps_montecarlo_simulated_temperatures_consumption.pdf", device = "pdf", dpi = "retina", 
         width = 10, height = 6)
  


# -- Prediction --
test



ggplot(results.all, aes(x = temp, y = C, col = year)) + geom_point()+facet_wrap(~age)
ggplot(results.all, aes(x = f, y = C)) + geom_point()+  geom_abline(intercept = 0,slope = 1)+
  theme_bw()+
  facet_wrap(~age)
names(results.all)[10:11] <- c("Empirical consum", "Model prediction")
ggplot(pivot_longer(results.all, cols = 10:11, names_to = "series", values_to = "value"),  
       aes(x = B, y = value,col = series)) +
  geom_point()+geom_line() + theme_bw() +
  facet_wrap(~age) +
  xlab("Capelin maturing biomass")+
  ylab("Consumption per hour") +
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
# ---------------------------------------------------------------
# --------------------------------------------------------------
