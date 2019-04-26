setwd("~/Documents/GitHub/EPP")

library(dplyr)
library(magrittr)

devtools::install_github("mrc-ide/eppasm@new-master")
library(eppasm)
install.packages("testthat")
devtools::load_all()

## Preparing EPP-ASM inputs
pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
bw <- prepare_spec_fit(pjnz, proj.end=2024.5)
bw_hhs <- prepare_spec_fit(pjnz, proj.end=2024.5)
bw_hhs_anc <- prepare_spec_fit(pjnz, proj.end=2024.5)



## Simulate model with direct incidence input (a la Spectrum)

fp <- prepare_directincid(pjnz)
fp_new <- prepare_directincid(pjnz)

x <- c(2010, 2020)
y <- c(fp[["incidinput"]][41], fp[["incidinput"]][41]*0.25)

test2 <- data.frame(approx(x, y, xout=2010:2020))

fp_new[["incidinput"]][41:51] <- test2$y

plot(fp_new[["incidinput"]])

mod <- simmod(fp)
mod_new <- simmod(fp_new)

outputs <- list(prevalence = prev(mod), incidence = incid(mod))
names(outputs[[1]]) <- 1970:2025
names(outputs[[2]]) <- 1970:2025

outputs_new <- list(prevalence = prev(mod_new), incidence = incid(mod_new))
names(outputs_new[[1]]) <- 1970:2025
names(outputs_new[[2]]) <- 1970:2025

outputs[["prevalence"]]["2020"]
outputs_new[["prevalence"]]["2020"]

#####################################

attr(bw_hhs_anc$Urban, "eppd")$hhs %<>%
  rbind(data.frame("year" = 2020, "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = 0.2, "se"=0.01,"deff"= 2, "deff_approx"= 2, "used"=TRUE))

attr(bw_hhs_anc$Rural, "eppd")$hhs %<>%
  rbind(data.frame("year" = 2020, "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = 0.21, "se"=0.01,"deff"= 2, "deff_approx"= 2, "used"=TRUE))
  
attr(bw_hhs_anc$Urban, "eppd")$ancsitedat %<>%
  rbind(data.frame("site" = "Francistown (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.32 , "n" = 547 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Gaborone (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.26 , "n" = 572 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Jwaneng (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.137 , "n" = 300 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Kweneng East (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.236 , "n" = 262 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Lobatse (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.214 , "n" = 169 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Mahalapye (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.305 , "n" = 463 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Selebi-Pikwe (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.382 , "n" = 282 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "SerowePalapye (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.349 , "n" = 458 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "SouthEast (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.128 , "n" = 242 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Southern (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.194 , "n" = 302 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35))

attr(bw_hhs_anc$Rural, "eppd")$ancsitedat %<>%
  rbind(data.frame("site" = "Bobirwa district (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.368 , "n" = 363 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Boteti (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.224 , "n" = 197 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Chobe (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.321 , "n" = 143 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Gantsi (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.123 , "n" = 135 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Goodhope (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.258 , "n" = 194 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Hukuntsi (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.255 , "n" = 118 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Kgalagadi (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.222 , "n" = 191 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Kgatleng (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.295 , "n" = 330 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Kweneng West (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.28 , "n" = 157 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Mabutsane (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.317 , "n" = 300 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Ngami (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.228 , "n" = 135 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "North East (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.246 , "n" = 198 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Okavango (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.212 , "n" = 327 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35)) %>%
  rbind(data.frame("site" = "Tutume (%)", "year" = 2020 , "used" = TRUE , "prev" = 0.309 , "n" = 403 , "type" = "ancss" , "agegr" = "15-49" , "age" = 15 , "agspan" = 35))

ancRural <- attr(bw_hhs_anc$Rural, "eppd")$ancsitedat
ancUrban <- attr(bw_hhs_anc$Urban, "eppd")$ancsitedat
## EPP model with single set of parameter inputs

# Prepare model parameters from parameter vector input

theta_ur <- c(-0.63758, -2.76655, -1.26204, 1996.65945, 0.00778, 0.05195,
              0.05103, 0.032, 0.01765, 0.01154, -0.00028, 0.01627, -0.00051,
              0.01439, -0.00937, -0.01135, 0.03692, 0.14959, 0.00803, 0.02424,
              -0.03548, 3.65223, -0.02515, -4.74563, 0.26259, -6.90124, 0.01583)

fp <- attr(bw$Urban, "specfp")
fp <- prepare_rhybrid(fp, rw_start = 2005, rw_dk = 1)

## Set some flags that are set in fitmod(), (later improve this code...)
fp$ancsitedata <- TRUE
fp$ancrt <- "both"
fp$logitiota <- TRUE
fp$rw_start <- 2005
fp$incidmod <- "eppspectrum"

param <- fnCreateParam(theta_ur, fp)
fp_par <- update(fp, list = param)

## Simulate the model once.

mod <- simmod(fp_par)

prev(mod)
incid(mod)

## R version of the model
modR <- simmod(fp_par, VERSION = "R")

prev(modR)
incid(modR)

## Prepare likelihood and calculate the likelihood once

## Prior 
lprior(theta_ur, fp)

## Prepare likelihood data
likdat <- prepare_likdat(attr(bw$Urban, "eppd"), fp)

## Calculate likelihood
ll(theta_ur, fp, likdat)

## Components of likelihood calculation
ll_hhsage(mod, likdat$hhs.dat)

ll_ancsite(mod, fp_par,
           coef = c(fp_par$ancbias, fp_par$ancrtsite.beta),
           vinfl = fp_par$v.infl,
           dat = likdat$ancsite.dat)

ll_ancrtcens(mod, likdat$ancrtcens.dat, fp_par)

## Fitting the EPP-ASM model
{
bwfit <- list()

bwfit$Urban <- fitmod(bw$Urban, eppmod = "rhybrid", rw_start = 2005,
                      B0=1e3, B=1e2, opt_iter = 1, number_k=50)
bwfit$Rural <- fitmod(bw$Rural, eppmod = "rhybrid", rw_start = 2005,
                      B0=1e3, B=1e2, opt_iter = 1, number_k=50)

#' When fitting, the random-walk based models only simulate through the end of the
#  data period. The `extend_projection()` function extends the random walk for $r(t)$
#  through the end of the projection period.

bwfit <- lapply(bwfit, extend_projection, proj_years = 52)


## Simulating model outptus

bwout <- Map(tidy_output, bwfit, "r-hybrid", "Botswana", names(bwfit))

old_urban  <- data.frame("year" = bwout[["Urban"]][["core"]][["year"]], "indicator" = bwout[["Urban"]][["core"]][["indicator"]], "mean" = bwout[["Urban"]][["core"]][["mean"]], "source" = "old")
old_rural  <- data.frame("year" = bwout[["Rural"]][["core"]][["year"]], "indicator" = bwout[["Rural"]][["core"]][["indicator"]], "mean" = bwout[["Rural"]][["core"]][["mean"]], "source" = "old")

}
##########################

bwfit_hhs <- list()

bwfit_hhs$Urban <- fitmod(bw_hhs$Urban, eppmod = "rhybrid", rw_start = 2005,
                          B0=1e3, B=1e2, opt_iter = 1, number_k=50)
bwfit_hhs$Rural <- fitmod(bw_hhs$Rural, eppmod = "rhybrid", rw_start = 2005,
                          B0=1e3, B=1e2, opt_iter = 1, number_k=50)

#' When fitting, the random-walk based models only simulate through the end of the
#  data period. The `extend_projection()` function extends the random walk for $r(t)$
#  through the end of the projection period.

bwfit_hhs <- lapply(bwfit_hhs, extend_projection, proj_years = 52)


## Simulating model outptus

bwout_hhs <- Map(tidy_output, bwfit_hhs, "r-hybrid", "Botswana", names(bwfit_hhs))

hhs_urban  <- data.frame("year" = bwout_hhs[["Urban"]][["core"]][["year"]], "indicator" = bwout_hhs[["Urban"]][["core"]][["indicator"]], "mean" = bwout_hhs[["Urban"]][["core"]][["mean"]], "source" = "hhs")
hhs_rural  <- data.frame("year" = bwout_hhs[["Rural"]][["core"]][["year"]], "indicator" = bwout_hhs[["Rural"]][["core"]][["indicator"]], "mean" = bwout_hhs[["Rural"]][["core"]][["mean"]], "source" = "hhs")


##############################

bwfit_hhs_anc <- list()

bwfit_hhs_anc$Urban <- fitmod(bw_hhs_anc$Urban, eppmod = "rhybrid", rw_start = 2005,
                              B0=1e3, B=1e2, opt_iter = 1, number_k=50)
bwfit_hhs_anc$Rural <- fitmod(bw_hhs_anc$Rural, eppmod = "rhybrid", rw_start = 2005,
                              B0=1e3, B=1e2, opt_iter = 1, number_k=50)

#' When fitting, the random-walk based models only simulate through the end of the
#  data period. The `extend_projection()` function extends the random walk for $r(t)$
#  through the end of the projection period.

bwfit_hhs_anc <- lapply(bwfit_hhs_anc, extend_projection, proj_years = 52)


## Simulating model outptus

bwout_hhs_anc <- Map(tidy_output, bwfit_hhs_anc, "r-hybrid", "Botswana", names(bwfit_hhs_anc))

hhs_anc_urban  <- data.frame("year" = bwout_hhs_anc[["Urban"]][["core"]][["year"]], "indicator" = bwout_hhs_anc[["Urban"]][["core"]][["indicator"]], "mean" = bwout_hhs_anc[["Urban"]][["core"]][["mean"]], "source" = "hhs_anc")
hhs_anc_rural  <- data.frame("year" = bwout_hhs_anc[["Rural"]][["core"]][["year"]], "indicator" = bwout_hhs_anc[["Rural"]][["core"]][["indicator"]], "mean" = bwout_hhs_anc[["Rural"]][["core"]][["mean"]], "source" = "hhs_anc")

#####################################

plot_hhs_urban <- attr(bw_hhs$Urban, "eppd")$hhs %>%
  select(year, prev, se) %>%
  rename(mean = prev)

urban_plot <- old_urban %>%
  rbind(hhs_urban) %>%
  rbind(hhs_anc_urban) %>%
  mutate(source = factor(source)) %>%
  filter(indicator == "prev") %>%
  filter(year<2021) %>%
  ggplot(aes(x=year, y=mean)) +
    geom_line(aes(color=source)) +
    geom_point(data=plot_hhs_urban %>% filter(year<2018)) +
    geom_errorbar(data=plot_hhs_urban %>% filter(year<2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), width=1)+
    geom_point(data=plot_hhs_urban %>% filter(year>2018), color="red") +
    geom_errorbar(data=plot_hhs_urban %>% filter(year>2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), color="red", linetype = 2, width=1)+
    geom_point(data=plot_hhs_urban %>% filter(year<2018)) +
    geom_errorbar(data=plot_hhs_urban %>% filter(year<2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), width=1)+
    geom_point(data=plot_hhs_urban %>% filter(year>2018), color="red") +
    geom_errorbar(data=plot_hhs_urban %>% filter(year>2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), color="red", linetype = 2, width=1)+
    geom_point(data=ancUrban, aes(x=year, y=prev), size=0.5, color="purple") +
    geom_line(data=ancUrban, aes(x=year, y=prev, group=site), linetype=3, color="purple")

urban_plot

plot_hhs_rural <- attr(bw_hhs$Rural, "eppd")$hhs %>%
  select(year, prev, se) %>%
  rename(mean = prev)

rural_plot <- old_rural %>%
  rbind(hhs_rural) %>%
  rbind(hhs_anc_rural) %>%
  mutate(source = factor(source)) %>%
  filter(indicator == "prev") %>%
  filter(year<2021) %>%
  ggplot(aes(x=year, y=mean)) +
  geom_line(aes(color=source)) +
  geom_point(data=plot_hhs_rural %>% filter(year<2018)) +
  geom_errorbar(data=plot_hhs_rural %>% filter(year<2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), width=1)+
  geom_point(data=plot_hhs_rural %>% filter(year>2018), color="red") +
  geom_errorbar(data=plot_hhs_rural %>% filter(year>2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), color="red", linetype = 2, width=1)+
  geom_point(data=plot_hhs_rural %>% filter(year<2018)) +
  geom_errorbar(data=plot_hhs_rural %>% filter(year<2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), width=1)+
  geom_point(data=plot_hhs_rural %>% filter(year>2018), color="red") +
  geom_errorbar(data=plot_hhs_rural %>% filter(year>2018), aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se)), color="red", linetype = 2, width=1)+
  geom_point(data=ancRural, aes(x=year, y=prev), size=0.5, color="purple") +
  geom_line(data=ancRural, aes(x=year, y=prev, group=site), linetype=3, color="purple")

rural_plot

grid.arrange(urban_plot, rural_plot, ncol=2)

## Pooling EPP subpopulation results

# The function `aggr_specfit()` 
# This involves simulating the model for all resamples in each subregion and summing the following `pop`, `hivpop`, and `artpop` arrays for each of the 3000 resamples to generate 3000 national outputs.

bwaggr <- aggr_specfit(bwfit)

fpp <- data.frame("incidence_red" = 1-sapply(1:3000, function(x) {
  attr(bwaggr[[x]], "incid15to49")[51]/attr(bwaggr[[x]], "incid15to49")[41]
})) 
fpp %>%
  ggplot(aes(x=incidence_red)) +
    geom_histogram()

prop_75 <- sum(fpp$incidence_red >= 0.75)/3000