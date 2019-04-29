setwd("~/Documents/GitHub/epp_2020_target")

devtools::install_github("mrc-ide/eppasm@new-master")
library(eppasm)
library(dplyr)
library(magrittr)
library(tidyverse)

############### FUNCTIONS TO RUN ##############

fp_old_incid <- prepare_directincid(pjnz)

fp_new_incid <- prepare_directincid(pjnz)

x <- c(2010, 2020)
y <- c(fp_old_incid[["incidinput"]][41], fp_old_incid[["incidinput"]][41]*0.25)

test2 <- data.frame(approx(x, y, xout=2010:2020))

fp_new_incid[["incidinput"]][41:51] <- test2$y
fp_new_incid[["incidinput"]][52:length(fp_new_incid[["incidinput"]])] <- fp_new_incid[["incidinput"]][51]

test_mod <- simmod(fp_new_incid)

#####################

anc_project <- function(obj, theta, ancsite=TRUE) {
  
  # idvars <- data.frame(country = "country",
  #                      eppregion = "eppregion",
  #                      modlab = "modlab")
  
  ss <- fp$ss
  
  param_list <- list(fnCreateParam(theta, fp))
  fp_list <- lapply(param_list, function(par) update(fp, list=par))
  # mod_list <- lapply(fp_list, simmod)
  mod_list <- list(test_mod) #Is this a heinous crime? Who knows
  
  ancsite.dat <- prepare_ancsite_likdat(attr(obj$Urban, "eppd")$ancsitedat, fp)
  
  b_site <- Map(sample_b_site, mod_list, fp_list, list(ancsite.dat), resid = FALSE)
  
  b_site_sigma <- sapply(b_site, anclik::sample.sigma2)
  
  ancsite_b <- data.frame(site = names(b_site[[1]]), mean = b_site[[1]])
  
  newdata <- expand.grid(site = unique(ancsite.dat$df$site),
                         year = 1985:2030,
                         type = "ancss",
                         age = 15,
                         agspan = 35,
                         n = 300)
  new_df <- ancsite_pred_df(newdata, fp)
  new_df$df$offset <- 0
  
  ancsite_pred <- mapply(sample_ancsite_pred, mod_list, fp_list,
                         b_site = b_site,
                         MoreArgs = list(newdata = new_df))
  
  ancsite_pred <- pnorm(ancsite_pred)
  
  
  # ancsite_pred <- data.frame(newdata, estci2(ancsite_pred))
  ancsite_pred <- data.frame(newdata, "prev" = ancsite_pred)
  #
  ancsite_pred <- merge(ancsite_pred,
                        ancsite.dat$df[c("site", "year", "type", "age", "agspan", "n", "prev")],
                        by = c("site", "year", "type", "age", "agspan"),
                        suffixes = c("_sim", "_obs"), all.x=TRUE)
  #
  # ancsite_pred <- data.frame(idvars, ancsite_pred)
  
  return(ancsite_pred)
  
} 

debugonce(anc_project)
anc_project(obj, theta_u)


merge_sim_anc <- function(obj, theta, clean_anc) {
  
  get_max_year <- anc_project(obj, theta) %>%
    group_by(site) %>%
    filter(!is.na(n_obs)) %>%
    summarise(max = max(year))
  
  years <- as.numeric(get_max_year$max)
  sites <- as.character(get_max_year$site)
  
  nest_sites <- anc_project(obj, theta) %>%
    group_by(site) %>%
    nest()
  
  filter_years <- function(a, b) {
    a%<>%
      filter(year>b) %>%
      mutate(prev=prev_sim) %>%
      mutate(n = n_sim)
  }
  
  filtered_max_year <- map2(nest_sites$data, years, filter_years)
  
  mapply(function(a, b){filtered_max_year[[a]]$site <<- b}, b=sites, a=1:10)
  
  # for (i in 1:length(filtered_max_year)) {
  #   filtered_max_year[[i]]$site <- sites[i]
  # }
  
  filtered_max_year <- data.frame(do.call(rbind, filtered_max_year)) %>%
    mutate(agegr = "15-49") %>%
    mutate(used = TRUE) %>%
    select(c("site", "year", "used",  "prev", "n", "type", "agegr", "age","agspan"))
  
  merge_anc <- clean_anc %>%
    rbind(filtered_max_year)
    

  return(merge_anc)
  
}

################ NEW FILES###############################


files <- c("~/Documents/Data/Spectrum files/2018 final/SSA/Botswana_ 2018 updated ART.PJNZ", "~/Documents/2018-12 2020 targets/Botswana_ 2030.PJNZ")

pjnz <- files[2]

obj <- prepare_spec_fit(pjnz, 2030.5)

## Don't understand what Jeff has done here ##
opt <- fitmod(obj$Urban, eppmod = "rhybrid", optfit = TRUE, B0=1e3, opthess = FALSE)
opt$resample <- matrix(opt$par, nrow = 1)
opt_full <- extend_projection(opt, opt$fp$ss$PROJ_YEARS)

theta_u <- as.numeric(opt_full$resample[1,])
fp <- opt_full$fp

# get_incid <- prepare_directincid(pjnz)

## Set some flags that are set in fitmod(), (later improve this code...)
fp <- prepare_anc_model(fp, attr(obj$Urban, "eppd"))
# fp$ancsitedata <- TRUE
# fp$ancrt <- "both"
fp$logitiota <- TRUE
fp$rw_start <- 2005
fp$incidmod <- "eppspectrum"

param <- fnCreateParam(theta_u, fp)
fp_par <- update(fp, list = param)

clean_anc <- attr(obj$Urban, "eppd")$ancsitedat

attr(obj$Urban, "eppd")$ancsitedat <- merge_sim_anc(obj, theta_u, clean_anc)

mod <- simmod(fp_par)

plot(prev(mod))
prev(mod)

## Prepare likelihood and calculate the likelihood once

## Prior
lprior(theta_u, fp)

## Prepare likelihood data
likdat <- prepare_likdat(attr(obj$Urban, "eppd"), fp)

## Calculate likelihood
ll(theta_u, fp, likdat)

## Components of likelihood calculation
ll_hhsage(mod, likdat$hhs.dat)

ll_ancsite(mod, fp_par,
           coef = c(fp_par$ancbias, fp_par$ancrtsite.beta),
           vinfl = fp_par$v.infl,
           dat = likdat$ancsite.dat)

ll_ancrtcens(mod, likdat$ancrtcens.dat, fp_par)

## Fitting the EPP-ASM model
bwfit <- list()

bwfit$Urban <- fitmod(obj$Urban, eppmod = "rhybrid", rw_start = 2005,
                      B0=1e4, B=1e2, opt_iter = 1, number_k=50)
# bwfit$Rural <- fitmod(bw$Rural, eppmod = "rhybrid", rw_start = 2005,
#                        B0=1e3, B=1e2, opt_iter = 1, number_k=50)

#' When fitting, the random-walk based models only simulate through the end of the
#  data period. The `extend_projection()` function extends the random walk for $r(t)$
#  through the end of the projection period.

bwfit <- lapply(bwfit, extend_projection, proj_years = 60)
# for 2030, extend to 60.


## Simulating model outptus
debugonce(tidy_output)
bwout <- Map(tidy_output, bwfit, "r-hybrid", attr(obj$Urban, "eppd")$country, names(bwfit))

bwaggr <- aggr_specfit(bwfit)



#########################################

ggplot()+  
  geom_line(data=min_year(bw) %>% filter(!is.na(prev)) %>% filter(year<2012), aes(x=year, y=prev, color=site)) +
  geom_line(data=min_year(bw) %>% filter(!is.na(prev)) %>% filter(year>2010), aes(x=year, y=prev, color=site), linetype = 2)

ggplot(data=attr(bw$Urban, "eppd")$hhs, aes(x=year, y=prev)) +
      geom_point()+
      geom_line(data=bwout$Urban$core %>% filter(indicator=="prev"), aes(x=year, y=mean))

