setwd("~/Documents/GitHub/EPP")

library(dplyr)
library(magrittr)

devtools::install_github("mrc-ide/eppasm@new-master")
library(eppasm)
install.packages("testthat")
devtools::load_all()

files <- c("~/Documents/Data/Spectrum files/2018 final/SSA/Ghana_2018_final_v5_65.PJNZ", "~/Documents/Data/Spectrum files/2018 final/SSA/Botswana_ 2018 updated ART.PJNZ")
# prop_75 <- vector("list", length(files))

prop_75 <- rep(list(rep(list(vector("list", 100)), length(files))), 2)
names(prop_75) <- c("Incidence", "Mortality")

for (j in 1:length(files)) {

  ## Preparing EPP-ASM inputs
  pjnz <- "~/Documents/Data/Spectrum files/2018 final/SSA/Rwanda _2018_final.PJNZ"
  # pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
  pjnz <- files[j]
  bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
  
  fp <- prepare_directincid(pjnz)
  
  ####### Direct incidence #########
  
  fp_new <- prepare_directincid(pjnz)
  
  x <- c(2010, 2020)
  y <- c(fp[["incidinput"]][41], fp[["incidinput"]][41]*0.25)
  
  test2 <- data.frame(approx(x, y, xout=2010:2020))
  
  fp_new[["incidinput"]][41:51] <- test2$y
  
  plot(fp_new[["incidinput"]])
  
  mod_new <- simmod(fp_new)
  
  ########## Simulate HHS #############
  
  sim_2020 <- data.frame(rbinom(100, 3000, prev(mod_new)[51])/3000) # 51st item in prev(mod_new) is the 2020 prevalence generated from 2020 incidence at 25% of 2010 incidence.
  sim_hhs <- vector("list", 100)
  
  # prop_75 <- data.frame("prop_above_75" = 0)
  
  for (i in 1:3) {
    
    sim_hhs[[i]] <- attr(bw$Urban, "eppd")$hhs %>%
        rbind(data.frame("year" = 2020, "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = sim_2020[i,], "se"=sqrt((sim_2020[i,]*(1-sim_2020[i,]))/3000),"deff"= 2, "deff_approx"= 2, "used"=TRUE))
    
    attr(bw$Urban, "eppd")$hhs <- sim_hhs[[i]]
    
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
    bwfit <- list()
      
    bwfit$Urban <- fitmod(bw$Urban, eppmod = "rhybrid", rw_start = 2005,
                            B0=1e4, B=1e2, opt_iter = 1, number_k=50)
    # bwfit$Rural <- fitmod(bw$Rural, eppmod = "rhybrid", rw_start = 2005,
    #                        B0=1e3, B=1e2, opt_iter = 1, number_k=50)
      
      #' When fitting, the random-walk based models only simulate through the end of the
      #  data period. The `extend_projection()` function extends the random walk for $r(t)$
      #  through the end of the projection period.
      
    bwfit <- lapply(bwfit, extend_projection, proj_years = 52)
      
      
      ## Simulating model outptus
      
    bwout <- Map(tidy_output, bwfit, "r-hybrid", attr(bw$Urban, "eppd")$country, names(bwfit))
      
    bwaggr <- aggr_specfit(bwfit)
    
    reduction <- data.frame("incidence_red" = 1-sapply(1:3000, function(x) {
      attr(bwaggr[[x]], "incid15to49")[51]/attr(bwaggr[[x]], "incid15to49")[41]
    }), "mort_red" = 1-sapply(1:3000, function(x) {
      attr(bwaggr[[x]], "hivdeaths")[51]/attr(bwaggr[[x]], "hivdeaths")[41]
    })) 
    
    prop_75[[1]][[j]][i] <- sum(reduction$incidence_red >= 0.75)/3000
    prop_75[[2]][[j]][i] <- sum(reduction$mort_red >= 0.75)/3000
  }
  
  names(prop_75)[[1]] <- attr(bw$Urban, "eppd")$country

}

