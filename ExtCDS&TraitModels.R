##################################################################
### Code to assemble data for analyses                         ###
### Seedclim climate difference spp responses by traits paper  ###
### Exctinction by species                                     ###
### Lynn et al.                                                ###
##################################################################

# load requiredd packages
library(tidyverse); library(car); library(rjags); library(R2jags); library(coda)
library(broom.mixed)

# inverse logit function
invlogit<-function(x){a <- exp(x)/(1+exp(x))
a[is.nan(a)]=1
return(a)
}

# load in data
dat <- read_csv("LynnETAL_Macrocontext_trait_data.csv")

### start with model of extinction and MAT CDS
# group level effects
site <- factor(dat$destSiteID)
block <- factor(dat$destBlockID)
spp <- factor(dat$species)

# independent variables
mat <- dat$mat_meandiff

#dependent variable
ext <- dat$extinct
N <- as.numeric(length(ext))

# assemble the data and required params
jags.data <- list("site", "block","spp", "mat", "ext", "N")
jags.param <- c("a", "b", "prec1","prec2", "asig1","asig2")

# Step 1 - CDS model
matdifmod <- function(){ 
  # group effects
  for (j in 1:12){nettstedet[j]~dnorm(0, prec1)}
  for (j in 1:60){blokkere[j]~dnorm(0, prec2)}
  #likelihood
  for (i in 1:N){
    ext[i]~dbern(mu[i])
    logit(mu[i]) <- a[spp[i]]+ b[spp[i]]*mat[i] + nettstedet[site[i]]+ 
      blokkere[block[i]]
  }
  # priors
  for(j in 1:75){a[j]~dnorm(0, 1E-6)}
  for(j in 1:75){b[j]~dnorm(0, 1E-6)}
  prec1~dgamma(0.001,0.001)
  prec2~dgamma(0.001,0.001)
  asig1 <- 1/sqrt(prec1)
  asig2 <- 1/sqrt(prec2)
}

# run the model
matdif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=500000 ,model.file=matdifmod, n.thin=5, n.chains=3)

# View results
matdif 

# put sims in a list
matdif.paramlist <- matdif$BUGSoutput$sims.list

# save poseteriors
matint <- data.frame(matdif.paramlist$a)
names(matint) <- c(levels(spp))
matint <- data.frame(apply(matint, 2, sort))
write.csv(matint, "SppExtIntMAT_post_pos.csv")
matint <- matint %>% gather(key="spp", value="int")

matslop <- as.data.frame(matdif.paramlist$b)
names(matslop) <- c(levels(spp))
matslop <- data.frame(apply(matslop, 2, sort))
write.csv(matslop,"SppExtSlopeMAT_post_pos.csv" )
matslop <- matslop %>% gather(key="spp", value="slope")

matint$slope <- matslop$slope

# summarize posteriors
matdifparm <- matint %>% group_by(spp) %>%
  summarize(m_int=median(int),
            max_int=max(int), 
            min_int=min(int), 
            sig_int=sd(int),
            m_slope=median(slope),
            max_slope=max(slope),
            min_slope=min(slope),
            sig_slope=sd(slope))

# prepare data to combine
spprangemat <- dat %>% group_by(species) %>%
  summarize(minmat = min(mat_meandiff),
            maxmat = max(mat_meandiff),
            meanmat = mean(mat_meandiff))%>%
  rename(spp=species)

matdifparm <- full_join(matdifparm, spprangemat, by="spp", copy=TRUE)

# add in the trait data
traitformod <-  dat %>% select(species, m_height_spp:seed_mass) %>% unique()
matdifparm <- left_join(matdifparm, traitformod, by=c("spp"="species"))

# drop species that didn't converge
matdifparm <- matdifparm %>% filter(!spp %in% c("Car.pan","Dia.del","Gal.ver", "Noc.cae", "Sax.aiz"))
dat2 <- dat %>% filter(!spp %in% c("Car.pan","Dia.del","Gal.ver", "Noc.cae", "Sax.aiz"))

# classify CDS
matdifparm <- matdifparm %>% 
  mutate(slope_cat = if_else(m_slope/sig_slope <  -1, "Negative",
                             if_else(m_slope/sig_slope < 1 & m_slope/sig_slope > -1,
                                     "Zero", "Positive")))

# construct models of trait predictions
# make the different data sets
sladat <- filter(matdifparm, !is.na(m_sla_spp))
slacat <- model.matrix(~factor(sladat$slope_cat, levels=c("Zero","Negative",  "Positive")))
sla <- log(sladat$m_sla_spp);hist(sla)
N_sla <- as.numeric(length(sla))

heightdat <- filter(matdifparm, !is.na(m_height_spp))
hcat <- model.matrix(~factor(heightdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
height <- log(heightdat$m_height_spp); hist(height)
N_h <- as.numeric(length(height))

seeddat <- matdifparm %>% filter(!is.na(seed_mass)) 
seedcat <- model.matrix(~factor(seeddat$slope_cat, levels=c("Zero","Negative",  "Positive")))
seed <- log(seeddat$seed_mass); hist(seed)
N_seed <- as.numeric(length(seed))

ldmcdat <- filter(matdifparm, !is.na(m_ldmc_spp))
ldmccat <- model.matrix(~factor(ldmcdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
ldmc <- log(ldmcdat$m_ldmc_spp); hist(ldmc)
N_ldmc <- as.numeric(length(ldmc))

thickdat <- filter(matdifparm, !is.na(m_leafT_spp))
thickcat <- model.matrix(~factor(thickdat$slope_cat, levels=c("Zero", "Negative", "Positive")))
thick <- log(thickdat$m_leafT_spp); hist(thick)
N_thick <- as.numeric(length(ldmc))

Ndat <- filter(matdifparm, !is.na(m_Nper_spp))
Ncat <- model.matrix(~factor(Ndat$slope_cat, levels=c("Zero", "Negative",  "Positive")))
leafN <- log(Ndat$m_Nper_spp); hist(leafN)
N_N <- as.numeric(length(leafN))

Cdat <- filter(matdifparm, !is.na(m_Cper_spp))
Ccat <- model.matrix(~factor(Cdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
leafC <- log(Cdat$m_Cper_spp); hist(leafC)
N_C <- as.numeric(length(leafC))

CNdat <- filter(matdifparm, !is.na(m_CN_spp))
CNcat <- model.matrix(~factor(CNdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
leafCN <- log(CNdat$m_CN_spp); hist(leafCN)
N_CN <- as.numeric(length(leafCN))

# assemble the dat for trait model
jags.data <- list("slacat", "sla", "N_sla", "hcat", "height",
                  "N_h", "seedcat","seed", "N_seed", "ldmccat", 
                  "ldmc", "N_ldmc", "thickcat", "thick", "N_thick", "Ncat",
                  "leafN", "N_N", "Ccat","leafC", 'N_C', "CNcat","leafCN", "N_CN")
jags.params <- c("asla", "aH", "aseed",  "aLDMC", "aT", "aN", "aC", "aCN", 
                 "prec_sla", "prec_H", "prec_seed", 
                 "prec_ldmc", "prec_T", "prec_N", 
                 "prec_C", "prec_CN", "sig_sla", 
                 "sig_H", "sig_seed", "sig_ldmc", 
                 "sig_T", "sig_N",  "sig_C", "sig_CN")

# Step 2 - Trait model 
traitmod <- function(){
  for(i in 1:N_sla){
    sla[i] ~ dnorm(mu[i], prec_sla)
    mu[i] <- inprod(asla, slacat[i,])
  }
  prec_sla~dgamma(0.001, 0.001)
  sig_sla <- 1/sqrt(prec_sla)
  for(i in 1:3){asla[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_h){
    height[i] ~dnorm(mu1[i], prec_H)
    mu1[i] <- inprod(aH, hcat[i,])
  }
  prec_H~dgamma(0.001, 0.001)
  sig_H <- 1/sqrt(prec_H)
  for(i in 1:3){aH[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_seed){
    seed[i] ~dnorm(mu2[i], prec_seed)
    mu2[i] <- inprod(aseed, seedcat[i,])
  }
  prec_seed~dgamma(0.001, 0.001)
  sig_seed <- 1/sqrt(prec_seed)
  for(i in 1:3){aseed[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_ldmc){
    ldmc[i] ~dnorm(mu3[i], prec_ldmc)
    mu3[i] <- inprod(aLDMC, ldmccat[i,])
  }
  prec_ldmc~dgamma(0.001, 0.001)
  sig_ldmc <- 1/sqrt(prec_ldmc)
  for(i in 1:3){aLDMC[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_thick){
    thick[i] ~dnorm(mu4[i], prec_T)
    mu4[i] <- inprod(aT , thickcat[i,])
  }
  prec_T~dgamma(0.001, 0.001)
  sig_T <- 1/sqrt(prec_T)
  for(i in 1:3){aT[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_N){
    leafN[i] ~dnorm(mu5[i], prec_N)
    mu5[i] <- inprod(aN, Ncat[i,])
  }
  prec_N~dgamma(0.001, 0.001)
  sig_N <- 1/sqrt(prec_N)
  for(i in 1:3){aN[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_C){
    leafC[i] ~dnorm(mu6[i], prec_C)
    mu6[i] <- inprod(aC , Ccat[i,])
  }
  prec_C~dgamma(0.001, 0.001)
  sig_C <- 1/sqrt(prec_C)
  for (i in 1:3){aC[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_CN){
    leafCN[i] ~dnorm(mu7[i], prec_CN)
    mu7[i] <- inprod(aCN ,CNcat[i,])
  }
  prec_CN~dgamma(0.001, 0.001)
  sig_CN <- 1/sqrt(prec_CN)
  for(i in 1:3){aCN[i]~dnorm(0, 1E-6)}
}

# run the trait model
traitextmat <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.params,
                             n.iter=50000 ,model.file=traitmod, n.thin=5, n.chains=3)

# triat model results
traitextmat

####################################################################################
## AP differences
# group level effects
site <- factor(dat$destSiteID)
block <- factor(dat$destBlockID)
spp <- factor(dat$species)

# independent variables
ap <- dat$ap_meandiff

#dependent variable
ext <- dat$extinct
N <- as.numeric(length(ext))

# set up the data and params
jags.data <- list("site", "block","spp", "ap", "ext", "N")
jags.param <- c("a", "b",  "prec1","prec2", "asig1","asig2")

# step 1 - CDS model
apdifmod <- function(){ 
  # group effects
  for (j in 1:12){nettstedet[j]~dnorm(0, prec1)}
  for (j in 1:60){blokkere[j]~dnorm(0, prec2)}
  #likelihood
  for (i in 1:N){
    ext[i]~dbern(mu[i])
    logit(mu[i]) <- a[spp[i]]+ b[spp[i]]*ap[i] + nettstedet[site[i]]+ 
      blokkere[block[i]]
  }
  # priors
  for(j in 1:75){a[j]~dnorm(0, 1E-6)}
  for(j in 1:75){b[j]~dnorm(0, 1E-6)}
  prec1~dgamma(0.001,0.001)
  prec2~dgamma(0.001,0.001)
  asig1 <- 1/sqrt(prec1)
  asig2 <- 1/sqrt(prec2)
}

# run the model
apdif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                       n.iter=500000 ,model.file=apdifmod, n.thin=5, n.chains=3)

# see results
apdif 

# pull posteriors in list 
apdif.paramlist <- apdif$BUGSoutput$sims.list

# save posteriors
apint <- data.frame(apdif.paramlist$a)
names(apint) <- c(levels(spp))
apint <- data.frame(apply(apint, 2, sort))
write.csv(apint, "SppExtIntAP_post.csv")
apint <- apint %>% gather(key="spp", value="int")

apslop <- as.data.frame(apdif.paramlist$b)
names(apslop) <- c(levels(spp))
apslop <- data.frame(apply(apslop, 2, sort))
write.csv(apslop,"SppExtSlopeAP_post.csv" )
apslop <- apslop %>% gather(key="spp", value="slope")

apint$slope <- apslop$slope

# summarize the posteriors
apdifparm <- apint %>% group_by(spp) %>%
  summarize(m_int=median(int),
            max_int=max(int), 
            min_int=min(int), 
            sig_int=sd(int),
            m_slope=median(slope),
            max_slope=max(slope),
            min_slope=min(slope),
            sig_slope=sd(slope))

spprangemat <- dat %>% group_by(species) %>%
  summarize(minap = min(ap_meandiff),
            maxap = max(ap_meandiff),
            meanap = mean(ap_meandiff))%>%
  rename(spp=species)

# combine the sets
apdifparm <- full_join(apdifparm, select(spprangemat, meanap, spp),by="spp", copy=TRUE)

# add in the trait data
traitformod <-  dat %>% select(species, m_height_spp:seed_mass) %>% unique()
apdifparm <- left_join(apdifparm, traitformod, by=c("spp"="species"))

# remove species that didn't converge
apdifparm <- apdifparm %>% 
  filter(!spp %in% c("Car.lep","Car.pan", "Gal.ver","Luz.pil", "Mel.pra", "Noc.cae"))

# Classify CDSs
apdifparm <- apdifparm %>% 
  mutate(slope_cat = if_else(m_slope/sig_slope <  -1, "Negative",
                             if_else(m_slope/sig_slope < 1 & m_slope/sig_slope > -1,
                                     "Zero", "Positive")))

# Construct models for trait predictions
# make the different datasets 
sladat <- filter(apdifparm, !is.na(m_sla_spp))
slacat <- model.matrix(~factor(sladat$slope_cat, levels=c("Zero","Negative",  "Positive")))
sla <- log(sladat$m_sla_spp);hist(sla)
N_sla <- as.numeric(length(sla))

heightdat <- filter(apdifparm, !is.na(m_height_spp))
hcat <- model.matrix(~factor(heightdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
height <- log(heightdat$m_height_spp); hist(height)
N_h <- as.numeric(length(height))

seeddat <- apdifparm %>% filter(!is.na(seed_mass)) 
seedcat <- model.matrix(~factor(seeddat$slope_cat, levels=c("Zero","Negative",  "Positive")))
seed <- log(seeddat$seed_mass); hist(seed)
N_seed <- as.numeric(length(seed))

ldmcdat <- filter(apdifparm, !is.na(m_ldmc_spp))
ldmccat <- model.matrix(~factor(ldmcdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
ldmc <- log(ldmcdat$m_ldmc_spp); hist(ldmc)
N_ldmc <- as.numeric(length(ldmc))

thickdat <- filter(apdifparm, !is.na(m_leafT_spp))
thickcat <- model.matrix(~factor(thickdat$slope_cat, levels=c("Zero", "Negative", "Positive")))
thick <- log(thickdat$m_leafT_spp); hist(thick)
N_thick <- as.numeric(length(ldmc))

Ndat <- filter(apdifparm, !is.na(m_Nper_spp))
Ncat <- model.matrix(~factor(Ndat$slope_cat, levels=c("Zero", "Negative",  "Positive")))
leafN <- log(Ndat$m_Nper_spp); hist(leafN)
N_N <- as.numeric(length(leafN))

Cdat <- filter(apdifparm, !is.na(m_Cper_spp))
Ccat <- model.matrix(~factor(Cdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
leafC <- log(Cdat$m_Cper_spp); hist(leafC)
N_C <- as.numeric(length(leafC))

CNdat <- filter(apdifparm, !is.na(m_CN_spp))
CNcat <- model.matrix(~factor(CNdat$slope_cat, levels=c("Zero","Negative",  "Positive")))
leafCN <- log(CNdat$m_CN_spp); hist(leafCN)
N_CN <- as.numeric(length(leafCN))

jags.data <- list("slacat", "sla", "N_sla", "hcat", "height",
                  "N_h", "seedcat","seed", "N_seed", "ldmccat", 
                  "ldmc", "N_ldmc", "thickcat", "thick", "N_thick", "Ncat",
                  "leafN", "N_N", "Ccat","leafC", 'N_C', "CNcat","leafCN", "N_CN")
jags.params <- c("asla", "aH", "aseed",  "aLDMC", "aT", "aN", "aC", "aCN", 
                 "prec_sla", "prec_H", "prec_seed", 
                 "prec_ldmc", "prec_T", "prec_N", 
                 "prec_C", "prec_CN", "sig_sla", 
                 "sig_H", "sig_seed", "sig_ldmc", 
                 "sig_T", "sig_N",  "sig_C", "sig_CN")

# step 2 - trait model
traitmod <- function(){
  for(i in 1:N_sla){
    sla[i] ~ dnorm(mu[i], prec_sla)
    mu[i] <- inprod(asla, slacat[i,])
  }
  prec_sla~dgamma(0.001, 0.001)
  sig_sla <- 1/sqrt(prec_sla)
  for(i in 1:3){asla[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_h){
    height[i] ~dnorm(mu1[i], prec_H)
    mu1[i] <- inprod(aH, hcat[i,])
  }
  prec_H~dgamma(0.001, 0.001)
  sig_H <- 1/sqrt(prec_H)
  for(i in 1:3){aH[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_seed){
    seed[i] ~dnorm(mu2[i], prec_seed)
    mu2[i] <- inprod(aseed, seedcat[i,])
  }
  prec_seed~dgamma(0.001, 0.001)
  sig_seed <- 1/sqrt(prec_seed)
  for(i in 1:3){aseed[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_ldmc){
    ldmc[i] ~dnorm(mu3[i], prec_ldmc)
    mu3[i] <- inprod(aLDMC, ldmccat[i,])
  }
  prec_ldmc~dgamma(0.001, 0.001)
  sig_ldmc <- 1/sqrt(prec_ldmc)
  for(i in 1:3){aLDMC[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_thick){
    thick[i] ~dnorm(mu4[i], prec_T)
    mu4[i] <- inprod(aT , thickcat[i,])
  }
  prec_T~dgamma(0.001, 0.001)
  sig_T <- 1/sqrt(prec_T)
  for(i in 1:3){aT[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_N){
    leafN[i] ~dnorm(mu5[i], prec_N)
    mu5[i] <- inprod(aN, Ncat[i,])
  }
  prec_N~dgamma(0.001, 0.001)
  sig_N <- 1/sqrt(prec_N)
  for(i in 1:3){aN[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_C){
    leafC[i] ~dnorm(mu6[i], prec_C)
    mu6[i] <- inprod(aC , Ccat[i,])
  }
  prec_C~dgamma(0.001, 0.001)
  sig_C <- 1/sqrt(prec_C)
  for (i in 1:3){aC[i]~dnorm(0, 1E-6)}
  
  for(i in 1:N_CN){
    leafCN[i] ~dnorm(mu7[i], prec_CN)
    mu7[i] <- inprod(aCN ,CNcat[i,])
  }
  prec_CN~dgamma(0.001, 0.001)
  sig_CN <- 1/sqrt(prec_CN)
  for(i in 1:3){aCN[i]~dnorm(0, 1E-6)}
}

# run the model
traitextap <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.params,
                            n.iter=50000 ,model.file=traitmod, n.thin=5, n.chains=3)

# see results
traitextap

