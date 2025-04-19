#Analysis of pombe clone data

#library(dplyr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(brms)
library(coda)
library(latex2exp)

### * Load the data

filtered <- read.csv("ODdata_filtered.csv", header = T, sep = ",", dec = ".")

#Kun kasvu on "flat" korvaa k = 0.01 ja r = 0.01

filtered$k <- replace(filtered$k, filtered$Kasvu == "flat", 0.01)
filtered$r <- replace(filtered$r, filtered$Kasvu == "flat", 0.01)

#Kun kasvu on "stat" korvaa K = maxOD

for(i in 1:nrow(filtered)) {
  if(filtered$Kasvu[i] == "stat") { filtered$k[i] <- filtered$maxOD[i] }
}

### * Analyysi alkaa tästä

### ** Piirretään kuvia kloonien kasvusta
#aineisto <- filtered

#Voi tarkastaa jakaumaa myös histogrammeilla
hist(filtered$r)

hist(filtered$k)

#Tutkittavia kysymyksiä:

#Miten reaktionormit ovat kehittyneet evoluution aikana?

#Tehdään ensin kuvia reaktionormeista keskiarvoilla

#jaetaan aineisto ancestorien mukaan, koska verrataan omaan klooniin

ancestors <- filter(filtered, Population.ID == "anc")

anc.means <- summarise(group_by(ancestors, Temperature, genotype), mean.r = mean(r), mean.K = mean(k))


#Evolved strains

#Extreme

evolved <- filter(filtered, Population.ID != "anc")

evolved.extreme <- filter(evolved, Evolution.treatment == "Extreme")

evol.extreme.clonemeans <- summarise(group_by(evolved.extreme, Temperature, Clone.ID, genotype), mean.r = mean(r), mean.K = mean(k))

ggplot(evol.extreme.clonemeans, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#Yleisesti parempi kasvunopeus kuin ancestorilla 38C

ggplot(evol.extreme.clonemeans, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#saman suuntaisia viitteitä K:lle

#Mean

evolved.mean <- filter(evolved, Evolution.treatment == "Mean")

evol.mean.clonemeans <- summarise(group_by(evolved.mean, Temperature, Clone.ID, genotype), mean.r = mean(r), mean.K = mean(k))

ggplot(evol.mean.clonemeans, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#34C kehittyminen huonontaa selkeiten äärilämmön sietoa

ggplot(evol.mean.clonemeans, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#Nähtävissä selkeämmin K:n arvoista

#Fast

evolved.fast <- filter(evolved, Evolution.treatment == "Fast")

evol.fast.clonemeans <- summarise(group_by(evolved.fast, Temperature, Clone.ID, genotype), mean.r = mean(r), mean.K = mean(k))

ggplot(evol.fast.clonemeans, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#nopea vaihtelu on valinnut nopeampaa kasvua optimaalisemmissa lämpötiloissa, tiputus 38C

ggplot(evol.fast.clonemeans, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#Kasvun määrä näyttää jäävän kaikissa käsittelyissä alle ancestorin

#Intermediate

evolved.intermediate <- filter(evolved, Evolution.treatment == "Intermediate")

evol.intermediate.clonemeans <- summarise(group_by(evolved.intermediate, Temperature, Clone.ID, genotype), mean.r = mean(r), mean.K = mean(k))

ggplot(evol.intermediate.clonemeans, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#Keskinopeassa vaihtelussa kasvunopeus yleisesti parantunut paitsi 38C

ggplot(evol.intermediate.clonemeans, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#Kssa selkeä tiputus 38C

#Slow

evolved.slow <- filter(evolved, Evolution.treatment == "Slow")

evol.slow.clonemeans <- summarise(group_by(evolved.slow, Temperature, Clone.ID, genotype), mean.r = mean(r), mean.K = mean(k))

ggplot(evol.slow.clonemeans, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#hitaassa vaihtelussa kasvu parantunut optimaalisissa lämpötiloissa, 38C enemmän vaihtelua

ggplot(evol.slow.clonemeans, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
  geom_point() +
  geom_line() +
  geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
  facet_wrap(~ genotype)

#Yllättäen kasvun määrä on pienentynyt melkein joka käsittelyssä, mutta ancestorien välillä on eroja


#Tehdään kuvat joissa kaikki evoluutiokäsittelyt yhdessä
evolved.means.all <- data.frame(rbind(evol.extreme.clonemeans, evol.mean.clonemeans, evol.intermediate.clonemeans, evol.fast.clonemeans, evol.slow.clonemeans), Treatment = factor( c( rep("Extreme", nrow(evol.extreme.clonemeans) ), rep("Mean", nrow(evol.mean.clonemeans)), rep("Intermediate", nrow(evol.intermediate.clonemeans)), rep("Fast", nrow(evol.fast.clonemeans)), rep("Slow", nrow(evol.slow.clonemeans)) )) )

#Save clone means for fast loading
#save(evolved.means.all, anc.means, file = "./pombe_evo/data/clonemeans.RData")


#Plot all together
plot.clones.r <- ggplot(evolved.means.all, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
    geom_point() +
    geom_line() +
    geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
    ylab("Maximum growth rate (r)") +
    xlab("Assay environment") +    
    facet_grid(Treatment ~ genotype)

plot.clones.K <- ggplot(evolved.means.all, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
    geom_point() +
    geom_line() +
    geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
    ylab("Carrying capacity (K)") +
    xlab("Assay environment") +    
    facet_grid(Treatment ~ genotype)

save_plot(filename = "./pombe_evo/fig/clones_r.pdf", plot.clones.r, base_height = 10, base_width=12)
save_plot(filename = "./pombe_evo/fig/clones_K.pdf", plot.clones.K, base_height = 10, base_width=12)

###################################################################################################
### ** Kloonien kasvu suhteessa ancestoriin eri lämpötiloissa. Onko trade-offeja havaittavissa? ###
###################################################################################################

evolved <- filter(filtered, Population.ID != "anc")

ancestors <- filter(filtered, Population.ID == "anc")

#Tähän käy samat mallit kuin variansseihin ja ne on laskettu valmiiksi.
load("~/modelfits/var.Mean.RData")
load("~/modelfits/var.Ext.RData")
load("~/modelfits/var.Fast.RData")
load("~/modelfits/var.Int.RData")
load("~/modelfits/var.Slow.RData")
load("~/modelfits/var.anc.RData")

##Tarvitsee ennustaa intercepti sekä populaatioille, että ancestorille pregrowth OD:ssä 1.3
#(Pregrowth OD = 1.3 on lähellä molempien keskiarvoa (30 astetta)

var.anc.30[[3]]
var.anc.30[[3]]

koe <- (var.Mean.30[[2]][,1]+10)  / (var.anc.30[[2]][,1]+10)

koe2 <- (var.Mean.30[[2]][,1])  / (var.anc.30[[2]][,1])

#Need to add intercept + 1.3*pregrowth
r.anc.30 <- var.anc.30[[2]][,1] + 1.0*var.anc.30[[2]][,2]
r.Mean.30 <- var.Mean.30[[2]][,1] + 1.0*var.Mean.30[[2]][,2]

#Loadin the variance calculations as needed
load("~/modelfits/var.Mean.K.RData")
load("~/modelfits/var.Ext.K.RData")
load("~/modelfits/var.Fast.K.RData")
load("~/modelfits/var.Int.K.RData")
load("~/modelfits/var.Slow.K.RData")
load("~/modelfits/var.anc.K.RData")

#Malli rmaxille

#This works! sensible growth rates
malli.r.anc <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD,
                   data = ancestors, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.anc <- posterior_samples(malli.r.anc)[,-c(1,18,19)]
head(post.r.anc)

#(1 || Population.ID/Clone.ID),
malli.r.Mean <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Mean"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.Mean <- posterior_samples(malli.r.Mean)[,2:17]

#Results for Mean evol treatment
res.r.Mean <- data.frame(estimate = apply(post.r.Mean/post.r.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.r.Mean/post.r.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Mean", 16))

ggplot(res.r.Mean, aes(x = assay, y = estimate, ymin = lower, ymax = upper, group = genotype)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    xlab("Assay environment") +
    ylab("Relative growth rate") +
    geom_hline(yintercept = 1, lty = "dashed")

### ** Varianssimallit

#Miten paljon populaatiossa on muuntelua? -> Varianssimalli brms
#Esim. ylläpitääkö vaihteleva ympäristö enemmän geneettistä muuntelua?

#miten pärjäävät lämpötilavaihtelussa?
#Onko kyse muuntelusta ja plastisuudesta

#Populaatioiden välinen ja sisäinen muuntelu, myös kloonin sisäinen

ancestors <- filter(filtered, Population.ID == "anc")

evolved <- filter(filtered, Population.ID != "anc")


#
#ancestors.30 <- filter(ancestors, Temperature == "30")
#Centering variables
#ancestors.30$r <- scale(as.numeric(ancestors.30$r, center = T, scale = F))
#ancestors.30$k <- scale(as.numeric(ancestors.30$k, center = T, scale = F))
#ancestors.30$Pregrowth.OD <- scale(as.numeric(ancestors.30$Pregrowth.OD, center = T, scale = F))

#Running the model
#anc.30.m1 <- brm(r ~ Pregrowth.OD + genotype,
#                 data = ancestors.30, family = gaussian(), iter = 4000, cores = 4)

calc.variances.ancestors <- function(ancestors, assay.temp, trait) {
    #### Filter the correct data and scale it to mean 0
    ancestors.filt <- filter(ancestors, Temperature == assay.temp)
    ancestors.filt$r <- scale(as.numeric(ancestors.filt$r, center = T, scale = F))
    ancestors.filt$k <- scale(as.numeric(ancestors.filt$k, center = T, scale = F))
    ancestors.filt$Pregrowth.OD <- scale(as.numeric(ancestors.filt$Pregrowth.OD, center = T, scale = F))
    #Fit the model
    #Check for which variances are to be calculated
    if(trait == "r") {
    anc.m1 <- brm(r ~ Pregrowth.OD + genotype,
                 data = ancestors.filt, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
    }

    if(trait == "k") {
    anc.m1 <- brm(k ~ Pregrowth.OD + genotype,
                 data = ancestors.filt, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
    }

    #Summary of model results
    yhteenveto <- summary(anc.m1)

    #Extract posteriors
    post <- posterior_samples(anc.m1)

    #Make the results matrix
    res.mat <- data.frame(matrix(rep(0, 6*1), ncol = 6))
    colnames(res.mat) <- c("type", "treatment", "assay", "estimate", "lower", "upper")
    res.mat$type  <- c("environmental")
    res.mat$treatment <- rep("ancestor",1)
    res.mat$assay <- rep(assay.temp, 1)

    #Store the results
    res.mat[1,4:6] <- c(median(post$sigma^2), HPDinterval(as.mcmc(post$sigma^2)))

    #Return the results
    mylist <- list(res.mat, post, yhteenveto)
    return(mylist)
}

###Calculate environmental variances for the ancestors
var.anc.30 <- calc.variances.ancestors(ancestors, assay.temp = "30", trait = "r")
var.anc.34 <- calc.variances.ancestors(ancestors, assay.temp = "34", trait = "r")
var.anc.38 <- calc.variances.ancestors(ancestors, assay.temp = "38", trait = "r")
var.anc.fast <- calc.variances.ancestors(ancestors, assay.temp = "Fast 30-38", trait = "r")

save(var.anc.30, var.anc.34, var.anc.38, var.anc.fast, file = "~/modelfits/var.anc.RData")

#For K
var.anc.30.K <- calc.variances.ancestors(ancestors, assay.temp = "30", trait = "k")
var.anc.34.K <- calc.variances.ancestors(ancestors, assay.temp = "34", trait = "k")
var.anc.38.K <- calc.variances.ancestors(ancestors, assay.temp = "38", trait = "k")
var.anc.fast.K <- calc.variances.ancestors(ancestors, assay.temp = "Fast 30-38", trait = "k")

save(var.anc.30.K, var.anc.34.K, var.anc.38.K, var.anc.fast.K, file = "~/modelfits/var.anc.K.RData")
####



######################################################################
### For one assay temperature and one evolution treatment
### Change this so this to a function that does this

calc.variances.evolved <- function(evolved, assay.temp, evo.treatment, trait) {
    ### Filter the correct data and scale it to mean 0
    evolved.filt <- filter(evolved, Temperature == assay.temp, Evolution.treatment == evo.treatment)
    evolved.filt$r <- scale(as.numeric(evolved.filt$r, center = T, scale = F))
    evolved.filt$k <- scale(as.numeric(evolved.filt$k, center = T, scale = F))
    evolved.filt$Pregrowth.OD <- scale(as.numeric(evolved.filt$Pregrowth.OD, center = T, scale = F))

    #Fit the model
    #Check for which variances are to be calculated
    if(trait == "r") {
    evo.m1 <- brm(r ~ Pregrowth.OD + genotype + (1 || Population.ID/Clone.ID),
                 data = evolved.filt, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
    }

    if(trait == "k") {
    evo.m1 <- brm(k ~ Pregrowth.OD + genotype + (1 || Population.ID/Clone.ID),
                 data = evolved.filt, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
    }        

    #Summary of model results
    yhteenveto <- summary(evo.m1)

    #Extract posteriors
    post <- posterior_samples(evo.m1)

    #Make the results matrix
    res.mat <- data.frame(matrix(rep(0, 6*3), ncol = 6))
    colnames(res.mat) <- c("type", "treatment", "assay", "estimate", "lower", "upper")
    res.mat$type  <- c("population", "clone", "environmental")
    res.mat$treatment <- rep(evo.treatment,3)
    res.mat$assay <- rep(assay.temp, 3)

    #Store the results
    res.mat[1,4:6] <- c(median(post$sd_Population.ID__Intercept^2), HPDinterval(as.mcmc(post$sd_Population.ID__Intercept^2)))
    res.mat[2,4:6] <- c(median(post$'sd_Population.ID:Clone.ID__Intercept'^2), HPDinterval(as.mcmc(post$'sd_Population.ID:Clone.ID__Intercept'^2)))
    res.mat[3,4:6] <- c(median(post$sigma^2), HPDinterval(as.mcmc(post$sigma^2)))

    #Return the results
    mylist <- list(res.mat, post, yhteenveto)
    return(mylist)
}

#Fit variance calculation models for different assay and evolution treatment combinations
### Evolved in Mean ###
var.Mean.30 <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Mean", trait = "r")
var.Mean.34 <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Mean", trait = "r")
var.Mean.38 <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Mean", trait = "r")
var.Mean.fast <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Mean", trait = "r")

#Save the results so no need to run the models again
save(var.Mean.30, var.Mean.34, var.Mean.38, var.Mean.fast, file = "~/modelfits/var.Mean.RData")

#Run for K
var.Mean.30.K <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Mean", trait = "k")
var.Mean.34.K <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Mean", trait = "k")
var.Mean.38.K <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Mean", trait = "k")
var.Mean.fast.K <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Mean", trait = "k")

#Save
save(var.Mean.30.K, var.Mean.34.K, var.Mean.38.K, var.Mean.fast.K, file = "~/modelfits/var.Mean.K.RData")

###########

### Evolved in extreme ####
var.Ext.30 <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Extreme", trait = "r")
var.Ext.34 <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Extreme", trait = "r")
var.Ext.38 <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Extreme", trait = "r")
var.Ext.fast <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Extreme", trait = "r")

save(var.Ext.30, var.Ext.34, var.Ext.38, var.Ext.fast, file = "~/modelfits/var.Ext.RData")

#For K
var.Ext.30.K <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Extreme", trait = "k")
var.Ext.34.K <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Extreme", trait = "k")
var.Ext.38.K <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Extreme", trait = "k")
var.Ext.fast.K <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Extreme", trait = "k")

save(var.Ext.30.K, var.Ext.34.K, var.Ext.38.K, var.Ext.fast.K, file = "~/modelfits/var.Ext.K.RData")

############################

### Evolved in fast #####
var.Fast.30 <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Fast", trait = "r")
var.Fast.34 <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Fast", trait = "r")
var.Fast.38 <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Fast", trait = "r")
var.Fast.fast <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Fast", trait = "r")

save(var.Fast.30, var.Fast.34, var.Fast.38, var.Fast.fast, file = "~/modelfits/var.Fast.RData")

#For K
var.Fast.30.K <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Fast", trait = "k")
var.Fast.34.K <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Fast", trait = "k")
var.Fast.38.K <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Fast", trait = "k")
var.Fast.fast.K <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Fast", trait = "k")

save(var.Fast.30.K, var.Fast.34.K, var.Fast.38.K, var.Fast.fast.K, file = "~/modelfits/var.Fast.K.RData")

##############################

### Evolved in slow ####
var.Slow.30 <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Slow", trait = "r")
var.Slow.34 <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Slow", trait = "r")
var.Slow.38 <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Slow", trait = "r")
var.Slow.fast <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Slow", trait = "r")

save(var.Slow.30, var.Slow.34, var.Slow.38, var.Slow.fast, file = "~/modelfits/var.Slow.RData")

#For K
var.Slow.30.K <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Slow", trait = "k")
var.Slow.34.K <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Slow", trait = "k")
var.Slow.38.K <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Slow", trait = "k")
var.Slow.fast.K <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Slow", trait = "k")

save(var.Slow.30.K, var.Slow.34.K, var.Slow.38.K, var.Slow.fast.K, file = "~/modelfits/var.Slow.K.RData")

##########################


### Evolved in intermediate ####
var.Int.30 <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Intermediate", trait = "r")
var.Int.34 <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Intermediate", trait = "r")
var.Int.38 <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Intermediate", trait = "r")
var.Int.fast <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Intermediate", trait = "r")

save(var.Int.30, var.Int.34, var.Int.38, var.Int.fast, file = "~/modelfits/var.Int.RData")

#For K
var.Int.30.K <- calc.variances.evolved(evolved, assay.temp = "30", evo.treatment = "Intermediate", trait = "k")
var.Int.34.K <- calc.variances.evolved(evolved, assay.temp = "34", evo.treatment = "Intermediate", trait = "k")
var.Int.38.K <- calc.variances.evolved(evolved, assay.temp = "38", evo.treatment = "Intermediate", trait = "k")
var.Int.fast.K <- calc.variances.evolved(evolved, assay.temp = "Fast 30-38", evo.treatment = "Intermediate", trait = "k")

save(var.Int.30.K, var.Int.34.K, var.Int.38.K, var.Int.fast.K, file = "~/modelfits/var.Int.K.RData")

################################

#Loadin the variance calculations as needed
load("~/modelfits/var.Mean.RData")
load("~/modelfits/var.Ext.RData")
load("~/modelfits/var.Fast.RData")
load("~/modelfits/var.Int.RData")
load("~/modelfits/var.Slow.RData")
load("~/modelfits/var.anc.RData")

#Make a dataframe of the results
var.results <- rbind(var.Mean.30[[1]], var.Mean.34[[1]], var.Mean.38[[1]], var.Mean.fast[[1]], var.Ext.30[[1]], var.Ext.34[[1]], var.Ext.38[[1]], var.Ext.fast[[1]], var.Fast.30[[1]], var.Fast.34[[1]], var.Fast.38[[1]], var.Fast.fast[[1]], var.Int.30[[1]], var.Int.34[[1]], var.Int.38[[1]], var.Int.fast[[1]], var.Slow.30[[1]], var.Slow.34[[1]], var.Slow.38[[1]], var.Slow.fast[[1]], var.anc.30[[1]], var.anc.34[[1]], var.anc.38[[1]], var.anc.fast[[1]] )


ggplot(var.results, aes(x = assay, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    ylab("Variance") +    
    facet_grid(treatment ~ type)


#Colors
#RColorBrewer::brewer.pal(5, "Set1")

variances.r.plot <- ggplot(var.results, aes(x = assay, y = estimate, ymin = lower, ymax = upper, colour = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    scale_colour_manual(name = "Evolutionary\ntreatment", values = c(Extreme = "#E41A1C", Fast = "#377EB8", Intermediate = "#4DAF4A", Mean = "#984EA3", Slow = "#FF7F00", ancestor = "black")) +
    ylab("Variance") +
    xlab("Assay environment") +    
    facet_grid( ~ type)

save_plot(filename = "./pombe_evo/fig/variances_r.pdf", variances.r.plot, base_height = 3, base_width = 12)

variances.r.plot2 <- ggplot(var.results, aes(x = treatment, y = estimate, ymin = lower, ymax = upper, colour = assay)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    #scale_colour_manual(name = "Assay\nenvironment") +
    ylab("Variance") +
    xlab("Evolutionary treatment") +
    labs(colour = "Assay\nenvironment") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    #theme(text=element_text(size=16)) +
    facet_grid( ~ type)

save_plot(filename = "./pombe_evo/fig/variances_r2.pdf", variances.r.plot2, base_height = 3, base_width = 12)

#Do the analysis for K

load("~/modelfits/var.Mean.K.RData")
load("~/modelfits/var.Ext.K.RData")
load("~/modelfits/var.Fast.K.RData")
load("~/modelfits/var.Int.K.RData")
load("~/modelfits/var.Slow.K.RData")
load("~/modelfits/var.anc.K.RData")


#Make a dataframe of the results
var.results.K <- rbind(var.Mean.30.K[[1]], var.Mean.34.K[[1]], var.Mean.38.K[[1]], var.Mean.fast.K[[1]], var.Ext.30.K[[1]], var.Ext.34.K[[1]], var.Ext.38.K[[1]], var.Ext.fast.K[[1]], var.Fast.30.K[[1]], var.Fast.34.K[[1]], var.Fast.38.K[[1]], var.Fast.fast.K[[1]], var.Int.30.K[[1]], var.Int.34.K[[1]], var.Int.38.K[[1]], var.Int.fast.K[[1]], var.Slow.30.K[[1]], var.Slow.34.K[[1]], var.Slow.38.K[[1]], var.Slow.fast.K[[1]], var.anc.30.K[[1]], var.anc.34.K[[1]], var.anc.38.K[[1]], var.anc.fast.K[[1]] )


variances.k.plot <- ggplot(var.results.K, aes(x = assay, y = estimate, ymin = lower, ymax = upper, colour = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    scale_colour_manual(name = "Evolutionary\ntreatment", values = c(Extreme = "#E41A1C", Fast = "#377EB8", Intermediate = "#4DAF4A", Mean = "#984EA3", Slow = "#FF7F00", ancestor = "black")) +
    ylab("Variance") +
    xlab("Assay environment") +    
    facet_grid( ~ type)

save_plot(filename = "./pombe_evo/fig/variances_k.pdf", variances.k.plot, base_height = 3, base_width = 12)


variances.k.plot2 <- ggplot(var.results.K, aes(x = treatment, y = estimate, ymin = lower, ymax = upper, colour = assay)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    ylab("Variance") +
    xlab("Evolutionary treatment") +
    labs(colour = "Assay environment") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +    
    facet_grid( ~ type)

save_plot(filename = "./pombe_evo/fig/variances_k2.pdf", variances.k.plot2, base_height = 3, base_width = 12)

#legend <- get_legend(variances.k.plot2)

legend_b <- get_legend(variances.k.plot2 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top", legend.justification = "center")
)


variances.r.plot2 <- variances.r.plot2 + theme(legend.position = "none")
variances.k.plot2 <- variances.k.plot2 + theme(legend.position = "none")

var.final.plot <- plot_grid(legend_b, variances.r.plot2, variances.k.plot2, labels = c("", "A", "B"), ncol = 1, rel_heights = c(.1,1,1))

save_plot(filename = "./pombe_evo/fig/variances_both.pdf", var.final.plot, base_height = 6.3, base_width = 12)

plot_grid(variances.r.plot2, variances.k.plot2, nrow = 2, labels = c("A", "B"))

### Posterior comparisons

#We need to compare variances over assay environments and evolutionary treatments
#Are there some environments that lead to canalization or some assay environments

#Make a dataframe, use model with measurement error
#clone and environmental variances, population variances have so large errors that they are not that interesting

treatment <- c( rep("Extreme", 4), rep("Fast", 4), rep("Intermediate", 4), rep("Mean", 4), rep("Slow", 4))
assay <- c( rep(c("30", "34", "38", "Fast30-38"), 5))
treatment2 <- c( rep("Extreme", 4), rep("Fast", 4), rep("Intermediate", 4), rep("Mean", 4), rep("Slow", 4), rep("Ancestor", 4))
assay2 <- c( rep(c("30", "34", "38", "Fast30-38"), 6))

#Clone variance, growth rate
clone.r.vardata <- data.frame(treatment = treatment , assay = assay, estimate = c( var.Ext.30[[1]][2,4], var.Ext.34[[1]][2,4], var.Ext.38[[1]][2,4], var.Ext.fast[[1]][2,4], var.Fast.30[[1]][2,4], var.Fast.34[[1]][2,4], var.Fast.38[[1]][2,4], var.Fast.fast[[1]][2,4], var.Int.30[[1]][2,4], var.Int.34[[1]][2,4], var.Int.38[[1]][2,4], var.Int.fast[[1]][2,4], var.Mean.30[[1]][2,4], var.Mean.34[[1]][2,4], var.Mean.38[[1]][2,4], var.Mean.fast[[1]][2,4], var.Slow.30[[1]][2,4], var.Slow.34[[1]][2,4], var.Slow.38[[1]][2,4], var.Slow.fast[[1]][2,4]), var.er = c( sd(var.Ext.30[[2]][,7]^2), sd(var.Ext.34[[2]][,7]^2), sd(var.Ext.38[[2]][,7]^2), sd(var.Ext.fast[[2]][,7]^2), sd(var.Fast.30[[2]][,7]^2), sd(var.Fast.34[[2]][,7]^2), sd(var.Fast.38[[2]][,7]^2), sd(var.Fast.fast[[2]][,7]^2), sd(var.Int.30[[2]][,7]^2), sd(var.Int.34[[2]][,7]^2), sd(var.Int.38[[2]][,7]^2), sd(var.Int.fast[[2]][,7]^2), sd(var.Mean.30[[2]][,7]^2), sd(var.Mean.34[[2]][,7]^2), sd(var.Mean.38[[2]][,7]^2), sd(var.Mean.fast[[2]][,7]^2), sd(var.Slow.30[[2]][,7]^2), sd(var.Slow.34[[2]][,7]^2), sd(var.Slow.38[[2]][,7]^2), sd(var.Slow.fast[[2]][,7]^2)  )   )

ilist <- list(estimate = clone.r.vardata$estimate)
inits_list <- list(ilist, ilist, ilist, ilist) #Needs to be a list of lists of size nchains

#Using a measurement error model for the variances
var.r.clone.model1 <- brm(data = clone.r.vardata, family = gaussian,
                          estimate | se(var.er, sigma = TRUE) ~ 0 + intercept + treatment + assay,
                          iter = 5000, warmup = 1000, cores = 4, chains = 4,
                          inits = inits_list)

pp_check(var.r.clone.model1, nsamples = 25) #Looks OK

#For treatment, Extreme is the intercept, and for assay it is 30 C
r.clone.results.treatment <- data.frame(fixef(var.r.clone.model1))[-c(1,6:8),] #Intercept can be dropped
r.clone.results.assay <- data.frame(fixef(var.r.clone.model1))[-c(1:5),] #Intercept can be dropped
r.clone.results.treatment$treatment <- c("Fast", "Intermediate", "Mean", "Slow")
r.clone.results.assay$assay <- c("34", "38", "Fast 30-38")


#Environmental variance, growth rate
env.r.vardata <- data.frame(treatment = treatment2, assay = assay2, estimate = c( var.Ext.30[[1]][3,4], var.Ext.34[[1]][3,4], var.Ext.38[[1]][3,4], var.Ext.fast[[1]][3,4], var.Fast.30[[1]][3,4], var.Fast.34[[1]][3,4], var.Fast.38[[1]][3,4], var.Fast.fast[[1]][3,4], var.Int.30[[1]][3,4], var.Int.34[[1]][3,4], var.Int.38[[1]][3,4], var.Int.fast[[1]][3,4], var.Mean.30[[1]][3,4], var.Mean.34[[1]][3,4], var.Mean.38[[1]][3,4], var.Mean.fast[[1]][3,4], var.Slow.30[[1]][3,4], var.Slow.34[[1]][3,4], var.Slow.38[[1]][3,4], var.Slow.fast[[1]][3,4], var.anc.30[[1]][1,4], var.anc.34[[1]][1,4], var.anc.38[[1]][1,4], var.anc.fast[[1]][1,4]), var.er = c( sd(var.Ext.30[[2]][,8]^2), sd(var.Ext.34[[2]][,8]^2), sd(var.Ext.38[[2]][,8]^2), sd(var.Ext.fast[[2]][,8]^2), sd(var.Fast.30[[2]][,8]^2), sd(var.Fast.34[[2]][,8]^2), sd(var.Fast.38[[2]][,8]^2), sd(var.Fast.fast[[2]][,8]^2), sd(var.Int.30[[2]][,8]^2), sd(var.Int.34[[2]][,8]^2), sd(var.Int.38[[2]][,8]^2), sd(var.Int.fast[[2]][,8]^2), sd(var.Mean.30[[2]][,8]^2), sd(var.Mean.34[[2]][,8]^2), sd(var.Mean.38[[2]][,8]^2), sd(var.Mean.fast[[2]][,8]^2), sd(var.Slow.30[[2]][,8]^2), sd(var.Slow.34[[2]][,8]^2), sd(var.Slow.38[[2]][,8]^2), sd(var.Slow.fast[[2]][,8]^2), sd(var.anc.30[[2]][,6]^2), sd(var.anc.34[[2]][,6]^2), sd(var.anc.38[[2]][,6]^2), sd(var.anc.fast[[2]][,6]^2) ) )

ilist <- list(estimate = env.r.vardata$estimate)
inits_list <- list(ilist, ilist, ilist, ilist) #Needs to be a list of lists of size nchains

#Using a measurement error model for the variances
var.r.env.model1 <- brm(data = env.r.vardata, family = gaussian,
                          estimate | se(var.er, sigma = TRUE) ~ 0 + intercept + treatment + assay,
                          iter = 5000, warmup = 1000, cores = 4, chains = 4,
                          inits = inits_list)

pp_check(var.r.env.model1, nsamples = 25) #Looks OK

#For treatment, Extreme is the intercept, and for assay it is 30 C
r.env.results.treatment <- data.frame(fixef(var.r.env.model1))[-c(1,7:9),] #Intercept can be dropped
r.env.results.assay <- data.frame(fixef(var.r.env.model1))[-c(1:6),] #Intercept can be dropped
r.env.results.treatment$treatment <- c("Extreme", "Fast", "Intermediate", "Mean", "Slow")
r.env.results.assay$assay <- c("34", "38", "Fast 30-38")


#### Saving the results so models don't need to be run again
save(r.clone.results.treatment, r.clone.results.assay, r.env.results.treatment, r.env.results.assay, file = "./pombe_evo/data/var_results_r.RData")

###Models for variance for carrying capacity, K

#Clone variance, growth rate
clone.K.vardata <- data.frame(treatment = treatment , assay = assay, estimate = c( var.Ext.30.K[[1]][2,4], var.Ext.34.K[[1]][2,4], var.Ext.38.K[[1]][2,4], var.Ext.fast.K[[1]][2,4], var.Fast.30.K[[1]][2,4], var.Fast.34.K[[1]][2,4], var.Fast.38.K[[1]][2,4], var.Fast.fast.K[[1]][2,4], var.Int.30.K[[1]][2,4], var.Int.34.K[[1]][2,4], var.Int.38.K[[1]][2,4], var.Int.fast.K[[1]][2,4], var.Mean.30.K[[1]][2,4], var.Mean.34.K[[1]][2,4], var.Mean.38.K[[1]][2,4], var.Mean.fast.K[[1]][2,4], var.Slow.30.K[[1]][2,4], var.Slow.34.K[[1]][2,4], var.Slow.38.K[[1]][2,4], var.Slow.fast.K[[1]][2,4]), var.er = c( sd(var.Ext.30.K[[2]][,7]^2), sd(var.Ext.34.K[[2]][,7]^2), sd(var.Ext.38.K[[2]][,7]^2), sd(var.Ext.fast.K[[2]][,7]^2), sd(var.Fast.30.K[[2]][,7]^2), sd(var.Fast.34.K[[2]][,7]^2), sd(var.Fast.38.K[[2]][,7]^2), sd(var.Fast.fast.K[[2]][,7]^2), sd(var.Int.30.K[[2]][,7]^2), sd(var.Int.34.K[[2]][,7]^2), sd(var.Int.38.K[[2]][,7]^2), sd(var.Int.fast.K[[2]][,7]^2), sd(var.Mean.30.K[[2]][,7]^2), sd(var.Mean.34.K[[2]][,7]^2), sd(var.Mean.38.K[[2]][,7]^2), sd(var.Mean.fast.K[[2]][,7]^2), sd(var.Slow.30.K[[2]][,7]^2), sd(var.Slow.34.K[[2]][,7]^2), sd(var.Slow.38.K[[2]][,7]^2), sd(var.Slow.fast.K[[2]][,7]^2)  )   )

ilist <- list(estimate = clone.K.vardata$estimate)
inits_list <- list(ilist, ilist, ilist, ilist) #Needs to be a list of lists of size nchains

#Using a measurement error model for the variances
var.K.clone.model1 <- brm(data = clone.K.vardata, family = gaussian,
                          estimate | se(var.er, sigma = TRUE) ~ 0 + intercept + treatment + assay,
                          iter = 5000, warmup = 1000, cores = 4, chains = 4,
                          inits = inits_list)

pp_check(var.K.clone.model1, nsamples = 25) #Looks OK

#For treatment, Extreme is the intercept, and for assay it is 30 C
K.clone.results.treatment <- data.frame(fixef(var.K.clone.model1))[-c(1,6:8),] #Intercept can be dropped
K.clone.results.assay <- data.frame(fixef(var.K.clone.model1))[-c(1:5),] #Intercept can be dropped
K.clone.results.treatment$treatment <- c("Fast", "Intermediate", "Mean", "Slow")
K.clone.results.assay$assay <- c("34", "38", "Fast 30-38")

#Environmental variance, growth rate
env.K.vardata <- data.frame(treatment = treatment2, assay = assay2, estimate = c( var.Ext.30.K[[1]][3,4], var.Ext.34.K[[1]][3,4], var.Ext.38.K[[1]][3,4], var.Ext.fast.K[[1]][3,4], var.Fast.30.K[[1]][3,4], var.Fast.34.K[[1]][3,4], var.Fast.38.K[[1]][3,4], var.Fast.fast.K[[1]][3,4], var.Int.30.K[[1]][3,4], var.Int.34.K[[1]][3,4], var.Int.38.K[[1]][3,4], var.Int.fast.K[[1]][3,4], var.Mean.30.K[[1]][3,4], var.Mean.34.K[[1]][3,4], var.Mean.38.K[[1]][3,4], var.Mean.fast.K[[1]][3,4], var.Slow.30.K[[1]][3,4], var.Slow.34.K[[1]][3,4], var.Slow.38.K[[1]][3,4], var.Slow.fast.K[[1]][3,4], var.anc.30.K[[1]][1,4], var.anc.34.K[[1]][1,4], var.anc.38.K[[1]][1,4], var.anc.fast.K[[1]][1,4]), var.er = c( sd(var.Ext.30.K[[2]][,8]^2), sd(var.Ext.34.K[[2]][,8]^2), sd(var.Ext.38.K[[2]][,8]^2), sd(var.Ext.fast.K[[2]][,8]^2), sd(var.Fast.30.K[[2]][,8]^2), sd(var.Fast.34.K[[2]][,8]^2), sd(var.Fast.38.K[[2]][,8]^2), sd(var.Fast.fast.K[[2]][,8]^2), sd(var.Int.30.K[[2]][,8]^2), sd(var.Int.34.K[[2]][,8]^2), sd(var.Int.38.K[[2]][,8]^2), sd(var.Int.fast.K[[2]][,8]^2), sd(var.Mean.30.K[[2]][,8]^2), sd(var.Mean.34.K[[2]][,8]^2), sd(var.Mean.38.K[[2]][,8]^2), sd(var.Mean.fast.K[[2]][,8]^2), sd(var.Slow.30.K[[2]][,8]^2), sd(var.Slow.34.K[[2]][,8]^2), sd(var.Slow.38.K[[2]][,8]^2), sd(var.Slow.fast.K[[2]][,8]^2), sd(var.anc.30.K[[2]][,6]^2), sd(var.anc.34.K[[2]][,6]^2), sd(var.anc.38.K[[2]][,6]^2), sd(var.anc.fast.K[[2]][,6]^2) ) )

ilist <- list(estimate = env.K.vardata$estimate)
inits_list <- list(ilist, ilist, ilist, ilist) #Needs to be a list of lists of size nchains

#Using a measurement error model for the variances
var.K.env.model1 <- brm(data = env.K.vardata, family = gaussian,
                          estimate | se(var.er, sigma = TRUE) ~ 0 + intercept + treatment + assay,
                          iter = 5000, warmup = 1000, cores = 4, chains = 4,
                          inits = inits_list)

pp_check(var.K.env.model1, nsamples = 25) #Looks OK

#For treatment, Extreme is the intercept, and for assay it is 30 C
K.env.results.treatment <- data.frame(fixef(var.K.env.model1))[-c(1,7:9),] #Intercept can be dropped
K.env.results.assay <- data.frame(fixef(var.K.env.model1))[-c(1:6),] #Intercept can be dropped
K.env.results.treatment$treatment <- c("Extreme", "Fast", "Intermediate", "Mean", "Slow")
K.env.results.assay$assay <- c("34", "38", "Fast 30-38")


##Plotting
my.xlab2 <- expression(paste("Relative to 30 ", degree,"C", " assay"))

##Clone variance, growth rate
clone.r.treat.plot <- ggplot(r.clone.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("Relative to Extreme treatment") +
    ylab("") +    
    ggtitle("Clone") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25)) +
    coord_flip()

#save_plot(filename = "clone_r_treatment.pdf", clone.r.trest.plot)


clone.r.assay.plot <- ggplot(r.clone.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab(my.xlab2) +
    ggtitle("Clone") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25)) +    
    coord_flip()

#clone.r.plot <- plot_grid(clone.r.trest.plot, clone.r.assay.plot)
#save_plot(filename = "./pombe_evo/fig/clone_r_treatment.pdf", clone.r.plot, base_height = 4, base_width = 7) 

##Environmental variance, growth rate
env.r.treat.plot <- ggplot(r.env.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("Relative to ancestor") +
    ggtitle("Environmental") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip()

env.r.assay.plot <- ggplot(r.env.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab(my.xlab2) +
    ggtitle("Environmental") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip()

clone.r.plot <- plot_grid(clone.r.treat.plot, clone.r.assay.plot, env.r.treat.plot, env.r.assay.plot, ncol = 4)

##Clone variance, K
clone.K.treat.plot <- ggplot(K.clone.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("Relative to Extreme treatment") +
    ggtitle("Clone") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip()

clone.K.assay.plot <- ggplot(K.clone.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab(my.xlab2) +
    ggtitle("Clone") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip()


##Environmental variance, K
env.K.treat.plot <- ggplot(K.env.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("Relative to ancestor") +
    ggtitle("Environmental") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.50, 0, 0.25)) +    
    coord_flip()

env.K.assay.plot <- ggplot(K.env.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab(my.xlab2) +
    ggtitle("Environmental") +
    ylab("") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip()

clone.K.plot <- plot_grid(clone.K.treat.plot, clone.K.assay.plot, env.K.treat.plot, env.K.assay.plot, ncol = 4)

#plot_grid(clone.r.plot, clone.K.plot, nrow = 2, labels)

var.final.plot <- plot_grid(legend_b, variances.r.plot2, variances.k.plot2, clone.r.plot, clone.K.plot, labels = c("", "A", "B", "C", "D"), ncol = 1, rel_heights = c(.1,1,1,1,1))

save_plot(filename = "./pombe_evo/fig/variances_both_2.pdf", var.final.plot, base_height = 14.3, base_width = 12)          
