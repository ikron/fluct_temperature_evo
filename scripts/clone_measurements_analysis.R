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
mylaby <- TeX("Maximum growth rate ($r_{max}$)")
plot.clones.r <- ggplot(evolved.means.all, aes(y = mean.r, x = Temperature, group = Clone.ID)) +
    geom_point() +
    geom_line() +
    geom_line(data = anc.means, aes(y = mean.r, x = Temperature, group = NA), col = "red", lwd = 1.5) +
    ylab(mylaby) +
    xlab("Assay environment") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +    
    facet_grid(Treatment ~ genotype) +
    theme(strip.text = element_text(size = 14)) +    
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))

plot.clones.K <- ggplot(evolved.means.all, aes(y = mean.K, x = Temperature, group = Clone.ID)) +
    geom_point() +
    geom_line() +
    geom_line(data = anc.means, aes(y = mean.K, x = Temperature, group = NA), col = "red", lwd = 1.5) +
    ylab("Carrying capacity (K)") +
    xlab("Assay environment") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +     
    facet_grid(Treatment ~ genotype) +
    theme(strip.text = element_text(size = 14)) +    
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))

save_plot(filename = "./fig/clones_r.pdf", plot.clones.r, base_height = 8, base_width=10)
save_plot(filename = "./fig/clones_K.pdf", plot.clones.K, base_height = 8, base_width=10)

###################################################################################################
### ** Kloonien kasvu suhteessa ancestoriin eri lämpötiloissa. Onko trade-offeja havaittavissa? ###
###################################################################################################

evolved <- filter(filtered, Population.ID != "anc")

ancestors <- filter(filtered, Population.ID == "anc")

### *** Kloonien kasvu suhteessa ancestoriin, r_max malli
#Malli rmaxille

#This works! sensible growth rates
malli.r.anc <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD,
                   data = ancestors, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.anc <- posterior_samples(malli.r.anc)[,-c(1,18,19)]
head(post.r.anc)

#(1 || Population.ID/Clone.ID),
#Mean evolution environment
malli.r.Mean <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Mean"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.Mean <- posterior_samples(malli.r.Mean)[,2:17]

#Extreme evolution environment
malli.r.Extreme <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Extreme"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.Extreme <- posterior_samples(malli.r.Extreme)[,2:17]

#Fast evolution environment
malli.r.Fast <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Fast"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.Fast <- posterior_samples(malli.r.Fast)[,2:17]

#Intermediate evolution environment
malli.r.Int <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Intermediate"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.Int <- posterior_samples(malli.r.Int)[,2:17]

#Intermediate evolution environment
malli.r.Slow <- brm(r ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Slow"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.r.Slow <- posterior_samples(malli.r.Slow)[,2:17]


#Results for Mean evol treatment
res.r.Mean <- data.frame(estimate = apply(post.r.Mean/post.r.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.r.Mean/post.r.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Mean", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Mean <- cbind( apply( (post.r.Mean/post.r.anc)[,1:4], MARGIN = 1, mean ), apply( (post.r.Mean/post.r.anc)[,5:8], MARGIN = 1, mean ), apply( (post.r.Mean/post.r.anc)[,9:12], MARGIN = 1, mean ), apply( (post.r.Mean/post.r.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Mean <- data.frame(estimate = apply(post.mean.across.anc.Mean, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Mean)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Mean", 4), genotype = rep(NA,4))


#Results for Extreme evol treatment
res.r.Extreme <- data.frame(estimate = apply(post.r.Extreme/post.r.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.r.Extreme/post.r.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Extreme", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Extreme <- cbind( apply( (post.r.Extreme/post.r.anc)[,1:4], MARGIN = 1, mean ), apply( (post.r.Extreme/post.r.anc)[,5:8], MARGIN = 1, mean ), apply( (post.r.Extreme/post.r.anc)[,9:12], MARGIN = 1, mean ), apply( (post.r.Extreme/post.r.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Extreme <- data.frame(estimate = apply(post.mean.across.anc.Extreme, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Extreme)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Extreme", 4), genotype = rep(NA,4))

#Results for Fast evol treatment
res.r.Fast <- data.frame(estimate = apply(post.r.Fast/post.r.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.r.Fast/post.r.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Fast", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Fast <- cbind( apply( (post.r.Fast/post.r.anc)[,1:4], MARGIN = 1, mean ), apply( (post.r.Fast/post.r.anc)[,5:8], MARGIN = 1, mean ), apply( (post.r.Fast/post.r.anc)[,9:12], MARGIN = 1, mean ), apply( (post.r.Fast/post.r.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Fast <- data.frame(estimate = apply(post.mean.across.anc.Fast, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Fast)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Fast", 4), genotype = rep(NA,4))

#Results for Intermediate evol treatment
res.r.Int <- data.frame(estimate = apply(post.r.Int/post.r.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.r.Int/post.r.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Intermediate", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Int <- cbind( apply( (post.r.Int/post.r.anc)[,1:4], MARGIN = 1, mean ), apply( (post.r.Int/post.r.anc)[,5:8], MARGIN = 1, mean ), apply( (post.r.Int/post.r.anc)[,9:12], MARGIN = 1, mean ), apply( (post.r.Int/post.r.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Int <- data.frame(estimate = apply(post.mean.across.anc.Int, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Int)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Intermediate", 4), genotype = rep(NA,4))

#Results for Slow evol treatment
res.r.Slow <- data.frame(estimate = apply(post.r.Slow/post.r.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.r.Slow/post.r.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Slow", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Slow <- cbind( apply( (post.r.Slow/post.r.anc)[,1:4], MARGIN = 1, mean ), apply( (post.r.Slow/post.r.anc)[,5:8], MARGIN = 1, mean ), apply( (post.r.Slow/post.r.anc)[,9:12], MARGIN = 1, mean ), apply( (post.r.Slow/post.r.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Slow <- data.frame(estimate = apply(post.mean.across.anc.Slow, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Slow)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Slow", 4), genotype = rep(NA,4))

#Combine all the results
res.r.clones <- rbind(res.r.Extreme, res.r.Fast, res.r.Int, res.r.Mean, res.r.Slow)
res.r.mean.across.anc <- rbind(res.mean.across.anc.Extreme, res.mean.across.anc.Fast, res.mean.across.anc.Int, res.mean.across.anc.Mean, res.mean.across.anc.Slow)

#Make a plot
ylab2 <- TeX("Relative $r_{max}$")
p.r.clones <- ggplot(res.r.clones, aes(x = assay, y = estimate, ymin = lower, ymax = upper, group = genotype)) +
    geom_pointrange(position = position_dodge(width = 0.5), alpha = 0.5, linewidth = 1) +
    geom_pointrange(data = res.r.mean.across.anc, aes(x = assay, y = estimate, ymin = lower, ymax = upper), colour = "blue", alpha = 0.8, linewidth = 1) +
    xlab("Assay environment") +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +         
    ylab(ylab2) +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap(. ~ evolution, scales = "free", nrow = 1) +
    theme(strip.text = element_text(size = 14)) +    
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1), axis.title.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))    

### *** Kloonien kasvu suhteessa ancestoriin K mallit

#Malli K:lle

#This works! sensible growth rates
malli.K.anc <- brm(k ~ -1 + genotype:Temperature + Pregrowth.OD,
                   data = ancestors, family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.K.anc <- posterior_samples(malli.K.anc)[,-c(1,18,19)]
head(post.K.anc)

#(1 || Population.ID/Clone.ID),
#Mean evolution environment
malli.K.Mean <- brm(k ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Mean"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.K.Mean <- posterior_samples(malli.K.Mean)[,2:17]

#Extreme evolution environment
malli.K.Extreme <- brm(k ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Extreme"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.K.Extreme <- posterior_samples(malli.K.Extreme)[,2:17]

#Fast evolution environment
malli.K.Fast <- brm(k ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Fast"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.K.Fast <- posterior_samples(malli.K.Fast)[,2:17]

#Intermediate evolution environment
malli.K.Int <- brm(k ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Intermediate"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.K.Int <- posterior_samples(malli.K.Int)[,2:17]

#Intermediate evolution environment
malli.K.Slow <- brm(k ~ -1 + genotype:Temperature + Pregrowth.OD + (1 | Clone.ID),
                    data = filter(evolved, Evolution.treatment == "Slow"), family = gaussian(), control = list(adapt_delta = 0.98), iter = 8000, thin = 2, cores = 4)
post.K.Slow <- posterior_samples(malli.K.Slow)[,2:17]


#Results for Mean evol treatment
res.K.Mean <- data.frame(estimate = apply(post.K.Mean/post.K.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.K.Mean/post.K.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Mean", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Mean.K <- cbind( apply( (post.K.Mean/post.K.anc)[,1:4], MARGIN = 1, mean ), apply( (post.K.Mean/post.K.anc)[,5:8], MARGIN = 1, mean ), apply( (post.K.Mean/post.K.anc)[,9:12], MARGIN = 1, mean ), apply( (post.K.Mean/post.K.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Mean.K <- data.frame(estimate = apply(post.mean.across.anc.Mean.K, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Mean.K)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Mean", 4), genotype = rep(NA,4))


#Results for Extreme evol treatment
res.K.Extreme <- data.frame(estimate = apply(post.K.Extreme/post.K.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.K.Extreme/post.K.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Extreme", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Extreme.K <- cbind( apply( (post.K.Extreme/post.K.anc)[,1:4], MARGIN = 1, mean ), apply( (post.K.Extreme/post.K.anc)[,5:8], MARGIN = 1, mean ), apply( (post.K.Extreme/post.K.anc)[,9:12], MARGIN = 1, mean ), apply( (post.K.Extreme/post.K.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Extreme.K <- data.frame(estimate = apply(post.mean.across.anc.Extreme.K, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Extreme.K)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Extreme", 4), genotype = rep(NA,4))

#Results for Fast evol treatment
res.K.Fast <- data.frame(estimate = apply(post.K.Fast/post.K.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.K.Fast/post.K.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Fast", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Fast.K <- cbind( apply( (post.K.Fast/post.K.anc)[,1:4], MARGIN = 1, mean ), apply( (post.K.Fast/post.K.anc)[,5:8], MARGIN = 1, mean ), apply( (post.K.Fast/post.K.anc)[,9:12], MARGIN = 1, mean ), apply( (post.K.Fast/post.K.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Fast.K <- data.frame(estimate = apply(post.mean.across.anc.Fast.K, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Fast.K)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Fast", 4), genotype = rep(NA,4))

#Results for Intermediate evol treatment
res.K.Int <- data.frame(estimate = apply(post.K.Int/post.K.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.K.Int/post.K.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Intermediate", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Int.K <- cbind( apply( (post.K.Int/post.K.anc)[,1:4], MARGIN = 1, mean ), apply( (post.K.Int/post.K.anc)[,5:8], MARGIN = 1, mean ), apply( (post.K.Int/post.K.anc)[,9:12], MARGIN = 1, mean ), apply( (post.K.Int/post.K.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Int.K <- data.frame(estimate = apply(post.mean.across.anc.Int.K, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Int.K)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Intermediate", 4), genotype = rep(NA,4))

#Results for Slow evol treatment
res.K.Slow <- data.frame(estimate = apply(post.K.Slow/post.K.anc, MARGIN = 2, median), HPDinterval(as.mcmc(post.K.Slow/post.K.anc)), assay = c(rep("30",4), rep("34", 4), rep("38",4), rep("Fast 30-38", 4)), genotype = rep(c("h+M210", "h+M216", "h-M210", "h-M216"),4), evolution = rep("Slow", 16))
#Calculate means across different ancestor clones
post.mean.across.anc.Slow.K <- cbind( apply( (post.K.Slow/post.K.anc)[,1:4], MARGIN = 1, mean ), apply( (post.K.Slow/post.K.anc)[,5:8], MARGIN = 1, mean ), apply( (post.K.Slow/post.K.anc)[,9:12], MARGIN = 1, mean ), apply( (post.K.Slow/post.K.anc)[,13:16], MARGIN = 1, mean ) )
res.mean.across.anc.Slow.K <- data.frame(estimate = apply(post.mean.across.anc.Slow.K, MARGIN = 2, median), HPDinterval(as.mcmc(post.mean.across.anc.Slow.K)), assay = c("30", "34", "38", "Fast 30-38"), evolution = rep("Slow", 4), genotype = rep(NA,4))

#Combine all the results
res.K.clones <- rbind(res.K.Extreme, res.K.Fast, res.K.Int, res.K.Mean, res.K.Slow)
res.K.mean.across.anc <- rbind(res.mean.across.anc.Extreme.K, res.mean.across.anc.Fast.K, res.mean.across.anc.Int.K, res.mean.across.anc.Mean.K, res.mean.across.anc.Slow.K)

save(res.r.clones, res.r.mean.across.anc, res.K.clones, res.K.mean.across.anc, file = "./data/clone_results.RData") #Saving the results for later

#Make a plot
p.K.clones <- ggplot(res.K.clones, aes(x = assay, y = estimate, ymin = lower, ymax = upper, group = genotype)) +
    geom_pointrange(position = position_dodge(width = 0.5), alpha = 0.5, linewidth = 1) +
    geom_pointrange(data = res.K.mean.across.anc, aes(x = assay, y = estimate, ymin = lower, ymax = upper), colour = "blue", alpha = 0.8, linewidth = 1) +
    xlab("Assay environment") +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +         
    ylab("Relative K") +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap(. ~ evolution, scales = "free", nrow = 1) +
    theme(strip.text = element_text(size = 14)) +    
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1), axis.title.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))    



### *** Piirretään lopulliset kuvat kloonien kasvusta

p.clones.final.r <- plot_grid(plot.clones.r, p.r.clones, rel_heights = c(2.5,1), nrow = 2, labels = c("A", "B"))

save_plot(filename = "./fig/clones_r2.pdf", p.clones.final.r, base_height = 11, base_width=10)

p.clones.final.K <- plot_grid(plot.clones.K, p.K.clones, rel_heights = c(2.5,1), nrow = 2, labels = c("A", "B"))

save_plot(filename = "./fig/clones_K2.pdf", p.clones.final.K, base_height = 11, base_width=10)

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

r.clone.results.treatment$variance <- "clone"
r.clone.results.assay$variance <- "clone"
r.clone.results.treatment$what <- "treatment effect"
r.clone.results.assay$what <- "assay effect"

#Calculating percentages for text


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

r.env.results.treatment$variance <- "environmental"
r.env.results.assay$variance <- "environmental"
r.env.results.treatment$what <- "treatment effect"
r.env.results.assay$what <- "assay effect"

#### Saving the results so models don't need to be run again
#save(r.clone.results.treatment, r.clone.results.assay, r.env.results.treatment, r.env.results.assay, file = "./data/var_results_r.RData")

load("./data/var_results_r.RData")

###Models for variance for carrying capacity, K

#Clone variance, K
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

K.clone.results.treatment$variance <- "clone"
K.clone.results.assay$variance <- "clone"
K.clone.results.treatment$what <- "treatment effect"
K.clone.results.assay$what <- "assay effect"

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

K.env.results.treatment$variance <- "environmental"
K.env.results.assay$variance <- "environmental"
K.env.results.treatment$what <- "treatment effect"
K.env.results.assay$what <- "assay effect"

#### Saving the results so models don't need to be run again
#save(K.clone.results.treatment, K.clone.results.assay, K.env.results.treatment, K.env.results.assay, file = "./data/var_results_K.RData")

### **** Making the final plot

#Load the data
load("./data/var_results_r.RData")
load("./data/var_results_K.RData")

variances.r.plot2 <- ggplot(var.results, aes(x = treatment, y = estimate, ymin = lower, ymax = upper, colour = assay)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    #scale_colour_manual(name = "Assay\nenvironment") +
    ylab("Variance") +
    xlab("Evolutionary treatment") +
    labs(colour = "Assay environment") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    guides(color = guide_legend(nrow = 1)) +    
    #theme(text=element_text(size=16)) +
    facet_grid( ~ type) +
    theme(legend.position = "top", legend.justification = "center") +
    theme(strip.text = element_text(size = 14)) +    
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))    

variances.k.plot2 <- ggplot(var.results.K, aes(x = treatment, y = estimate, ymin = lower, ymax = upper, colour = assay)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    ylab("Variance") +
    xlab("Evolutionary treatment") +
    labs(colour = "Assay environment") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +    
    facet_grid( ~ type) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))

#save_plot(filename = "./pombe_evo/fig/variances_k2.pdf", variances.k.plot2, base_height = 3, base_width = 12)


#Get legend function has some problems with this particular version of ggplot2
#The function below is now working
get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)

  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

#Now this works with the new function
legend_b <- get_legend(variances.r.plot2)

#variances.r.plot2 <- variances.r.plot2 + theme(legend.position = "none")
#variances.k.plot2 <- variances.k.plot2 + theme(legend.position = "none")

var.final.plot <- plot_grid(legend_b, variances.r.plot2 + theme(legend.position = "none"), variances.k.plot2 + theme(legend.position = "none"), labels = c("", "A", "B"), ncol = 1, rel_heights = c(.1,1,1))


my.xlab2 <- expression(paste("Difference to 30 ", degree,"C"))

##Clone variance, growth rate
clone.r.treat.plot <- ggplot(r.clone.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("") +
    ylab("Difference to Extreme") +    
    #ggtitle("Clone") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25)) +
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))

#save_plot(filename = "clone_r_treatment.pdf", clone.r.trest.plot)


clone.r.assay.plot <- ggplot(r.clone.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Clone") +
    ylab(my.xlab2) +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25)) +    
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))    

#clone.r.plot <- plot_grid(clone.r.trest.plot, clone.r.assay.plot)
#save_plot(filename = "./pombe_evo/fig/clone_r_treatment.pdf", clone.r.plot, base_height = 4, base_width = 7)
##Plotting


##Environmental variance, growth rate
env.r.treat.plot <- ggplot(r.env.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Environmental") +
    ylab("Difference to ancestor") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))    

env.r.assay.plot <- ggplot(r.env.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Environmental") +
    ylab(my.xlab2) +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))

#clone.r.plot <- plot_grid(clone.r.treat.plot, clone.r.assay.plot, env.r.treat.plot, env.r.assay.plot, ncol = 4)


clone.r.plot <- plot_grid(NULL, clone.r.treat.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, clone.r.assay.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, env.r.treat.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, env.r.assay.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), nrow = 1, rel_widths = c(-0.05, 1, -0.05, 1, -0.05, 1, -0.05, 1 ))



##Clone variance, K
clone.K.treat.plot <- ggplot(K.clone.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Clone") +
    ylab("Difference to Extreme") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.4, 0, 0.4)) +     
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))    

clone.K.assay.plot <- ggplot(K.clone.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Clone") +
    ylab(my.xlab2) +    
    geom_hline(yintercept = 0, lty = "dashed") +
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12)) 


##Environmental variance, K
env.K.treat.plot <- ggplot(K.env.results.treatment, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = treatment)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Environmental") +
    ylab("Difference to ancestor") +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.50, 0, 0.25)) +    
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12)) 

env.K.assay.plot <- ggplot(K.env.results.assay, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5, x = assay)) +
    geom_pointrange() +
    xlab("") +
    #ggtitle("Environmental") +
    ylab(my.xlab2) +    
    geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2)) +     
    coord_flip() +
    facet_wrap(what ~ variance) +
    theme(strip.text = element_text(size = 14)) +        
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12)) 

#clone.K.plot <- plot_grid(clone.K.treat.plot, clone.K.assay.plot, env.K.treat.plot, env.K.assay.plot, ncol = 4)

clone.K.plot <- plot_grid(NULL, clone.K.treat.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, clone.K.assay.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, env.K.treat.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, env.K.assay.plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), nrow = 1, rel_widths = c(-0.05, 1, -0.05, 1, -0.05, 1, -0.05, 1 ))

#plot_grid(clone.r.plot, clone.K.plot, nrow = 2, labels)

#var.final.plot <- plot_grid(legend_b, variances.r.plot2 + theme(legend.position = "none"), variances.k.plot2 + theme(legend.position = "none"), clone.r.plot, clone.K.plot, labels = c("", "A", "B", "C", "D"), ncol = 1, rel_heights = c(.1,1,1,1,1))

var.final.plot <- plot_grid(legend_b, variances.r.plot2 + theme(legend.position = "none"),  clone.r.plot, variances.k.plot2 + theme(legend.position = "none"), clone.K.plot, labels = c("", "A", "", "B", ""), ncol = 1, rel_heights = c(.1,1,0.7,1,0.7))

save_plot(filename = "./fig/variances_both_final.pdf", var.final.plot, base_height = 12, base_width = 13)


var.final.plot <- plot_grid(legend_b, variances.r.plot2 + theme(legend.position = "none"),  clone.r.plot + theme(axis.title.x = element_text(size = 12)), variances.k.plot2 + theme(legend.position = "none"), clone.K.plot, labels = c("", "A", "", "B", ""), ncol = 1, rel_heights = c(.1,1,0.7,1,0.7))

save_plot(filename = "./fig/variances_both_final.pdf", var.final.plot, base_height = 12, base_width = 10)      
