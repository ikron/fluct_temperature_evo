#Data analysis for competition experiments and clone measurements made for 
#S. pombe strains in  Räsänen et al. 2024.

#Load packages

library(dplyr)
library(forcats)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(brms)
library(coda)
library(signs)
library(scales)
library(grid)

#Load the data

data <- read.csv("Pombe_competition_experiment_data.csv", header = T, dec = ".", sep = ",")

#Add number of generations during competition experiment, for intermediate 
#fluctuations 2.5 generations*13 days and other treatments 2.5 generations*4 days

data$generation <- ifelse(data$Competition.environment == "Intermediate", 32.5, 10)

#Function to calculate relative fitness from competition assays

calc.w <- function(data) {
    ncomp <- nrow(data)
    res <- rep(0, ncomp) #Store results here
    #When evolved strain had zero frequency
    check <- data$Evolved.strain.final.CFU == 0
    data$Evolved.strain.final.CFU[check] <- 0.5 #Adding a small number, when evolved strain lost completely
    check2 <- data$Ancestor.final.CFU == 0
    data$Ancestor.final.CFU[check2] <- 0.5
    #Calculate log-ratios
    data$initf <- log(data$Evolved.strain.initial.CFU/data$Ancestor.initial.CFU)
    data$finalf <- log(data$Evolved.strain.final.CFU/data$Ancestor.final.CFU)
    #
    for(i in 1:ncomp) { #Loop over all competitions
        dat <- data[i,] #Take a row
        #ymat is the data for fitting the model
        ymat <- data.frame( ylog = c(dat$initf, dat$finalf), gen = c(0, dat$generation))
        #Check if values exist, TEMPORARY!
        if(is.infinite(ymat[2,1]) == T | is.nan(ymat[2,1]) == T) { res[i] <- NA } else {
        #Fit the model
        modelfit <- lm(ylog ~ gen,  data = ymat)
        wlog <- coef(modelfit)[2] #Extract slope
        res[i] <- exp(wlog) #Store results and transform to normal scale
        }
    }
#
    return(res)
}

#Use function and add relative fitness to data

data$relW <- calc.w(data) 

#Code strain ID's so that they are unique by combining evolution environment
#and evolved strain ID (well on a plate)

data$StrainID <- paste(data$Evolution.environment, data$Evolved.strain.ID, sep = "")

#Code treatment factor by combining evolution environment and competition environment

data$treatment <- paste("EVO_",data$Evolution.environment, "COMP_", data$Competition.environment, sep = "")

#Code evolved strain genotype factor by combining mating type and ade6 marker

data$evol.strain.genotype <- paste(data$Evolved.strain.mating.type, data$Evolved.strain.ade6.marker, sep = "/")

#Check the structure of the data

str(data)

#Population size needs to be a factor

data$Population.size <- factor(data$Population.size)

#Filter small and large populations as separate data

data.smallpop <- filter(data, Population.size == 96)
data.largepop <- filter(data, Population.size == 24)

#Running the model with a Bayesian approach by brms package

#Models for small and large populations separately to get posterior estimates
#for calculations and figures

model_small <- brm(relW ~ -1 + treatment + Ancestor.genotype + (1 | StrainID), data = data.smallpop, family = "gaussian", cores = 4,iter = 4000, warmup = 2000, thin = 2)

model_large <- brm(relW ~ -1 + treatment + Ancestor.genotype + (1 | StrainID), data = data.largepop, family = "gaussian", cores = 4, iter = 4000, warmup = 2000, thin = 2)

#Model for complete data to test population size effect and to get main results

model_all <- brm(relW ~ -1 + treatment + treatment*Population.size + Ancestor.genotype + (1 | StrainID), data = data, family = "gaussian", cores = 4, iter = 4000, warmup = 2000, thin = 2)
              
#Check results

summary(model_small)
plot(model_small)

summary(model_large)
plot(model_large)

summary(model_all)
plot(model_all)

#Save posteriors so that there is no need for running brms every time
#save(model_small, model_large, model_all, file = "Pombe_competition_experiment_models_original.RData")

#Load posteriors saved in Rdata
load("Pombe_competition_experiment_models_original.RData")

#Figure 2: The relative fitness of the evolved strain compared to the ancestor over
#all populations, shown as a posterior mean with 95 % HPDI. The dashed line is
#the ancestor fitness set to one.

post.all <- data.frame(fixef(model_all)[1:12,])
post.all$Evolution <- c("Extreme", "Extreme", "Fast", "Fast", "Intermediate", "Intermediate", rep("Mean",4), rep("Slow", 2))
post.all$Competition <- c("Extreme", "Mean", "Fast", "Mean", "Intermediate", "Mean", "Extreme", "Fast", "Intermediate", "Mean", "Extreme", "Mean")

relw_all <- ggplot(data = post.all,  aes(x = Evolution, y = Estimate, ymin = Q2.5, ymax = Q97.5, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  theme_classic()+
  xlab("Evolution environment") +
  ylab("Relative fitness") +  
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_color_manual(values=c("#f50020","#99cf00","#c415b7","#0a74ea"))+
  labs(colour = "Competition environment")+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=16))

save_plot("Figure_2.png", relw_all, base_height=3.71, base_width = 3.71*2*1.618)

#Figure 3: The relative fitness of the evolved populations compared 
#to their ancestors for data divided by the population size. The relative fitness is
#shown as a posterior mean with 95 % HPDI. The dashed line is ancestor fitness set to one.

post.small <- data.frame(fixef(model_small)[1:12,])
post.small$Evolution <- c("Extreme", "Extreme", "Fast", "Fast", "Intermediate", "Intermediate", rep("Mean",4), rep("Slow", 2))
post.small$Competition <- c("Extreme", "Mean", "Fast", "Mean", "Intermediate", "Mean", "Extreme", "Fast", "Intermediate", "Mean", "Extreme", "Mean")

post.large <- data.frame(fixef(model_large)[1:12,])
post.large$Evolution <- c("Extreme", "Extreme", "Fast", "Fast", "Intermediate", "Intermediate", rep("Mean",4), rep("Slow", 2))
post.large$Competition <- c("Extreme", "Mean", "Fast", "Mean", "Intermediate", "Mean", "Extreme", "Fast", "Intermediate", "Mean", "Extreme", "Mean")

post.popsize <- rbind(post.large, post.small)
post.popsize$popsize <- c(rep("Large", 12), rep("Small", 12))

relw_popsize <- ggplot(data = post.popsize,  aes(x = Evolution, y = Estimate, ymin = Q2.5, ymax = Q97.5, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  xlab("Evolution environment") +
  ylab("Relative fitness") +  
  theme_classic()+
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_color_manual(values=c("#f50020","#99cf00","#c415b7","#0a74ea"))+
  labs(colour = "Competition environment")+
  facet_grid(popsize ~ .)+
  theme(strip.text = element_text(
    size = 14))+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=16))

save_plot("Figure_3.png", relw_popsize, base_height=3.71, base_width = 3.71*2*1.618)


#Figure 4: The difference in relative fitness between population sizes at 
#different temperature treatments. The differences are statistically significant
#when the 95 % HPDI does not overlap with zero.

#Matrix for posterior results, note that estimate is mean

differences.matrix.pop <- data.frame(matrix(rep(0, 12*6), ncol = 6))
colnames(differences.matrix.pop) <- c("Evolution", "Competition", "Comparison", "mean", "lower", "upper")

#Evolution environment is given first, then competition environments

differences.matrix.pop$Evolution <- c(rep("Extreme", 2), rep("Fast",2), rep("Intermediate",2), rep("Mean", 4),
                                      rep("Slow",2))

differences.matrix.pop$Competition <- c("Extreme", "Mean", "Fast", "Mean", "Intermediate", "Mean", "Extreme", 
                                        "Fast", "Intermediate", "Mean", "Mean", "Extreme")

differences.matrix.pop$Comparison <- rep("Large-Small", 12)

differences.matrix.pop

#Select posteriors for fixed factors

differences.post.small <- posterior_samples(model_small)[,1:12]
differences.post.large <- posterior_samples(model_large)[,1:12]

colnames(differences.post.small)
colnames(differences.post.large)

#Differences between population sizes within competition environment

ext.ext.pop <- differences.post.large[,1] - differences.post.small[,1]
ext.mean.pop <- differences.post.large[,2] - differences.post.small[,2]

fast.fast.pop <- differences.post.large[,3] - differences.post.small[,3]
fast.mean.pop <- differences.post.large[,4] - differences.post.small[,4]

int.int.pop <- differences.post.large[,5] - differences.post.small[,5]
int.mean.pop <- differences.post.large[,6] - differences.post.small[,6]

mean.ext.pop <- differences.post.large[,7] - differences.post.small[,7]
mean.fast.pop <- differences.post.large[,8] - differences.post.small[,8]
mean.int.pop <- differences.post.large[,9] - differences.post.small[,9]
mean.mean.pop <- differences.post.large[,10] - differences.post.small[,10]

slow.ext.pop <- differences.post.large[,11] - differences.post.small[,11]
slow.mean.pop <- differences.post.large[,12] - differences.post.small[,12]

#Fill matrix with posterior estimates of the differences

differences.matrix.pop[1,4:6] <- c(mean(ext.ext.pop), HPDinterval(as.mcmc(ext.ext.pop)))
differences.matrix.pop[2,4:6] <- c(mean(ext.mean.pop), HPDinterval(as.mcmc(ext.mean.pop)))

differences.matrix.pop[3,4:6] <- c(mean(fast.fast.pop), HPDinterval(as.mcmc(fast.fast.pop)))
differences.matrix.pop[4,4:6] <- c(mean(fast.mean.pop), HPDinterval(as.mcmc(fast.mean.pop)))

differences.matrix.pop[5,4:6] <- c(mean(int.int.pop), HPDinterval(as.mcmc(int.int.pop)))
differences.matrix.pop[6,4:6] <- c(mean(int.mean.pop), HPDinterval(as.mcmc(int.mean.pop)))

differences.matrix.pop[7,4:6] <- c(mean(mean.ext.pop), HPDinterval(as.mcmc(mean.ext.pop)))
differences.matrix.pop[8,4:6] <- c(mean(mean.fast.pop), HPDinterval(as.mcmc(mean.fast.pop)))
differences.matrix.pop[9,4:6] <- c(mean(mean.int.pop), HPDinterval(as.mcmc(mean.int.pop)))
differences.matrix.pop[10,4:6] <- c(mean(mean.mean.pop), HPDinterval(as.mcmc(mean.mean.pop)))

differences.matrix.pop[11,4:6] <- c(mean(slow.ext.pop), HPDinterval(as.mcmc(slow.ext.pop)))
differences.matrix.pop[12,4:6] <- c(mean(slow.mean.pop), HPDinterval(as.mcmc(slow.mean.pop)))

#Check matrix and change to data frame

differences.matrix.pop
differences.matrix.pop <- data.frame(differences.matrix.pop)

#Make separate plots by evolution environment

#Figure 4 A) Evolution Extreme, competition Extreme - Mean

plot_evoext_pop <- ggplot(data = differences.matrix.pop[c(1,2),], aes(x = Comparison, y = mean, ymin = lower, ymax = upper, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n population sizes") +
  xlab(" ")+
  labs(title = "Evolution environment Extreme")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#f50020","#0a74ea"))+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

##theme(axis.text.x = element_text(face="bold", size = 15))

save_plot("Figure_4A.png", plot_evoext_pop)

#Figure 4 B) Evolution Mean, competition Mean - Extreme, Mean - Fast, Mean - Intermediate

plot_evomean_pop <- ggplot(data = differences.matrix.pop[c(7:10),], aes(x = Comparison, y = mean, ymin = lower, ymax = upper, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n population sizes") +
  xlab(" ")+
  labs(title = "Evolution environment Mean")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#f50020","#99cf00","#c415b7","#0a74ea"))+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_4B.png", plot_evomean_pop)

#Figure 4 C) Evolution Fast, competition Fast - Mean 

plot_evofast_pop <- ggplot(data = differences.matrix.pop[c(3,4),], aes(x = Comparison, y = mean, ymin = lower, ymax = upper, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n population sizes") +
  xlab(" ")+
  labs(title = "Evolution environment Fast")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#99cf00","#0a74ea"))+
  scale_y_continuous(limits=c(-0.12,0.4), breaks=seq(-0.12,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_4C.png", plot_evofast_pop)

#Figure 4 D) Evolution Intermediate, competition Intermediate - Mean 

plot_evoint_pop <- ggplot(data = differences.matrix.pop[c(5,6),], aes(x = Comparison, y = mean, ymin = lower, ymax = upper, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n population sizes") +
  xlab(" ")+
  labs(title = "Evolution environment Intermediate")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#c415b7","#0a74ea"))+
  scale_y_continuous(limits=c(-0.12,0.4), breaks=seq(-0.12,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_4D.png", plot_evoint_pop)

#Figure 4 E) Evolution Slow, competition Mean - Extreme 

plot_evoslow_pop <- ggplot(data = differences.matrix.pop[c(11,12),], aes(x = Comparison, y = mean, ymin = lower, ymax = upper, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n population sizes") +
  xlab(" ")+
  labs(title = "Evolution environment Slow")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#f50020","#0a74ea"))+
  scale_y_continuous(limits=c(-0.12,0.4), breaks=seq(-0.12,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_4E.png", plot_evoslow_pop)

#Separate plot for legend

plot_legend_comp <- ggplot(data = differences.matrix.pop[c(7:10),], aes(x = Comparison, y = mean, ymin = lower, ymax = upper, colour = Competition)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n population sizes") +
  xlab(" ")+
  labs(title = "Evolution environment Mean")+
  theme_classic()+
  scale_color_manual(values=c("#f50020","#99cf00","#c415b7","#0a74ea"))+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=16, face = "bold"))

legend_comp <- get_legend(plot_legend_comp)
grid.newpage()
grid.draw(legend_comp)

save_plot("Figure_4_legend.png", legend_comp)

#Combine separate plots and legend into one figure

differences_pop <-plot_grid(plot_evoext_pop, plot_evomean_pop, plot_evofast_pop, plot_evoint_pop, plot_evoslow_pop, legend_comp,
                            nrow = 3,
                            ncol = 2,
                            rel_heights = c(2,2),
                            rel_widths = c(1,1,2),
                            labels = c("(A)", "(B)", "(C)", "(D)", "(E)"))

save_plot("Figure_4.png", differences_pop, base_height=7.71, base_width = 3.71*2*1.618)


#Figure 5: The difference in relative fitness between competition environments for
#large (green) and small (yellow) populations. For the statistical significance, the
#95 % HPDI does not overlap with zero.

#Matrix for posterior results, note that estimate is mean

differences.matrix.comp <- data.frame(matrix(rep(0, 14*6), ncol = 6))
colnames(differences.matrix.comp) <- c("Evolution", "Competition", "populationsize", "mean", "lower", "upper")

#Evolution environment is given first, then competition environments between which relative fitness is compared 

differences.matrix.comp$Evolution <- rep(c("Extreme", "Fast", "Intermediate", "Mean", "Mean",
                                       "Mean", "Slow"), 2)

differences.matrix.comp$Competition <- rep(c("Extreme - Mean", "Fast - Mean", "Intermediate - Mean", "Mean - Extreme", 
                                        "Mean - Fast", "Mean - Intermediate", "Mean - Extreme"), 2)

differences.matrix.comp$populationsize <- c(rep("Small", 7), rep("Large", 7))

differences.matrix.comp

#To remember the columns

colnames(differences.post.small)
colnames(differences.post.large)

#Differences between competition environments within population size

#Small populations

ext.ext_mean.small <- differences.post.small[,1] - differences.post.small[,2]
fast.fast_mean.small <- differences.post.small[,3] - differences.post.small[,4] 
int.int_mean.small <- differences.post.small[,5] - differences.post.small[,6]
mean.mean_ext.small <- differences.post.small[,10] - differences.post.small[,7]  
mean.mean_fast.small <- differences.post.small[,10] - differences.post.small[,8]
mean.mean_int.small <- differences.post.small[,10] - differences.post.small[,9]
slow.mean_ext.small <- differences.post.small[,12] - differences.post.small[,11]

#Large populations

ext.ext_mean.large <- differences.post.large[,1] - differences.post.large[,2]
fast.fast_mean.large <- differences.post.large[,3] - differences.post.large[,4] 
int.int_mean.large <- differences.post.large[,5] - differences.post.large[,6]
mean.mean_ext.large <- differences.post.large[,10] - differences.post.large[,7]
mean.mean_fast.large <- differences.post.large[,10] - differences.post.large[,8]
mean.mean_int.large <- differences.post.large[,10] - differences.post.large[,9]
slow.mean_ext.large <- differences.post.large[,12] - differences.post.large[,11]

#Fill matrix with posterior estimates of the differences

#Small populations

differences.matrix.comp[1,4:6] <- c(mean(ext.ext_mean.small), HPDinterval(as.mcmc(ext.ext_mean.small)))
differences.matrix.comp[2,4:6] <- c(mean(fast.fast_mean.small), HPDinterval(as.mcmc(fast.fast_mean.small)))
differences.matrix.comp[3,4:6] <- c(mean(int.int_mean.small), HPDinterval(as.mcmc(int.int_mean.small)))
differences.matrix.comp[4,4:6] <- c(mean(mean.mean_ext.small), HPDinterval(as.mcmc(mean.mean_ext.small)))
differences.matrix.comp[5,4:6] <- c(mean(mean.mean_fast.small), HPDinterval(as.mcmc(mean.mean_fast.small)))
differences.matrix.comp[6,4:6] <- c(mean(mean.mean_int.small), HPDinterval(as.mcmc(mean.mean_int.small)))
differences.matrix.comp[7,4:6] <- c(mean(slow.mean_ext.small), HPDinterval(as.mcmc(slow.mean_ext.small)))

#Large populations

differences.matrix.comp[8,4:6] <- c(mean(ext.ext_mean.large), HPDinterval(as.mcmc(ext.ext_mean.large)))
differences.matrix.comp[9,4:6] <- c(mean(fast.fast_mean.large), HPDinterval(as.mcmc(fast.fast_mean.large)))
differences.matrix.comp[10,4:6] <- c(mean(int.int_mean.large), HPDinterval(as.mcmc(int.int_mean.large)))
differences.matrix.comp[11,4:6] <- c(mean(mean.mean_ext.large), HPDinterval(as.mcmc(mean.mean_ext.large)))
differences.matrix.comp[12,4:6] <- c(mean(mean.mean_fast.large), HPDinterval(as.mcmc(mean.mean_fast.large)))
differences.matrix.comp[13,4:6] <- c(mean(mean.mean_int.large), HPDinterval(as.mcmc(mean.mean_int.large)))
differences.matrix.comp[14,4:6] <- c(mean(slow.mean_ext.large), HPDinterval(as.mcmc(slow.mean_ext.large)))

#Check matrix and change to data frame

differences.matrix.comp
differences.matrix.comp <- data.frame(differences.matrix.comp)

#Make separate plots by evolution environment

#Figure 5 A) Evolution Extreme, competition Extreme - Mean

plot_evoext_comp <- ggplot(data = differences.matrix.comp[c(1,8),], aes(x = Competition, y = mean, ymin = lower, ymax = upper, colour = populationsize)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n competition environments") +
  xlab(" ")+
  labs(title = "Evolution environment Extreme")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#0cd7c5","#efae2e"))+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_5A.png", plot_evoext_comp)

#Figure 5 B) Evolution Mean, competition Mean - Extreme, Mean - Fast, Mean - Intermediate 

plot_evomean_comp <- ggplot(data = differences.matrix.comp[c(4:6,11:13),], aes(x = Competition, y = mean, ymin = lower, ymax = upper, colour = populationsize)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n competition environments") +
  xlab(" ")+
  labs(title = "Evolution environment Mean")+
  theme_classic()+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#0cd7c5","#efae2e"))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_5B.png", plot_evomean_comp)

#Figure 5 C) Evolution Fast, competition Fast - Mean 

plot_evofast_comp <- ggplot(data = differences.matrix.comp[c(2,9),], aes(x = Competition, y = mean, ymin = lower, ymax = upper, colour = populationsize)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n competition environments") +
  xlab(" ")+
  labs(title = "Evolution environment Fast")+
  theme_classic()+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#0cd7c5","#efae2e"))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_5C.png", plot_evofast_comp)


#Figure 5 D) Evolution Intermediate, competition Intermediate - Mean 

plot_evoint_comp <- ggplot(data = differences.matrix.comp[c(3,10),], aes(x = Competition, y = mean, ymin = lower, ymax = upper, colour = populationsize)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n competition environments") +
  xlab(" ")+
  labs(title = "Evolution environment Intermediate")+
  theme_classic()+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#0cd7c5","#efae2e"))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_5D.png", plot_evoint_comp)

#Figure 5 E) Evolution Slow, competition Mean - Extreme 

plot_evoslow_comp <- ggplot(data = differences.matrix.comp[c(7,14),], aes(x = Competition, y = mean, ymin = lower, ymax = upper, colour = populationsize)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n competition environments") +
  xlab(" ")+
  labs(title = "Evolution environment Slow")+
  theme_classic()+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#0cd7c5","#efae2e"))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))

save_plot("Figure_5E.png", plot_evoslow_comp)

#Separate plot for legend

legend_plot_pop <- ggplot(data = differences.matrix.comp[c(7,14),], aes(x = Competition, y = mean, ymin = lower, ymax = upper, colour = populationsize)) +
  geom_pointrange(position= position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ylab("Difference in relative fitness between \n competition environments") +
  xlab(" ")+
  labs(title = "Evolution environment Slow")+
  theme_classic()+
  scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4, 0.1), labels=signs_format(accuracy=0.1))+
  guides(color= guide_legend(title= "Population size"))+
  scale_color_manual(values=c("#0cd7c5","#efae2e"))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(plot.title=element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=16, face = "bold"))

legend_pop <- get_legend(legend_plot_pop)
grid.newpage()
grid.draw(legend_pop)

save_plot("Figure_5_legend.png", legend_pop)

#Combine separate plots and legend into one figure

differences_comp <-plot_grid(plot_evoext_comp,plot_evomean_comp,plot_evofast_comp,plot_evoint_comp,plot_evoslow_comp, legend_pop,
                        nrow = 3,
                        ncol = 2,
                        rel_heights = c(2,2),
                        rel_widths = c(1,1,2),
                        labels = c("(A)", "(B)", "(C)", "(D)", "(E)"))

save_plot("Figure_5.png", differences_comp, base_height=7.71, base_width = 3.71*2*1.618)
