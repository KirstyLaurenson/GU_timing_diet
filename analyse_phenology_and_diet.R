#### packages ####
library(tidyverse)
library(brms)
library(loo)
library(jagsUI)
library(cowplot)
library(bayestestR)
library(ggeffects)

#### data ####
breed_df <- read.csv("new breeding/breed_dates_success.csv") # breeding success and lay dates
cov_df <- read.csv("new breeding/enviro_cov_breeding.csv")   # environmental covariates
cov_df$medianlay <- breed_df$medianlay
diet_df <- read.csv("new breeding/diet_data.csv")

# scale covariates
cov_z <- cov_df %>% mutate_at(c("spring_air", "spring_sst", "winter_air", "winter_sst",
                               "winter_nao", "chickrear_sst"), ~ (scale(.) %>% as.vector))

# create new dataframe for analysis of timing of breeding, combining covariates
# with median lay date
phen_df <- cov_z
phen_df$medianlay <- breed_df$medianlay

#### temporal trend in lay date ####
# estimate whether guillemots have advanced median lay date over time

# model
brm_time <- brm(
  medianlay ~ year,
  data = phen_df, 
  iter = 6000, 
  warmup = 2000
)

summary(brm_time)

# use function from bayestestR to obtain probability of direction (PD)
describe_posterior(brm_time, centrality = "mean")

# get model predicted estimates of median lay date to plot model results
predicted_time <- predict_response(brm_time, terms = c("year"))

# plot trend in median lay date over time
(laytrend <- ggplot(predicted_time, aes(x = x, y = predicted)) +
    geom_line(data = phen_df, aes(x = year, y = medianlay), col = "grey40", lwd = 0.3) +
    geom_line(col = "#56B4E9") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#56B4E9", col = NA, alpha = 0.3) +
    labs(x = "Year", y = "Median lay date", title = "") +
    theme_bw() +
    scale_y_continuous(limits = c(119, 143), breaks = c(115, 120, 125, 130, 135, 140, 145)) +
    scale_x_continuous(limits = c(1991, 2026), breaks = c(1990, 2000, 2010, 2020)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8)))

#### estimate inclusion probability ####
# estimate inclusion probability of spring and winter sea and air temperatures and 
# winter NAO to identify strongest predictor of timing of breeding
# use jags

# write model file in jags
cat(file = "laymodelinc.txt", "
model {
  
  # priors
  alpha ~ dnorm(0, 0.001)
  tau ~ dgamma(1, 0.001)
  omega ~ dbeta(1,1) # prior distribution for the number of included variables (let model estimate)

  # likelihood
  for (i in 1:N) { # for every row (occasion)
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + inprod(X[i,], beta)
  } 
  
  for (j in 1:p) { # for every covariate (column)
    gamma[j] ~ dbern(omega) # estimate inclusion probability
    betaT[j] ~ dnorm(0, 0.001) # prior for slope
    beta[j] <- gamma[j] * betaT[j]
  }

}
")

# covariate matrix
X <- cov_z
X <- X %>% select(2:6)
X <- as.matrix(X)

# create jags data object
jags.data <- list(y = phen_df$medianlay, X = X, N = nrow(X), p = ncol(X))
str(jags.data)

# specify monitors, initial values, iterations, burn-in and chains
params <- c("alpha", "beta", "gamma", "omega", "tau")
inits <- function() {list(alpha = rnorm(1))}
ni <- 30000; nb <- 10000; nc <- 4; nt <-5

# run
inclusionprob <- jags(jags.data, inits, params, "laymodelinc.txt", n.iter = ni, n.burnin = nb, n.chains = nc, 
               n.thin = nt, parallel = T)

print(inclusionprob)
# spring SST has the highest inclusion probability - use this in models.

#### trend in spring SST ####

# model
brm_spring <- brm(
  spring_sst ~ year,
  data = phen_df, 
  iter = 6000, 
  warmup = 2000
)

summary(brm_spring)
describe_posterior(brm_spring, centrality = "mean")

# get model predicted estimates of spring SST to plot model results
predicted_spring <- predict_response(brm_spring, terms = c("year"))

# convert back from standardised SST to plot on normal scale
mu <- mean(cov_df$spring_sst)
stdev <- sd(cov_df$spring_sst)
predicted_spring$sst_raw <- (predicted_spring$predicted*stdev) + mu
predicted_spring$low_raw <- (predicted_spring$conf.low*stdev) + mu
predicted_spring$high_raw <- (predicted_spring$conf.high*stdev) + mu

# plot results of year model with non-standardised x axis
(springtrend <- ggplot(predicted_spring, aes(x = x, y = sst_raw)) +
    geom_line(data = cov_df, aes(x = year, y = spring_sst), col = "grey40", lwd = 0.3) +
    geom_line(col = "#56B4E9") +
    geom_ribbon(aes(ymin = low_raw, ymax = high_raw), fill = "#56B4E9", col = NA, alpha = 0.3) +
    labs(x = "Year", y = "Spring SST (Celsius)", title = "") +
    theme_bw() +
    scale_x_continuous(limits = c(1991, 2026), breaks = c(1995, 2005, 2015, 2025)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8)))

#### impact of spring SST on median lay date ####
# has median lay date advanced in response to warming spring SST?

# model
brm_sst_lay <- brm(
  medianlay ~ spring_sst,
  data = phen_df, 
  iter = 6000, 
  warmup = 2000)

summary(brm_sst_lay)
describe_posterior(brm_sst_lay, centrality = "mean")

# get model predicted estimates of spring SST to plot model results
predicted_sst_lay <- predict_response(brm_sst_lay, terms = c("spring_sst"))

# convert back from standardised SST to plot on normal scale
predicted_sst_lay$sst_raw <- (predicted_sst_lay$x*stdev) + mu

# plot results of spring sst model with non-standardised x axis
(ssttrend <- ggplot(predicted_sst_lay, aes(x = sst_raw, y = predicted)) +
    geom_point(data = cov_df, aes(x = spring_sst, y = medianlay), col = "grey40", alpha = 0.5, size = 1) +
    geom_line(col = "#56B4E9") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#56B4E9", col = NA, alpha = 0.3) +
    labs(x = "Spring SST (celsius)", y = "Median lay date", title = "") +
    theme_bw() +
    scale_y_continuous(limits = c(115, 145), breaks = c(115, 120, 125, 130, 135, 140, 145)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8)))

# test relationship between lay date and detrended SST
phen_df$spring_d <- resid(lm(spring_sst ~ year, data = phen_df))

brm_sstd <- brm(
  medianlay ~ spring_d,
  data = phen_df, 
  iter = 6000, 
  warmup = 2000)

summary(brm_sstd)
describe_posterior(brm_sstd) # no support for any effect of detrended spring SST

#### investigate relationship between breeding success and median lay date ####
# has change in median lay date affected mean breeding success?

# limit success analysis to 1995 onwards
success_df <- breed_df %>% slice(4:35)

# model
brm_success_lay <- brm(success | trials(pairs) ~ medianlay, 
             family = binomial("logit"), 
             data = success_df
)

summary(brm_success_lay)
describe_posterior(brm_success_lay)

# plot
# Create a sequence of predictor values
newdata <- data.frame(medianlay = seq(min(success_df$medianlay), max(success_df$medianlay), length.out = 100), pairs = 1)

# Get predicted values
fits <- fitted(brm_success_lay, newdata = newdata, scale = "response")

# combine
predicted_lay_success <- cbind(newdata, fits)

# plot
(successdate <- ggplot(success_df, aes(x = medianlay, y = success / pairs)) +
    geom_point(col = "grey40", size = 1, alpha = 0.5) +
    geom_line(data = predicted_lay_success, aes(x = medianlay, y = Estimate), color = "#56B4E9", lty = 2) +
    labs(y = "Breeding success probability", x = "Median lay date") +
    ylim(0.6, 0.92) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8)))

#### trend in breeding success over time ####

# model
brm_success_time <- brm(success | trials(pairs) ~ year, 
             family = binomial("logit"), 
             data = success_df
)

summary(brm_success_time)
describe_posterior(brm_success_time)

# plot
# Create a sequence of predictor values
newdata2 <- data.frame(year = seq(min(success_df$year),  max(success_df$year), length.out = 100), pairs = 1)

# Get predicted values
fits2 <- fitted(brm_success_time, newdata = newdata2, scale = "response")

# combine
predicted_success_time <- cbind(newdata2, fits2)

# plot
(successtrend <- ggplot(success_df, aes(x = year, y = success / pairs)) +
    geom_line(col = "grey40", lwd = 0.3) +
    geom_line(data = predicted_success_time, aes(x = year, y = Estimate), color = "#56B4E9") +
    geom_ribbon(data = predicted_success_time, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), alpha = 0.3, fill = "#56B4E9") +
    labs(y = "Breeding success probability", x = "Year") +
    ylim(0.6, 0.92) +
    scale_x_continuous(labels = c(1995, 2005, 2015, 2025), 
                       breaks = c(1995, 2005, 2015, 2025)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8)))

#### create figure 2 ####
(comb1 <- plot_grid(laytrend, springtrend, ssttrend, nrow = 1, labels = c("A", "B", "C"), label_size = 8, 
                    label_x = 0.05, label_y = 0.9, align = "hv"))

(comb1.1 <- plot_grid(NA, comb1, nrow = 2, rel_heights = c(-0.05, 1)))

(comb2.1 <- plot_grid(NA, successdate, successtrend, NA, nrow = 1, labels = c("", "D", "E", ""), label_size = 8, 
                      label_x = 0.05, label_y = 0.97, rel_widths = c(0.5, 1, 1, 0.5), align = "hv"))

(all <- plot_grid(comb1.1, comb2.1, ncol = 1))

(all_final <- ggdraw() +
    draw_plot(all) +
    theme(panel.background = element_rect(fill = "white")) +
    panel_border(colour = "white"))

#save_plot("breeding success/all_breeding_results.pdf", all_final, base_height = 5, base_width = 7, dpi = 600)
#save_plot("breeding success/all_breeding_results.png", all_final, base_height = 5, base_width = 7, dpi = 600)

#### diet analysis ####

# create stacked bar chart of diet composition over time
diet_red <- diet_df %>% select(1:5) # leave out total feeds column
diet_long <- pivot_longer(diet_red, !year, names_to = "taxa", values_to = "count")

# re-order factor levels - smallest proportions at top of chart
diet_long$taxa <- as.factor(diet_long$taxa)
levels(diet_long$taxa)
diet_long$taxa <- factor(diet_long$taxa, levels = c("squid", "gadids", "sandeels", "clupeids"))

str(diet_long)

# stacked bar chart of proportion of species per year
(bar <- ggplot(diet_long, aes(fill=taxa, y=count, x=year)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=c("#F0E442", "#009E73", "#56B4E9", "#0072B2")) +
    ylab("Proportion feeds") + xlab("Year") +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1), expand=c(0,0)) +
    scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020, 2025), expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 7), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.key.spacing.x = unit(0.1, "cm"), 
          legend.key.size = unit(0.2, "cm"), 
          legend.position = "bottom", 
          legend.margin = margin(0,0,0,0), 
          legend.box.margin = margin(-6,-2,-2,-2)))

#### model chick diet composition ####

# reshape diet data for use in brms models - one feed per row
prey <- diet_red %>% select(1:4) # remove squid - only observed towards end of study

diet_reshape <- prey %>% pivot_longer(cols = -year, names_to = "species", values_to = "count") %>%
  uncount(count)
diet_reshape$year <- as.numeric(diet_reshape$year)
diet_reshape$species <- as.factor(diet_reshape$species)

# scale year - will be combined with SST in models
diet_reshape$year_z <- scale(diet_reshape$year)

# add chick-rearing SST
phen_df$chickrear_d <- resid(lm(chickrear_sst ~ year, data = phen_df)) # detrended
chickrearsst <- phen_df %>% select(1, 7, 10) # select just chick rearing SST to join by year

diet_reshape <- diet_reshape |> left_join(chickrearsst)

### models

# diet ~ year
brm_diet_time <- brm(species ~ year_z,
            data = diet_reshape, 
            family = categorical(), 
            chains = 4, 
            cores = 4,
            iter = 6000, 
            warmup = 2000)

summary(brm_diet_time)
describe_posterior(brm_diet_time)

# diet ~ year + sst
brm_diet_time_sst <- brm(species ~ year_z + chickrear_sst,
                     data = diet_reshape, 
                     family = categorical(), 
                     chains = 4, 
                     cores = 4,
                     iter = 6000, 
                     warmup = 2000)

summary(brm_diet_time_sst)
describe_posterior(brm_diet_time_sst)

# diet ~ year + sst_d
brm_diet_time_sstd <- brm(species ~ year_z + chickrear_d,
                         data = diet_reshape, 
                         family = categorical(), 
                         chains = 4, 
                         cores = 4,
                         iter = 6000, 
                         warmup = 2000)

summary(brm_diet_time_sstd)
describe_posterior(brm_diet_time_sstd)

# diet ~ year + sst + year*sst
brm_diet_timexsst <- brm(species ~ year_z + chickrear_sst + year_z*chickrear_sst,
                          data = diet_reshape, 
                          family = categorical(), 
                          chains = 4, 
                          cores = 4,
                          iter = 6000, 
                          warmup = 2000)

summary(brm_diet_timexsst)
describe_posterior(brm_diet_timexsst)

# diet ~ year + sst_d + year*sst_d
brm_diet_timexsstd <- brm(species ~ year_z + chickrear_d + year_z*chickrear_d,
                         data = diet_reshape, 
                         family = categorical(), 
                         chains = 4, 
                         cores = 4,
                         iter = 6000, 
                         warmup = 2000)

summary(brm_diet_timexsstd)
describe_posterior(brm_diet_timexsstd)

# compare models using leave-one-out cross validation
loo1 <- loo(brm_diet_time)
loo2 <- loo(brm_diet_time_sst)
loo3 <- loo(brm_diet_time_sstd)
loo4 <- loo(brm_diet_timexsst)
loo5 <- loo(brm_diet_timexsstd)

loo_compare(loo1, loo2, loo3, loo4, loo5)
# interaction between year and detrended chick-rearing is best supported model

# plot results of year*detrended SST model

# get predicted proportions of fish in diet
predicted_diet <- predict_response(brm_diet_timexsstd, terms = c("year_z", "chickrear_d"))

# create new group column to name interaction levels with helpful labels
predicted_diet$group2 <- predicted_diet$group
levels(predicted_diet$group2) <- c("-1 SD", "mean", "+1 SD")
colnames(predicted_diet)[7] <- "SST (detrended)"

# get year back on normal scale
mu <- mean(diet_reshape$year)
stdev <- sd(diet_reshape$year)

predicted_diet$year <- (predicted_diet$x*stdev) + mu

# get proportions of each species to plot points
diet_df$p_clupeids <- diet_df$clupeids / diet_df$feeds
diet_df$p_gadids <- diet_df$gadids / diet_df$feeds
diet_df$p_sandeels <- diet_df$sandeels / diet_df$feeds

diet_prop <- diet_df %>% select(1,7:9)
colnames(diet_prop) <- c("year", "clupeids", "gadids", "sandeels")

prop_long <- pivot_longer(diet_prop, !year, names_to = "taxa", values_to = "prop")

colnames(prop_long)[2] <- "response.level"

(intplot <- ggplot(predicted_diet, aes(x = year, y = predicted)) +
    geom_point(data = prop_long, aes(x = year, y = prop, group = response.level), size = 1, col = "grey40", alpha = 0.5) +
    geom_line(aes(x = year, y = predicted, colour = `SST (detrended)`)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = `SST (detrended)`), alpha = 0.3, color = NA) +
    ylab("Probability present in diet") + xlab("Year") +
    scale_colour_manual(values = c("#0072B2", "#009E73", "#56B4E9")) +
    scale_fill_manual(values = c("#0072B2", "#009E73", "#56B4E9")) +
    facet_wrap(vars(response.level)) +
    scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020, 2025)) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8), 
          axis.text = element_text(size = 7), 
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8), 
          legend.key.spacing.y = unit(0.0001, "cm"), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "bottom", 
          legend.margin = margin(0,0,0,0), 
          legend.box.margin = margin(-6,-2,-2,-2)))

#### investigate relationships between diet and breeding success ####

# combine breeding and diet data
diet_success <- success_df %>% left_join(diet_prop, by = "year")

# model success as a function of proportion clupeids in diet
# use clupeids as they form largest proportion of diet, and proportions of other 
# species correlated 
brm_diet_success <- brm(success | trials(pairs) ~ clupeids, 
             family = binomial("logit"), 
             data = diet_success
)

summary(brm_diet_success)
describe_posterior(brm_diet_success, centrality = "mean")
# clupeids have positive impact on breeding success

# plot
# Create a sequence of predictor values
newdata3 <- data.frame(clupeids = seq(min(diet_success$clupeids), max(diet_success$clupeids),  length.out = 100), pairs = 1)

# Get fitted values (posterior mean + intervals)
fits3 <- fitted(brm_diet_success, newdata = newdata3, scale = "response")

# combine
predicted_diet_success <- cbind(newdata3, fits3)

# plot
(clupeid_success <- ggplot(diet_success, aes(x = clupeids, y = success / pairs)) +
    geom_point(col = "grey40", size = 1, alpha = 0.5) +
    geom_line(data = predicted_diet_success, aes(x = clupeids, y = Estimate), color = "#56B4E9") +
    geom_ribbon(data = predicted_diet_success, aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), alpha = 0.3, fill = "#56B4E9") +
    labs(y = "Breeding success", x = "Proportion clupeids") +
    ylim(0.6, 0.93) +
    theme_bw() +
    theme(axis.text = element_text(size = 7), 
          axis.title = element_text(size = 9)))

# plot all together
(p1 <- plot_grid(bar, clupeid_success, nrow = 1, rel_widths = c(1.2, 0.9), labels = c("A", "B"), label_size = 9))
(p2 <- plot_grid(intplot, labels = c("C"), label_size = 9))

(p12 <- plot_grid(p1, p2, nrow = 2))

(p12_final <- ggdraw() +
    draw_plot(p12) +
    theme(panel.background = element_rect(fill = "white")) +
    panel_border(colour = "white"))

#save_plot("feeding/all_results_2026.pdf", p12_final, base_height = 6, base_width = 7, dpi = 600)
#save_plot("feeding/all_results_2026.png", p12_final, base_height = 6, base_width = 7, dpi = 600)