rm(list = ls())
library(tidyverse); library(data.table); library(lme4); library(lmerTest); library(hausekeep); library(sjPlot)

# ggplot theme
mytheme <- theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text.x = element_text(size = 18, colour = 'black'),
          axis.text.y = element_text(size = 18, colour = 'black'),
          plot.title = element_text(size = 18, hjust = 0.5),
          text = element_text(size = 18), panel.spacing.x = unit(1, "lines"),
          legend.key = element_blank(), legend.key.height = unit(2, 'line'))

# read data
dt <- fread("./Data_analysis/pupil_coefs.csv")

####  fit model to cue ####

m_cue1 <- lmer(pupil_responseZ ~ Efficacy * Reward + Trial + (1 | subject), dt[event == "cue"])
summaryh(m_cue1)

m_cue2 <- lmer(pupil_responseZ ~ Efficacy * Reward + Trial + (1 + Efficacy + Reward + Trial | subject),
               control = lmerControl(optimizer ='bobyqa', optCtrl = list(maxfun = 2e4)),
               dt[event == "cue"])
summaryh(m_cue2)
tab_model(m_cue2)

#### fit model to Stroop target ####

# congruency as dummy variable
m_target1 <- lmer(pupil_responseZ ~ congruency + Trial + (1 | subject), dt[event == "stimulus"])
summaryh(m_target1)

# congruency as linear variable
m_target2 <- lmer(pupil_responseZ ~ congruencyC + Trial + (1 | subject), dt[event == "stimulus"])
summaryh(m_target2)

m_target3 <- lmer(pupil_responseZ ~ congruencyC + Trial + (1 + congruency + Trial | subject),
               control = lmerControl(optimizer ='bobyqa', optCtrl = list(maxfun = 2e4)),
               dt[event == "stimulus"])
summaryh(m_target3)

#### fit model to feedback ####
m_feedback1 <- lmer(pupil_responseZ ~ isRewarded + Trial + (1 | subject), dt[event == "feedback"])
summaryh(m_feedback1)

m_feedback2 <- lmer(pupil_responseZ ~ isRewarded + Trial + (1 + isRewarded + Trial | subject),
               control = lmerControl(optimizer ='bobyqa', optCtrl = list(maxfun = 2e4)),
               dt[event == "feedback"])
summaryh(m_feedback2)

#### plot scaled erlang gamma function ####

erlang <- function(t, tmax, n) {
    return(t ^ n * exp(-(n * t) / tmax))
}

tmax <- unique(dt$tmax) / 1000 # tmax in seconds
times <- seq(0, 5, by = 0.05)
plot(times, erlang(times, tmax, 10.1))

# compute subject average cue coef
avg <- dt[event == "cue",
          .(pupil_responseZ = mean(pupil_responseZ, na.rm = T)),
          by = .(subject, Efficacy, Reward)]

# compute mean and CI
ci <- seWithin(data = avg, measurevar = c("pupil_responseZ"), withinvars = c("Efficacy", "Reward"), idvar = "subject")
ci

labels <- c("Low efficacy, low reward",
            "Low efficacy, high reward",
            "High efficacy, low reward",
            "High efficacy, high reward")

# create scaled pupil responses for plotting
fx <- erlang(times, tmax, 10.1)

fx <- fx / max(fx) # rescale for plotting purposes
dt1 <- lapply(1:length(labels), function(x) { # for each condition
    temp_dt <- data.table(times, response = fx * ci$pupil_responseZ[x], condition = labels[x]) # scaled response
    temp_dt[times == tmax, ci := ci$ci[x]] # errorbar at peak
    return(temp_dt)
})
dt1 <- rbindlist(dt1)
dt1$Condition <- factor(dt1$condition, labels)

ggplot(dt1, aes(times, response, col = Condition)) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = response - ci, ymax = response + ci, col = condition), size = 0.8, width = 0.05) +
    mytheme +
    theme(legend.title = element_blank(), legend.position = c(.80, .80)) +
    labs(y = "Erlang gamma pupil response function", x = "Time")
ggsave("pupil_response.pdf", width = 8, height = 5)
