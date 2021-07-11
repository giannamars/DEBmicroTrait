library(lme4)
library(lmerTest)
library(MuMIn)
library(sjstats)
library(dplyr)
library(broom)
library(AICcmodavg)


df <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/batch_model_isolates.csv")
df_pn <- filter(df, response == 'positive' | response == 'negative')
df_isolates <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv")
df_isolates_pn <- filter(df_isolates, Rhizosphere_response == 'positive' | Rhizosphere_response == 'negative')
df_isolates_pn <- rename(df_isolates_pn, species = Isolate)

df_out <- left_join(df_pn, df_isolates_pn)


mlm <- lmer(BGE ~ as.factor(response) + (1|species), data=df_pn)
anova(mlm)
effectsize::eta_squared(mlm, partial=TRUE)
r.squaredGLMM(mlm)
summary(mlm)
VarCorr(mlm)
AICc(mlm)


out.stats <- lmer(BGE ~ as.factor(species)|Genome_size, data = df_out)
summary(out.stats)


