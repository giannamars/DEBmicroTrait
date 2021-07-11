library(dplyr)
library(ggplot2)

df_plt <- read.csv("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/batch_model_isolates_minerals.csv")
df_molli <- filter(df_plt, soil == 'Mollisols' & name == 'glucose' | name == 'alanine' | name == 'salicylic acid')


ggplot(df_molli, aes(x = BGE, y = Dads)) + geom_point(aes(color = factor(name))) 
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
  geom_smooth(method="lm", se= F, aes(colour = name, group = name)) +
  labs(x = "Bacterial growth efficiency [-]",y = "MAOM [\u03bcM]",color = "") + 
  theme_bw()

ggsave(filename="/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/plots/batch_model_isolates_minerals_1.png", plot = p)
