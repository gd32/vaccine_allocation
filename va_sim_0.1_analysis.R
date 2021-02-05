library(ggpubr)
library(tidyverse)
library(igraph)
library(Matrix)
library(data.table)
library(readr)

rm(list=ls())
setwd("~/vaccine_alloc/Results")

result_VHC = read_csv("sim_result_post_vac_VHC_0205.csv")
result_VCD = read_csv("sim_result_post_vac_VCD_0205.csv")
result_VCN = read_csv("sim_result_post_vac_VCN_0205.csv")

func_75CI = function(x) { quantile(x, probs = c(0.25,0.5,0.75), na.rm = TRUE) }

dt 

q_VCD = result_VCD %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                     'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                     'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                          'q50' = quantile(., p = 0.5,  na.rm = T),
                                                          'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))


t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))

t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))

t_VHC %>% filter(round == 300) #0.865, (0.848, 0.88)
t_VCD %>% filter(round == 300) #0.881, (0.857, 0.896)
t_VCN %>% filter(round == 300) #0.87, (0.849, 0.896)

t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = 'VHC')) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.5)+
  geom_line(data= t_VCD, aes(y=q50, color = 'VCD')) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.5)+
  geom_line(data= t_VCN, aes(y=q50, color = "VCN")) +
  scale_color_manual(values = c("VHC" = "black", "VCD" = "red", "VCN" = "blue")) +
  xlab(" ") + scale_y_continuous("Cumulative incidence", limits = c(0, 1))
