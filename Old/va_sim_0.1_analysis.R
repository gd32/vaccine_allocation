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

#avoid local social network of previously infected?
#goal @ day 300 - CI = 0.4
# how many vaccines for VHC = 0.6?

t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = 'VHC')) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.5)+
  geom_line(data= t_VCD, aes(y=q50, color = 'VCD')) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.5)+
  geom_line(data= t_VCN, aes(y=q50, color = "VCN")) +
  scale_color_manual(values = c("VHC" = "black", "VCD" = "red", "VCN" = "blue")) +
  xlab(" ") + scale_y_continuous("Cumulative incidence", limits = c(0, 1))

# sim_result_total_pv[, order(colnames(sim_result_total_pv))]

## Changing vaccination level
setwd("~/vaccine_alloc/Results")

result_VHC100 = read_csv("sim_result_post_vac_sym_only_VHC_v100_0205.csv")
result_VHC200 = read_csv("sim_result_post_vac_sym_only_VHC_v200_0205.csv")
result_VHC500 = read_csv("sim_result_post_vac_sym_only_VHC_v500_0205.csv")

q_VHC100 = result_VHC100[, order(colnames(result_VHC100))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                               'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                               'q75' = quantile(., p = 0.75, na.rm = T)))


t_VHC100 = as.data.frame(t(q_VHC100))
colnames(t_VHC100) = c("q25", "q50", "q75")
t_VHC100$round = parse_number(row.names(t_VHC100))
t_VHC100 = t_VHC100[order(t_VHC100$round),]

q_VHC200 = result_VHC200[, order(colnames(result_VHC200))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                          'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                          'q75' = quantile(., p = 0.75, na.rm = T)))


t_VHC200 = as.data.frame(t(q_VHC200))
colnames(t_VHC200) = c("q25", "q50", "q75")
t_VHC200$round = parse_number(row.names(t_VHC200))
t_VHC200 = t_VHC200[order(t_VHC200$round),]

q_VHC500 = result_VHC500[, order(colnames(result_VHC500))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                          'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                          'q75' = quantile(., p = 0.75, na.rm = T)))

t_VHC500 = as.data.frame(t(q_VHC500))
colnames(t_VHC500) = c("q25", "q50", "q75")
t_VHC500$round = parse_number(row.names(t_VHC500))
t_VHC500 = t_VHC500[order(t_VHC500$round),]


t_VHC100 %>% ggplot(aes(x = round)) + 
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgreen", alpha = 0.4)+
  geom_line(aes(y=q50, color = '100')) +
  geom_ribbon(data = t_VHC200, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.4)+
  geom_line(data= t_VHC200, aes(y=q50, color = '200')) +
  geom_ribbon(data = t_VHC500, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.4)+
  geom_line(data= t_VHC500, aes(y=q50, color = "500")) +
  scale_color_manual(name = "Number of Vaccinations", values = c("100" = "purple", "200" = "red", "500" = "blue")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 1)) + 
  theme_minimal()
  
## 1 Case per week, 2 per vaccination block
result_VHC_1cb2 = read_csv("sim_result_post_vac_1cperweek_block2_0205.csv")

q_VHC_1cb2 = result_VHC_1cb2[, order(colnames(result_VHC_1cb2))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                          'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                          'q75' = quantile(., p = 0.75, na.rm = T)))


t_VHC_1cb2 = as.data.frame(t(q_VHC_1cb2))
colnames(t_VHC_1cb2) = c("q25", "q50", "q75")
t_VHC_1cb2$round = parse_number(row.names(t_VHC_1cb2))
t_VHC_1cb2 = t_VHC_1cb2[order(t_VHC_1cb2$round),]

t_VHC_1cb2 %>% ggplot(aes(x = round)) + 
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightblue", alpha = 0.4)+
  geom_line(aes(y=q50, color = "CI")) +
  scale_color_manual(values = c("CI" = "blue")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 1)) + 
  theme_minimal() + ggtitle("1 case per week, centrality calculated in blocks of 2")

tail(t_VHC_1cb2)

## 1 Case per week, 1 per vaccination block
result_VHC_1cb1 = read_csv("sim_result_post_vac_1cperweek_block1_0205.csv")

q_VHC_1cb1 = result_VHC_1cb1[, order(colnames(result_VHC_1cb1))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                                'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                                'q75' = quantile(., p = 0.75, na.rm = T)))
t_VHC_1cb1 = as.data.frame(t(q_VHC_1cb1))
colnames(t_VHC_1cb1) = c("q25", "q50", "q75")
t_VHC_1cb1$round = parse_number(row.names(t_VHC_1cb1))
t_VHC_1cb1 = t_VHC_1cb1[order(t_VHC_1cb1$round),]

t_VHC_1cb1 %>% ggplot(aes(x = round)) + 
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lavender", alpha = 0.6)+
  geom_line(aes(y=q50, color = "Blocks of 1")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 1)) + 
  theme_minimal() +
  geom_ribbon(data = t_VHC_1cb2, aes(ymin = q25, ymax = q75), fill = "lightgreen", alpha = 0.6)+
  geom_line(data = t_VHC_1cb2, aes(y=q50, color = "Blocks of 2")) +
  scale_color_manual(values = c("Blocks of 1" = "purple", "Blocks of 2" = "green")) + 
  ggtitle("1 case per week, centrality calculation modified")

#Compare 3 strategies - 1 case per week
result_VHC_1cb1 = read_csv("sim_result_post_vac_1cperweek_block1_0205.csv")

q_VHC_1cb1 = result_VHC_1cb1[, order(colnames(result_VHC_1cb1))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                                'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                                'q75' = quantile(., p = 0.75, na.rm = T)))
t_VHC_1cb1 = as.data.frame(t(q_VHC_1cb1))
colnames(t_VHC_1cb1) = c("q25", "q50", "q75")
t_VHC_1cb1$round = parse_number(row.names(t_VHC_1cb1))
t_VHC_1cb1 = t_VHC_1cb1[order(t_VHC_1cb1$round),]

result_VCN_1cpw = read_csv("sim_result_post_vac_VCN_1cpw_0205.csv")

q_VCN_1cpw = result_VCN_1cpw[, order(colnames(result_VCN_1cpw))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                                'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                                'q75' = quantile(., p = 0.75, na.rm = T)))
t_VCN_1cpw = as.data.frame(t(q_VCN_1cpw))
colnames(t_VCN_1cpw) = c("q25", "q50", "q75")
t_VCN_1cpw$round = parse_number(row.names(t_VCN_1cpw))
t_VCN_1cpw = t_VCN_1cpw[order(t_VCN_1cpw$round),]

result_VCD_1cpw = read_csv("sim_result_post_vac_VCN_cpw1_0205.csv")

q_VCD_1cpw = result_VCD_1cpw[, order(colnames(result_VCD_1cpw))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                                'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                              'q75' = quantile(., p = 0.75, na.rm = T)))
t_VCD_1cpw = as.data.frame(t(q_VCD_1cpw))
colnames(t_VCD_1cpw) = c("q25", "q50", "q75")
t_VCD_1cpw$round = parse_number(row.names(t_VCD_1cpw))
t_VCD_1cpw = t_VCD_1cpw[order(t_VCD_1cpw$round),]


t_VHC_1cb1 %>% ggplot(aes(x = round)) + 
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "red", alpha = 0.2)+
  geom_line(aes(y=q50, color = 'VHC')) +
  geom_ribbon(data = t_VCN_1cpw, aes(ymin = q25, ymax=q75), fill = "lavender", alpha = 0.6)+
  geom_line(data= t_VCN_1cpw, aes(y=q50, color = 'VCN')) +
  geom_ribbon(data = t_VCD_1cpw, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.6)+
  geom_line(data= t_VCD_1cpw, aes(y=q50, color = "VCD")) +
  scale_color_manual(name = "Number of Vaccinations", values = c("VHC" = "red", "VCN" = "black", "VCD" = "blue")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 1)) + 
  theme_minimal() + ggtitle("50 Vaccinations, 1 case per week")

t_VCN_1cpw %>% filter(round == 50)
t_VCD_1cpw %>% filter(round == 50)


##10 starting infections, no cases per day
result_VCN_10seed = read_csv("sim_result_post_vac_sym_only_VCN_10seed_0205.csv")

q_VCN_10seed = result_VCN_10seed[, order(colnames(result_VCN_10seed))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                                'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                                'q75' = quantile(., p = 0.75, na.rm = T)))
t_VCN_10seed = as.data.frame(t(q_VCN_10seed))
colnames(t_VCN_10seed) = c("q25", "q50", "q75")
t_VCN_10seed$round = parse_number(row.names(t_VCN_10seed))
t_VCN_10seed = t_VCN_10seed[order(t_VCN_10seed$round),]

result_VCD_10seed = read_csv("sim_result_post_vac_VCD_10seed_0205.csv")

q_VCD_10seed = result_VCD_10seed[, order(colnames(result_VCD_10seed))] %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                                                                'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                                                                'q75' = quantile(., p = 0.75, na.rm = T)))
t_VCD_10seed = as.data.frame(t(q_VCD_10seed))
colnames(t_VCD_10seed) = c("q25", "q50", "q75")
t_VCD_10seed$round = parse_number(row.names(t_VCD_10seed))
t_VCD_10seed = t_VCD_10seed[order(t_VCD_10seed$round),]

tail(t_VCD_10seed)


t_VCN_10seed %>% ggplot(aes(x = round)) + 
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "red", alpha = 0.2)+
  geom_line(aes(y=q50, color = 'VCN')) +
  geom_ribbon(data = t_VCD_10seed, aes(ymin = q25, ymax=q75), fill = "lavender", alpha = 0.6)+
  geom_line(data= t_VCD_10seed, aes(y=q50, color = 'VCD')) +
  scale_color_manual(name = "Strategy", values = c("VCD" = "magenta", "VCN" = "red")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 1)) + 
  theme_minimal() + ggtitle("50 Vaccinations, 10 initial infections, no cases per week")
