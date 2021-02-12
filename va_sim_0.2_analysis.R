### Analysis for 2-Panel Figure


setwd("~/vaccine_alloc/Results")


#100 vaccinations
rm(list=ls())
result_VHC_100 = read_csv("sim_result_post_vac_VHC_1cperweek_block1_v100_0205.csv")
result_VCD_100 = read_csv("sim_result_post_vac_VCD_v100_1cpw_0205.csv")
result_VCN_100 = read_csv("sim_result_post_vac_VCN_v100_cpw1_0205.csv")

q_VCD = result_VCD_100 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC_100 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN_100 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))

t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))
t_VCD = t_VCD[order(t_VCD$round),]

t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),]

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),]

t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
  scale_color_manual(name = "Strategy", values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                                    "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
                                                    "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 0.5)) +
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies (10% vaccination rate)") +
  theme(legend.position = "bottom", legend.direction = "vertical")

#200 vaccinations
rm(list=ls())

result_VHC_200 = read_csv("sim_result_post_vac_VHC_1cperweek_block1_v200_0205.csv")
result_VCD_200 = read_csv("sim_result_post_vac_VCD_v200_1cpw_0205.csv")
result_VCN_200 = read_csv("sim_result_post_vac_VCN_v200_cpw1_0205.csv")

q_VCD = result_VCD_200 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC_200 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN_200 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))
t_VCD = t_VCD[order(t_VCD$round),]

t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),]

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),]

t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
  scale_color_manual(name = "Strategy", values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                                   "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
                                                   "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
  xlab(" ") + scale_y_continuous("Cumulative incidence", limits = c(0, 0.5)) +
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies (20% vaccination rate)") +
  theme(legend.position = "bottom", legend.direction = "vertical")

#500 vaccinations
setwd("~/vaccine_alloc/Results")

library(ggpubr)
library(tidyverse)
library(igraph)
library(Matrix)
library(data.table)
library(readr)

result_VCD_500 = read_csv("sim_result_post_vac_VCD_v500_0209.csv")
result_VCN_500 = read_csv("sim_result_post_vac_VCN_v500_cpw1_0205.csv")
result_VHC_500 = read_csv("sim_result_post_vac_VHC_v500_0209.csv")


q_VCD = result_VCD_500 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC_500 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN_500 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))



t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))
t_VCD = t_VCD[order(t_VCD$round),]
t_VCD = t_VCD %>% filter(between(t_VCD$round, 21, 300))
head(t_VCD)
t_VCD %>% ggplot(aes(x = round, y = q50)) + geom_line()


t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),]
t_VHC = t_VHC %>% filter(between(t_VHC$round, 30, 300))
t_VHC %>% ggplot(aes(x = round, y = q50)) + geom_line()

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),]
t_VCN = t_VCN %>% filter(between(t_VCN$round, 25, 300))



g_500_1 = t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
  scale_color_manual(name = "Strategy", values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                                   "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
                                                   "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
  xlab(" ") + scale_y_continuous("Cumulative incidence \n", limits = c(0, 0.5)) +
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies: N = 1000, iterations = 100, \n vaccinations = 500, 1 new case per week") +
  theme(legend.position = "none", legend.direction = "vertical")


q_VCD_n = result_VCD_500 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC_n = result_VHC_500 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN_n = result_VCN_500 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

t_VCD_n = as.data.frame(t(q_VCD_n))
colnames(t_VCD_n) = c("q25", "q50", "q75")
t_VCD_n$round = parse_number(row.names(t_VCD_n))
t_VCD_n = t_VCD_n[order(t_VCD_n$round),]
t_VCD_n = t_VCD_n %>% filter(between(t_VCD_n$round, 21, 300))

t_VHC_n = as.data.frame(t(q_VHC_n))
colnames(t_VHC_n) = c("q25", "q50", "q75")
t_VHC_n$round = parse_number(row.names(t_VHC_n))
t_VHC_n = t_VHC_n[order(t_VHC_n$round),]
t_VHC_n = t_VHC_n %>% filter(between(t_VHC_n$round, 30, 300))

t_VCN_n = as.data.frame(t(q_VCN_n))
colnames(t_VCN_n) = c("q25", "q50", "q75")
t_VCN_n$round = parse_number(row.names(t_VCN_n))
t_VCN_n = t_VCN_n[order(t_VCN_n$round),]
t_VCN_n = t_VCN_n %>% filter(between(t_VCN_n$round, 28, 300))


g_500_2 = t_VHC_n %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
  geom_ribbon(data = t_VCD_n, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD_n, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
  geom_ribbon(data = t_VCN_n, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN_n, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
  scale_color_manual(name = "Strategy", values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                                   "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
                                                   "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
  xlab("Day") + scale_y_continuous("Proportion of new cases \n", limits = c(0, 0.005)) +
  theme_minimal() +
  theme(legend.position = "none")


ggarrange(g_500_1, g_500_2, nrow = 1, ncol = 2, align = "hv", legend = "bottom", common.legend = TRUE)

### Calculate CIRs

head(t_VHC)
head(t_VCN)

CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))

