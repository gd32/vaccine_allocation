### Analysis for 2-Panel 
library(ggpubr)
library(tidyverse)
library(Matrix)
library(data.table)
library(readr)


setwd("~/vaccine_alloc/Results")

#50 vaccinations
rm(list=ls())
result_VHC_50 = read_csv("sim_result_VHC_v50_0209.csv")
result_VCD_50 = read_csv("sim_result_VCD_v50_0209.csv")
result_VCN_50 = read_csv("sim_result_VCN_v50_0210.csv")

q_VCD = result_VCD_50 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC_50 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN_50 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))


CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))



t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))
t_VCD = t_VCD[order(t_VCD$round),] %>% filter(between(t_VCD$round, 30, 300))

t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),] %>% filter(between(t_VHC$round, 30, 300))

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),] %>% filter(between(t_VCN$round, 32, 300))

g_50_1 = t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
  scale_color_manual(name = "Strategy", values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                                   "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
                                                   "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 0.55)) +
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies: N = 1000, iterations = 100, \n vaccinations = 50, 1 new case per week") +
  theme(legend.position = "bottom", legend.direction = "vertical")



q_VCD_n = result_VCD_50 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC_n = result_VHC_50 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN_n = result_VCN_50 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))
t_VCD_n = as.data.frame(t(q_VCD_n))
colnames(t_VCD_n) = c("q25", "q50", "q75")
t_VCD_n$round = parse_number(row.names(t_VCD_n))
t_VCD_n = t_VCD_n[order(t_VCD_n$round),]
t_VCD_n = t_VCD_n %>% filter(between(t_VCD_n$round, 30, 300))

t_VHC_n = as.data.frame(t(q_VHC_n))
colnames(t_VHC_n) = c("q25", "q50", "q75")
t_VHC_n$round = parse_number(row.names(t_VHC_n))
t_VHC_n = t_VHC_n[order(t_VHC_n$round),]
t_VHC_n = t_VHC_n %>% filter(between(t_VHC_n$round, 30, 300))

t_VCN_n = as.data.frame(t(q_VCN_n))
colnames(t_VCN_n) = c("q25", "q50", "q75")
t_VCN_n$round = parse_number(row.names(t_VCN_n))
t_VCN_n = t_VCN_n[order(t_VCN_n$round),]
t_VCN_n = t_VCN_n %>% filter(between(t_VCN_n$round, 32, 300))


g_50_2 = t_VHC_n %>% ggplot(aes(x = round)) +
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

final_50 = ggarrange(g_50_1, g_50_2, nrow = 1, ncol = 2, align = "hv", legend = "bottom", common.legend = TRUE)
setwd("~/vaccine_alloc/Results/Images")
ggsave("v50.png")

#100 vaccinations
rm(list=ls())
setwd("~/vaccine_alloc/Results/")
result_VHC_100 = read_csv("sim_result_VHC_v100_0209.csv")
result_VCD_100 = read_csv("sim_result_VCD_v100_0209.csv")
result_VCN_100 = read_csv("sim_result_VCN_v100_0210.csv")

q_VCD = result_VCD_100 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC_100 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN_100 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                  'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                  'q75' = quantile(., p = 0.75, na.rm = T)))


CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))



t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))
t_VCD = t_VCD[order(t_VCD$round),] %>% filter(between(t_VCD$round, 30, 300))

t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),] %>% filter(between(t_VHC$round, 30, 300))

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),] %>% filter(between(t_VCN$round, 32, 300))


t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
  scale_color_manual(values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                                    "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
                                                    "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
  scale_x_continuous("Day after CI reaches 10%", limits = c(30, 293)) + scale_y_continuous("Cumulative incidence", limits = c(0, 0.5)) +
  annotate(geom = "text", x = 204.5, y = 0.155, label = "CIR for Empirical vs. reference: 1.0104, (1.0091, 1.0118)") +
  annotate(geom = "text", x = 204.5, y = 0.13, label = "CIR for Proactive vs. reference: 0.9865, (0.9813, 0.9917)") +
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies: N = 1000, iterations = 100, \n vaccinations = 100, 1 new case per week") +
  theme(legend.position = c(0.66, 0.15), legend.direction = "vertical", legend.title = element_blank())


CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))


setwd("~/vaccine_alloc/Results/Images")
ggsave("v100.png", width = 8, height = 5, dpi=400)



# ----
# 
# q_VCD_n = result_VCD_100 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
#                                                                         'q50' = quantile(., p = 0.5,  na.rm = T),
#                                                                         'q75' = quantile(., p = 0.75, na.rm = T)))
# 
# q_VHC_n = result_VHC_100 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
#                                                                         'q50' = quantile(., p = 0.5,  na.rm = T),
#                                                                         'q75' = quantile(., p = 0.75, na.rm = T)))
# 
# q_VCN_n = result_VCN_100 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
#                                                                         'q50' = quantile(., p = 0.5,  na.rm = T),
#                                                                         'q75' = quantile(., p = 0.75, na.rm = T)))
# t_VCD_n = as.data.frame(t(q_VCD_n))
# colnames(t_VCD_n) = c("q25", "q50", "q75")
# t_VCD_n$round = parse_number(row.names(t_VCD_n))
# t_VCD_n = t_VCD_n[order(t_VCD_n$round),]
# t_VCD_n = t_VCD_n %>% filter(between(t_VCD_n$round, 30, 300))
# 
# t_VHC_n = as.data.frame(t(q_VHC_n))
# colnames(t_VHC_n) = c("q25", "q50", "q75")
# t_VHC_n$round = parse_number(row.names(t_VHC_n))
# t_VHC_n = t_VHC_n[order(t_VHC_n$round),]
# t_VHC_n = t_VHC_n %>% filter(between(t_VHC_n$round, 30, 300))
# 
# t_VCN_n = as.data.frame(t(q_VCN_n))
# colnames(t_VCN_n) = c("q25", "q50", "q75")
# t_VCN_n$round = parse_number(row.names(t_VCN_n))
# t_VCN_n = t_VCN_n[order(t_VCN_n$round),]
# t_VCN_n = t_VCN_n %>% filter(between(t_VCN_n$round, 32, 300))
# 
# 
# g_100_2 = t_VHC_n %>% ggplot(aes(x = round)) +
#   geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2)+
#   geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)")) +
#   geom_ribbon(data = t_VCD_n, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
#   geom_line(data= t_VCD_n, aes(y=q50, color = "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)")) +
#   geom_ribbon(data = t_VCN_n, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
#   geom_line(data= t_VCN_n, aes(y=q50, color = "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)", )) +
#   scale_color_manual(name = "Strategy", values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
#                                                    "Vaccinating Close Contacts of Previously Diagnosed Individuals (Empirical)" = "magenta", 
#                                                    "Vaccinating Close Contacts of Previously Undiagnosed Individuals (Proactive)" = "cyan")) +
#   xlab("Day") + scale_y_continuous("Proportion of new cases \n", limits = c(0, 0.005)) +
#   theme_minimal() +
#   theme(legend.position = "none")




#200 vaccinations
rm(list=ls())
setwd("~/vaccine_alloc/Results")

result_VHC_250 = read_csv("sim_result_VHC_v250_0209.csv")
result_VCD_250 = read_csv("sim_result_VCD_v250_0209.csv")
result_VCN_250 = read_csv("sim_result_VCN_v250_0210.csv")

q_VCD = result_VCD_250 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC = result_VHC_250 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN = result_VCN_250 %>% select(contains("c")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                      'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                      'q75' = quantile(., p = 0.75, na.rm = T)))



t_VCD = as.data.frame(t(q_VCD))
colnames(t_VCD) = c("q25", "q50", "q75")
t_VCD$round = parse_number(row.names(t_VCD))
t_VCD = t_VCD[order(t_VCD$round),] %>% filter(between(t_VCD$round, 30, 300))

t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),] %>% filter(between(t_VHC$round, 30, 300))

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),] %>% filter(between(t_VCN$round, 32, 300))

CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))

g_250_1 = t_VHC %>% ggplot(aes(x = round)) +
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
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies: N = 1000, iterations = 100, \n vaccinations = 250, 1 new case per week") +
  theme(legend.position = "bottom", legend.direction = "vertical")

q_VCD_n = result_VCD_250 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))

q_VHC_n = result_VHC_250 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))

q_VCN_n = result_VCN_250 %>% select(contains("n")) %>% summarize_all(~c('q25' = quantile(., p = 0.25, na.rm = T),
                                                                        'q50' = quantile(., p = 0.5,  na.rm = T),
                                                                        'q75' = quantile(., p = 0.75, na.rm = T)))
t_VCD_n = as.data.frame(t(q_VCD_n))
colnames(t_VCD_n) = c("q25", "q50", "q75")
t_VCD_n$round = parse_number(row.names(t_VCD_n))
t_VCD_n = t_VCD_n[order(t_VCD_n$round),]
t_VCD_n = t_VCD_n %>% filter(between(t_VCD_n$round, 30, 300))

t_VHC_n = as.data.frame(t(q_VHC_n))
colnames(t_VHC_n) = c("q25", "q50", "q75")
t_VHC_n$round = parse_number(row.names(t_VHC_n))
t_VHC_n = t_VHC_n[order(t_VHC_n$round),]
t_VHC_n = t_VHC_n %>% filter(between(t_VHC_n$round, 30, 300))

t_VCN_n = as.data.frame(t(q_VCN_n))
colnames(t_VCN_n) = c("q25", "q50", "q75")
t_VCN_n$round = parse_number(row.names(t_VCN_n))
t_VCN_n = t_VCN_n[order(t_VCN_n$round),]
t_VCN_n = t_VCN_n %>% filter(between(t_VCN_n$round, 32, 300))


g_250_2 = t_VHC_n %>% ggplot(aes(x = round)) +
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

final_250 = ggarrange(g_250_1, g_250_2, nrow = 1, ncol = 2, align = "hv", legend = "bottom", common.legend = TRUE)
final_250
setwd("~/vaccine_alloc/Results/Images")
ggsave("v250.png")


#500 vaccinations
rm(list = ls())
setwd("~/vaccine_alloc/Results")

result_VHC_500 = read_csv("sim_result_post_vac_VHC_v500_0209.csv")
result_VCD_500 = read_csv("sim_result_post_vac_VCD_v500_0209.csv")
result_VCN_500 = read_csv("sim_result_VCN_v500_0210.csv")

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
t_VCD = t_VCD[order(t_VCD$round),] %>% filter(between(t_VCD$round, 30, 300))

t_VHC = as.data.frame(t(q_VHC))
colnames(t_VHC) = c("q25", "q50", "q75")
t_VHC$round = parse_number(row.names(t_VHC))
t_VHC = t_VHC[order(t_VHC$round),] %>% filter(between(t_VHC$round, 30, 300))

t_VCN = as.data.frame(t(q_VCN))
colnames(t_VCN) = c("q25", "q50", "q75")
t_VCN$round = parse_number(row.names(t_VCN))
t_VCN = t_VCN[order(t_VCN$round),] %>% filter(between(t_VCN$round, 32, 300))

CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))

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
  xlab("Day") + scale_y_continuous("Cumulative incidence", limits = c(0, 0.25)) +
  theme_minimal() + ggtitle("Comparison of Vaccination Strategies: N = 1000, iterations = 100, \n vaccinations = 500, 1 new case per week") +
  theme(legend.position = "bottom", legend.direction = "vertical")

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
t_VCD_n = t_VCD_n %>% filter(between(t_VCD_n$round, 30, 300))

t_VHC_n = as.data.frame(t(q_VHC_n))
colnames(t_VHC_n) = c("q25", "q50", "q75")
t_VHC_n$round = parse_number(row.names(t_VHC_n))
t_VHC_n = t_VHC_n[order(t_VHC_n$round),]
t_VHC_n = t_VHC_n %>% filter(between(t_VHC_n$round, 30, 300))

t_VCN_n = as.data.frame(t(q_VCN_n))
colnames(t_VCN_n) = c("q25", "q50", "q75")
t_VCN_n$round = parse_number(row.names(t_VCN_n))
t_VCN_n = t_VCN_n[order(t_VCN_n$round),]
t_VCN_n = t_VCN_n %>% filter(between(t_VCN_n$round, 32, 300))


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

final_500 = ggarrange(g_500_1, g_500_2, nrow = 1, ncol = 2, align = "hv", legend = "bottom", common.legend = TRUE)
final_500
setwd("~/vaccine_alloc/Results/Images")
ggsave("v500.png")

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

