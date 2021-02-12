library(tidyverse)
library(ggpubr)

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


t_VHC %>% ggplot(aes(x = round)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "lightgray", alpha = 0.2,)+
  geom_line(aes(y=q50, color = "Vaccinating the Highest Centrality Individuals (Reference)"), size = 1.4) +
  geom_ribbon(data = t_VCD, aes(ymin = q25, ymax=q75), fill = "lightpink", alpha = 0.2)+
  geom_line(data= t_VCD, aes(y=q50, color = "Vaccinating Close Contacts of \n Previously Diagnosed Individuals (Empirical) \n"), size = 1.4) +
  geom_ribbon(data = t_VCN, aes(ymin = q25, ymax=q75), fill = "lightblue", alpha = 0.2)+
  geom_line(data= t_VCN, aes(y=q50, color = "Vaccinating Close Contacts of \n Previously Undiagnosed Individuals (Proactive) \n", ), size = 1.4) +
  scale_color_manual(values = c("Vaccinating the Highest Centrality Individuals (Reference)" = "grey20", 
                                "Vaccinating Close Contacts of \n Previously Diagnosed Individuals (Empirical) \n" = "magenta", 
                                "Vaccinating Close Contacts of \n Previously Undiagnosed Individuals (Proactive) \n" = "cyan"),
                     guide = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous("\n Day after Cumulative Incidence Reaches 10%", limits = c(30, 300)) + 
  scale_y_continuous("Cumulative Incidence \n", limits = c(0, 0.5)) +
  theme_minimal() +
  theme(legend.position = c(0.55, 0.25), 
        legend.direction = "vertical", 
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "in"),
        axis.text = element_text(size = 40),
        axis.title = element_text(size = 40))


# annotate(geom = "text", x = 204.5, y = 0.155, label = "CIR for Empirical vs. reference: 1.0104, (1.0091, 1.0118)") +
#   annotate(geom = "text", x = 204.5, y = 0.13, label = "CIR for Proactive vs. reference: 0.9865, (0.9813, 0.9917)") +

CI_table = t_VHC %>% full_join(t_VCN, by = 'round', suffix = c("_VHC", "_VCN")) %>% 
  right_join(t_VCD) %>% column_to_rownames(var = "round") %>% rename(q25_VCD = q25, 
                                                                     q50_VCD = q50, 
                                                                     q75_VCD = q75)

CI_table %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC) %>% 
  summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
            ci_vcd = list(mean_cl_normal(CIR_VCD)))

CI_table %>% rownames_to_column(var = "round") %>% filter(round >= 30) %>% 
  mutate(CIR_VCN = q50_VCN/q50_VHC, CIR_VCD = q50_VCD/q50_VHC)
  # summarize(ci_vcn = list(mean_cl_normal(CIR_VCN)), 
  #           ci_vcd = list(mean_cl_normal(CIR_VCD)))

CI_table %>% rownames_to_column(var = "round") 

setwd("~/vaccine_alloc/Results/Images")
ggsave("v100_v1.1.png", width = 16, height = 10, dpi=400)
