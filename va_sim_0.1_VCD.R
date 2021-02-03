###############################################################################
# Part A: Network Construction 

#Clean environment
rm(list = ls())

#Load required packages
library(Matrix)
library(igraph)
library(tidyverse)

#Create sample convenience function
sample1 = function(x){
  sample(c(0, 1), size = 1, replace = FALSE, prob = c(1-x, x))
}

#Set network parameters
n = 1000 #n = 10000 For later
p = 0.2 #rewiring rate for Watts-Strogatz small world model

#Number of iterations
iters = 1 #n = 1000 for later

#Create network result vector
net_result = data.frame(h=1:iters,
                    mean = NA,
                    sd = NA,
                    networkN = NA,
                    Reff = NA)

setwd("C:/Users/George/Documents/vaccine_alloc/Log")

## Do it this way?
for(h in 1:iters){
  set.seed(h)
  g0 = sample_smallworld(dim = 1, size = n, nei = 7, p = p) #How to set transitivity = 0.18?
  xdata = as_adjacency_matrix(g0)
  diag(xdata) = 0
  writeMM(obj = xdata, file = paste0("xdata_h", h, "_0203.mtx"))
  
#Checking network appearance
# plot(g0, layout = layout.kamada.kawai, vertex.size = 5, vertex.label = NA)

  ndata = data.frame(ID = 1:n, 
                   infection = 0,
                   vaccinated = 0,
                   state = NA)
  ndata$round = 0
  ndata$new_e = 0 #new infection at the round: this sudo applies only at round=0
  ndata$new_i = 0 #new infection at the round: this sudo applies only at round=0
  ndata$new_r = 0 #new recovery at the round: this sudo applies only at round=0
  ndata$state_end = ifelse(ndata$state==1,ndata$round+e_period,999) #this sudo applies only at round=0, 999 is default
  ndata$contacts = NA
  head(ndata)
  write_csv(ndata, file = paste0("ndata_R0_h", h, "_0203.csv"))

#Calculate network metrics
  degree = rowSums(xdata, na.rm = T, dims = 1)
  mean = mean(degree)
  sd = sd(degree)
  networkN = mean + ((sd^2)/mean)
  Reff = 0.041*3*networkN

#Create result table

  net_result[net_result$h == h,] = c(h, mean, sd, networkN, Reff) 
  
  print(paste0("h = ", h, " is done at ", Sys.time()))
}

write_csv(net_result, file = paste0("net_results_0203.csv"))
################################################################################
#Part B: Simulation of infection dynamics

#Clean environment
rm(list = ls())

#Loading packages
library(tidyverse)

#Create sample convenience function
sample1 = function(x){
  sample(c(0, 1), size = 1, replace = FALSE, prob = c(1-x, x))
}

#Set infection parameters
n = 1000
infection_rate = 0.045 #transmissibility is 0.045 per day per individual (from Nishi et. al, PNAS 2020)
new_infections = 1 # 1 new infection per round
e_period = 3
i_period = 3
r_period = 300
period = 300
iters = 1
# infection_rate = 0.005

sim_result = c(0, rep(0, period+1), rep(0, period+1), rep(0, period+1))
names(sim_result) = c("h",
                      paste0("n_",0:period),#n_x = x-th round new exposure (infection)
                      paste0("p_",0:period),#p_x = x-th round prevalence
                      paste0("c_",0:period))#c_x = x-th round cumulative incidence,

#interpret diagnosed and symptomatic as the same - these are the people who can be vaccinated 

for (h in 1:iters){
    set.seed(h)
  
    #Simulate infection using SEIR model
    ndata1 = read.csv(paste0("ndata_R0_h", h, "_0203.csv"))
    xdata0 = readMM(paste0("xdata_h", h, "_0203.mtx"))

    #Setting initial infection states
    ndata1$state = sample(c(rep(1, new_infections), rep(0, n - new_infections)), size = n)
    ndata1$state_end = ifelse(ndata1$state==1, ndata1$round + e_period, 999)

    #Set initial parameters
    new_exp = new_infections/n #new infections vector (first component)
    cis = new_infections/n #cumulative incidences vector (first component)
    prevs = 0 #prevalence (proportion of the latent/infectious people)

    ndata1$inf_length = rgeom(n,prob=(1/i_period)) + 1 #mean=3 (Infectious period is determined by geometric distribution (for R, requires + 1))
    #NOTE: If infected, they stay in I(4) or I(5)+I(6) for "inf_length" days (after 3 days of E). 
    ndata1$symptomatic = sample(c(0,1),n,replace=T,prob=c(0.45,0.55))
    #NOTE: If infected, each individual has 45% probability for I(4) and 55% probability for I(5)+I(6).
    ndata1$when_symptomatic = rbinom(n,size=ndata1$inf_length,prob=0.5) #mean=1.5 (aka duration of presymptomatic period)
      #NOTE: If selected for symptomatic, when does I(6) start? 
      
      #NOTE: I(2) will go, and I-asymptomatic(4), I-presymptomatic(5), and I-postsymptomatic(6) will be used.
      #NOTE: I-asymptomatic(4) lasts for all the I periods once determined.
      #NOTE: I-presymptomatic(5) lasts for "when_symptomatic" days, and the state moves to I(6) 
      #NOTE: This is fine because individuals would have at most 1 infection in lifetime here in the simulations.
      
      for (m in 1:period) {

          ndata1$round = m
          
          #Step A. Implement first all the automatic changes when the day/round changes
          #NOTE: The tag for the state change is "state_end", which turns to be the same as the present "round"
          
          #A0 - Add cases per round from outside the network
          selected_ID = sample(ndata1[ndata1$state == 0, "ID"],size = new_infections,replace=F)
          ndata1[ndata1$ID %in% selected_ID,"state"] = 1 
          ndata1[ndata1$ID %in% selected_ID,"new_e"] = 1 
          ndata1[ndata1$ID %in% selected_ID,"state_end"] = m + e_period  
          
          #A1. E(1) to I(4 or 5)
          ndata1[ndata1$state == 1 & (ndata1$round == ndata1$state_end),"new_i"] = 1
          #I(4)
          ndata1[ndata1$state==1 & ndata1$new_i==1 & ndata1$symptomatic==0,"state"] = 4     
          ndata1[ndata1$state==4 & ndata1$new_i==1 & ndata1$symptomatic==0,"state_end"] = ndata1[ndata1$state==4 & ndata1$new_i==1 & ndata1$symptomatic==0,"round"] + ndata1[ndata1$state==4 & ndata1$new_i==1 & ndata1$symptomatic==0,"inf_length"]
          #NOTE: When state changes, state_end date/round needs to be updated.
          #I(5)
          ndata1[ndata1$state==1 & ndata1$new_i==1 & ndata1$symptomatic==1,"state"] = 5     
          ndata1[ndata1$state==5 & ndata1$new_i==1 & ndata1$symptomatic==1,"state_end"] = ndata1[ndata1$state==5 & ndata1$new_i==1 & ndata1$symptomatic==1,"round"] + ndata1[ndata1$state==5 & ndata1$new_i==1 & ndata1$symptomatic==1,"when_symptomatic"]
          
          #A1x. I(5) to I(6)
          ndata1[ndata1$state == 5 & (ndata1$round == ndata1$state_end),"new_i"] = 2 #2 is new flag for symptom acquired
          ndata1[ndata1$state == 5 & (ndata1$round == ndata1$state_end),"state"] = 6    
          ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"state_end"] = ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"round"] + ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"inf_length"] - ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"when_symptomatic"]
          
          #A2. I(4,6) to R(3)
          ndata1[ndata1$state %in% c(4,6) & (ndata1$round == ndata1$state_end),"new_r"] = 1
          ndata1[ndata1$state %in% c(4,6) & ndata1$new_r==1,"state"] = 3     
          ndata1[ndata1$state==3 & ndata1$new_r==1,"state_end"] = ndata1[ndata1$state==3 & ndata1$new_r==1,"state_end"] + r_period
          
          #Step B. Implement second all the infection events after all the statuses are updated
          #NOTE: Here, making secondary infections from primary cases
          #B1. For possible state change from S(0) to E(1)
          ndata1$contacts = as.matrix(xdata0 %*% ifelse(ndata1$state %in% c(4,5),1,0)) #contact with infected (state== 4 or 5)
          # ndata1$family_contacts = as.matrix(xdata0_family %*% ifelse(ndata1$state==6,1,0)) #family member only
          #NOTE: people with S(0) gets infected from individuals who are not self-isolated = I(4) or I(5)
          #NOTE: or from individual who are self-isolated but are their family members
          #NOTE: sum them all (with weights) are the number of contacts in the present round
          ndata1$new_e = as.numeric(lapply(1-(1-infection_rate)^(ndata1$contacts), sample1)) #possible infections 
          ndata1$new_e = ifelse(ndata1$state==0,ndata1$new_e,0) #only state=0 can get a new infection (latent period)
          ndata1[ndata1$state==0 & ndata1$new_e==1,"state"] = 1 
          ndata1[ndata1$state==1 & ndata1$new_e==1,"state_end"] = ndata1[ndata1$state==1 & ndata1$new_e==1,"round"] + e_period
          
          #If people are vaccinated (regardless of what state they were in), they go directly to R(3)
          # ndata1[ndata1$state %in% 1:6 & (ndata1$vaccinated == 1), "state"] = 3
          
          #Step C. Record the results
          if(tail(cis, n = 1) <= 0.1){
            new_exp = c(new_exp,sum(ndata1$new_e==1)/n)
            cis = c(cis,sum(ndata1$state %in% c(1,3,4,5,6))/n)
            prevs = c(prevs,sum(ndata1$state %in% c(4,5,6))/n)
          }
          else {
            new_exp = c(new_exp, 999)
            cis = c(cis, 999)
            prevs = c(prevs, 999)
            print(paste0("Simulation stopped when reaching CI = 0.1 at round ", m))
            
            #Save ndata at this point to implement vaccination
            if (h < 10){
              write_csv(ndata1, file = paste0("ndata1_h_0",h,"_ci_limit_reached_0201_rd_",m,".csv"))
            }
            else {
              write_csv(ndata1, file = paste0("ndata1_h_",h,"_ci_limit_reached_0201_rd_",m,".csv"))
            }
            break
          }
          

          #Step D. Erasing "new" states because they are no longer new after all the above actions
          ndata1$new_e = 0
          ndata1$new_i = 0
          ndata1$new_r = 0
          
        }
          
        result_vector = c(h,new_exp,prevs,cis) 
        names(result_vector) = c("h",
                                 paste0("n_",1:length(new_exp)),#n_x = x-th round new exposure (infection)
                                 paste0("p_",1:length(prevs)),#p_x = x-th round prevalence
                                 paste0("c_",1:length(cis)))#c_x = x-th round cumulative incidence,
        sim_result = bind_rows(sim_result, result_vector)
    }

#Save results
sim_result
to_save = sim_result %>% filter(h != 0)
write_csv(to_save, file = paste0("sim_result_pre_vac_0203.csv"))

###############################################################################
# Part C: Vaccination and continue simulation

# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse)

#Create sample convenience function
sample1 = function(x){
  sample(c(0, 1), size = 1, replace = FALSE, prob = c(1-x, x))
}

#Set infection parameters
n = 1000
infection_rate = 0.045 #transmissibility is 0.04 per day per individual (from Nishi et. al, PNAS 2020)
new_infections = 1
e_period = 3
i_period = 3
r_period = 300
period = 300
available_vaccines = 500

sim2_result = NULL

#interpret diagnosed and symptomatic as the same - these are the people who can be vaccinated 
# 
# for (h in 1:1){
h = 1
  # Load data from when simulation stopped
  # and continue simulation
  if(h < 10){
    file_to_use = list.files(path = getwd(), pattern = paste0("ndata1_h_0",h))[1]
    } else {
    file_to_use = list.files(path = getwd(), pattern = paste0("ndata1_h_",h))[1]
    }
  ndata1 = read_csv(file_to_use) 
  xdata0 = readMM(paste0("xdata_h", h, "_0203.mtx"))
  pre_vac_results = read_csv("sim_result_pre_vac_0203.csv")
  last_day = as.numeric(str_extract(file_to_use, pattern = regex("\\d{1,2}\\."))) 
  new_start_day = last_day + 1
  
  #Set first value of result vectors as value from last round of pre-vacc simulation
  new_exp = pre_vac_results %>% filter(h == h) %>% 
    select(one_of(paste0("n_",as.character(last_day)))) %>% as.numeric()
  cis = pre_vac_results %>% filter(h == h) %>% 
    select(one_of(paste0("c_",as.character(last_day)))) %>% as.numeric()
  prevs = pre_vac_results %>% filter(h == h) %>% 
    select(one_of(paste0("p_",as.character(last_day)))) %>% as.numeric()

  #Strategy 1 - Vaccinate Close Contacts of Previously Diagnosed/Symptomatic Individuals
  contact_matrix = as.matrix(xdata0 %*% ifelse(ndata1$symptomatic == 1,1,0))
  contacted_ids = which(contact_matrix > 0)
  
  if(length(contacted_ids) <= available_vaccines){
      ndata1[ndata1$ID %in% contacted_ids, "vaccinated"] = 1
    } else if (length(contacted_ids) > available_vaccines) {
      to_vacc = sample(contacted_ids, size = available_vaccines, replace = FALSE)
      ndata1[ndata1$ID %in% to_vacc, "vaccinated"] = 1
    }

  
  
#Continue the simulation on the day after we stopped before. Simulation behavior is the same
#except for people who were vaccinated (immediately -> R(3))
  
  for (m in new_start_day:period) {
    
    ndata1$round = m

    #Step A. Implement first all the automatic changes when the day/round changes
    #NOTE: The tag for the state change is "state_end", which turns to be the same as the present "round"
      
      #If people are vaccinated (regardless of what state they were in), they go directly to R(3)
      ndata1[ndata1$vaccinated == 1, "state"] = 7
      
      #A1. E(1) to I(4 or 5)
      ndata1[ndata1$state == 1 & (ndata1$round == ndata1$state_end),"new_i"] = 1
      #I(4)
      ndata1[ndata1$state==1 & ndata1$new_i==1 & ndata1$symptomatic==0,"state"] = 4     
      ndata1[ndata1$state==4 & ndata1$new_i==1 & ndata1$symptomatic==0,"state_end"] = ndata1[ndata1$state==4 & ndata1$new_i==1 & ndata1$symptomatic==0,"round"] + ndata1[ndata1$state==4 & ndata1$new_i==1 & ndata1$symptomatic==0,"inf_length"]
      #NOTE: When state changes, state_end date/round needs to be updated.
      #I(5)
      ndata1[ndata1$state==1 & ndata1$new_i==1 & ndata1$symptomatic==1,"state"] = 5     
      ndata1[ndata1$state==5 & ndata1$new_i==1 & ndata1$symptomatic==1,"state_end"] = ndata1[ndata1$state==5 & ndata1$new_i==1 & ndata1$symptomatic==1,"round"] + ndata1[ndata1$state==5 & ndata1$new_i==1 & ndata1$symptomatic==1,"when_symptomatic"]
      
      #A1x. I(5) to I(6)
      ndata1[ndata1$state == 5 & (ndata1$round == ndata1$state_end),"new_i"] = 2 #2 is new flag for symptom acquired
      ndata1[ndata1$state == 5 & (ndata1$round == ndata1$state_end),"state"] = 6    
      ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"state_end"] = ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"round"] + ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"inf_length"] - ndata1[ndata1$state==6 & ndata1$new_i==2 & ndata1$symptomatic==1,"when_symptomatic"]
      
      #A2. I(4,6) to R(3)
      ndata1[ndata1$state %in% c(4,6) & (ndata1$round == ndata1$state_end),"new_r"] = 1
      ndata1[ndata1$state %in% c(4,6) & ndata1$new_r==1,"state"] = 3     
      ndata1[ndata1$state==3 & ndata1$new_r==1,"state_end"] = ndata1[ndata1$state==3 & ndata1$new_r==1,"state_end"] + r_period
           
      #Step B. Implement second all the infection events after all the statuses are updated
      #NOTE: Here, making secondary infections from primary cases
      
      #B1. For possible state change from S(0) to E(1)
     
      ndata1$contacts = as.matrix(xdata0 %*% ifelse(ndata1$state %in% c(4,5),1,0)) #contact with infected (state== 4 or 5)
      
      #NOTE: sum them all (with weights) are the number of contacts in the present round
      ndata1$new_e = as.numeric(lapply(1-(1-infection_rate)^(ndata1$contacts), sample1)) #possible infections 
      ndata1$new_e = ifelse(ndata1$state==0,ndata1$new_e,0) #only state=0 can get a new infection (latent period)
      ndata1[ndata1$state==0 & ndata1$new_e==1,"state"] = 1 
      ndata1[ndata1$state==1 & ndata1$new_e==1,"state_end"] = ndata1[ndata1$state==1 & ndata1$new_e==1,"round"] + e_period
      
      #Step C. Record the results
      new_exp = c(new_exp,sum(ndata1$new_e==1)/n)
      cis = c(cis,sum(ndata1$state %in% c(1,3,4,5,6))/n)
      prevs = c(prevs,sum(ndata1$state %in% c(4,5,6))/n)
      
      #Step D. Erasing "new" states because they are no longer new after all the above actions
      ndata1$new_e = 0
      ndata1$new_i = 0
      ndata1$new_r = 0
      
    }
  }  

  # Save the results    
  result_vector = c(h,new_exp,prevs,cis)
  result_vector
  names(result_vector) = c("h",
                           paste0("n_",new_start_day:period),#n_x = x-th round new exposure (infection)
                           paste0("p_",new_start_day:period),#p_x = x-th round prevalence
                           paste0("c_",new_start_day:period))#c_x = x-th round cumulative incidence
 
  sim2_result = as.data.frame(rbind(sim2_result, result_vector))
  
  setwd("C:/Users/George/Documents/vaccine_alloc/Log")
  write_csv(sim2_result, file = "sim_result_post_vac_0203.csv")
  
}
sim2_result
