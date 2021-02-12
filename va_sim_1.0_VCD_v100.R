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
iters = 100 #n = 1000 for later

#Create network result vector
net_result = data.frame(h=1:iters,
                    mean = NA,
                    sd = NA,
                    networkN = NA,
                    Reff = NA)

setwd("C:/Users/George/Documents/vaccine_alloc/Log")

## Create network
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
people_n = 1000
infection_rate = 0.045 #transmissibility is 0.045 per day per individual (from Nishi et. al, PNAS 2020)
seed_number = 1 # 1 new infection per round
e_period = 3
i_period = 3
r_period = 300
period = 300
available_vaccines = 100 #500 for 10000 people
iters = 100

sim_result_total = NULL
sim_result_sym = NULL
sim_result_total_pv = NULL 
sim_result_sym_pv = NULL


for (h in 1:iters) {
  # h = 1
  set.seed(h)
  
  ndata1 = read.csv(paste0("ndata_R0_h",h,"_0203.csv"))
  xdata0 = readMM(paste0("xdata_h",h,"_0203.mtx"))
  
  #4.2. Setting in ndata1
  ndata1$state  = sample(c(rep(1, seed_number), #state2=#infectious, #state1=E (#seed=1),#state3=R
                           rep(0, people_n-seed_number)), #state0=#susceptible
                         size = people_n)
  ndata1$state_end = ifelse(ndata1$state==1,ndata1$round+e_period,999)
  
  #4.3. Simulations
    new_exp = seed_number/people_n #new infections vector (first component)
    cis = seed_number/people_n #cumulative incidences vector (first component)
    prevs = 0 #prevalence (proportion of the latent/infectious people)
    
    ndata1$inf_length = rgeom(people_n,prob=(1/i_period)) + 1 #mean=3 (Infectious period is determined by geometric distribution (for R, requires + 1))
    #NOTE: If infected, they stay in I(4) or I(5)+I(6) for "inf_length" days (after 3 days of E). 
    ndata1$symptomatic = sample(c(0,1),people_n,replace=T,prob=c(0.45,0.55))
    #NOTE: If infected, each individual has 45% probability for I(4) and 55% probability for I(5)+I(6).
    ndata1$when_symptomatic = rbinom(people_n,size=ndata1$inf_length,prob=0.5) #mean=1.5 (aka duration of presymptomatic period)
    #NOTE: If selected for symptomatic, when does I(6) start? 
    
    n_sym = unlist((ndata1 %>% count(symptomatic == 1))[2, 2]) #number of symptomatics 
    ci_sym = seed_number/n_sym #CI among symptomatics
    prev_sym = 0 #Prevalence among symptomatics
    
    halt_round = NULL
    
    for (m in 1:period){ 
      if(tail(ci_sym, n = 1) <0.1){
            
        ndata1$round = m
        
        #Step A. Implement first all the automatic changes when the day/round changes
        #NOTE: The taf for the state chage is "state_end", which turns to be the same as the present "round"
        
        #A0 - Add cases per round from outside the network
        if(m %in% 1:42*7){
          selected_ID = sample(ndata1[ndata1$state == 0, "ID"],size = seed_number, replace=F) #Only susceptibles can be new infections
          ndata1[ndata1$ID %in% selected_ID,"state"] = 1 
          ndata1[ndata1$ID %in% selected_ID,"new_e"] = 1 
          ndata1[ndata1$ID %in% selected_ID,"state_end"] = m + e_period  
        }
        
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
        ndata1$contacts = as.matrix(xdata0 %*% ifelse(ndata1$state %in% c(4,5),1,0)) #contact with I (state== 4 or 5)
        #NOTE: people with S(0) gets infected from individuals who are not self-isolated = I(4) or I(5)
        #NOTE: sum them all (with weights) are the number of contacts in the present round
        ndata1$new_e = as.numeric(lapply(1-(1-infection_rate)^(ndata1$contacts), sample1)) #possible infections 
        ndata1$new_e = ifelse(ndata1$state==0,ndata1$new_e,0) #only state=0 can get a new infection (latent period)
        ndata1[ndata1$state==0 & ndata1$new_e==1,"state"] = 1 
        ndata1[ndata1$state==1 & ndata1$new_e==1,"state_end"] = ndata1[ndata1$state==1 & ndata1$new_e==1,"round"] + e_period
        
        #Step C1. Record the results for everyone
        new_exp = c(new_exp,sum(ndata1$new_e==1)/people_n)
        cis = c(cis,sum(ndata1$state %in% c(1,3,4,5,6))/people_n)
        prevs = c(prevs,sum(ndata1$state %in% c(4,5,6))/people_n)
        
        #Step C2. Record the results among symptomatic
        new_exp_sym = new_exp
        ci_sym = c(ci_sym, sum((ndata1 %>% filter(symptomatic == 1))$state %in% c(1, 3, 4, 5, 6))/n_sym)
        prev_sym = c(prev_sym, sum((ndata1 %>% filter(symptomatic == 1))$state %in% 4:6)/n_sym)
        
        #Step D. Erasing "new" states because they are no longer new after all the above actions
        ndata1$new_e = 0
        ndata1$new_i = 0
        ndata1$new_r = 0
      } else {
      print(paste0("Sim ",h," stopped when CI among symptomatic individuals reached 10% at round ",m))
      halt_round = m
      break
      }  
    }
    
    #Once CI among symptomatics is 10%, initiate vaccination strategy
    remaining_period = period - halt_round
    
    #Strategy 1 - Vaccinate Close Contacts of Previously Diagnosed/Symptomatic Individuals
    contact_matrix = as.matrix(xdata0 %*% ifelse(ndata1$symptomatic == 1,1,0))
    contacted_ids = which(contact_matrix > 0)

    sum(ndata1$state %in% c(1,3,4,5,6))
    
    if(length(contacted_ids) <= available_vaccines){
        ndata1[ndata1$ID %in% contacted_ids, "vaccinated"] = 1
        } else if (length(contacted_ids) > available_vaccines) {
          to_vacc = sample(contacted_ids, size = available_vaccines, replace = FALSE)
          ndata1[ndata1$ID %in% to_vacc, "vaccinated"] = 1}
    
    #Continue simulation as before until full period is reached
    
    new_exp_pv = tail(new_exp, 1) #new infections vector (first component)
    cis_pv = tail(cis, 1) #cumulative incidences vector (first component)
    prevs_pv = tail(prevs, 1) #prevalence (proportion of the latent/infectious people)
    
    # n_sym_pv = n_sym #number of symptomatic 
    ci_sym_pv = tail(ci_sym, 1) #CI among symptomatic
    prev_sym_pv = tail(prev_sym, 1) #Prevalence among symptomatic
    
    
    
    for (m in (halt_round+1):period){
        # m = halt_round + 1
        ndata1$round = m
        
        #Step A. Implement first all the automatic changes when the day/round changes
        #NOTE: The taf for the state change is "state_end", which turns to be the same as the present "round"

        #A0 - Add cases per round from outside the network
        if(m %in% 1:42*7){
          selected_ID = sample(ndata1[ndata1$state == 0, "ID"],size = seed_number, replace=F) #Only susceptibles can be new infections
          ndata1[ndata1$ID %in% selected_ID,"state"] = 1 
          ndata1[ndata1$ID %in% selected_ID,"new_e"] = 1 
          ndata1[ndata1$ID %in% selected_ID,"state_end"] = m + e_period  
        }
        sum(ndata1$state %in% c(1,3,4,5,6))
        
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
        ndata1$contacts = as.matrix(xdata0 %*% ifelse(ndata1$state %in% c(4,5),1,0)) #contact with I (state== 4 or 5)
        #NOTE: people with S(0) gets infected from individuals who are not self-isolated = I(4) or I(5)
        #NOTE: sum them all (with weights) are the number of contacts in the present round
        ndata1$new_e = as.numeric(lapply(1-(1-infection_rate)^(ndata1$contacts), sample1)) #possible infections 
        ndata1$new_e = ifelse(ndata1$state==0 & ndata1$vaccinated == 0,ndata1$new_e,0) #only state=0 can get a new infection (latent period)
        ndata1[ndata1$state==0 & ndata1$new_e==1,"state"] = 1 
        ndata1[ndata1$state==1 & ndata1$new_e==1,"state_end"] = ndata1[ndata1$state==1 & ndata1$new_e==1,"round"] + e_period
        sum(ndata1$state %in% c(1,3,4,5,6))
        
        #Step C1. Record the post-vaccination results for everyone
        new_exp_pv = c(new_exp_pv,sum(ndata1$new_e==1)/people_n)
        cis_pv = c(cis_pv,sum(ndata1$state %in% c(1,3,4,5,6))/people_n)
        prevs_pv = c(prevs_pv,sum(ndata1$state %in% c(4,5,6))/people_n)
        
        #Step C2. Record the results among symptomatic
        new_exp_sym = new_exp_pv
        ci_sym_pv = c(ci_sym_pv, sum((ndata1 %>% filter(symptomatic == 1))$state %in% c(1, 3, 4, 5, 6))/n_sym)
        prev_sym_pv = c(prev_sym_pv, sum((ndata1 %>% filter(symptomatic == 1))$state %in% 4:6)/n_sym)
        
        cis_pv
        
        #Step D. Erasing "new" states because they are no longer new after all the above actions
        ndata1$new_e = 0
        ndata1$new_i = 0
        ndata1$new_r = 0
        
      } 
    

        result_vector_total = c(h,new_exp,prevs,cis)

        result_vector_sym = c(h, new_exp, prev_sym, ci_sym)
        
        result_vector_total_pv = c(h, new_exp_pv, prevs_pv, cis_pv)
        
        result_vector_sym_pv = c(h, new_exp_pv, prev_sym_pv, ci_sym_pv)
        
        names(result_vector_total) = c("h",
                                     paste0("n_",0:(length(new_exp)-1)),#n_x = x-th round new exposure (infection)
                                     paste0("p_",0:(length(prevs)-1)),#p_x = x-th round prevalence
                                     paste0("c_",0:(length(cis)-1)))#c_x = x-th round cumulative incidence,
        

        names(result_vector_sym) = c("h", 
                                     paste0("n_",0:(length(new_exp)-1)), 
                                     paste0("psym_",0:(length(prev_sym)-1)), 
                                     paste0("csym_",0:(length(ci_sym)-1)))
        
        names(result_vector_total_pv) = c("h",
                                       paste0("n_",(halt_round):period),#n_x = x-th round new exposure (infection)
                                       paste0("p_",(halt_round):period),#p_x = x-th round prevalence
                                       paste0("c_",(halt_round):period))#c_x = x-th round cumulative incidence,
                                       
        names(result_vector_sym_pv) = c("h", 
                                        paste0("n_",(halt_round):period),#n_x = x-th round new exposure (infection)
                                        paste0("p_sym",(halt_round):period),#p_x = x-th round prevalence
                                        paste0("c_sym",(halt_round):period))#c_x = x-th round cumulative incidence,
                                                                                 
        sim_result_total = bind_rows(sim_result_total, result_vector_total)
        sim_result_sym = bind_rows(sim_result_sym, result_vector_sym)
        sim_result_total_pv = bind_rows(sim_result_total_pv, result_vector_total_pv)
        sim_result_sym_pv = bind_rows(sim_result_sym_pv, result_vector_sym_pv)
        
        print(paste0("Sim ",h, " is done at ", Sys.time()))
}
#Save results

setwd("C:/Users/George/Documents/vaccine_alloc/Results")

write_csv(sim_result_total_pv, "sim_result_VCD_v100_0209.csv")
write_csv(sim_result_sym_pv, "sim_result_sym_only_VCD_v100_0209.csv")
