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
#n = 10000 For later
n = 100
p = 0.2 #connectivity of the random graph

#Number of iterations
iters = 10

#Creating baseline network
# ndata = data.frame(ID = 1:n, 
#                    infection = 0,
#                    vaccinated = 0,
#                    state = NA)
# ndata$round = 0
# ndata$new_e = 0 #new infection at the round: this sudo applies only at round=0
# ndata$new_i = 0 #new infection at the round: this sudo applies only at round=0
# ndata$new_r = 0 #new recovery at the round: this sudo applies only at round=0
# ndata$state_end = ifelse(ndata$state==1,ndata$round+e_period,999) #this sudo applies only at round=0, 999 is default
# ndata$contacts = NA
# head(ndata)

# xdata = Matrix(data = 0, nrow = n, ncol = n, sparse = TRUE)
# for(i in 1:2) {
#     ndata[ndata[,names(ndata)[i+3]] == i, "ID"]
#     xdata_temp = Matrix(data = 0, nrow = n, ncol = n, sparse = T)
#     ID_list = ndata[ndata[,names(ndata)[i+2]] == i, "ID"]
#     xdata_temp[ID_list, ID_list] = xdata_temp[ID_list, ID_list] + 1
#     xdata = xdata + xdata_temp
# }
# 
# head(xdata)  
# diag(xdata) = 0
# xdata_sum = Reduce("+", xdata)

## Do it this way?
for(h in 1:iters){
  set.seed(h)
  g0 = sample_gnp(n, p, directed = FALSE, loops = FALSE)
  xdata = as_adjacency_matrix(g0)
  diag(xdata) = 0

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

#Calculate network metrics
  degree = rowSums(xdata, na.rm = T, dims = 1)
  mean = mean(degree)
  sd = sd(degree)
  networkN = mean + ((sd^2)/mean)
  Reff = 0.041*3*networkN

#Create result table
  result = data.frame(h=1:iters,
                    mean = NA,
                    sd = NA,
                    networkN = NA,
                    Reff = NA)

  result[result$h == h,] = c(h, mean, sd, networkN, Reff) 
  
  print(paste0("h= ", h, " is done at ", Sys.time()))
}

################################################################################
#Part B: Simulation of infection dynamics

#Clean environment
# rm(list = ls())

#Set infection parameters
infection_rate = 0.05
beta = 0.04 #transmissibility is 0.04 per day per individual (from Nishi et. al, PNAS 2020)
initial_infections = 5
e_period = 3
i_period = 3
r_period = 300
period = 30
# infection_rate = 0.005

#Simulate infection using SEIR model
ndata1 = ndata
xdata0 = xdata

#Setting initial infection states

ndata1$state = sample(c(rep(1, initial_infections), rep(0, n - initial_infections)), size = n)

ndata1$state_end = ifelse(ndata1$state==1, ndata1$round + e_period, 999)

ndata1$infection = sample(c(rep(1, initial_infections), rep(0, n - initial_infections)), size = n)

head(ndata1)

new_exp = initial_infections/n #new infections vector (first component)
cis = initial_infections/n #cumulative incidences vector (first component)
prevs = 0 #prevalence (proportion of the latent/infectious people)

c(new_exp, cis, prevs)

#interpret diagnosed and symptomatic as the same - these are the people who can be vaccinated 

for (h in 1:10){
  for (m in 1:period) {
    # NOT USING FAMILY TIES FOR NOW  
    # xdata0_family = xdata0
    # xdata0_family[as.numeric(xdata0_family)<1] = 0 #Make all the non-family ties disappear 
    # ndata1$family_contacts = NA  #not used
  
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
      new_exp = c(new_exp,sum(ndata1$new_e==1)/n)
      cis = c(cis,sum(ndata1$state %in% c(1,3,4,5,6))/n)
      prevs = c(prevs,sum(ndata1$state %in% c(4,5,6))/n)
      
      #Step D. Erasing "new" states because they are no longer new after all the above actions
      ndata1$new_e = 0
      ndata1$new_i = 0
      ndata1$new_r = 0
      
      # print(m)
    }
  }

    result_vector = c(h,new_exp,prevs,cis) #add h for random seed
    names(result_vector) = c("h",
                             paste0("n_",0:period),#n_x = x-th round new exposure (infection)
                             paste0("p_",0:period),#p_x = x-th round prevalence
                             paste0("c_",0:period))#c_x = x-th round cumulative incidence
    result_vector
    
    result = as.data.frame(rbind(result, result_vector))
}      

result

ndata1
