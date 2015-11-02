#######################
##' This function takes in a WAIFW and returns a WAIFW that corresponds to a specified R_0
##' Waifw is the relative mixing between different age groups
##' R_0 is the desired R_0
##' state is the number of individuals in each age group.
#######################

output.waifw<- function(waifw, R_0, state){ 
    denom <- sum(state)
    next.gen <- as.numeric(state)*(1-exp(-waifw/denom))
    
    #get the first eigen value
    cur.R0 <- Re(eigen(next.gen)$value[1])
    
    #More correct transform
    R.ratio <- R_0/cur.R0 #print(R0); #print(cur.R0); #print(R.ratio)
    waifw <- -log(1-R.ratio*(1-exp(-waifw/denom)))*denom
    return(waifw)
}



#######################
##' This function calculates the force of infection on each age group - method as described in Petra Klepac paper.
##' Input the waifw, the cases at time point in question and the size of the population.
#######################

calc.phi <- function(waifw, x, denom, i.power){
    
    phi <- (waifw %*% t(x[, ] ) ^ i.power) /denom
    
    phi <- 1 - exp(-phi)
    
    phi <- t(phi)
    return(phi)
}

