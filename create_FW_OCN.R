rm(list=ls())

library(deSolve)
library(R.matlab)
library(OCNet)
# Niche model
# bodymass: a numeric vector of species body masses
# NF.Prop: a 0-1 value or "free": the proportion of nutrient feeders, which include basal species (producers) and consumers feeding also on nutrient. "free" means no constraint and take whatever the model generated.
# Connectance: a 0-1 value for the desired connectance of the generated web (not counting links to the nutrient block).
# Cannibalism: TRUE/FALSE, whether cannibalistic links are allowed in the generated web.
# Mutual.Feeding: TRUE/FALSE, whether bi-directional consumption is allowed in the generated web.
# Precise.Control: TRUE/FALSE, whether to force the generated connectance to be closer to the desired value by removing one layer of stochasticity.
# Beta.Scale: a scaling value that influence the beta distribution from with the foraging centre and radius are drawn. In general, a large value avoids super small radius.
# Checkpoint: TRUE/FALSE, whether to active a check point to disallow isolated nodes.

Nic = function(bodymass,
               NF.Prop = "free",  # when "free" usually it's the smallest one(s) being the NF.
               Connectance = 0.25,
               Cannibalism = TRUE,
               Mutual.Feeding = TRUE,
               Precise.Control = FALSE,
               Beta.Scale = 100,
               Checkpoint = FALSE){
  
  richness = length(bodymass)
  if(NF.Prop != "free") no.NF = round(richness * NF.Prop)
  nutri.feeder = rep("N", richness) # a vector to later store identities of nutrient feeders.
  
  iter = 0
  repeat{ # The checkpoint repeat loop.
    iter = iter + 1
    # Assign niche values
    if(Precise.Control){
      n = rank(bodymass)/richness
    } else {
      n = runif(richness)
      n = n[order(n)][rank(bodymass)] # regenerated, uniform distributed niche values according to bodymass ranking
    }
    
    # An empty diet matrix to be filled
    D = matrix(0, richness, richness)
    
    # Corrected connectance depending on the selected criteria
    if(Cannibalism){
      corrected.connectance = Connectance
    } else {
      corrected.connectance = (Connectance * richness^2)/(richness * (richness - 1))
    }
    
    # Assign diets
    alpha = (2 * corrected.connectance) * Beta.Scale
    beta = (1 - 2 * corrected.connectance) * Beta.Scale
    radius = n * rbeta(richness, alpha, beta)
    centre = unlist(lapply(c(1:richness), function(x){runif(1, radius[x]/2, n[x])})) # due to the setting here, none of the diets will cover n=0 or lower.
    for(i in 1:richness){
      diet = which(n <= centre[i] + radius[i]/2 & n >= centre[i] - radius[i]/2)
      D[diet, i] = 1
      if(length(setdiff(diet, i)) == 0) nutri.feeder[i] = "Y"  # Those eat nothing (other than itself) are assigned to feed on nutrient.
      if(Cannibalism == F) diag(D) = 0 # just wire-up as no constraint then remove cannibalistic links, as in Cannibalism == F case connectances is already corrected.
    } # End of i for-loop.
    no.diet = which(nutri.feeder == "Y") # those without feeding on any others, and already be assigned as NF.
    
    if(iter > 100) { message("The given community and niche model settings are unlikely to generate a valid food web."); D = NA; break}
    if(NF.Prop != "free"){
      
      # IF TOO MANY, OK; IF TOO FEW, ADD THE LIGHTEST
      if(length(which(nutri.feeder == "Y")) < no.NF){
        no.new.NF = no.NF - length(which(nutri.feeder == "Y"))
        nutri.feeder[which(nutri.feeder == "N")[(length(which(nutri.feeder == "N"))-no.new.NF+1):length(which(nutri.feeder == "N"))]] <- "Y"                           
      }
    }
    if(Mutual.Feeding == FALSE) {if(any((D + t(D))[upper.tri(D)] > 1)) next} # if mutual feeding i.e., two-node loop detected, regenerate a new web.
    if(Checkpoint == FALSE) {break} # If the checkpoint isn't activated, accept whatever web generated with the first iteration.
    if(Checkpoint == TRUE & prod(apply(D, 2, sum)[-no.diet]) > 0 & prod(apply(D, 1, sum)[no.diet]) > 0) {break} # If no isolated nodes, pass. Otherwise regenerate a new web.
    
  } # End of the checkpoint repeat loop.
  
  return(list(D = D, NF = which(nutri.feeder == "Y")))
}

nSp <- 99 # makes 100 with nutrients
MassMean <- 10^-2
MassSD <-  10
a0 <- 1

nFW <- 100
bodymass_list <- dietMatrix_list <- nutrientFeeders_list <- r_list <- A_list <- vector("list",nFW)
nms <- NULL
for (i in 1:nFW){
  nms <- c(nms, paste0('w',i))
}
names(bodymass_list) <- names(dietMatrix_list) <- names(nutrientFeeders_list) <- nms
names(r_list) <- names(A_list) <- nms

for (ind in 1:nFW){
  cat(sprintf('ind: %d \n',ind))
  set.seed(1000+floor((ind-1)/10) + 1) # first 10 FW share the same bodymass vector
  bodymass <- rlnorm(nSp, log(MassMean), log(MassSD))
  bodymass <- sort(bodymass,decreasing = T) # indexing sorted by mass (1st is largest)
  set.seed(ind)
  Web <-  Nic(bodymass, NF.Prop=0.05, Connectance=0.1,Cannibalism=T)
  D <- Web$D # the generated food web as a diet matrix, cols eat rows.
  NF <- Web$NF # nutrient feeders in the web.
  bodymass_list[[ind]] <- bodymass
  dietMatrix_list[[ind]] <- D
  nutrientFeeders_list[[ind]] <- NF

  r_list[[ind]] <- -4.15e-8*bodymass^-0.25
  D_N <- cbind(rbind(D, numeric(nSp)),numeric(nSp+1)) # Diet matrix with nutrients added
  D_N[nSp+1,NF] <- 1
  bodymass_N <- c(bodymass,2e-5) # nutrient mass is set to 2e-5
  
  A <- matrix(0,nSp+1,nSp+1)
  for (i in 1:(nSp+1)){
    for (j in 1:(nSp+1)){
      A[i,j] <- -2.72*D_N[i,j]*bodymass_N[i]^0.63*bodymass_N[j]^0.42 +
        0.5*2.72*D_N[j,i]*bodymass_N[i]^(-0.37)*bodymass_N[j]^1.42
    }
  }
  diag(A) <- diag(A) + c(-a0*bodymass^0.5, 0)
  A_list[[ind]] <- A
}

writeMat('utilities/100FW_nSp100_c01.mat',bodymass_list=bodymass_list, dietMatrix_list=dietMatrix_list,
         nutrientFeeders_list=nutrientFeeders_list, r_list=r_list, A_list=A_list)

# create OCN

if (!file.exists("utilities/OCN.rda")){
  library(OCNet)
  set.seed(8)
  OCN <- create_OCN(400,400,outletPos=20,typeInitialState = "V",
                    coolingRate = 0.5, initialNoCoolingPhase = 0.1, nIter=50*400*400,cellsize=100,
                    displayUpdates = 2, nUpdates=200)
  
  OCN <- landscape_OCN(OCN, zMin=200, slope0 = 0.005)
  OCN <- aggregate_OCN(OCN, thrA=5e6, maxReachLength = 4000)
  OCN <- paths_OCN(OCN)
  OCN <- rivergeometry_OCN(OCN, widthMax = 24, depthMax = 1, velocityMax = 1.7) #??? Qoutlet=40 m3/s (=0.025*Acatchment [km2])
  # with slope_outlet=0.005, it is Ks = 25 m^1/3 s^-1 (average value for main channels)
  # slope = 0.005 is for mountainous catchments (~2500 m altitude range)
  save(file="OCN.rda",OCN)
} else {load('utilities/OCN.rda')}

distToOutlet <- OCN$AG$downstreamPathLength[,OCN$AG$outlet]
writeMat('utilities/OCN.mat',W=as.matrix(OCN$AG$W), downNode=OCN$AG$downNode, A=OCN$AG$AReach, width=OCN$AG$width, depth=OCN$AG$depth, 
         velocity=OCN$AG$velocity, As=OCN$SC$ALocal, L=OCN$AG$leng, streamOrder=OCN$AG$streamOrder, distToOutlet=distToOutlet)