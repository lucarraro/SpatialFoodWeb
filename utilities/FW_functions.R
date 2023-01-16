# Calculate unweighted TL (prey-averaged, see Williams & Martinez 2004) and Omnivory index & coherence index based on this TL.
Unw.TL.Omn.q = function(DietMatrix){
  if(nrow(DietMatrix) != ncol(DietMatrix)) stop("nrow != ncol, is this a unipartite web?")
  Spe.Richness. = nrow(DietMatrix)
  basal. = which(apply(DietMatrix, 2, sum) == 0)
  
  if(length(basal.) != 0 ){
    TL1 = -t(DietMatrix)/apply(DietMatrix, 2, sum)
    TL1[basal.,] = 0
    diag(TL1) = 1
    TL2 = c(rep(1, Spe.Richness.))
    
    if(det(TL1) < 1e-14) {
      TL = NA; Omn = NA; q = NA # Avoid error from computational singular matrix.
    } else {
      TL = solve(TL1, TL2)
      
      Omn = rep(NA, Spe.Richness.)
      for(s in 1:Spe.Richness.){
        if(s %in% basal.) next
        if(apply(DietMatrix, 2, sum)[s] == 1){ Omn[s] = 0
        } else {Omn[s] = sd(TL[which(DietMatrix[, s]==1)])}
      } # end of s for-loop.
      
      Distance1 = matrix(TL, Spe.Richness., Spe.Richness., byrow = T)
      Distance2 = t(Distance1)
      Distance = DietMatrix * (Distance2 - Distance1)
      q = sd(Distance[which(Distance != 0)])
    } # end of else.  
  } else {TL = NA; Omn = NA; q = NA } # end of basal. if statement.
  
  return(list(TL = TL, Omn = Omn, q = q))
} # end of the function.

# A function repeats the calculation of modularity for given times, then derive the max value to avoid local maximum. (igraph needed)
Max.modu = function(DietMatrix, times = 10){
  all.modu = c()
  for (t in 1:times) {
    all.modu = append(all.modu, modularity(multilevel.community(graph.adjacency(DietMatrix, mode = "undirected"))))
  }
  return(max(all.modu)) 
} # end of the function.

