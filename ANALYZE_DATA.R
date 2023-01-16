rm(list=ls())

library(R.matlab)
library(OCNet)
library(UNODF) # for nested calculation.
library(vegan)
library(igraph) # for modularity calculation.
library(bipartite) # for networklevel().
library(car) # for Anova (type III)

source('utilities/FW_functions.R')

load('utilities/OCN.rda')
distToOutlet <- OCN$AG$downstreamPathLength[,OCN$AG$outlet]

createPlots <- T # want to create the main figures?

# LOAD RESULTS ####
# define 4 regions of the river networks 
low = OCN$AG$A < median(OCN$AG$A) & distToOutlet < median(distToOutlet)
high = OCN$AG$A < median(OCN$AG$A) & distToOutlet > median(distToOutlet)
down = OCN$AG$A > median(OCN$AG$A) & distToOutlet < median(distToOutlet)
mid = OCN$AG$A > median(OCN$AG$A) & distToOutlet > median(distToOutlet)

position = 1*high + 2*mid + 3*low + 4*down
kol <- hcl.colors(4, "Temps") # define colors

foodWebs <- readMat("utilities/100FW_nSp100.mat") # load food webs
mat <- readMat("utilities/results_for_R.mat") # load simulation results

if (!file.exists("utilities/FW_indices.rda")){
  df <- data.frame(alphaDiv = mat$alphaDiv.vec,  nutFeed=mat$nNutFeed.vec,
                   nL = mat$nL.vec, nT = mat$nT.vec, dR = mat$dR.vec, dB = mat$dB.vec, vUp = mat$vUp.vec, 
                   pos = mat$pos.vec, FW = mat$FW.vec)
  
  # evaluate FW metrics
  nNodes <- OCN$AG$nNodes
  nSp <- length(foodWebs$bodymass.list)
  nSim <- length(mat$FW.raw)
  thr <- 5e-7
  Connectance <- LinkDensity <- Nestedness <- numeric(length(df$alphaDiv))
  Modularity <- NicheOverlap <- Omnivory <-  numeric(length(df$alphaDiv))
  ind <- 1
  for (indSim in 1:nSim){
    cat(sprintf("indSim: %d \r",indSim))
    y_mat <- matrix(mat$y.raw[,indSim],nSp, nNodes)
    PA_mat <- y_mat > thr
    FW_ind <- mat$FW.raw[indSim]
    NF <- foodWebs$nutrientFeeders.list[[FW_ind]]
    D.withN <- matrix(0, nSp, nSp)
    D.withN[1:(nSp-1), 1:(nSp-1)] <- foodWebs$dietMatrix.list[[FW_ind]]
    D.withN[nSp, NF] <- 1
    bodymass <- c(foodWebs$bodymass.list[[FW_ind]], 2e-5)
    for (node in 1:nNodes){
      PA_vec <- which(PA_mat[,node])
      Local.D <- D.withN[PA_vec,PA_vec]
      Connectance[ind] <- sum(Local.D)/nrow(Local.D)^2
      LinkDensity[ind] <- sum(Local.D)/nrow(Local.D)
      Modularity[ind] <-  Max.modu(Local.D)
      Nestedness[ind] <- unname(unlist(unodf(Local.D, selfloop = T)[2]))
      NicheOverlap[ind] <- tryCatch(
        unname(networklevel(Local.D, index = "niche overlap")[2]),
        error = function(e) {NA}) 
      Local.TL.info = Unw.TL.Omn.q(Local.D)
      Omnivory[ind] <- mean(Local.TL.info$Omn, na.rm = T) # In terms of TL, not primary producer
      ind <- ind + 1
    }
  }
  cat("\n")
  df[["Connectance"]] <- Connectance
  df[["LinkDensity"]] <- LinkDensity
  df[["Modularity"]] <- Modularity
  df[["Nestedness"]] <- Nestedness
  df[["NicheOverlap"]] <- NicheOverlap
  df[["Omnivory"]] <- Omnivory
  save(df,file="FW_indices.rda")
  
} else {load("utilities/FW_indices.rda")}


# Fig. 1c - OCN ####
if (createPlots){
  pdf(file="Fig_1c.pdf",width=8/2.54,height=8/2.54)
  draw_thematic_OCN(position,OCN,discreteLevels = T, colPalette = kol)
  dev.off()
}

# Fig. 2 - PLOT OF SP. RICHNESS & 6 FW METRICS ####
if (createPlots){
  pdf(file="Fig2.pdf",width=20/2.54,height=12/2.54)
  layout(matrix(c(1,2,3,4,1,5,6,7), 2, 4, byrow = TRUE))
  boxplot(alphaDiv ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
          ylim=c(0,80),bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0); axis(2,pos=0.35)
  boxplot(Connectance ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
          ylim=c(0.1,0.4), bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0.1); axis(2,pos=0.35)
  boxplot(LinkDensity ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
           ylim=c(0,16),bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0); axis(2,pos=0.35,at=seq(0,16,2))
  boxplot(Modularity ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
          ylim=c(0,0.4),bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0); axis(2,pos=0.35,at=seq(0,0.4,0.1))
  boxplot(Nestedness ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
          ylim=c(0,0.8),bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0); axis(2,pos=0.35,at=seq(0,0.8,0.2))
  boxplot(NicheOverlap ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
          ylim=c(0.1,0.7),bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0.1); axis(2,pos=0.35,at=seq(0.1,0.7,0.2))
  boxplot(Omnivory ~ pos,df,subset= nT==2 & nL==2 & dR==2 & dB==2 & vUp==2,outline=F,col=kol,
          ylim=c(0.2,0.7),bty="n",xaxt="n",yaxt="n")
  axis(1,pos=0.2); axis(2,pos=0.35,at=seq(0.2,0.7,0.1))
  dev.off()
}

# Fig. 3 - EXPLAIN PATTERNS OF SPECIES RICHNESS ####
subset <- df$nL==2 & df$nT==2 & df$dB==2 & df$dR==2 & df$vUp==2
DF <- data.frame(alphaDiv = df$alphaDiv[subset], Connectance = df$Connectance[subset],
                 LinkDensity = df$LinkDensity[subset], Modularity = df$Modularity[subset],
                 Nestedness = df$Nestedness[subset], NicheOverlap = df$NicheOverlap[subset],
                Omnivory = df$Omnivory[subset],
                 FW=df$FW[subset], distToOutlet = rep(distToOutlet,100), A = rep(OCN$AG$A,100))
meanAlphaDiv <- meanConnectance <- meanLinkDensity <- meanModularity <-  numeric(100)
meanNestedness <- meanNicheOverlap <- meanOmnivory <-  numeric(100)
for(i in 1:100){
  meanAlphaDiv[i] <- mean(DF$alphaDiv[DF$FW==i])
  meanConnectance[i] <- mean(DF$Connectance[DF$FW==i])
  meanLinkDensity[i] <- mean(DF$LinkDensity[DF$FW==i])
  meanModularity[i] <- mean(DF$Modularity[DF$FW==i])
  meanNestedness[i] <- mean(DF$Nestedness[DF$FW==i])
  meanNicheOverlap[i] <- mean(DF$NicheOverlap[DF$FW==i])
  meanOmnivory[i] <- mean(DF$Omnivory[DF$FW==i])
    }
DF[["meanAlphaDiv"]] <- meanAlphaDiv[DF$FW]
DF[["meanConnectance"]] <- meanConnectance[DF$FW]
DF[["meanLinkDensity"]] <- meanLinkDensity[DF$FW]
DF[["meanModularity"]] <- meanModularity[DF$FW]
DF[["meanNestedness"]] <- meanNestedness[DF$FW]
DF[["meanNicheOverlap"]] <- meanNicheOverlap[DF$FW]
DF[["meanOmnivory"]] <- meanOmnivory[DF$FW]
DF[["FW_char"]] <- as.character(DF$FW)

evalFracVar <- function(field,DF){
  eval(parse(text=paste0("aa <- lm(", field, " ~ FW_char + distToOutlet + log10(A), DF)")))
  #aa <- lm(field ~ FW_char + distToOutlet + log10(A), DF)
  AA <- Anova(aa,type="III")
  fracVar <- AA$`Sum Sq`[2:5]/sum(AA$`Sum Sq`[2:5])
  p_values <- AA$`Pr(>F)`[2:4]
  coef <- aa$coefficients[101:102]
  out <- list(fracVar=fracVar, p_values=p_values, coef=coef)
}

fracVar <- data.frame(matrix(0,4,7), row.names=c("FW","dO","logA", "res"))
pVal <- data.frame(matrix(0,3,7), row.names=c("FW","dO","logA"))
coefMat <- data.frame(matrix(0,2,7), row.names=c("dO","logA"))
CV <- data.frame(matrix(0,1,7), row.names=c("CV"))
names(fracVar) <- names(pVal) <- names(coefMat) <- names(CV) <- 
  c("alphaDiv","Connectance","LinkDensity","Modularity","Nestedness","NicheOverlap","Omnivory")
for (i in 1:7){
  nam <- names(fracVar)[i]
  eval(parse(text=paste0("out <- evalFracVar('",nam,"',DF)")))
  fracVar[nam] <- out$fracVar*100
  pVal[nam] <- out$p_values
  coefMat[nam] <- out$coef
  CV[nam] <- sd(DF[[nam]],na.rm=T)/mean(DF[[nam]],na.rm=T)
}

singleMFW_trends <- vector("list",7)
names(singleMFW_trends) <- names(fracVar)
for (j in 1:7){
  singleMFW_trends[[j]] <- data.frame(matrix(0,100,4))
  names(singleMFW_trends[[j]]) <- c("slope_dO","p_dO","slope_logA","p_logA")
}
#look at trends within meta-food webs
for (i in 1:100){
  for (j in 1:7){
    nam <- names(fracVar)[j]
  eval(parse(text=paste0("ss <- summary(lm(",nam," ~ distToOutlet + log(A), DF, subset = which(DF$FW==i)))")))
  singleMFW_trends[[j]][i,] <- c(ss$coefficients[2,1], ss$coefficients[2,4],  ss$coefficients[3,1], ss$coefficients[3,4])
  }
}

singleTrends_dO <- singleTrends_logA <-data.frame(matrix(0,4,7))
names(singleTrends_dO) <- names(singleTrends_logA) <- names(fracVar)
for (j in 1:7){
  nam <- names(fracVar)[j]
  singleTrends_dO[1,j] <-  sum(singleMFW_trends[[nam]][,1] > 0 &  singleMFW_trends[[nam]][,2] < 0.05) # significant pos
  singleTrends_dO[2,j] <-  sum(singleMFW_trends[[nam]][,1] > 0 &  singleMFW_trends[[nam]][,2] > 0.05) # non-sig pos
  singleTrends_dO[3,j] <-  sum(singleMFW_trends[[nam]][,1] < 0 &  singleMFW_trends[[nam]][,2] > 0.05) # non-sig neg
  singleTrends_dO[4,j] <-  sum(singleMFW_trends[[nam]][,1] < 0 &  singleMFW_trends[[nam]][,2] < 0.05) # significant neg
  singleTrends_logA[1,j] <-  sum(singleMFW_trends[[nam]][,3] > 0 &  singleMFW_trends[[nam]][,4] < 0.05) # significant pos
  singleTrends_logA[2,j] <-  sum(singleMFW_trends[[nam]][,3] > 0 &  singleMFW_trends[[nam]][,4] > 0.05) # non-sig pos
  singleTrends_logA[3,j] <-  sum(singleMFW_trends[[nam]][,3] < 0 &  singleMFW_trends[[nam]][,4] > 0.05) # non-sig neg
  singleTrends_logA[4,j] <-  sum(singleMFW_trends[[nam]][,3] < 0 &  singleMFW_trends[[nam]][,4] < 0.05) # significant neg
}

if(createPlots){
pdf(file="Fig3.pdf",width=20/2.54,height=8/2.54)
par(mfrow=c(1,2))
barplot(as.matrix(singleTrends_dO), col=hcl.colors(4,"Blue-Red"), 
        main="Effect of distance to outlet",ylab="Fraction of meta-food web realizations")
abline(h=c(5,95),col="#909090")
barplot(as.matrix(singleTrends_logA), col=hcl.colors(4,"Blue-Red"), main="Effect of drainage area")
abline(h=c(5,95),col="#909090")
dev.off()
}

# FIG. 4 - SPECIES RICHNESS SENSITIVITY ####
if(createPlots){
  df[["alphaDiv_norm"]] <- df$alphaDiv/meanAlphaDiv[df$FW]
  pdf(file="Fig4.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(alphaDiv_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F,ylim=c(0,3), 
          xlab = "", ylab="alphaDiv_norm",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35)
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(alphaDiv_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F,ylim=c(0,2), 
          xlab = "", ylab="alphaDiv_norm",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35)
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(alphaDiv_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F,  ylim=c(0,2),
          xlab = "", ylab="alphaDiv_norm",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35)
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(alphaDiv_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="alphaDiv_norm",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35)
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(alphaDiv_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35)
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}

# FIG. 5 - CONNECTANCE SENSITIVITY ####
if(createPlots){
  df[["Connectance_norm"]] <- df$Connectance/meanConnectance[df$FW]
  pdf(file="Fig5.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(Connectance_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0.5,2), 
          xlab = "", ylab="Connectance (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0.5, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0.5,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Connectance_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0.5,1.5),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0.5, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35,at=seq(0.5,1.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Connectance_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F, ylim=c(0.5,1.5),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0.5, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0.5,1.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Connectance_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0.5,1.5),
          xlab = "", ylab="Connectance (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0.5, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0.5,1.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Connectance_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0.5,1.5),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0.5, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0.5,1.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}

# run RMW model ####
# (this takes ~1 day to run!)
if(!file.exists("utilities/results_RMW.rda")){
  subset <- df$nL==2 & df$nT==2 & df$dB==2 & df$dR==2 & df$vUp==2
  df.Resh <- df[subset,] # only consider default simulation
  
  nResample <- 100
  Connectance.Resh.Mean <- Connectance.Resh.SD <- numeric(length(df.Resh$alphaDiv))
  LinkDensity.Resh.Mean <- LinkDensity.Resh.SD <- numeric(length(df.Resh$alphaDiv))
  Omnivory.Resh.Mean <- Omnivory.Resh.SD <- numeric(length(df.Resh$alphaDiv))
  Nestedness.Resh.Mean <- Nestedness.Resh.SD <- numeric(length(df.Resh$alphaDiv))
  Modularity.Resh.Mean <- Modularity.Resh.SD <- numeric(length(df.Resh$alphaDiv))
  NicheOverlap.Resh.Mean <- NicheOverlap.Resh.SD <- numeric(length(df.Resh$alphaDiv))
  
  t0 <- Sys.time()
  set.seed(22)
  for (i in 1:length(df.Resh$alphaDiv)){
    alphaDiv_val <- df.Resh$alphaDiv[i] 
    Connectance.Resh <- LinkDensity.Resh <- Omnivory.Resh <- Nestedness.Resh <- Modularity.Resh <- NicheOverlap.Resh <- numeric(nResample)
    indFW <- df.Resh$FW[i]
    nutFeed <- df.Resh$nutFeed[i]
    for (j in 1:nResample){
      bodymass <- c(foodWebs$bodymass.list[[indFW]], 2e-5)
      nF <- sample(foodWebs$nutrientFeeders.list[[indFW]], nutFeed) # sample nut feeders
      ss <- sample(setdiff(1:99,nF),alphaDiv_val-length(nF)-1) # sample non-nut feeders (-1 to exclude nutrients, which are not in the diet matrix)
      ss <- c(nF,ss) # merge subsets - nut feeders are now first, but it doesn't matter as long as the order is kept in rows and columns
      PA <- numeric(length(bodymass)); PA[ss] <- 1
      
      D <- foodWebs$dietMatrix.list[[indFW]][ss,ss]
      D_withN <- rbind(cbind(D, 0), 0)
      D_withN[alphaDiv_val, 1:length(nF)] <- 1 # add nutrient link
      Local.TL.info = Unw.TL.Omn.q(D_withN)
      
      Connectance.Resh[j] <- sum(D_withN)/alphaDiv_val^2
      LinkDensity.Resh[j] <- sum(D_withN)/alphaDiv_val
      Omnivory.Resh[j] <-  mean(Local.TL.info$Omn, na.rm = T)
      Nestedness.Resh[j]<-  unname(unlist(unodf(D_withN, selfloop = T)[2]))
      Modularity.Resh[j] <- Max.modu(D_withN)
      NicheOverlap.Resh[j] <- tryCatch(unname(networklevel(D_withN, index = "niche overlap")[2]), error = function(e) {NA}) 
    }
    Connectance.Resh.Mean[i] <- mean(Connectance.Resh);  Connectance.Resh.SD[i] <- sd(Connectance.Resh)
    LinkDensity.Resh.Mean[i] <- mean(LinkDensity.Resh);  LinkDensity.Resh.SD[i] <- sd(LinkDensity.Resh)
    Omnivory.Resh.Mean[i] <- mean(Omnivory.Resh);  Omnivory.Resh.SD[i] <- sd(Omnivory.Resh)
    Nestedness.Resh.Mean[i] <- mean(Nestedness.Resh);  Nestedness.Resh.SD[i] <- sd(Nestedness.Resh)
    Modularity.Resh.Mean[i] <- mean(Modularity.Resh);  Modularity.Resh.SD[i] <- sd(Modularity.Resh)
    NicheOverlap.Resh.Mean[i] <- mean(NicheOverlap.Resh);  NicheOverlap.Resh.SD[i] <- sd(NicheOverlap.Resh) 
    
    cat(sprintf('Percent done: %.4f%%   -  Elapsed time: %.1f s \r',i/length(df.Resh$alphaDiv)*100,difftime(Sys.time(),t0, units="secs")))
  }
  
  df.Resh$Connectance.Resh.Mean <- Connectance.Resh.Mean;   df.Resh$Connectance.Resh.SD <- Connectance.Resh.SD         
  df.Resh$LinkDensity.Resh.Mean <- LinkDensity.Resh.Mean;   df.Resh$LinkDensity.Resh.SD <- LinkDensity.Resh.SD
  df.Resh$Omnivory.Resh.Mean <- Omnivory.Resh.Mean;         df.Resh$Omnivory.Resh.SD <- Omnivory.Resh.SD
  df.Resh$Nestedness.Resh.Mean <- Nestedness.Resh.Mean;     df.Resh$Nestedness.Resh.SD <- Nestedness.Resh.SD
  df.Resh$Modularity.Resh.Mean <- Modularity.Resh.Mean;     df.Resh$Modularity.Resh.SD <- Modularity.Resh.SD
  df.Resh$NicheOverlap.Resh.Mean <- NicheOverlap.Resh.Mean; df.Resh$NicheOverlap.Resh.SD <- NicheOverlap.Resh.SD
  save(df.Resh,file="utilities/results_RMW.rda")
  
} else {load("utilities/results_RMW.rda")}

# read results from UWB model  ####
matLoc <- readMat("utilities/resultsUWB_for_R.mat") # load simulation results
dfLoc <- data.frame(alphaDiv = matLoc$alphaDiv.vec,  alphaDiv_norm = matLoc$alphaDivNorm.vec, nL = matLoc$nL.vec, 
                    FW = matLoc$FW.vec)

Connectance.Loc <- LinkDensity.Loc <- Omnivory.Loc <- Nestedness.Loc <- Modularity.Loc <- NicheOverlap.Loc <- MeanBodyMass.Loc <- numeric(length(dfLoc$alphaDiv))
for (i in 1:length(dfLoc$alphaDiv)){
  indFW <- dfLoc$FW[i]
  bodymass <- c(foodWebs$bodymass.list[[indFW]], 2e-5)
  alphaDiv_val <- dfLoc$alphaDiv[i] 
  nutFeed <- foodWebs$nutrientFeeders.list[[indFW]]
  species <- which(matLoc$y.raw[1:99,i]>0) # do not sample nutrients, but add them later
  PA <- numeric(length(bodymass)); PA[species] <- 1
  
  D <- foodWebs$dietMatrix.list[[indFW]][species,species]
  D_withN <- rbind(cbind(D, 0), 0)
  D_withN[alphaDiv_val, match(nutFeed, species)] <- 1 # add nutrient link
  Local.TL.info = Unw.TL.Omn.q(D_withN)
  
  Connectance.Loc[i] <- sum(D_withN)/alphaDiv_val^2
  LinkDensity.Loc[i] <- sum(D_withN)/alphaDiv_val
  Omnivory.Loc[i] <-  mean(Local.TL.info$Omn, na.rm = T)
  Nestedness.Loc[i] <-  unname(unlist(unodf(D_withN, selfloop = T)[2]))
  Modularity.Loc[i] <- Max.modu(D_withN)
  NicheOverlap.Loc[i] <- tryCatch(unname(networklevel(D_withN, index = "niche overlap")[2]), error = function(e) {NA}) 
}

dfLoc$Connectance <- Connectance.Loc
dfLoc$LinkDensity <- LinkDensity.Loc 
dfLoc$Nestedness <- Nestedness.Loc 
dfLoc$Modularity <- Modularity.Loc
dfLoc$Omnivory <- Omnivory.Loc
dfLoc$NicheOverlap <- NicheOverlap.Loc


# Fig. 6: COMPARE DELTA SP.RICH & DELTA FW INDICES WITH UWB MODEL ####
if (!file.exists('utilities/delta_SFW_UWB.mat')){
deltaFW <- deltaPos <- deltaSpeciesRichness <- deltaConnectance <- deltaLinkDensity <- numeric(0)
deltaModularity <- deltaNestedness <- deltaNicheOverlap <- deltaOmnivory <- deltaMeanBodyMass <- numeric(0)

for (ind_FW in 1:100){
  subset <- which(df$FW==ind_FW & df$nL==2 & df$nT==2 & df$dR==2 & df$dB==2 & df$vU==2)
  deltaFW <- c(deltaFW, df$FW[subset])
  deltaSpeciesRichness <- c(deltaSpeciesRichness, df$alphaDiv[subset]-dfLoc$alphaDiv[dfLoc$FW==ind_FW])
  deltaConnectance <- c(deltaConnectance, df$Connectance[subset]-dfLoc$Connectance[dfLoc$FW==ind_FW])
  deltaLinkDensity <- c(deltaLinkDensity, df$LinkDensity[subset]-dfLoc$LinkDensity[dfLoc$FW==ind_FW])
  deltaModularity <- c(deltaModularity, df$Modularity[subset]-dfLoc$Modularity[dfLoc$FW==ind_FW])
  deltaNestedness <- c(deltaNestedness, df$Nestedness[subset]-dfLoc$Nestedness[dfLoc$FW==ind_FW])
  deltaNicheOverlap <- c(deltaNicheOverlap, df$NicheOverlap[subset]-dfLoc$NicheOverlap[dfLoc$FW==ind_FW])
  deltaOmnivory <- c(deltaOmnivory, df$Omnivory[subset]-dfLoc$Omnivory[dfLoc$FW==ind_FW])
  deltaMeanBodyMass <- c(deltaMeanBodyMass, df$MeanBodyMass[subset]-dfLoc$MeanBodyMass[dfLoc$FW==ind_FW])
  deltaPos <- c(deltaPos, df$pos[subset])
}
dfDelta <- data.frame(FW = deltaFW, pos = deltaPos,
  alphaDiv = deltaSpeciesRichness,
  Connectance = deltaConnectance,
  LinkDensity = deltaLinkDensity,
  Modularity = deltaModularity,
  Nestedness = deltaNestedness,
  NicheOverlap = deltaNicheOverlap,
  Omnivory = deltaOmnivory,
  MeanBodyMass = deltaMeanBodyMass)

writeMat("utilities/delta_SFW_UWB.mat",dfDelta=dfDelta)} 
# then run Matlab for the figure!


# Run RND model ####
# pick random species and re-create random diet matrix with same no. links
if(!file.exists("utilities/results_RND.rda")){
  subset <- df$nL==2 & df$nT==2 & df$dB==2 & df$dR==2 & df$vUp==2
  df.ReshAlt <- df[subset,] # only consider default simulation
  
  nResample <- 100
  Connectance.ReshAlt.Mean <- Connectance.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  LinkDensity.ReshAlt.Mean <- LinkDensity.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  Omnivory.ReshAlt.Mean <- Omnivory.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  Nestedness.ReshAlt.Mean <- Nestedness.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  Modularity.ReshAlt.Mean <- Modularity.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  NicheOverlap.ReshAlt.Mean <- NicheOverlap.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  MeanBodyMass.ReshAlt.Mean <- MeanBodyMass.ReshAlt.SD <- numeric(length(df.ReshAlt$alphaDiv))
  
  t0 <- Sys.time()
  set.seed(45)
  for (i in 1:length(df.ReshAlt$alphaDiv)){
    alphaDiv_val <- df.ReshAlt$alphaDiv[i] 
    Connectance.ReshAlt <- LinkDensity.ReshAlt <- Omnivory.ReshAlt <- Nestedness.ReshAlt <- Modularity.ReshAlt <- NicheOverlap.ReshAlt <- MeanBodyMass.ReshAlt <- numeric(nResample)
    indFW <- df.ReshAlt$FW[i]
    nutFeed <- df.ReshAlt$nutFeed[i]
    bodymass <- c(foodWebs$bodymass.list[[indFW]], 2e-5)
    nLinks <- df.ReshAlt$LinkDensity[i]*alphaDiv_val
    for (j in 1:nResample){
      
      ss <- sample(100,alphaDiv_val) # sample non-nut feeders (-1 to exclude nutrients, which are not in the diet matrix)
      PA <- numeric(length(bodymass)); PA[ss] <- 1
      
      Dnew <- matrix(0,alphaDiv_val,alphaDiv_val)
      while (all(Dnew[alphaDiv_val,1:(alphaDiv_val-1)]==0)  ){ # check that at least one species (different than the basal resource) consumes the basal resource
        newLinks <- sample(alphaDiv_val*(alphaDiv_val-1), nLinks)
        tmp <- numeric(alphaDiv_val*(alphaDiv_val-1)); tmp[newLinks] <- 1
        Dnew <- matrix(tmp,alphaDiv_val,alphaDiv_val-1)
        Dnew <- cbind(Dnew, numeric(alphaDiv_val)) # impose last species as basal resource
        if (any(colSums(Dnew)[1:(alphaDiv_val-1)]==0)){ # if there is a second basal resource, reject 
          Dnew <- matrix(0,alphaDiv_val,alphaDiv_val)
        }
      }
      
      Local.TL.info = Unw.TL.Omn.q(Dnew)
      
      Connectance.ReshAlt[j] <- sum(Dnew)/alphaDiv_val^2
      LinkDensity.ReshAlt[j] <- sum(Dnew)/alphaDiv_val
      Omnivory.ReshAlt[j] <-  mean(Local.TL.info$Omn, na.rm = T)
      Nestedness.ReshAlt[j]<-  unname(unlist(unodf(Dnew, selfloop = T)[2]))
      Modularity.ReshAlt[j] <- Max.modu(Dnew)
      NicheOverlap.ReshAlt[j] <- tryCatch(unname(networklevel(Dnew, index = "niche overlap")[2]), error = function(e) {NA}) 
      MeanBodyMass.ReshAlt[j] <- PA %*% bodymass # unweighted body mass (for comparison with reshuffled model - does the comparison make sense?)
      
    }
    Connectance.ReshAlt.Mean[i] <- mean(Connectance.ReshAlt);  Connectance.ReshAlt.SD[i] <- sd(Connectance.ReshAlt)
    LinkDensity.ReshAlt.Mean[i] <- mean(LinkDensity.ReshAlt);  LinkDensity.ReshAlt.SD[i] <- sd(LinkDensity.ReshAlt)
    Omnivory.ReshAlt.Mean[i] <- mean(Omnivory.ReshAlt);  Omnivory.ReshAlt.SD[i] <- sd(Omnivory.ReshAlt)
    Nestedness.ReshAlt.Mean[i] <- mean(Nestedness.ReshAlt);  Nestedness.ReshAlt.SD[i] <- sd(Nestedness.ReshAlt)
    Modularity.ReshAlt.Mean[i] <- mean(Modularity.ReshAlt);  Modularity.ReshAlt.SD[i] <- sd(Modularity.ReshAlt)
    NicheOverlap.ReshAlt.Mean[i] <- mean(NicheOverlap.ReshAlt);  NicheOverlap.ReshAlt.SD[i] <- sd(NicheOverlap.ReshAlt) 
    MeanBodyMass.ReshAlt.Mean[i] <- mean(MeanBodyMass.ReshAlt);  MeanBodyMass.ReshAlt.SD[i] <- sd(MeanBodyMass.ReshAlt) 
    
    cat(sprintf('Percent done: %.4f%%   -  Elapsed time: %.1f s \r',i/length(df.ReshAlt$alphaDiv)*100,difftime(Sys.time(),t0, units="secs")))
  }
  
  df.ReshAlt$Connectance.ReshAlt.Mean <- Connectance.ReshAlt.Mean;   df.ReshAlt$Connectance.ReshAlt.SD <- Connectance.ReshAlt.SD         
  df.ReshAlt$LinkDensity.ReshAlt.Mean <- LinkDensity.ReshAlt.Mean;   df.ReshAlt$LinkDensity.ReshAlt.SD <- LinkDensity.ReshAlt.SD
  df.ReshAlt$Omnivory.ReshAlt.Mean <- Omnivory.ReshAlt.Mean;         df.ReshAlt$Omnivory.ReshAlt.SD <- Omnivory.ReshAlt.SD
  df.ReshAlt$Nestedness.ReshAlt.Mean <- Nestedness.ReshAlt.Mean;     df.ReshAlt$Nestedness.ReshAlt.SD <- Nestedness.ReshAlt.SD
  df.ReshAlt$Modularity.ReshAlt.Mean <- Modularity.ReshAlt.Mean;     df.ReshAlt$Modularity.ReshAlt.SD <- Modularity.ReshAlt.SD
  df.ReshAlt$NicheOverlap.ReshAlt.Mean <- NicheOverlap.ReshAlt.Mean; df.ReshAlt$NicheOverlap.ReshAlt.SD <- NicheOverlap.ReshAlt.SD
  df.ReshAlt$MeanBodyMass.ReshAlt.Mean <- MeanBodyMass.ReshAlt.Mean; df.ReshAlt$MeanBodyMass.ReshAlt.SD <- MeanBodyMass.ReshAlt.SD
  save(df.ReshAlt,file="utilities/results_RND.rda")
  
} else {load("utilities/results_RND.rda")}

# FIG. 7 - MULTICOMPARISON AMONG NULL MODELS ####
df.Long <- data.frame(Connectance = c(df.Resh$Connectance, df.Resh$Connectance.Resh.Mean, df.ReshAlt$Connectance.ReshAlt.Mean),
                      LinkDensity = c(df.Resh$LinkDensity, df.Resh$LinkDensity.Resh.Mean, df.ReshAlt$LinkDensity.ReshAlt.Mean),
                      Modularity = c(df.Resh$Modularity, df.Resh$Modularity.Resh.Mean, df.ReshAlt$Modularity.ReshAlt.Mean),
                      Nestedness = c(df.Resh$Nestedness, df.Resh$Nestedness.Resh.Mean, df.ReshAlt$Nestedness.ReshAlt.Mean),
                      NicheOverlap = c(df.Resh$NicheOverlap, df.Resh$NicheOverlap.Resh.Mean, df.ReshAlt$NicheOverlap.ReshAlt.Mean),
                      Omnivory = c(df.Resh$Omnivory, df.Resh$Omnivory.Resh.Mean, df.ReshAlt$Omnivory.ReshAlt.Mean),
                     pos= c(df.Resh$pos,df.Resh$pos,df.ReshAlt$pos),
                     model=c(rep("a_SFW",length(df.Resh$pos)), rep("b_RMW",length(df.Resh$pos)), rep("c_RND",length(df.Resh$pos))))

if (createPlots){
  pdf(file="Fig7.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(Connectance ~ pos + model,df.Long,outline=F,col=kol,xaxt="n",yaxt="n",ylim=c(0.1,0.35),
          ylab="",xlab="",bty="n",main="Connectance")
  axis(1,pos=0.1,at=seq(2.5,10.5,4)); axis(2,pos=0.35,at=seq(0.1,0.35,0.05)); abline(v=c(4.5,8.5))
  boxplot(LinkDensity ~ pos + model,df.Long,outline=F,col=kol,xaxt="n",yaxt="n",ylim=c(0,15),
          ylab="",xlab="",bty="n",main="Link density")
  axis(1,pos=0,at=seq(2.5,10.5,4)); axis(2,pos=0.35,at=seq(0,15,5)); abline(v=c(4.5,8.5))
  boxplot(Modularity ~ pos + model,df.Long,outline=F,col=kol,xaxt="n",yaxt="n",ylim=c(0,0.6),
          ylab="",xlab="",bty="n",main="Modularity")
  axis(1,pos=0,at=seq(2.5,10.5,4)); axis(2,pos=0.35,at=seq(0,0.6,0.2)); abline(v=c(4.5,8.5))
  boxplot(Nestedness ~ pos + model,df.Long,outline=F,col=kol,xaxt="n",yaxt="n",ylim=c(0,0.8),
          ylab="",xlab="",bty="n",main="Nestedness")
  axis(1,pos=0,at=seq(2.5,10.5,4)); axis(2,pos=0.35,at=seq(0,0.8,0.2)); abline(v=c(4.5,8.5))
  boxplot(NicheOverlap ~ pos + model,df.Long,outline=F,col=kol,xaxt="n",yaxt="n",ylim=c(0.1,0.7),
          ylab="",xlab="",bty="n",main="Niche overlap") 
  axis(1,pos=0.1,at=seq(2.5,10.5,4)); axis(2,pos=0.35,at=seq(0.1,0.7,0.2)); abline(v=c(4.5,8.5))
  boxplot(Omnivory ~ pos + model,df.Long,outline=F,col=kol,log="y",xaxt="n",yaxt="n",ylim=c(0.25,4),
          ylab="",xlab="",bty="n",main="Omnivory") 
  axis(1,pos=0.25,at=seq(2.5,10.5,4)); axis(2,pos=0.35,at=c(0.25,0.5,1,2,4)); abline(v=c(4.5,8.5))
  dev.off()
}

# FIG. S1 - SCATTERPLOT DIST OUTLET vs DRAINAGE AREA ####
if (createPlots){
  pdf(file="FigS1.pdf",width=12/2.54,height=12/2.54)
  plot(OCN$AG$A[high], distToOutlet[high], log="x", pch=19, col=kol[1], xlim=c(5e6, 5e9), ylim=c(0,80000), cex=0.5)
  points(OCN$AG$A[mid], distToOutlet[mid], pch=19, col=kol[2], cex=0.5)
  points(OCN$AG$A[low], distToOutlet[low], pch=19, col=kol[3], cex=0.5) 
  points(OCN$AG$A[down], distToOutlet[down], pch=19, col=kol[4], cex=0.5)
  abline(h=median(distToOutlet),col="#909090")
  abline(v=median(OCN$AG$A),col="#909090")
  dev.off()
}

# FIG. S2 - LINK DENSITY SENSITIVITY ####
if(createPlots){
  df[["LinkDensity_norm"]] <- df$LinkDensity/meanLinkDensity[df$FW]
  pdf(file="FigS2.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(LinkDensity_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2.5), 
          xlab = "", ylab="Link Density (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(LinkDensity_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(LinkDensity_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(LinkDensity_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="Link Density (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(LinkDensity_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}

# FIG. S3 - MODULARITY SENSITIVITY ####
if(createPlots){
  df[["Modularity_norm"]] <- df$Modularity/meanModularity[df$FW]
  pdf(file="FigS3.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(Modularity_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2), 
          xlab = "", ylab="Modularity (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Modularity_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Modularity_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Modularity_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="Modularity (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Modularity_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}

# FIG. S4 - NESTEDNESS SENSITIVITY ####
if(createPlots){
  df[["Nestedness_norm"]] <- df$Nestedness/meanNestedness[df$FW]
  pdf(file="FigS4.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(Nestedness_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2.5), 
          xlab = "", ylab="Nestedness (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Nestedness_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Nestedness_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F, ylim=c(0,2.5),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Nestedness_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2.5),
          xlab = "", ylab="Nestedness (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Nestedness_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}

# FIG. S5 - NICHE OVERLAP SENSITIVITY ####
if(createPlots){
  df[["NicheOverlap_norm"]] <- df$NicheOverlap/meanNicheOverlap[df$FW]
  pdf(file="FigS5.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(NicheOverlap_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,3), 
          xlab = "", ylab="Niche Overlap (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,3,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(NicheOverlap_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(NicheOverlap_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F, ylim=c(0,2.5),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(NicheOverlap_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0,2.5),
          xlab = "", ylab="Niche Overlap (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(NicheOverlap_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0,2),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}

# FIG. S6 - OMNIVORY SENSITIVITY ####
if(createPlots){
  df[["Omnivory_norm"]] <- df$Omnivory/meanOmnivory[df$FW]
  pdf(file="FigS6.pdf",width=20/2.54,height=15/2.54)
  par(mfrow=c(2,3))
  boxplot(Omnivory_norm ~ pos + nL, df, col=kol,main = "Effect of nutrient load",
          subset = nT==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0,1.5), 
          xlab = "", ylab="Omnivory (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,3,0.5))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Omnivory_norm ~ pos + nT, df, col=kol, main = "Effect of nutrient distribution",
          subset = nL==2 & dR==2 & dB==2 & vUp==2, outline=F, ylim=c(0.75,1.25),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0.75, at = seq(2.5,10.5,4), labels=c("Flat","Downstream","Random")); axis(2,pos=0.35,at=seq(0,2,0.25))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Omnivory_norm ~ pos + vUp, df, col=kol, main = "Effect of uptake velocity",
          subset = nL==2 & nT==2 & dR==2 & dB==2, outline=F, ylim=c(0.5,1.25),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0.5, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.25))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Omnivory_norm ~ pos + dR, df, col=kol, main = "Effect of dispersal rate",
          subset = nL==2 & nT==2 & dB==2 & vUp==2, outline=F, ylim=c(0.75,1.25),
          xlab = "", ylab="Omnivory (normalized)",xaxt="n",yaxt="n")
  axis(1, pos=0.75, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2.5,0.25))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  boxplot(Omnivory_norm ~ pos + dB, df, col=kol, main = "Effect of downstream bias",
          subset = nL==2 & nT==2 & dR==2 & vUp==2, outline=F, ylim=c(0.75,1.25),
          xlab = "", ylab="",xaxt="n",yaxt="n")
  axis(1, pos=0.75, at = seq(2.5,10.5,4), labels=c("Low","Medium","High")); axis(2,pos=0.35,at=seq(0,2,0.25))
  abline(v=c(4.5,8.5),col="#909090"); abline(h=1)
  dev.off()
}
