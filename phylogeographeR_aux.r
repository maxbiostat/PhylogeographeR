###Companion of main script 'phylogeosim.R'###
# Copyleft (or the one to blame): Carvalho, LMF (2012)
# Last update: 31/10/2012
#######################################
#Loading packages
get.packs<- function(packs,repos){
  for (i in 1:length(packs)){  
    if (!packs[i] %in% installed.packages()) {
      install.packages(packs[i],repos=repos)
    }
    library(packs[i],character.only = TRUE)
  } 
}
packs<-c("geiger","ape","phybase","phyclust","fields","spdep","RColorBrewer")
get.packs(packs,repos="http://cran.us.r-project.org")
bigpal<-c(brewer.pal(8,"Accent"),brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))
models<-c("HKY", "F84", "GTR", "JTT", "WAG", "PAM", "BLOSUM", "MTREV", "CPREV", "GENERAL")
#______________________________________________________#
#                 Auxiliary functions
#______________________________________________________#
#Building the CTMC's infinitesimal generator (Q)
#Transform into rate matrix
reg.matrix<-function(m){
  for (i in 1:ncol(m)){
    m[i,i]<-sum(m[i,][-i])*-1    
  }
  return(m)
}
#Migration
build.ss.mat<-function(neigh.mat,sources,sinks,rho,contiguity=FALSE,bidirectional=TRUE){
  K<-ncol(neigh.mat)
  ss.mat<-matrix(NA,ncol=K,nrow=K)
  ss.pos<-expand.grid(sources,sinks)
   for (i in 1:nrow(ss.pos)){
     if(bidirectional=="TRUE"){
     ss.mat[ss.pos[i,1],ss.pos[i,2]]<- 1/rho*K
     ss.mat[ss.pos[i,2],ss.pos[i,1]]<- rho/K
     } else{ss.mat[ss.pos[i,2],ss.pos[i,1]]<- rho/K}
   }  
  if(contiguity=="TRUE"){ss.mat[which(is.na(ss.mat))]<- (1/K)*neigh.mat[which(is.na(ss.mat))]
  }else ss.mat[which(is.na(ss.mat))]<- 1/K
return(ss.mat)}
#Build Q
mount.rate.matrix<-function(K,rho,r,map,sources=NULL,sinks=NULL,bidirectional=FALSE,contiguity=FALSE,geo.model,variables=NULL,method=NULL){  
  if(!is.null(map)){neigh.mat<-nb2mat(poly2nb(map),style="B");K<-ncol(neigh.mat)}
  if(!is.null(geo.model) && geo.model %in% c("patchy","distance","source-sink","gravity")){
    d.mat<-rdist.earth(coordinates(map));d.mat<-(1/d.mat);diag(d.mat)<-0 
    if(geo.model=="distance"){
      M<-reg.matrix(d.mat);cat("Distance-informed model","\n")
    } else if(geo.model=="gravity"){
      full.mat<-matrix(rep(1,K^2),ncol=K);diag(full.mat)<-0
      full.nb<-mat2listw(full.mat)$neighbours
      if(method=="mahalanobis"){
        cov<-var(variables)
        g.mat<-as.matrix(as.spam.listw(nb2listw(full.nb, nbcosts(full.nb, variables,method=method,cov=cov), style="B")))
      } else g.mat<-as.matrix(as.spam.listw(nb2listw(full.nb, nbcosts(full.nb, variables,method=method), style="B")))
      if(contiguity=="TRUE"){M<-reg.matrix(g.mat*neigh.mat)} else {M<-reg.matrix(g.mat)}
      cat("Gravity model","\n")
    } else if(!is.null(sources) && !is.null(sinks) && geo.model=="source-sink"){
      M<-reg.matrix(build.ss.mat(neigh.mat=neigh.mat,sources=sources,sinks=sinks,bidirectional=bidirectional,rho=rho,contiguity=contiguity))
      cat("Source-sink model","\n")
    }else if(geo.model=="patchy"){      
      if(!is.null(r) && r==1){
        M<-reg.matrix(neigh.mat*d.mat)
      }else{ knn.r <- knearneigh(as.matrix(coordinates(map)),r)
             neigh.mat.r<-nb2mat(knn2nb(knn.r),style="B")
             M<-reg.matrix(neigh.mat.r*d.mat)
      }
      cat(paste(" 'Patchy' model r=",r,sep=""),"\n")
    }
  }else{ if(contiguity=="TRUE"){
    M<-reg.matrix(matrix(rep(1/K,K^2),ncol=K)*neigh.mat);cat("Contiguity-constrained homogeneous model","\n")
  }else { M<-reg.matrix(matrix(rep(1/K,K^2),ncol=K));cat("Not using geography -- homogeneous model","\n")}
  }
  return(list(M))}
######
#Modified sim.char from 'geiger' package (sim_char) function. The original one crashes with large datasets (under Windows only -- Lapack error).
# sim_char<-function (phy, model.matrix, nsims = 1, model = "brownian", root.state = 1) 
# {
#   phy <- new2old.phylo(phy)
#   nchar <- nrow(model.matrix)
#   m <- get.simulation.matrix(phy)
#   nbranches <- ncol(m)
#   nspecies <- nrow(m)
#   if (model == "brownian" | model == "speciational") {
#     if (model == "speciational") {
#       m[m > 0] <- 1
#     }
#     rnd <- t(mvrnorm(nsims * nbranches, mu = rep(0, nchar), 
#                      Sigma = model.matrix))
#     rnd <- array(rnd, dim = c(nchar, nbranches, nsims))
#     simulate <- function(v) m %*% as.matrix(v)
#     result <- apply(rnd, 1, simulate)
#     result <- aperm(array(result, dim = c(nspecies, nsims, 
#                                           nchar)), c(1, 3, 2))
#     rownames(result) <- phy$tip.label
#   }
#   if (model == "discrete") {
#     nchar <- length(model.matrix)
#     node.value <- numeric(nbranches)
#     result <- array(0, dim = c(nspecies, nchar, nsims))
#     for (j in 1:nchar) {
#       m <- model.matrix[[j]]
#       for (k in 1:nsims) {
#         for (i in 1:nbranches) {
#           if (as.numeric(phy$edge[i, 1]) == min(phy$edge[,1])) 
#             s <- root.state
#           else {
#             parent <- which(phy$edge[, 2] == phy$edge[i,1])
#             s <- node.value[parent]
#           }
#           p <- matexpo(m * phy$edge.length[i]) #here I substituted the old MatrixExp() by ape's matexpo()
#           probs <- cumsum(p[s,])
#           r <- runif(1)
#           node.value[i] <- min(which(r < probs))
#         }
#         result[, j, k] <- node.value[as.numeric(phy$edge[,2]) > 0]
#       }
#     }
#     rownames(result) <- phy$tip.label
#   }
#   return(result)
# }
# environment(sim_char)<-environment(geiger::sim.char)