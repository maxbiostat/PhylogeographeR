# Code for simulating discrete character states and generating alignments
# Intended for assessing performance of phylogeographic inference
# Copyleft (or the one to blame): Carvalho, LMF (2012)
# Last update: 05/08/2013 
######################################
#Loading
#setwd("~/Dropbox/SAMPLING_BIAS/PHYLOGEOSIM/")
source("~/Dropbox/CTMC/CODE/phylogeographeR_alt_aux.R")
setwd("~/TEST/")
#______________________________________________________________________________________________________#
##################################Simulation definitions################################################
#__General
seed<-930982
set.seed(seed)
simu.no<-"N50K5"
dir.create(paste("simu_",simu.no,sep=""))
setwd(paste("simu_",simu.no,sep=""))
ntax<-50
#__Sequence-generating options
model<-"HKY"
seqlength<-1000
alpha<-0.92
rcat<-4
pinv<-0.15
kappa<-7.61
#__Tree settings
birth_rate = 0.45
death_rate = 0.45
growth_rate=birth_rate-death_rate
################################################################
#__Phylogeography settings
#readShapePoly("~/Dropbox/DATA/SHPs/orcounty.shp")
# Mapp<-readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))
# area.names<-Mapp$names
Mapp<-NULL
root_state<-4 #match("Donegal",area.names) # if a non-homogeneous geo.model is activated, use the number Id of the area (region) of interest
geo.model<-"homogeneous"
K<-15
r<-NULL
rho<-NULL#15 # remember 1<rho<=K
sources<-NULL#c(22)
sinks<-NULL#c(9)
bidirectional<-FALSE
variables<-NULL
method<-NULL
contiguity<-TRUE
plot<-TRUE
###############################################################
# A glimpse on the neighborhood structure
if(!is.null(Mapp)){
  summary(poly2nb(Mapp))  
  layout(matrix(c(0,1,0,2,2,2),2,3,byrow=T),TRUE)
  hist(colSums(nb2mat(poly2nb(Mapp),style="B")),xlab="Number of neighbours",main="Neighborhood structure")
  plot(Mapp,col=bigpal[1:length(Mapp)])
  plot(poly2nb(Mapp),coordinates(Mapp),lwd=2,add=T)
}
###############################################################
###############################################################
#__Writing your specs to a .log file
simu.log<-as.data.frame(rbind(Simu_Number=simu.no,Random_Seed=seed,Taxa_number=ntax,Seq_Length=seqlength,model=model,Alpha=alpha,Kappa=kappa,
                              Rate_Categories=rcat,
                              Proportion_Invariant=pinv,
                              Birth_Rate=birth_rate,
                              Death_rate=death_rate,
                              Growth_Rate=growth_rate,
                              Number_of_states=K,
                              Root_state=LETTERS[root_state]),
                        GEO.MODEL=geo.model);names(simu.log)<-"Value"
write.table(simu.log,paste("simu","_",simu.no,".log",sep=""),sep="\t")
#_______________________________________________________________________________________________________#
###Simulating the tree
tree<-sim.bdtree(b=birth_rate, d=death_rate, stop="taxa", n=ntax, seed=seed, extinct=FALSE)
####Now Simulating the discrete traits for each tip
if(!is.null(Mapp)){trait.names<-area.names}else{trait.names<-paste(LETTERS[1:K])}
q<-reg.matrix(matrix(rep(1/K,K^2),ncol=K))
#q<-mount.rate.matrix(K=K,r=r,rho=rho,sources=sources,sinks=sinks,
#bidirectional=bidirectional,map=Mapp,contiguity=contiguity,
#variables=variables,method=method,geo.model=geo.model)
q<-list(q)
summary(as.vector(q[[1]]))
traits <- sim.char(tree, q, model = "discrete",root=root_state, n =1)
traits_ <- as.data.frame(traits[,,1])
traitsvec <- trait.names[traits_[,1]]
tree$tip.label<-paste(tree$tip.label,"_",traitsvec,sep="")
pdf(file=paste("complete_tree","_",simu.no,".pdf",sep=""))
plot(tree,edge.width=3,show.tip.label=F);title("Complete History")#;windows()
tiplabels(pch=16,col=bigpal[match(traitsvec,trait.names)],cex=1.5)
dev.off()
tree_p <- drop.fossil(tree)
#
library(phytools)
Qest<-ace(traitsvec,tree,type="discrete",model="ER");Qest
#Q_p<-ace(country.vec,tree_p,type="discrete");Q
write.nexus(tree,file=paste("simu","_",simu.no,"_complete",".tree",sep=""))
###Simulate the sequences
seqgen.opts<-paste(  #building the options for seqgen
  paste("-m",models[match(model,models)],sep=""),
  paste("-l",seqlength,sep=""),
  paste("-a",alpha,sep=""),
  paste("-g",rcat,sep=""),
  paste("-i",pinv,sep=""),
  paste("-t",kappa,sep=""),
  "-on",
  sep=" ")
dna_np<-seqgen(opts =seqgen.opts, newick.tree =write.tree(tree))
write(dna_np,file = paste("simu","_",simu.no,"_all_seqs",".nex",sep=""))
if(!is.null(tree_p)){
  pdf(file=paste("extant_only","_",simu.no,".pdf",sep=""))
  plot(tree_p,edge.width=3,cex=.7,show.tip.label=F);title("Extant Only")
  tiplabels(pch=16,col=bigpal[match(traitsvec[match(tree_p$tip.label,tree$tip.label)],trait.names)],cex=1.5);title("Extant Only") 
  dev.off()
  write.nexus(tree_p,file=paste("simu","_",simu.no,".tree",sep=""))
  dna<-seqgen(opts =seqgen.opts, newick.tree =write.tree(tree_p))
  write(dna,file = paste("simu","_",simu.no,".nex",sep=""))
} 
#______________________________________#
#######Computing Summary Statistics#####
simu.stats<-list(
  LOG=simu.log,
  INI_FREQ=states_freqs_total<-table(traitsvec)/length(tree$tip.label),
  EXT_FREQ=states_freqs_extant<-table(traitsvec[match(tree_p$tip.label,tree$tip.label)])/ntax,
  Q=Q<-q[[1]],
  Q_exp=Q_stoch<-MatrixExp(q[[1]])
)
save(simu.stats,file=paste("statistics_simu","_",simu.no,".RData",sep=""))
simu.stats
########################################
if(!is.null(Mapp) && plot=="TRUE"){
  pdf(file=paste("map_and_tree","_",simu.no,".pdf",sep=""))  
  layout(matrix(c(1,1,2,2),2,2,byrow=T),TRUE)
  plot(Mapp,col=bigpal[1:length(Mapp)])
  title(paste("Origin =",trait.names[root_state]))
  plot(poly2nb(Mapp),coordinates(Mapp),pch=16,lty=1,add=T)
  legend(x="bottomleft",col=bigpal[match(names(which(simu.stats$INI_FREQ>0)),trait.names)],legend=paste(names(which(simu.stats$INI_FREQ>0))),pch=16,cex=.7,pt.cex=1.5)
  plot(tree,edge.width=3,show.tip.label=F)
  tiplabels(pch=16,col=bigpal[match(traitsvec,trait.names)],cex=1.5)
  dev.off()
}
########################################
graphics.off()
