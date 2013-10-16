###Code for analysing the output of 'phylogeosim.R'###
# Copyleft (or the one to blame): Carvalho, LMF (2012)
#Last update: 13/09/2012
#Alignment stats

#Tree stats
balance(true_tree)
balance(inferred_tree)

#Tree Distances (inferred  versus simulated)

#
barplot(table(traitsvec),)
Q<-exp(q[[1]])/sum(exp(q[[1]]))
#chi-square to assess deviance from the expected
#
#mcmc.popsize 