# Fonction pour calculer les termes d'IBD dans la covariance entre scores pour diff�rentes cat�gories,
# � l'int�rieur d'une seule famille

# par Jordie Croteau
# 1er avril 2011

# correspond au fichier ibd_terms_v7.R dans le dossier "programmes"


# modifi�e le 21 mai pour permettre que les paires "un sujet avec lui-m�me" soit pr�sentes ou non dans les 
# donn�es d'IBD (�a peut m�me �tre une partie seulement de ces paires "un sujet avec lui-m�me" qui soit pr�sentes).  
# Dans tous les cas, on attribue � toutes ces paires une proportion d'IBD inf�r�e de 1.0 (peu importe si la paire 
# est pr�sente ou non dans les donn�es d'IBD).

# modifi�e le 14 juin 2012 pour lister les paires manquantes (lorsqu'il y en a).

# modifi� le 19 octobre: pour toute paire manquante, poser pim �gale � 0. 

ibd.terms=function(y,subject.ids,l1,l2,pim)
{
################################################################################################
# y : vecteur des cat�gories des sujets, valeur entre 1 et nlevels(y). 
# subject.ids: ID des sujets (doit correspondre aux ID dans l1 et l2)
# l1 : liste des 1ers sujets des paires 
# l2 : liste des 2e sujets des paires
# pim : matrice des proportion d'IBD inf�r�e pi entre sujets l1 et l2 pour chaque locus dans x.
#
# Attention, on suppose aussi que tous les individus sont g�notyp�s.
################################################################################################

#nombre de sujets dans la famille
n=length(y)

#################################### V�rifications pour s'assurer que les donn�es sont dans le bon format (avec modifs appliqu�es si n�cessaire et possible) ##############################################
if(!is.factor(y)) stop("y should be of type factor")

########## ajouter automatiquement les paires (i,i) manquantes (s'il y a lieu). ##########
# enlever d'abord les paires (i,i) d�j� l�, pour ne pas les d�doubler
paires.i.i.boo=apply(cbind(l1,l2),1,function(x) x[1]==x[2])
l1=l1[!paires.i.i.boo]
l2=l2[!paires.i.i.boo]
pim=array(pim[!paires.i.i.boo,],c(sum(!paires.i.i.boo),ncol(pim)))

l1=c(subject.ids,l1)
l2=c(subject.ids,l2)
pim=rbind(array(1,c(n,ncol(pim))),pim)
###########################################################################################

pairs.table=as.matrix(table(factor(l1,levels=subject.ids),factor(l2,levels=subject.ids)))
if(any(((pairs.table+t(pairs.table))[lower.tri(pairs.table)])>1)) stop("all pairs of subjects must be present only once in the ibd file")

# Pour les paires manquantes, poser la valeur de pim �gale � 0 (c'est d'ailleurs la bonne chose � faire m�me pour les colonnes correspondant aux produits des paires).
# Poser la valeur de pim �gale � 0 est la bonne chose � faire pour les paires de fondateurs par exemple (simwalk2 ne les met pas dans ses sorties, 
# et en plus simwalk2 omet toute paire dont P1=P2=0 � toutes les positions)
# Il se pourrait cependant que des paires soient absentes par erreur (par exemple si on oublie de mettre un statut "atteint" � tout le monde dans Simwalk2
# (simwalk2 fait les ibd seulement entre les sujets atteints).  On met donc un avertissement donnant la liste des paires manquantes.

################## Faire d'abord s�par�ment pour chaque colonne de pim car il pourrait y avoir des diff�rences d'un locus � l'autre #####################
for(j in 1:ncol(pim))
 {
  pi.tmp=pim[,j]
  paires.pi.NA=is.na(pi.tmp)
  if(sum(paires.pi.NA)>0)
   {
    l1.pi.NA=l1[paires.pi.NA]
    l2.pi.NA=l2[paires.pi.NA]
    pi.tmp[paires.pi.NA]=0
	pim[,j]=pi.tmp
    missing.pairs=paste(l1.pi.NA,l2.pi.NA)
    cat("IBD data missing for column",j,"of pim, for the following pairs of subjects \n",paste(missing.pairs,"\n"))
    cat("Posterior probabilities of IBD P1 and P2 have been set to 0 for these pairs ! \n \n")
   }
 }
#########################################################################################################################################################
 
pairs.table.mod=pairs.table+t(pairs.table)
pairs.table.mod[lower.tri(pairs.table)]=1
which.missing=which(pairs.table.mod==0,arr.ind=TRUE)
if(nrow(which.missing)!=0)
 {
  missing1=subject.ids[which.missing[,1]]
  missing2=subject.ids[which.missing[,2]]
  missing.pairs=paste(missing1,missing2)

  l1=c(missing1,l1)
  l2=c(missing2,l2)
  pim=rbind(array(0,c(length(missing1),ncol(pim))),pim)

  cat("IBD data missing, from all of the IBD files, for the following",length(missing.pairs),"pairs of subjects (out of",n*(n-1)/2,"possible pairs): \n",paste(missing.pairs,"\n"))
  cat("Posterior probabilities of IBD P1 and P2 have been set to 0 for these pairs ! \n")
 }
############################################################################################################################################################################################################

## D�termination des listes de sujets dans chaque cat�gorie
liste.par.cat=tapply(1:n,y,function (vec) subject.ids[vec],simplify=FALSE)
 
## Calcul du nombre de sujets par cat�gorie
ny=table(y)

n.levels=nlevels(y)
# Les dimensions de l'array des r�sultats sont (nombre de locus) * K-1 * K-1
res=array(NA,c(ncol(pim),n.levels-1,n.levels-1))

# afin de simplifier les calculs suivants, d�doubler les paires (l1[j],l2[j]) (et leur pim[j,] correspondant) telles que l1[j]!=l2[j]
# en mettant (l1[j],l2[j]) dans l'ordre inverse: (l2[j],l1[j]). Ceci est dans le but de ne jamais manquer rien des pi_i.j et pi_j.i
# car seulement un des 2 n'appara�t dans le fichier des ibd.
l1.double=c(l1,l2[l1!=l2])
l2.double=c(l2,l1[l1!=l2])
pim=rbind(pim,array(pim[l1!=l2,],c(sum(l1!=l2),ncol(pim))))

for (k in 1:(n.levels-1))
  {
	# On fait les calculs seulement s'il y a au moins un sujet dans la cat�gorie k
	if (length(liste.par.cat[[k]])>0)
	 {
	   #premier terme de somme de pi's
	   dims=c(sum(l1.double %in% liste.par.cat[[k]] & l2.double %in% liste.par.cat[[k]]),ncol(pim))
	   term1=((n-ny[k])^2)*apply(array(pim[l1.double %in% liste.par.cat[[k]] & l2.double %in% liste.par.cat[[k]],],dims),2,sum) 

	   #deuxi�me terme de somme de pi's
	   dims=c(sum(!(l1.double %in% liste.par.cat[[k]]) & !(l2.double %in% liste.par.cat[[k]])),ncol(pim))
	   term2=(ny[k]^2)*apply(array(pim[!(l1.double %in% liste.par.cat[[k]]) & !(l2.double %in% liste.par.cat[[k]]),],dims),2,sum)

	   #troisi�me terme de somme de pi's
	   dims=c(sum(l1.double %in% liste.par.cat[[k]] & !(l2.double %in% liste.par.cat[[k]])),ncol(pim))
	   term3=2*ny[k]*(n-ny[k])*apply(array(pim[l1.double %in% liste.par.cat[[k]] & !(l2.double %in% liste.par.cat[[k]]),],dims),2,sum)
	   
	   res[,k,k]=(term1+term2-term3)/(n^2)
     }	 
  }

if(n.levels>2)
  {
   for (k in 2:(n.levels-1))
    {
      for (l in 1:(k-1))
       {
  	     # On fait les calculs seulement s'il y a au moins un sujet dans chacune des cat�gories k et l
	     if (length(liste.par.cat[[k]])>0 & length(liste.par.cat[[l]])>0)
          {
		   dims=c(sum(l1.double %in% liste.par.cat[[k]] & l2.double %in% liste.par.cat[[l]]),ncol(pim))
           term1=(n-ny[k])*(n-ny[l])*apply(array(pim[l1.double %in% liste.par.cat[[k]] & l2.double %in% liste.par.cat[[l]],],dims),2,sum)
		   dims=c(sum(!(l1.double %in% liste.par.cat[[k]]) & !(l2.double %in% liste.par.cat[[l]])),ncol(pim))		   
           term2=ny[k]*ny[l]*apply(array(pim[!(l1.double %in% liste.par.cat[[k]]) & !(l2.double %in% liste.par.cat[[l]]),],dims),2,sum)
		   dims=c(sum(l1.double %in% liste.par.cat[[k]] & !(l2.double %in% liste.par.cat[[l]])),ncol(pim))		   
           term3=ny[l]*(n-ny[k])*apply(array(pim[l1.double %in% liste.par.cat[[k]] & !(l2.double %in% liste.par.cat[[l]]),],dims),2,sum)
		   dims=c(sum(!(l1.double %in% liste.par.cat[[k]]) & l2.double %in% liste.par.cat[[l]]),ncol(pim))		   
           term4=ny[k]*(n-ny[l])*apply(array(pim[!(l1.double %in% liste.par.cat[[k]]) & l2.double %in% liste.par.cat[[l]],],dims),2,sum)
           res[,k,l]=res[,l,k]=(term1+term2-term3-term4)/(n^2)
          }
	   }
    }
  }
  
res
}

