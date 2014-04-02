# Fonction globale qui calcule les scores et les covariances pour toutes les familles et tous les locus.
# ATTENTION, M�ME SI CERTAINES PARTIES DE CE PROGRAMME FONCTIONNENT POUR UN NOMBRE ARBITRAIRE DE LOCUS,
# D'AUTRES PARTIES NE SONT VALIDES QUE POUR UN OU 2 LOCUS.  PAR EXEMPLE, LES LIGNES 255, 258 ET 283 DANS LA VERSION 12, 
# POUR LE MOMENT (EN DATE DU 18 MAI 2012) NE FONCTIONNENT BIEN QUE POUR LES CAS DE 1 OU 2 LOCUS.

# correspond au fichier fonction_scores_covs_v2.R dans le dossier "programmes"


# Jordie Croteau
# 11 mars 2011

# Modifi�e par Alexandre Bureau
# 24 mars 2011

# Ajout de l'argument design.constraint pour cr�er les matrices de design � l'INT�RIEUR chaque cat�gorie 
# pour des contraintes sp�cifiques � chaque cat�gorie, et incluant possiblement des termes d'interaction.
# Retrait de l'argument x.loc. Les listes de locus impliqu�s dans chaque terme sont maintenant produites
# par la fonction design.constraint.

# Re-modifi�e par Jordie Croteau
# 25 mars 2011

# Re-modifi�e par Jordie Croteau (seulement les commentaires)
# 31 mars 2011

# Re-modifi�e par Jordie Croteau
# 1er avril 2011

# Ajout de la gestion des analyses quand certains sujets ont un ph�notype inconnu.
# Ajout aussi d'un traitement s�par� des familles dont tous les sujets sont dans 
# la m�me cat�gorie ph�notypique. Ces familles ne contribuent ni au score, ni � sa variance.
# Cependant, pour �viter les divisions par 0, une petite variance de 1e-10 est comptabilis�e
# (Il faudra v�rifier quand ces divisions par 0 se produisent)

# Re-modifi� par Alexandre Bureau
# 4 avril 2011

# Correction d'une erreur qui faisait que les calculs de variance pour toutes les cat�gories se faisaient  
# tous avec les termes de variance par locus pour la cat�gorie 1
# Ajout de la gestion des contraintes inter-cat�gories sp�cifi�es par par.constrained et constraints.

# Re-modifi� par Alexandre Bureau
# 29 avril 2011

# Correction d'une erreur dans dans le recodage des y1 et y2 en cat�gories 1 � 4

# Re-modifi� par Alexandre Bureau
# 9 juin 2011

# Correction d'une erreur dans le calcul de la variance des termes d'interaction: il faut calculer la somme des variances du produit 
# des scores de chaque locus par paires et non le produit des variances de la somme des scores de chaque locus.
# Implantation du calcul des covariances entre scores pour le m�me locus dans des fonctions de r�gression distincte

# Re-modifi� par Alexandre Bureau
# 20 janvier 2012
# Correction des longueurs des vecteurs ind.par et ind.cat 

# Re-modifi� par Jordie les 17 et 18 mai 2012 pour corriger quelques "degenerate cases".

# Re-modifi� par Jordie le 21 mai 2012 pour permettre d�avoir des fichiers de pedigree avec des ensembles de sujets diff�rents 
# (ou seulement dans un ordre diff�rent) d�un fichier de pedigree � l�autre.  Un avertissement est toutefois affich� si 
# certains ensembles de sujets diff�rent.  Dans ce cas, seulement les sujets en commun � tous les fichiers de pedigree
# sont conserv�s.  De plus, on exclut les sujets dont le g�notype d'un des SNPs � analyser est manquant.  On exclut 
# aussi les sujets dont le ph�notype ou l'endoph�notype est manquant (valeur 0).

# Re-modifi� par Alexandre Bureau
# 23 mai 2012
# ajout de conditions "if (n.loc>1)" pour �viter des op�rations qui ne s'appliquent pas avec 1 seul locus
# S'il y a des termes d'IBD qui sont NA parce qu'une cat�gorie n'est pas repr�sent�e, une correction a �t� apport�e 
# pour que la contribution � la covariance soit 0, pas NA

# Re-modifi� par Alexandre Bureau
# 25 mai 2012
# Changement � l'appel de cov.score.interfunction pour passer de la version 2 � la version 3 de cette fonction

# Re-modifi� par Alexandre Bureau
# 5 juin 2012

# Implantation d'une nouvelle fa�on d'estimer la variance des termes de produit

# Re-modifi� le 25 juin 2012 par Jordie:
# partie du d�but (lecture des donn�es) rel�gu�e � une fonction s�par�e: read.merlin.files,
# o� on lit les donn�es de format merlin.

# petites modifs par Jordie le 17 ao�t 2012.

# modifi� par Jordie le 22 ao�t pour enlever des arguments superflus.

scores.covs=function(subject.ids,fam.id,y,n.levels,ibd.dat,n.loc,xp,xp.loc,xl,il,xibd.loc,ind.par,rep.par,ind.catl,ind.cat,contingency.file,descrip.file,lc=NULL,alpha.vec=rep(0,n.levels-1))
{
###################### D�finition des arguments #####################################################################################
# subject.ids: vecteur des id de chaque sujet
# fam.id: vecteur des id de la famille de chaque sujet
# y: vecteur de type "factor" donnant le statut, qui est la combinaison des 6e colonne (endoph�notype) et 7e colonne (ph�notype) des fichiers ped
# n.levels : nombre de cat�gories du statut (y)
# ibd.dat : data frame des donn�es d'IBD combin�es
# n.loc : nombre de locus
# xp : matrice de design pour le calcul du score
# xp.loc : Liste du vecteur de locus impliqu�s dans chaque param�tre
# xl : matrice de design pour le calcul des covariances (contient seulement les effets principaux des locus)
# il : indices apr�s conversion des termes de produits en indices des variables dans le produit
# ind.par : donne les indices des locus pour la cat�gorie � laquelle chaque terme appartient
# lc : locus sur lequel on conditionne le test du score
# alpha.vec : vecteur de log rapports de cote entre ph�notype et compte d'all�le au locus pr�cis� par lc, un �l�ment pour chaque cat�gorie k vs. les autres cat�gories
# autres param�tres � commenter..............
#####################################################################################################################################


######################################## donn�es d'IBD ##########################################################
# Prendre le produit des IBD des paires de locus pour le calcul de la variance des termes d'interaction
# s'il y a plus d'un locus et il y a pr�sence de termes d'interaction.
# Attention! Pour l'instant, soit tous les termes d'interaction sont pr�sents, soit il y en a aucun
l1=ibd.dat$ID1
l2=ibd.dat$ID2
fam.id.ibd=ibd.dat$FAMILY
pim=as.matrix(ibd.dat[,4:ncol(ibd.dat)])

if (n.loc>1 & max(xibd.loc) > n.loc)
{
pim2 <- apply(pim,1,produits.paires)

# Combiner les termes originaux et les produits en une seule matrice
# Jordie, 17 mai 2012: ajout de transposes pour avoir la matrice dans le bon sens lorsqu'il y a plus de 2 locus.
pim <- t(rbind(t(pim),pim2))
}
#################################################################################################################

#################### calculs des scores et covariances pour toutes les familles #####################################################
fam.id.u=unique(fam.id)
nb.fam=length(fam.id.u)

ibd.terms.mat=sigma2i.mat=array(NA,c(nb.fam,ncol(pim),n.levels-1,n.levels-1))
sigma2.mat=array(NA,c(nb.fam,dim(xl)[2],n.levels-1,n.levels-1))
scores.mat=array(NA,c(nb.fam,dim(xp)[2]))

for(i in 1:nb.fam)
 {
  indices.y=fam.id==fam.id.u[i]
  if(contingency.file)
   {
    cat("family",fam.id.u[i],"\n",file=descrip.file,append=TRUE)
    cat(table(y[indices.y]),"\n",file=descrip.file,append=TRUE)
    cat("\n",file=descrip.file,append=TRUE)
   }
  # Si plus d'une cat�gorie de y est repr�sent�e dans la famille
  if(length(unique(y[indices.y]))>1)
   {
    # nombre de sujets dans la famille
    ni = sum(indices.y)
 
    # dans le cas o� l'analyse est faite avec un sous-ensemble des sujets de la famille,
    # limiter les donn�es d'IBD � ce sous-ensemble de sujets.
    sub.tmp=subject.ids[indices.y]
    indices.ibd=fam.id.ibd==fam.id.u[i] & l1%in%sub.tmp & l2%in%sub.tmp
  
	  sigma2.mat[i,,,]=cov.score.poly(array(xl[indices.y,,],c(ni,dim(xl)[2],dim(xl)[3])),y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))),x.loc=xibd.loc)
	  # Les covariances entre scores pour le m�me locus dans diff�rente fonctions de r�gression impliquent des estimations de la variance des scores
	  # On les copie de sigma2.mat dans les �l�ments appropri�s de sigma2i.mat
    # Appel de la version 4
    sigma2i.mat[i,,,]=cov.score.interfunction(xibd.loc,ind.catl,array(sigma2.mat[i,,,],dim(sigma2.mat)[2:4]))
    
    # Si on ne conditionne pas sur un locus  
    if (is.null(lc))
      {
      ibd.terms.mat[i,,,]=ibd.terms(y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))))
      scores.mat[i,]=score.poly(array(xp[indices.y,,],c(ni,dim(xp)[2],dim(xp)[3])),y[indices.y])
      }
    # Sinon on conditionne sur un locus en pond�rant
    else
      {
      if (!(lc %in% 1:dim(xl)[2])) stop("The index lc does not correspond to a valid locus index in",1:dim(xl)[2])
      if (length(alpha.vec)!=dim(xl)[3]) stop("The number of alpha coefficients for the computation of weights does not equal the 3rd dimension of xl (",dim(xl)[3],")")
      # Calcul des poids pour la famille
      # On suppose que la matrice de design pour le niveau 1 contient le locus lc
      w = calcule.poids(array(xl[indices.y,,],c(ni,dim(xl)[2],dim(xl)[3])),y[indices.y],ind.par,rep.par,alpha.vec,lc,klc=1)      
      ibd.terms.mat[i,,,]=ibd.terms.w(y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))),w)
      scores.mat[i,]=score.poly.w(array(xp[indices.y,,],c(ni,dim(xp)[2],dim(xp)[3])),y[indices.y],w)      
      }
    }
  # s'il y a une seule cat�gorie de y repr�sent�e dans la famille, les programmes ne s'appliquent pas 
  # et on donne une valeur de 0 au score et une valeur tr�s petite � ibd.terms.mat et sigma2.mat.
  else 
   {
	ibd.terms.mat[i,,,]=1e-10
    scores.mat[i,]=0
   }
 }

# la matrice (array) des covariances des effets principaux est la somme sur les deux dimensions de cat�gories
# des produits des termes d'ibd et des sigma2.
#petite modif Jordie, 24 mars (il y avait une erreur si une ou des dimensions sont de longueur 1):
#ibd.terms.xl.mat=array(ibd.terms.mat[,xl.loc,,],c(dim(ibd.terms.mat)[1],length(xl.loc),dim(ibd.terms.mat)[3],dim(ibd.terms.mat)[4]))
ibd.terms.xl.mat=array(ibd.terms.mat[,xibd.loc,,],c(dim(ibd.terms.mat)[1],length(xibd.loc),dim(ibd.terms.mat)[3],dim(ibd.terms.mat)[4]))
cov.mat.l=apply(ibd.terms.xl.mat*sigma2.mat,1:2,sum,na.rm=TRUE)

# Calcul de covariances entre les diff�rentes cat�gories pour les m�mes locus
#     Hypoth�se: les matrices de design sont les m�mes pour diff�rentes cat�gories (ce qu'on obtient 
#     s'il n'y a pas de contrainte)
# Il faudra v�rifier si on a besoin de forcer un array s'il y a un seul locus
cov.inter.mat <- ibd.terms.mat*sigma2i.mat
# S'il y a des termes d'IBD qui sont NA parce qu'une cat�gorie n'est pas repr�sent�e, la contribution 
# � la covariance doit �tre 0, pas NA
cov.inter.mat[is.na(cov.inter.mat)] = 0
# V�rification que tous les termes de produits sont pr�sents dans sigm2i.mat
if (dim(sigma2i.mat)[2] > n.loc & dim(sigma2i.mat)[2] != n.loc*(n.loc+1)/2 ) stop ("sigma2i.mat does not contain variances for all product terms.")
####################################################################################################################################################

# On obtient la matrice des covariances de tous les effets en prenant la combinaison des termes appropri�s
cov.mat <- array(0,c(nb.fam,dim(xp)[2],dim(xp)[2]))
# Boucle sur les param�tres
for (j in 1:length(xp.loc))
  {
  # Boucle sur les termes associ�s � chaque param�tre (ordinairement un seul)
  for (h in 1:length(xp.loc[[j]]))
    {
	# Extraction des locus impliqu�s dans le terme h
    il <- as.numeric(unlist(strsplit(xp.loc[[j]][h],split="")))
    # On va chercher la variance du produit des locus impliqu�s dans l'�l�ment h du terme xp.loc[[j]]
	# Les indices des locus doivent s'interpr�ter dans la portion de la matrice de covariance
	# qui s'applique aux param�tres du terme xp.loc[[j]]
	# Pour l'instant, on suppose qu'il y a juste 2 locus
    if (length(il) > 1) 
	{   
    # AB: Voici ma solution pour le cas d'un seul terme d'interaction qui est consid�r� dans le programme.
    # indice o� se trouve le terme d'IBD
    # calcul� selon la formule n.loc + somme_i=1^(il[1]-1) (n.loc - i) + il[2] - il[1]
    # (on additionne n.loc pour sauter les termes d'IBD pour les effets principaux des locus)
    ii = n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1]
    cov.tmp <- cov.mat.l[,ind.par[[j]][ii]]
	}
    else cov.tmp <- cov.mat.l[,ind.par[[j]][il]]
	# On additionne le terme courant � la valeur de la covariance
    cov.mat[,j,j] <- cov.mat[,j,j] + cov.tmp
    }
  # s'il y a plus d'un param�tre
  if (length(xp.loc) > 1 & j > 1)
	{
	for (jj in 1:(j-1))
	  {
	  # On trouve l'intersection entre les termes des 2 listes
	  inter <- intersect(xp.loc[[j]],xp.loc[[jj]])
	  
	  if (length(inter) > 0) 
		{
        for (h in 1:length(inter))
          {
		  # On laisse tomber compl�tement pour l'instant la covariance entre termes d'une m�me fonction de r�gression
		  		  # Extraction des locus impliqu�s dans le terme h
		  il <- as.numeric(unlist(strsplit(inter[h],split="")))

	      # Si les indices de param�tres r�f�rent � des cat�gories diff�rentes 
          if (ind.cat[j] != ind.cat[jj])
            {
	        if (length(il) > 1) 
	        {
          # AB: Voici ma solution pour le cas d'un seul terme d'interaction qui est consid�r� dans le programme.
          # indice o� se trouve le terme d'IBD
          # calcul� selon la formule n.loc + somme_i=1^(il[1]-1) (n.loc - i) + il[2] - il[1]
          # (on additionne n.loc pour sauter les termes d'IBD pour les effets principaux des locus)
          ii = n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1]
          
          cov.tmp <- cov.inter.mat[,ii,ind.cat[j],ind.cat[jj]]
          }
	        # Sinon il y a un seul locus dans l'intersection
	        else cov.tmp <- cov.inter.mat[,il,ind.cat[j],ind.cat[jj]]
	        # On additionne le terme courant � la valeur de la covariance
            cov.mat[,j,jj] <- cov.mat[,j,jj] + cov.tmp
		    }		  
	      }
		}
	  cov.mat[,jj,j] <- cov.mat[,j,jj]
	  }
	}
  }
list(scores.mat=scores.mat,cov.mat=cov.mat)
}