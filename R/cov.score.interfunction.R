# Fonction pour calculer la covariance entre scores pour le m�me locus mais 
# des fonctions logit associ�es � diff�rentes cat�gories � l'int�rieur d'une famille
# Ceci est une r��criture compl�te de la fonction pour utiliser les valeurs 
# d�j� calcul�es par cov.score.poly

# correspond au fichier covariance_score_inter_v3.R dans le dossier "programmes"


# Attention! Cette fonction n'est pas en utilisation dans la fonction globale!
# Il semble qu'on est rest� � covariance_score_inter_v2.R

# par Alexandre Bureau
# version 3
# mai 2012

# xl.loc : Liste des locus impliqu�s dans chaque effet
# ind.cat : donne la cat�gorie � laquelle appartient chaque param�tre
# sigma2 : matrice de variance-covariance intra fonction pour une famille

cov.score.interfunction <- function(xl.loc,ind.catl,sigma2)
{
n.loc <- max(xl.loc)
nl <- dim(sigma2)[3] + 1
if(nl>2)
{

# Les dimensions de sigmai sont nombre de locus * K-1 * K-1
sigmai <- array(NA,c(n.loc,nl-1,nl-1))

for (k in 2:(nl-1))
  {
  for (l in 1:(k-1))
    {
      # On trouve l'intersection des locus pr�sents dans les cat�gories k et l
      ll = intersect(xl.loc[ind.catl==k],xl.loc[ind.catl==l])
      # On prend la moyenne des variances pour les deux cat�gories
      # Note: si une variance �gale 0 parce qu'il n'y a personne dans la cat�gorie en question, on ne devrait pas
      # l'utiliser, mais on laisse faire parce que la covariance sera mise � 0 avec les termes d'IBD
      # [xl.loc[ind.catl==k] sert � aller chercher les estimations de variance pour les locus dans l'intersection 
	    sigmai[ll,l,k] <- sigmai[ll,k,l] <- (sigma2[ind.catl==k,k,k][xl.loc[ind.catl==k]%in%ll]+sigma2[ind.catl==l,l,l][xl.loc[ind.catl==l]%in%ll])/2
    }
  }
sigmai
}
else NA
}