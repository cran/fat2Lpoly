# Fonction pour mod�le simple � 1 locus avec 4 cat�gories

# correspond au fichier design_1locus_v3.R dans le dossier "programmes"


# par Alexandre Bureau

# modifi� par Jordie le 22 ao�t dans le but d'ajouter n.levels (nombre de cat�gories de la variable r�ponse) comme attribut de la fonction

design.1locus <- function(x,par.constrained,constraints)
{

# Matrice avec tous les effets

x.e <- list(x,x,x)
# Liste des locus impliqu�s dans chaque effet
x.loc.e <- list("1","1","1")

# Matrice avec les effets principaux seulement
x.l <- x.e
# Liste des locus impliqu�s dans chaque effet
x.loc.l <- x.loc.e

# Liste des sorties
li <- list(x.e=x.e,x.loc.e=x.loc.e,x.l=x.l,x.loc.l=x.loc.l)
attributes(li)$n.levels <- 4
li
}