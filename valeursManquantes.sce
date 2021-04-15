/*-------------------------------------------------
 Auteur : Manon Cassagne & Valentin Labat
 Vous trouverez ci-dessous les fonctions de reconstruction du vecteur avec valeurs manquantes

---------------------------------------------------*/

////////////////// Partie sur zr/////////////////


// On importe les fichiers scilab nécessaires
exec("COSAMP.sce")
exec("Procede.sce")
exec("irls.sce")

// On lit les valeurs dans les csv associés
z = read_csv("Donnees/z.csv",";")
// On enlève la première ligne, qui est la légende des colonnes
z = z(2:86,:)
z = strtod(z,",") // On remplace les valeurs de X qui sont considérés comme des chaines de caractère par des nombres, et on définit le séparateur "," qui est le séparateur du fichier csv
D = read_csv('Resultats/Dico.csv',";")
D = strtod(D,".")
[lignesD,colonnesD] = size(D)

[lignes, colonnes] = size(z)

// On calcules les phis nécessaires
phi2 = phi2(lignes,lignesD)
phi4 = phi4(lignes,lignesD)

// Calcul des dictionnaires
Dico2 = phi2 * D
Dico4 = phi4 * D

signal = z
eps = 1e-4
parcimonie = 108
kmax = 150

// Calcul des alphas par cosamp et irls
[alpha_COSAMP2, iterations_COSAMP2, residuel_COSAMP2] = COSAMP(signal, Dico2, parcimonie, eps, kmax)
[alpha_COSAMP4, iterations_COSAMP4, residuel_COSAMP4] = COSAMP(signal, Dico4, parcimonie, eps, kmax)


// Reconstruction des vecteurs
zr2cosamp = Dico2 * alpha_COSAMP2
zr4cosamp = Dico4 * alpha_COSAMP4

// Ecriture dans csv
csvWrite(zr2cosamp, "Resultats/z/zr2cosamp.csv", ";")
csvWrite(zr4cosamp, "Resultats/z/zr4cosamp.csv", ";")
