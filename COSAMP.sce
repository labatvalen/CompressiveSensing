/*-------------------------------------------------
 Auteur : Manon Cassagne & Valentin Labat
 Vous trouverez ci-dessous la fonction COSAMP et ses sous-fonctions d'exécution

---------------------------------------------------*/

/* 
 ENTRÉES :
  X : signal donné
  D : dictionnaire 
  s : un entier 
  eps : seuil
  kMax : nombre maximal d'itérations
  
 SORTIES : 
  alpha : représentation parcimonieuse de x 
  nb_it : nombre d'itérations
  norme_residuel : la norme du résiduel de sortie
*/
function [alpha, nb_it, norme_residuel] = COSAMP(X,D,s,eps,kMax)
    // initialisation des variables
    [M,K] = size(D) // M lignes, K colonnes
    residuel = X // initialisation du résiduel
    support = [] // Initialisation du support
    k = 0
    nb_atomes_selectionnes = 2*s
    
    while (k < kMax) & (norm(residuel) > eps)
        //disp(size(D))
        //disp(size(residuel))
        contribution = abs(D'*residuel)
        
        // Sélection
        nb_atomes_selectionnes = min([nb_atomes_selectionnes, size(contribution)(1)]) // Si nb_atomes_selectionnes est plus grand que le nombre de lignes de contribution (= plus grand que le nombre de colonnes de D), nb_atomes_selectionnes prend la valeur du nombre de colonnes de D
        [vecteurTrie,positions] = gsort(contribution)// On trie le vecteur dans l'ordre décroissant
        support1 = positions(1:nb_atomes_selectionnes,:)' // On selectionne les premiers (les plus grands)
        support = union(support,support1) // On fait une union entre les deux supports
        AS = D(:,support) // On sélectionne les atomes
        
        z = pinv(AS) * X // Méthode des moindres carrés
        zAbsolu = abs(z)

        alpha = zeros(K,1)
        for i=1:s // Rejet
            [valeur_max,index_max] = max(zAbsolu)
            alpha(support(index_max)) = z(index_max)
            zAbsolu(index_max) = -1 // Comme c'est des valeurs absolues, en le mettant à -1, le nombre devient le dernier (en terme de valeur)du vecteur
        end
        
        [ligneAlpha,colonneAlpha] = size(alpha)
        support = []
        for i=1:ligneAlpha
            if ~(alpha(i) == 0)
                support = [support,i] // On reinitialise le support (qui va servir de support précédent à la prochaine itération)
            end
        end
        
        residuel = X - D * alpha // On met a jour le résiduel
        k = k+1
    end
    nb_it = k // On met à jour le nombre d'itérations qu'a produit l'algorithme
    norme_residuel = norm(residuel)
endfunction



// Fonction qui permet l'initialisation des variables X, D,s,eps,kMax
function [x,D,s,eps,kMax]=initialisationVariables()
    s = 39
    eps = 1e-4
    kMax = 1000 // Nombre max d'itérations

    // On lit les valeurs dans les csv associés
    Xtemp = read_csv("Donnees/xVal.csv",";")
    D = read_csv("Resultats/Dico.csv",";")
    
    // On enlève la première ligne, qui est la légende des colonnes, et on sépare en trois vecteurs
    x1 = Xtemp(2:99,1)
    x2 = Xtemp(2:99,2)
    x3 = Xtemp(2:99,3)
    
    // On remplace les valeurs des X qui sont considérés comme des chaines de caractère par des nombres, et on définit le séparateur du fichier csv comme étant la virgule
    // Pareil pour D, mais le séparateur est un point
    x1 = strtod(x1,",")
    x2 = strtod(x2,",")
    x3 = strtod(x3,",")
    D = strtod(D,".")
    
    // On choisit un des trois vecteurs à considérer
    x = x1
endfunction


// Permet l'execution des fonctions du fichier
// Decommentez la ligne qui fait appel à la fonction si vous voulez executer ce fichier seulement
// Recommentez la ensuite pour que, lors de l'appel de ce fichier par d'autres fichiers, le COSAMP ne s'execute pas dès le départ
function executionDuFichier()
    // On initialise les variables
    [x,D,s,eps,kMax]=initialisationVariables()
    // On execute le COSAMP avec les variables initialisées
    [alpha, nb_it, norme_residuel] = COSAMP(x,D,s,eps,kMax)
    disp("alpha = ",alpha)
    disp("nombre ditérations : ", nb_it)
    disp("norme du résiduel :",norme_residuel)
endfunction



//executionDuFichier()
