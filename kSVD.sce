/*-------------------------------------------------
 Auteur : Manon Cassagne & Valentin Labat

 L'apprentissage d'un dictionnaire par la méthode K-SVD nécessite
 l'utilisation de la méthode de codage parcimonieux OMP
 Ce fichier s'articule en 2 parties : 
    Partie 1 - OMP
    Partie 2 - K-SVD
    Veuillez décommenter la ligne 245 pour exécuter le fichier seul
---------------------------------------------------*/

/*--------------------------------------------------
                    Partie 1 - OMP
---------------------------------------------------*/

/*------------------------------
    FONCTIONS PRÉ REQUISES
-------------------------------*/

/* 
 Choix du m_k : 
 D = dictionnaire
 Rk = résiduel
 mk = poisition de l'atome de D de plus grande contribution
 m_k = argmax |<dj,R_k-1>| / ||dj||
*/
function mk = choix_mk(D,Rk)
    [M,N] = size(D) // D : M lignes, N colonnes
    X = zeros(N,1)  // X : vecteur utilisé pour trouver le max
    for j = 1:N     // Parcours des colonnes de D
        X(j) = abs(D(:,j)'*Rk)/norm(D(:,j)) // Pour chaque colonne, on enregistre la contribution de l'atome j
    end
    for j = 1:N // Parcours des colonnes de D
        if max(X) == X(j) then  // Recherche du maximum 
            mk = j // Enregistrement de la position j du maximum de X
        end
    end
endfunction


//---------------------------------------------------------------------------------------


/*----------------------------------------
      ORTHOGONAL MATCHING PURSUIT
----------------------------------------*/

/* 
 ENTRÉES : 
 x : vecteur à représenter parcimonieusement tq x=D*alpha
 D : dictionnaire
 eps : seuil 
 kMax : nombre maximum d'itérations

 SORTIES :
 nbIter = nombre d'itérations
 alpha = représentation parcimonieuse de x
 residu_final = résidu final 
*/
function [nbIter, residu_final, alpha] = OMP(x,D,eps,kMax)
    [M,N] = size(D) // M lignes, N colonnes
    alpha= zeros(N,1) // initialisation de la représentation parcimonieuse de x à 0 
    Rk = x // Initialisation du résiduel
    mk = choix_mk(D,Rk)
    phi = D(:,mk) // matrice contenant les colonnes retenues
    P=[mk] // Position de l'atome choisi 
    zk = (Rk'*D(:,mk))/norm(D(:,mk))^2
    k = 1
    
    while (k < kMax) & (norm(Rk) > eps)
        mk = choix_mk(D,Rk) //position de la kième colonne retenue
        P=[P,mk]
        phi = [phi,D(:,mk)]
        zk=phi\x
        alpha(P) = zk //Ajout des composantes retenues dans le vetceur de représentation parcimonieuse de x
        Rk = x-phi*zk // mise à jour du résidu
        k = k+1 
    end
    nbIter = k
    residu_final = norm(Rk)
endfunction



//---------------------------------------------------------------------------------------

/*--------------------------------------------------
                    Partie 2 - K-SVD
---------------------------------------------------*/

/*------------------------------
    FONCTIONS PRÉ REQUISES
-------------------------------*/


// Calcul de la matrice A contenant les représentations parcimonieuse par la méthode OMP 
// Entrée : x = vecteur d'apprentissage 
// D = dictionnaire
// Sortie : matrice A contenant les représentations parcimonieuse de X
function A = OMP_A(x,D)
    [M,N] = size(x) // M lignes, N colonnes
    A = []          // Initialisation matrice A
    //Données pour OMP
    eps = 1e-6;
    kMax = 40;
    for i = 1:N // Parcours des colonnes
        // On applique l'algo OMP à chaque colonne de x
        [nbIter, residu_final, alpha] = OMP(x(:,i),D,eps,kMax) 
        // On met à jour la matrice A en ajoutant la représentation parcimonieuse alpha de chaque colonne de x
        A = [A, alpha]
    end
endfunction

// Calcul de Wi = support de la ième ligne de A
// Contient les positions des composantes non nulles de la ième ligne de A
function wi = calcul_wi(A,l,i)
    wi = []
    for j=1:l //parcours des colonnes de A
        if ~(A(i,j)==0) then // Si A(i,j) non nul
            wi = [wi,j] // mise à jour de wi avec la jième position
        end
    end
endfunction


// Calcul de Ei : l'erreur entre les vecteurs d'apprentissage 
// et la contribution des atomes en enlevant le ième atome
function Ei = calcul_Ei(X,D,A,i,K)
    somme=0 
    for j=1:K // K colonnes 
        if ~(j==i) then // j différent de i
            somme=somme + D(:,j)*A(j,:) //Somme du produit dj*Aj
        end
    end
    Ei = X - somme 
endfunction


// Calcul de Omega_i = 1 pour i=1,...,card(wi)
// Omega_i = 0 sinon 
function Omega_i = calcul_Omega_i(l,wk)
    [valeurs,card]=size(wk)
    Omega_i = zeros(l,card)
    for i = 1:card
        Omega_i(wk(i),i)=1 // Omega vaut 1 si i est entre 1 et card
    end
endfunction

// Fonction qui met à jour la matrice A 
// qui contient les représentations parcimonieuse de X
// A : une matrice 
// i : position de la ligne de A que l'on met à jour
// li : ligne par laquelle on remplace la ième ligne 
function A_maj = maj_A(A,i,li)
    [M,N] = size(A)
    A_maj = A 
    j = 1 
    for k=1:N 
        // Si la valeur de la composante en ième ligne, kième colonne est non nulle
        if ~(A(i,k)==0) then 
            A_maj(i,k) = li(j) // On met A à jour
        end
    end
endfunction



//---------------------------------------------------------------------------------------


/*----------------------------------------
                   KSVD
----------------------------------------*/
// ENTRÉE :
// X = vecteurs d'apprentissage
// K = nombre de colonnes de D 
// L nb de répétitions
function [Dico,A]=KSVD_apprentissage(X,K,L)
    // Variables
    [N,l] = size(X) // N lignes, l colonnes
    Dico = []         // Initialisation D0
    A = zeros(K,l) // Initialisation de A
    
    // Initialisation du dictionnaire
    for i=1:K // Les K premières colonnes de X que l'on retient pour D0
        // Mise à jour de Dico avec la ième colonne de X normalisée
        Dico = [Dico, X(:,i)/norm(X(:,i))] 
    end
    
    // On répète le processus L fois, c'est là que s'effectue l'apprentissage
    for j=1:L 
       
        A = OMP_A(X,Dico)  // Calcul de A par la méthode OMP
       
        for i=1:K // Parcours de colonnes de D
            wi = calcul_wi(A,l,i) // Calcul du support
            card_wi = length(wi) // Cardinal de wi 
            
            if  ~(card_wi==0) then // Si wi n'est pas vide
                Ei = calcul_Ei(X,Dico,A,i,K) // Calcul de l'erreur
                Omega_i = calcul_Omega_i(l,wi) // Calcul de Omega_i
                Ei_r = Ei * Omega_i // Calcul de l'erreur rectifiée 
                [U Delta V] = svd(Ei_r) // Décomposition en valeurs singulières de l'erreur rectifiée 
                Dico(:,i)=U(:,1)
                A = maj_A(A,i,Delta(1,1)*V(1,:)) // Mise à jour de A pour que X soit parcimonieux dans D 
            end
        end
    end
endfunction


function [minParcimonie,maxParcimonie]=parcimonie(D,A)
    // On parcourt les colonnes de A et on regarde quel est l'ordre de parcimonie des vecteurs
    // Pour se donner un intervalle, on chosit de renvoyer le max et le min de la parcimonie
    [lignes,colonnes] = size(A)
    minParcimonie = lignes // On initialise le min au nombre de lignes comme ça si il y a un zéro, le min va changer
    maxParcimonie = 0 // On initialise le max a 0 comme ça si il y a un nombre qui n'est pas 0, le max va changer
    
    for i=1:colonnes
        parci = 0
        for j=1:lignes
            if ~(A(j,i) == 0)
                parci = parci + 1
            end
        end
        minParcimonie = min([parci,minParcimonie])
        maxParcimonie = max([parci,maxParcimonie])
    end
    disp('Ordre de Parcimonie min :', minParcimonie)
    disp('Ordre de Parcimonie max :', maxParcimonie)
endfunction


// Fonction qui permet l'initialisation des variables nécessaires au fichier
function [X,K,L]=initialisationVariables()
    K = 100 // Nombre de colonnes du dictionnaire
    L = 10 // Nombre de répétitions du KSVD

    // On lit les valeurs dans les csv associés
    X = read_csv("Donnees/x.csv",";")
    // On enlève la première ligne, qui est la légende des colonnes
    X = X(2:99,:)
    X = strtod(X,",") // On remplace les valeurs de X qui sont considérés comme des chaines de caractère par des nombres, et on définit le séparateur "," qui est le séparateur du fichier csv
endfunction


// Fonction qui permet d'écrire les résultats Dico et A dans leurs Fichiers
function []=ecritureDansFichier(Dico,A)
    csvWrite(Dico, "Resultats/Dico.csv", ";")
    csvWrite(A, "Resultats/A.csv", ";")
endfunction

// Permet l'execution des fonctions du fichier
// Decommentez la ligne qui fait appel à la fonction si vous voulez executer ce fichier seulement
// Recommentez la ensuite pour que, lors de l'appel de ce fichier par d'autres fichiers, le kSVD ne s'execute pas dès le départ
function executionDuFichier()
    // On initialise les variables
    [X,K,L]=initialisationVariables()
    // On execute le KSVD avec les variables initialisées
    [Dico,A] = KSVD_apprentissage(X,K,L)
    // Ordre de parcimonie
    [minParcimonie,maxParcimonie]=parcimonie(Dico,A)
    // On écrit les résultats obtenus dans les fichiers associés
    ecritureDansFichier(Dico,A)
endfunction



//executionDuFichier()



