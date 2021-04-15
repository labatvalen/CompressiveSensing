/*-------------------------------------------------
 Auteur : Manon Cassagne & Valentin Labat
 Vous trouverez ci-dessous les fonctions phi et les focntions associées

---------------------------------------------------*/


// Fonction de création de la matrice de mesure phi1
// Matrice aléatoire générée à partir d’un processus uniformément distribué
function phi = phi1(lignes, colonnes)
    phi = rand(lignes, colonnes)
endfunction


// Fonction de création de la matrice de mesure phi2
// Matrice aléatoire générée à partir d’un processus bernoullien {−1,1}. (Si rand()<p alors X=-1 sinon X=1)
function phi=phi2(lignes, colonnes)
    p = 0.5
    phi = zeros(lignes, colonnes)

    for i=1:lignes
        for j=1:colonnes
            if rand() < p
                phi(i, j) = -1
            else
                phi(i, j) = 1
            end
        end
    end
endfunction


// Fonction de création de la matrice de mesure phi3
// Matrice aléatoire générée à partir d’un processus bernoullien {0,1}
function phi=phi3(lignes, colonnes)
    p = 0.5
    phi = zeros(lignes, colonnes)

    for i=1:lignes
        for j=1:colonnes
            if rand() < p
                phi(i, j) = 0
            else
                phi(i, j) = 1
            end
        end
    end
endfunction


// Fonction de création de la matrice de mesure phi4
// Matrice aléatoire générée à partir d’un processus gaussien identique et indépendamment distribué
function phi = phi4(lignes, colonnes)
    phi = rand(lignes, colonnes,"normal")
endfunction


// Fonction de création de la matrice de mesure phi5
// Matrice creuse ou parcimonieuse générée de façon aléatoire.
function phi=phi5(lignes, colonnes)
    p = 0.5
    phi = zeros(lignes, colonnes)

    for i=1:lignes
        for j=1:colonnes
            if rand() < p
                phi(i,j) = 0
            else
                phi(i,j) = rand()
            end
        end
    end
endfunction


// Avec un x donné, on effectue la fonction Phi_numeroDuPhi et on print les y dans des fichiers
function calculDesY(x,M, N, numeroDuPhi, nomDuX)
    [lignesM, colonnesM] = size(M)
    for i=1:colonnesM
        texte = 'Resultats/y/' + nomDuX + '_phi' + string(numeroDuPhi) + '_' + "y" + string(i)
        textePhi = 'Resultats/phi/' + nomDuX + '_' + "y" + string(i) + '_phi' + string(numeroDuPhi)
        phi = []
        if (numeroDuPhi == 1) then
            phi = phi1(M(i), N)
        elseif (numeroDuPhi == 2) then
            phi = phi2(M(i), N)
        elseif (numeroDuPhi == 3) then
            phi = phi3(M(i), N)
        elseif (numeroDuPhi == 4) then
            phi = phi4(M(i), N)
        elseif (numeroDuPhi == 5) then
            phi = phi5(M(i), N)
        end
        y = phi*x
        texte = texte + ".csv"
        textePhi = textePhi + ".csv"
        csvWrite(y, texte, ";")
        csvWrite(phi, textePhi, ";")
    end
endfunction

// Fonction intermédiaire qui permet de calculer la cohérence mutuelle entre phi et D donnés
function CM = calculCoherenceMutuelle(phi,D)
    [M, N] = size(phi)
    matCM = zeros(M, N)
    
    for i=1:M
        for j=1:N
            matCM(i, j) = sqrt(N) * abs(phi(i,:)*D(:, j)) / (norm(phi(i, :), 2) * norm(D(:, j), 2))
        end
    end
    CM = max(matCM)
endfunction


// La matrice des cohérences mutuelles des différents phi
// En ligne 1, Vous avez la cohérence des phi1 en fonction des 6 pourcentages, etc.
function matriceCM = matriceDesCoherencesMutuelles(D,M,N)
    [lignes,colonnes] = size(M)
    for j=1:colonnes
            matriceCM(1,j) = calculCoherenceMutuelle(phi1(M(j), N),D)
            matriceCM(2,j) = calculCoherenceMutuelle(phi2(M(j), N),D)
            matriceCM(3,j) = calculCoherenceMutuelle(phi3(M(j), N),D)
            matriceCM(4,j) = calculCoherenceMutuelle(phi4(M(j), N),D)
            matriceCM(5,j) = calculCoherenceMutuelle(phi5(M(j), N),D)
    end

endfunction


// Fonction qui permet l'initialisation des variables P, X, lignesX, colonnesX, M,X1,X2,X3,D,nombreDePhi
function [P, X, lignesX, colonnesX, M,X1,X2,X3,D]=initialisationVariables()
    // Vecteur des pourcentages
    P = [15, 20, 25, 30, 50, 75]
    
    // On lit les valeurs dans les csv associés
    X = read_csv("Donnees/x.csv",";")
    // On enlève la première ligne, qui est la légende des colonnes
    X = X(2:99,:)
    X = strtod(X,",") // On remplace les valeurs de X qui sont considérés comme des chaines de caractère par des nombres, et on définit le séparateur "," qui est le séparateur du fichier csv
    [lignesX,colonnesX] = size(X)// lignesX est aussi appelé N
    
    [lignesP,colonnesP] = size(P)
    // M est la matrice des dimensions des lignes de phi
    M = zeros(lignesP,colonnesP)
    for i=1:colonnesP
        M(i) = ceil((P(i) * lignesX) / 100)
    end
    
    // On lit les valeurs dans les csv associés
    Xtemp = read_csv("Donnees/xVal.csv",";")
    
    // On enlève la première ligne, qui est la légende des colonnes, et on sépare en trois vecteurs
    X1 = Xtemp(2:99,1)
    X2 = Xtemp(2:99,2)
    X3 = Xtemp(2:99,3)
    
    // On remplace les valeurs des X qui sont considérés comme des chaines de caractère par des nombres, et on définit le séparateur du fichier csv comme étant la virgule
    X1 = strtod(X1,",")
    X2 = strtod(X2,",")
    X3 = strtod(X3,",")
    
    D = read_csv("Resultats/Dico.csv",";")
    D = strtod(D,".")
endfunction


// Permet l'execution des fonctions du fichier
// Decommentez la ligne qui fait appel à la fonction si vous voulez executer ce fichier seulement
// Recommentez la ensuite pour que, lors de l'appel de ce fichier par d'autres fichiers, le kSVD ne s'execute pas dès le départ
function executionDuFichier()
    [P, X, lignesX, colonnesX, M,X1,X2,X3,D]=initialisationVariables()
    // Calcule tous les y et les mets dans le fichier associé
    nombreDePhi = 5
    for i=1:nombreDePhi
        calculDesY(X1,M, lignesX, i, 'x1')
        calculDesY(X2,M, lignesX, i, 'x2')
        calculDesY(X3,M, lignesX, i, 'x3')
    end
    
    matriceCM = matriceDesCoherencesMutuelles(D,M, lignesX)
endfunction


//executionDuFichier()
