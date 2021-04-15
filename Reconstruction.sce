/*-------------------------------------------------
 Auteur : Manon Cassagne & Valentin Labat
 Vous trouverez ci-dessous les fonctions de reconstruction et les fonctions nécessaires à
 cette reconstruction

---------------------------------------------------*/


// Permet de pouvoir appeler les fonctions qui se trouvent dans les fichiers ci-dessous
exec("kSVD.sce")
exec("COSAMP.sce")
exec("irls.sce")
exec("Procede.sce")

// Fonction argmax mk
// Renvoie la position du max de rk
// Utilisée dans MP
function m_k = argmax_mk(D, R_k)
    [N, K] = size(D)
    X = zeros(K)

    for k=1:K
        X(k) = abs(D(:, k) * R_k) / norm(D(:, k), 2)
    end

    for k=1:K
        if max(X) == X(k) then
            m_k = k
        end
    end
endfunction

function mk=choix_mk(D, Rk)
    [M,N] = size(D)
    X = zeros(N,1)
    for j = 1:N
        X(j) = abs(D(:,j)'*Rk)/norm(D(:,j))
    end
    for j = 1:N
        if max(X) == X(j) then
            mk = j
        end
    end
endfunction



function [niter, residu_final, alpha]=MP(D, x0, kmax, eps)
    [M,N] = size(D)
    alpha=zeros(N,1)
    Rk = x0;
    x = 0
    k = 0
    while (k < kmax) & (norm(Rk) > eps)
        mk = choix_mk(D,Rk)
        x = x + (Rk'*D(:,mk))/norm(D(:,mk))^2*D(:,mk)
        alpha(mk)=alpha(mk)+(Rk'*D(:,mk))/norm(D(:,mk))^2;
        Rk = Rk - (Rk'*D(:,mk))/norm(D(:,mk))^2*D(:,mk)
        k = k+1
    end
    niter = k
    residu_final = norm(Rk);
endfunction


// Algorithme du MP
//function [k,normeRk,alpha]=MP(D, x, k_max, eps)
//    [N, K] = size(D)
//    alpha = zeros(K, 1)
//    R_k = x
//    k = 0
//
//    while (k < k_max) & (norm(R_k, 2) > eps)
//        //Cj = abs(D'*R_k)/norm()
//        m_k = choix_mk(D, R_k)
//        z_k = (R_k*D(:, m_k)) / (norm(D(:, m_k)) ** 2)
//        alpha(m_k) = alpha(m_k) + z_k
//        R_k = R_k - z_k * D(:, m_k)
//        k = k + 1
//    end
//    normeRk = norm(R_k)
//endfunction

// Calcule la contribution
// Appelée dans STOMP
function [Cj] = contribution(D, R_k)
    [N, K] = size(D)
    Cj = zeros(K, 1)
    for j=1:K
        Cj(j) = abs(D(:, j)*R_k) / norm(D(:, j), 2)
    end
endfunction

// Algorithme STOMP
function [k,normeRk,alpha] = STOMP(D, x, t, k_max, eps)
    [N, K] = size(D)
    alpha = zeros(K, 1)
    R_k = x
    P_k = []
    k = 0

    while (k < k_max) & (norm(R_k, 2) > eps)
        Cj = abs(D'*R_k)
        //Cj = contribution(D, R_k)  // On calcule la contribution de chaque atome
        S_k = t * norm(R_k, 2) / sqrt(K)  // Calcul de Sk

        for j=1:K
            if (Cj(j) > S_k) then
                P_k = [P_k,j]
            end
        end
        
        phi = D(:, P_k)

        z = pinv(phi)*x
        
        [lignes,colonnes] = size(P_k)
        for i=1:colonnes
            if (z(i) == 0) then
                disp(z(i))
            end
            alpha(P_k(i)) = z(i)
        end

        R_k = x - D*alpha
        k = k + 1
    end
    normeRk = norm(R_k)
endfunction

function ER = erreurRelative(reel,mesure)
    ER = norm(mesure - reel) / norm(reel)
endfunction


P = [15, 20, 25, 30, 50, 75]
[lignesP,colonnesP] = size(P)
// On choisit un des trois vecteurs à considérer, ici on prend x3

xConsidereLettre = 'x3'
phiConsidereLettre = 'phi4'
erreurRelativeMP = []
erreurRelativeOMP = []
erreurRelativeSTOMP = []
erreurRelativeCOSAMP = []
erreurRelativeIRLS = []

for p = 1:colonnesP
    chainey = 'Resultats/y/' + xConsidereLettre + '_' + phiConsidereLettre + '_y' + string(p) + '.csv'
    signal = read_csv(chainey,";")
    signal = strtod(signal,".")
    
    chainephi = 'Resultats/phi/' + xConsidereLettre + '_y' + string(p) + '_' + phiConsidereLettre + '.csv'
    phi = read_csv(chainephi,";")
    phi = strtod(phi,".")
    
    D = read_csv('Resultats/Dico.csv',";")
    D = strtod(D,".")
    //disp(size(phi))
    //disp(size(signal))
    //disp(size(signal))
//    signalmodif = signal'
    //disp(size(signalmodif))
//    signal = signalmodif * phi
//    signal = signal'
    

    Dico = phi * D
    eps = 1e-4
    parcimonie = 108
    kmax = 150
    [alpha_COSAMP, iterations_COSAMP, residuel_COSAMP] = COSAMP(signal, Dico, parcimonie, eps, kmax)
    [iterations_STOMP,residuel_STOMP,alpha_STOMP] = STOMP(Dico, signal, 2, kmax, eps)
    [iterations_MP,residuel_MP,alpha_MP]=MP(Dico, signal, kmax, eps)
    [iterations_OMP, residuel_OMP, alpha_OMP] = OMP(signal, Dico, eps, kmax)
    alpha_IRLS = IRLS(signal, Dico, 0.5, kmax, eps)
    //disp(size(signal))
    
    erreurRelativeMP = [erreurRelativeMP, erreurRelative(signal, Dico*alpha_MP)]
    erreurRelativeOMP = [erreurRelativeOMP, erreurRelative(signal, Dico*alpha_OMP)]
    erreurRelativeSTOMP = [erreurRelativeSTOMP, erreurRelative(signal, Dico*alpha_STOMP)]
    erreurRelativeCOSAMP = [erreurRelativeCOSAMP, erreurRelative(signal, Dico*alpha_COSAMP)]
    erreurRelativeIRLS = [erreurRelativeIRLS, erreurRelative(signal, Dico*alpha_IRLS)]
    
end
