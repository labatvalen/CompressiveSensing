/*-------------------------------------------------
 Auteur : Manon Cassagne & Valentin Labat
 Vous trouverez ci-dessous la fonction IRLS

---------------------------------------------------*/


/* 
 ENTRÉES :
  x : signal donné
  D : dictionnaire 
  p : un entier
  nb_it_max : nombre maximal d'itérations
  
 SORTIES : 
  alpha : représentation parcimonieuse de x
*/
function alpha = IRLS(x,D,p,nb_it_max,epsilon)
    alpha = pinv(D) * x // initilisation de alpha
    k = 1
    while (k < nb_it_max)
        // On décompose le calcul de X pour que cela soit plus clair
        absolu = abs(alpha.^2 + epsilon)
        exposant = p / 2 - 1
        W = absolu.^exposant
        Q = diag(1./W)
        // alphaTemp est le alpha0, le alpha "précédent"
        alphaTemp = alpha
        // (inv(D * Q * D') * x = ((D * Q * D')\x)
        //disp(Q * D' * (pinv(D * Q * D') * x) == Q * D' * ((D * Q * D')\x))
        alpha = Q * D' * (pinv(D * Q * D') * x)
        if (abs(norm(alphaTemp) - norm(alpha)) < sqrt(epsilon) / 100) then
            if epsilon > 0.00000001 then
                epsilon = epsilon / 10
            end
        end
        k = k + 1
    end
endfunction
