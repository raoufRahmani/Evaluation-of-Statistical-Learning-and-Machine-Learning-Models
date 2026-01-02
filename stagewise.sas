options source2;

proc iml;
   /* Paramètres globaux */
   n = 500;
   p = 5;

   /* ====================================
      FONCTIONS DE GÉNÉRATION DE DONNÉES
      ==================================== */

   /* CAS H0 : indépendance */
   start DGP_H0(X);
      n_h0 = nrow(X);
      p_h0 = ncol(X);
      beta_H0 = j(p_h0, 1, 0);
      eps = normal(j(n_h0, 1, 0)) * 0.1;
      Y_H0 = X * beta_H0 + eps;
      return(Y_H0);
   finish DGP_H0;

   /* CAS H1 : dépendance */
   start DGP_H1(X);
      n_h1 = nrow(X);
      eps = normal(j(n_h1, 1, 0)) * 0.1;
      beta_H1 = {0.8, -0.7, 0.3, 0, -3};
      Y_H1 = X * beta_H1 + eps;
      return(Y_H1);
   finish DGP_H1;

   /* Séries avec ruptures */
   start DGP_ruptures(X);
      n_r = nrow(X);
      p_r = ncol(X);
      t = floor(n_r / 3);
      beta1 = {-1, 0.2, 1, 0, 0};
      beta2 = {-0.8, 0, 0, -1, 0.5};
      eps = normal(j(n_r, 1, 0)) * 0.1;
      Y_rupt = j(n_r, 1, 0);
      do i = 1 to n_r;
         if i <= t then
            Y_rupt[i] = X[i,] * beta1;
         else
            Y_rupt[i] = X[i,] * beta2;
      end;
      Y_rupt = Y_rupt + eps;
      return(Y_rupt);
   finish DGP_ruptures;

   /* Colinéarité */
   start DGP_colinearite(n, p);
      X_cor = j(n, p, 0);
      X_cor[,1] = normal(j(n, 1, 0));
      X_cor[,2] = X_cor[,1] + normal(j(n, 1, 0)) * 0.01;
      X_cor[,3] = normal(j(n, 1, 0));
      X_cor[,4] = X_cor[,3] + normal(j(n, 1, 0)) * 0.02;
      X_cor[,5] = normal(j(n, 1, 0));
      /* Centrer / réduire par colonne */
      do j = 1 to ncol(X_cor);
         X_cor[,j] = (X_cor[,j] - mean(X_cor[,j])) / std(X_cor[,j]);
      end;
      beta_cor = {1, -1, 0.5, 0, 0};
      Y_cor = X_cor * beta_cor + normal(j(n, 1, 0)) * 0.1;
      return(X_cor || Y_cor);
   finish DGP_colinearite;

   /* Outliers simples */
   start DGP_outliers(X);
      n_out = nrow(X);
      eps2 = normal(j(n_out, 1, 0));
      U = uniform(j(n_out, 1, 0));
      outlier_index = loc(U > 0.9);
      if ncol(outlier_index) > 0 then do;
         eps2[outlier_index] = normal(j(ncol(outlier_index), 1, 0)) * 5;
      end;
      beta_out = {0.5, 0.5, 0, 0, 0};
      Y_out = X * beta_out + eps2;
      return(Y_out);
   finish DGP_outliers;

   /* ====================================
      ALGORITHME STAGEWISE
      ==================================== */
   start STAGEWISE_MC(X, Y, beta_true, tol, maxSteps);
      n = nrow(X);
      p = ncol(X);
      
      beta_hat = j(p, 1, 0);
      active = j(p, 1, 0);
      residual = Y;
      
      mse_old = ssq(residual) / n;
      
      /* Boucle Stagewise Forward */
      do step = 1 to maxSteps;
         best_j = 0;
         best_bj = 0;
         best_mse = mse_old;
         
         do j = 1 to p;
            if active[j] = 0 then do;
               denom = X[,j]` * X[,j];
               if denom > 1e-10 then do;
                  bj = (X[,j]` * residual) / denom;
                  res_j = residual - X[,j] * bj;
                  mse_j = ssq(res_j) / n;
                  
                  if mse_j < best_mse then do;
                     best_mse = mse_j;
                     best_j = j;
                     best_bj = bj;
                  end;
               end;
            end;
         end;
         
         /* Critère d'arrêt */
         if best_j = 0 then leave;
         if (mse_old - best_mse) < tol then leave;
         
         /* Mise à jour */
         beta_hat[best_j] = beta_hat[best_j] + best_bj;
         active[best_j] = 1;
         residual = residual - X[,best_j] * best_bj;
         mse_old = best_mse;
      end;
      
      /* ====================================
         CALCUL DES MÉTRIQUES
         ==================================== */
      
      /* Biais au carré */
      bias = beta_hat - beta_true;
      bias2 = bias` * bias;
      
      /* Variance */
      var_hat = var(beta_hat);
      
      /* Underfitting */
      idx_true = loc(beta_true ^= 0);
      underfit = 0;
      if ncol(idx_true) > 0 then do;
         underfit_check = loc(active[idx_true] = 0);
         if ncol(underfit_check) > 0 then underfit = 1;
      end;
      
      /* Overfitting */
      idx_null = loc(beta_true = 0);
      overfit = 0;
      if ncol(idx_null) > 0 then do;
         overfit_check = loc(active[idx_null] = 1);
         if ncol(overfit_check) > 0 then overfit = 1;
      end;
      
      /* MSE de prédiction */
      Y_hat = X * beta_hat;
      mse = ssq(Y - Y_hat) / n;
      
      metrics = bias2 || var_hat || mse || underfit || overfit;
      return(metrics);
   finish STAGEWISE_MC;

   /* ====================================
      SIMULATION MONTE CARLO
      ==================================== */
   
   nScen = 5;
   MC = 500;
   tol = 1e-6;
   maxSteps = p;
   
   Res = j(MC * nScen, 5, .);
   Scenario = j(MC * nScen, 1, .);
   
   /* Définition des vrais bêta pour chaque scénario */
   beta_H0 = j(p, 1, 0);
   beta_H1 = {0.8, -0.7, 0.3, 0, -3};
   beta_out = {0.5, 0.5, 0, 0, 0};
   beta_col = {1, -1, 0.5, 0, 0};
   
   /* Beta moyen pour ruptures */
   beta_r = {-0.9, 0.1, 0.5, -0.5, 0.25};
   
   row = 1;
   
   do mc = 1 to MC;
      /* Génération de X commun */
      X = normal(j(n, p, 0));
      do j = 1 to p;
         X[,j] = (X[,j] - mean(X[,j])) / std(X[,j]);
      end;
      
      /* Scénario 1 : H0 */
      Y = DGP_H0(X);
      Res[row,] = STAGEWISE_MC(X, Y, beta_H0, tol, maxSteps);
      Scenario[row] = 1;
      row = row + 1;
      
      /* Scénario 2 : H1 */
      Y = DGP_H1(X);
      Res[row,] = STAGEWISE_MC(X, Y, beta_H1, tol, maxSteps);
      Scenario[row] = 2;
      row = row + 1;
      
      /* Scénario 3 : H1 + outliers */
      Y = DGP_outliers(X);
      Res[row,] = STAGEWISE_MC(X, Y, beta_out, tol, maxSteps);
      Scenario[row] = 3;
      row = row + 1;
      
      /* Scénario 4 : Ruptures */
      Y = DGP_ruptures(X);
      Res[row,] = STAGEWISE_MC(X, Y, beta_r, tol, maxSteps);
      Scenario[row] = 4;
      row = row + 1;
      
      /* Scénario 5 : Colinéarité */
      Mat = DGP_colinearite(n, p);
      Xc = Mat[, 1:p];
      Yc = Mat[, p+1];
      Res[row,] = STAGEWISE_MC(Xc, Yc, beta_col, tol, maxSteps);
      Scenario[row] = 5;
      row = row + 1;
   end;
   
   /* ====================================
      CRÉATION DU RÉSUMÉ
      ==================================== */
   
   Summary = j(nScen, 5, .);
   
   do s = 1 to nScen;
      idx = loc(Scenario = s);
      if ncol(idx) > 0 then
         Summary[s,] = mean(Res[idx,]);
      else
         Summary[s,] = j(1, 5, .);
   end;
   
   /* Export du résumé */
   create STAGEWISE_SUMMARY from Summary[
      rowname={"H0" "H1" "H1_Outliers" "Ruptures" "Colinearite"}
      colname={"Bias2" "Variance" "MSE" "UnderfitRate" "OverfitRate"}
   ];
   append from Summary;
   close STAGEWISE_SUMMARY;
   
   /* Affichage des résultats */
   print "=== RÉSULTATS DE LA SIMULATION MONTE CARLO ===";
   print Summary[
      rowname={"H0" "H1" "H1_Outliers" "Ruptures" "Colinearite"}
      colname={"Bias2" "Variance" "MSE" "UnderfitRate" "OverfitRate"}
      format=12.6
   ];

quit;
