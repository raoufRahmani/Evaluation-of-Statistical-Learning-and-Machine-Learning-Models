options source2;

proc iml;
   /* Paramètres */
   n = 500;
   p = 5;

   /* Génération de X*/
   X = normal(j(n, p, 0));
   /* essayer centrer et réduire col par col */
   do j = 1 to ncol(X);
      X[,j] = (X[,j] - mean(X[,j])) / std(X[,j]);
   end;

   /* Bruit de base (eps global utilisable) */
   eps_global = normal(j(n,1,0)) * 0.1;

   /*CAS H0 : indépendance*/
   start DGP_H0(X);
      n_h0 = nrow(X);
      p_h0 = ncol(X);
      beta_H0 = j(p_h0,1,0); /* vecteur de zéros */
      eps = normal(j(n_h0,1,0)) * 0.1;
      Y_H0 = X * beta_H0 + eps;
      return(Y_H0);
   finish DGP_H0;

   /* ==========================
      CAS H1 : dépendance
      ========================== */
   start DGP_H1(X);
      n_h1 = nrow(X);
      /* eps spécifique */
      eps = normal(j(n_h1,1,0)) * 0.1;
      beta_H1 = {0.8, -0.7, 0.3, 0, -3};
      Y_H1 = X * beta_H1 + eps;
      return(Y_H1);
   finish DGP_H1;

   /*Séries avec ruptures */

   start DGP_ruptures(X);
      n_r = nrow(X);
      p_r = ncol(X);
      t = floor(n_r / 3);
      beta1 = {-1, 0.2, 1, 0, 0};
      beta2 = {-0.8, 0, 0, -1, 0.5};
      eps = normal(j(n_r,1,0)) * 0.1;
      Y_rupt = j(n_r,1,0);
      do i = 1 to n_r;
          if i <= t then
            Y_rupt[i] = X[i,] * beta1;
          else
            Y_rupt[i] = X[i,] * beta2;
      end;
      Y_rupt = Y_rupt + eps;
      return(Y_rupt);
   finish DGP_ruptures;

   /*Colinéarité*/
   start DGP_colinearite(n, p);
      X_cor = j(n,p,0);
      X_cor[,1] = normal(j(n,1,0));
      X_cor[,2] = X_cor[,1] + normal(j(n,1,0)) * 0.01;
      X_cor[,3] = normal(j(n,1,0));
      X_cor[,4] = X_cor[,3] + normal(j(n,1,0)) * 0.02;
      X_cor[,5] = normal(j(n,1,0));
      /* centrer / réduire par colonne */
      do j = 1 to ncol(X_cor);
          X_cor[,j] = (X_cor[,j] - mean(X_cor[,j])) / std(X_cor[,j]);
      end;
      beta_cor = {1, -1, 0.5, 0, 0};
      Y_cor = X_cor * beta_cor + normal(j(n,1,0)) * 0.1;
      return( X_cor || Y_cor ); /* concaténé pour usage si nécessaire */
   finish DGP_colinearite;

   /*Outliers simples*/
   start DGP_outliers(X);
      n_out = nrow(X);
      eps2 = normal(j(n_out,1,0));
      U = uniform(j(n_out,1,0));
      outlier_index = loc(U > 0.9); /* ~10% */
      if ncol(outlier_index) > 0 then do;
          eps2[outlier_index] = normal(j(ncol(outlier_index),1,0)) * 5;
      end;
      beta_out = {0.5, 0.5, 0, 0, 0};
      Y_out = X * beta_out + eps2;
      return(Y_out);
   finish DGP_outliers;

   /*Stagewise qui marche*/
   start ForwardStagewise(X, Y, epsStep, Titer);
      n_fs = nrow(X);
      p_fs = ncol(X);
      beta = j(p_fs,1,0);
      r = Y; /* résidu initial = Y (beta=0) */

      do t = 1 to Titer;
          c = X` * r; /* corrélations */
          maxabs = max(abs(c));
          if maxabs = 0 then leave;
          idx_list = loc(abs(c) = maxabs); /* peut retourner plusieurs indices */
          idx = idx_list[1]; /* prendre le premier si égalité */
          gamma = epsStep * sign(c[idx]);
          beta[idx] = beta[idx] + gamma;
          r = r - gamma * X[,idx];
      end;

      return(beta);
   finish ForwardStagewise;


   /*LARS*/
   /*ajoute a petits pas jusqu'a avant une variable ne devienne plus importante que la premiere*/
   start LARS(dataName, critere_arret); /*on peut mettre dans critere d'arret soit AIC, BIC, CV (pour validation croisée)*/
      submit dataName critere_arret;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=lars(stop=&critere_arret);
            output out=LarsOut p=predicted;
         run;
      endsubmit;

      use LarsOut;
      read all var {Y predicted};
      close LarsOut;

      n_lars = nrow(Y);
      mse = ((Y - predicted)` * (Y - predicted)) / n_lars;
      print "--- Résultats LARS ---", "Critère d'arrêt:" critere_arret, "MSE:" mse;
   finish LARS;

/*LARS ne fait que ajouter des vars, alots que lasso peut eliminier des vars*/

/*Lasso*/ 
 /*Lasso avec Lambda fixe */

      /*Lasso par minimiser l'erruer(logvraisemblance + norme L1 de nombre de coeffs), le lasso cherche a equilibrer entre l'erreur et
   la taille du modele*/

   start Lasso(dataName, lambdaValue);
      submit dataName lambdaValue;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=lasso(steps=120 lambda=&lambdaValue) stats=all;
            output out=LassoOut p=predicted;
         run;
      endsubmit;

      use LassoOut;
      read all var {Y predicted};
      close LassoOut;

      n_lasso = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n_lasso;
      print "--- Résultats Lasso Fixe ---", "Lambda:" lambdaValue, "MSE:" mse;
   finish Lasso;

   /*Lasso avec validation croisée car on peut pas le mettre dans une fonction ou on donne la valeur de lambda */
   start LassoCV(dataName, cvType);
      submit dataName cvType;
         proc glmselect data=&dataName;
            partition fraction(validate=0.3); /* Sépare les données */
            model Y = X1-X5 / selection=lasso(stop=cv) cvmethod=&cvType;
            output out=LassoCVOut p=predicted copyvar=_ROLE_; /* _ROLE_ permet de savoir qui est test/train */
         run;
      endsubmit;

      use LassoCVOut;
      /* On ne calcule le MSE que sur les données de validation (le groupe de test) */
      read all var {Y predicted} where(_ROLE_='VALIDATE'); 
      close LassoCVOut;

      n_lcv = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n_lcv;
      print "--- LASSO CV (Validation uniquement) ---", "CV type:" cvType, "MSE:" mse;
   finish LassoCV;

   /*Stepwise : Selecton forward et backward*/

   start Stepwise(dataName, selectionMethod, stopCriterion);
      submit dataName selectionMethod stopCriterion;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=&selectionMethod(stop=&stopCriterion) stats=all;
            output out=StepOut p=predicted;
         run;
      endsubmit;

      use StepOut;
      read all var {Y predicted};
      close StepOut;

      n_step = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n_step;
      print "--- Stepwise ---", "Méthode:" selectionMethod, "Stop:" stopCriterion, "MSE:" mse;
   finish Stepwise;

   /*Elastic Net*/
/*mélange les deux techniques que nous avons vues : le Lasso et la Ridge.*/

   start ElasticNet(dataName, stopCriterion);
      submit dataName stopCriterion;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=elasticnet(stop=&stopCriterion) stats=all;
            output out=ENetOut p=predicted;
         run;
      endsubmit;

      use ENetOut;
      read all var {Y predicted};
      close ENetOut;

      n_enet = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n_enet;
      print "--- Elastic Net ---", "Stop:" stopCriterion, "MSE:" mse;
   finish ElasticNet;

   /* ==========================
      Générer les jeux de données (H0, H1, ruptures, outliers)
      ========================== */

   Y_H0 = DGP_H0(X);
   Y_H1 = DGP_H1(X);
   Y_rupt = DGP_ruptures(X);
   Y_out = DGP_outliers(X);

   /* Exporter Data_H1 (X + Y_H1) pour PROC GLMSELECT */
   Z = X || Y_H1;
   create Data_H1 from Z[colname={"X1" "X2" "X3" "X4" "X5" "Y"}];
   append from Z;
   close Data_H1;

   /* Afficher 5 premières obs via PROC PRINT (hors IML) */
   submit;
      proc print data=Data_H1(obs=5);
      run;
   endsubmit;

   /*Tester Forward Stagewise sur Y_H1, Y_H0, Y_rupt   epsStep = 0.001;
   Titer = 5000;

   /* Sur Y_H1 */
   beta_FS_H1 = ForwardStagewise(X, Y_H1, epsStep, Titer);
   Y_hat_FS_H1 = X * beta_FS_H1;
   MSE_FS_H1 = (Y_H1 - Y_hat_FS_H1)` * (Y_H1 - Y_hat_FS_H1) / n;
   print "Forward Stagewise - Y_H1 - MSE:", MSE_FS_H1;
   print "Beta (Forward Stagewise) - Y_H1:", beta_FS_H1[colname={"b1" "b2" "b3" "b4" "b5"}];

   /* Sur Y_H0 */
   beta_FS_H0 = ForwardStagewise(X, Y_H0, epsStep, Titer);
   Y_hat_FS_H0 = X * beta_FS_H0;
   MSE_FS_H0 = (Y_H0 - Y_hat_FS_H0)` * (Y_H0 - Y_hat_FS_H0) / n;
   print "Forward Stagewise - Y_H0 - MSE:", MSE_FS_H0;
   print "Beta (Forward Stagewise) - Y_H0:", beta_FS_H0[colname={"b1" "b2" "b3" "b4" "b5"}];

   /* Sur Y_rupt */
   beta_FS_rupt = ForwardStagewise(X, Y_rupt, epsStep, Titer);
   Y_hat_FS_rupt = X * beta_FS_rupt;
   MSE_FS_rupt = (Y_rupt - Y_hat_FS_rupt)` * (Y_rupt - Y_hat_FS_rupt) / n;
   print "Forward Stagewise - Y_rupt - MSE:", MSE_FS_rupt;
   print "Beta (Forward Stagewise) - Y_rupt:", beta_FS_rupt[colname={"b1" "b2" "b3" "b4" "b5"}];

   /*Appel des autres méthodes*/
   run LARS("Data_H1", "AIC");
   run Lasso("Data_H1", 0.1);
   run Stepwise("Data_H1", "forward", "SL");
   run Stepwise("Data_H1", "backward", "AIC");
   run LassoCV("Data_H1", "random(5)");
   run ElasticNet("Data_H1", "AIC");

quit;
