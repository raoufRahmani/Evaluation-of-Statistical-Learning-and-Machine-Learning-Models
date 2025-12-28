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

   /* ==========================
      CAS H0 : indépendance
      ========================== */
   start DGP_H0(X);
      n = nrow(X);
      p = ncol(X);
      beta_H0 = j(p,1,0); /* vecteur de zéros */
      eps = normal(j(n,1,0)) * 0.1;
      Y_H0 = X * beta_H0 + eps;
      return(Y_H0);
   finish DGP_H0;

   /* ==========================
      CAS H1 : dépendance
      ========================== */
   start DGP_H1(X);
      n = nrow(X);
      /* eps spécifique */
      eps = normal(j(n,1,0)) * 0.1;
      beta_H1 = {0.8, -0.7, 0.3, 0, -3};
      Y_H1 = X * beta_H1 + eps;
      return(Y_H1);
   finish DGP_H1;

   /*Séries avec ruptures */

   start DGP_ruptures(X);
      n = nrow(X);
      p = ncol(X);
      t = floor(n / 3);
      beta1 = {-1, 0.2, 1, 0, 0};
      beta2 = {-0.8, 0, 0, -1, 0.5};
      eps = normal(j(n,1,0)) * 0.1;
      Y_rupt = j(n,1,0);
      do i = 1 to n;
         if i <= t then
            Y_rupt[i] = X[i,] * beta1;
         else
            Y_rupt[i] = X[i,] * beta2;
      end;
      Y_rupt = Y_rupt + eps;
      return(Y_rupt);
   finish DGP_ruptures;

   /* ==========================
      Colinéarité
      ========================== */
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
      n = nrow(X);
      eps2 = normal(j(n,1,0));
      U = uniform(j(n,1,0));
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
      n = nrow(X);
      p = ncol(X);
      beta = j(p,1,0);
      r = Y; /* résidu initial = Y (beta=0) */

      do t = 1 to Titer;
         c = X` * r; /* corrélations */
         maxabs = max(abs(c));
         if maxabs = 0 then leave;
         idx = loc(abs(c) = maxabs); /* peut retourner plusieurs indices */
         idx = idx[1]; /* prendre le premier si égalité */
         gamma = epsStep * sign(c[idx]);
         beta[idx] = beta[idx] + gamma;
         r = r - gamma * X[,idx];
      end;

      return(beta);
   finish ForwardStagewise;

   /* ==========================
      LARS, Lasso, LassoCV, Stepwise, ElasticNet (proc glmselect)
      ========================== */
   start LARS(dataName, stopCriterion);
      submit dataName stopCriterion;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=lars(stop=&stopCriterion);
            output out=LarsOut p=predicted;
         run;
      endsubmit;

      use LarsOut;
      read all var {Y predicted};
      close LarsOut;

      n = nrow(Y);
      mse = ((Y - predicted)` * (Y - predicted)) / n;
      print "--- Résultats LARS ---", "Critère d'arrêt:" stopCriterion, "MSE:" mse;
   finish LARS;

   start Lasso(dataName, lambdaValue, stopCriterion);
      submit dataName lambdaValue stopCriterion;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=lasso(stop=&stopCriterion  lambda=&lambdaValue) stats=all;
            output out=LassoOut p=predicted;
         run;
      endsubmit;

      use LassoOut;
      read all var {Y predicted};
      close LassoOut;

      n = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n;
      print "--- Résultats Lasso ---", "Lambda:" lambdaValue, "Stop:" stopCriterion, "MSE:" mse;
   finish Lasso;


   start LassoCV(dataName, cvType);
      submit dataName cvType;
         proc glmselect data=&dataName;
            partition fraction(validate=0.3);
            model Y = X1-X5 / selection=lasso(stop=cv) cvmethod=&cvType;
            output out=LassoCVOut p=predicted;
         run;
      endsubmit;

      use LassoCVOut;
      read all var {Y predicted};
      close LassoCVOut;

      n = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n;
      print "--- LASSO CV ---", "CV type:" cvType, "MSE:" mse;
   finish LassoCV;

   start Stepwise(dataName, selectionMethod, stopCriterion, selectCriterion);
      submit dataName selectionMethod stopCriterion selectCriterion;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=&selectionMethod(select=&selectCriterion stop=&stopCriterion) stats=all;
            output out=StepOut p=predicted;
         run;
      endsubmit;

      use StepOut;
      read all var {Y predicted};
      close StepOut;

      n = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n;
      print "--- Stepwise ---", "Méthode:" selectionMethod, "Critère sélection:" selectCriterion, "Stop:" stopCriterion, "MSE:" mse;
   finish Stepwise;

   start ElasticNet(dataName, alphaValue, lambdaValue, stopCriterion);
      submit dataName alphaValue lambdaValue stopCriterion;
         proc glmselect data=&dataName;
            model Y = X1-X5 / selection=elasticnet(L1=0.5 L2=0.5 stop=sl) stats=all;
            output out=ENetOut p=predicted;
         run;
      endsubmit;

      use ENetOut;
      read all var {Y predicted};
      close ENetOut;

      n = nrow(Y);
      residus = Y - predicted;
      mse = (residus`*residus)/n;
      print "--- Elastic Net ---", "Alpha:" alphaValue, "Lambda:" lambdaValue, "Stop:" stopCriterion, "MSE:" mse;
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

   /* ==========================
      Tester Forward Stagewise sur Y_H1, Y_H0, Y_rupt
      ========================== */
   epsStep = 0.001;
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
   beta_FS_rupt = ForwardStagewise(X, Y_, epsStep, Titer);
   Y_hat_FS_rupt = X * beta_FS_rupt;
   MSE_FS_rupt = (Y_rupt - Y_hat_FS_rupt)` * (Y_rupt - Y_hat_FS_rupt) / n;
   print "Forward Stagewise - Y_rupt - MSE:", MSE_FS_rupt;
   print "Beta (Forward Stagewise) - Y_rupt:", beta_FS_rupt[colname={"b1" "b2" "b3" "b4" "b5"}];



   /*Appel des autres méthodes*/
   call LARS("Data_H1", "SL");
   call Lasso("Data_H1", 0.1, "SL");
   call Stepwise("Data_H1", "forward", "SL", "SL");
   call Stepwise("Data_H1", "backward", "AIC", "AIC");
   call LassoCV("Data_H1", "random(5)");
   call ElasticNet("Data_H1", 0.5, 0.1, "SL");




