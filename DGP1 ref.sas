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

   /*CAS H1 : dépendance */
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




   /* --- Génération des vecteurs Y pour chaque cas --- */
   Y_H0    = DGP_H0(X);
   Y_H1    = DGP_H1(X);
   Y_rupt  = DGP_ruptures(X);
   Mat_col = DGP_colinearite(n, p); /* Renvoie X et Y déjà concaténés */
   Y_out   = DGP_outliers(X);

   /* 1. Création de DATA_H0 */
   Z_H0 = X || Y_H0;
   create DATA_H0 from Z_H0[colname={"X1" "X2" "X3" "X4" "X5" "Y"}];
   append from Z_H0;
   close DATA_H0;

   /* 2. Création de DATA_H1 */
   Z_H1 = X || Y_H1;
   create DATA_H1 from Z_H1[colname={"X1" "X2" "X3" "X4" "X5" "Y"}];
   append from Z_H1;
   close DATA_H1;

   /* 3. Création de DATA_rupt */
   Z_rupt = X || Y_rupt;
   create DATA_rupt from Z_rupt[colname={"X1" "X2" "X3" "X4" "X5" "Y"}];
   append from Z_rupt;
   close DATA_rupt;

   /* 4. Création de DATA_colin */
   /* Note: Ici on utilise directement Mat_col car votre fonction fait déjà la fusion */
   create DATA_colin from Mat_col[colname={"X1" "X2" "X3" "X4" "X5" "Y"}];
   append from Mat_col;
   close DATA_colin;

   /* 5. Création de DATA_out */
   Z_out = X || Y_out;
   create DATA_out from Z_out[colname={"X1" "X2" "X3" "X4" "X5" "Y"}];
   append from Z_out;
   close DATA_out;
quit;
   



proc print data=work.DATA_H1(obs=5);
run;


proc contents data=work.DATA_H1;
run;


proc means data=work.DATA_H1 mean std;
   var X1 X2 X3 X4 X5;
run;
