/*On nomme les DGP comme suit :
DGP : 1=H0(Ind�pendance), 2=H1 (REALTION lin�aire respectant les hypotheses), 3=Ruptures, 4=Colinéarité, 5=Outliers */
%let DGP_TYPE = 3; 

/* Méthode : on insert soit : FORWARD | BACKWARD | LARS | LASSO | ELASTICNET */
%let METHOD   = ELASTICNET;

/* Paramètres généraux */
%let n  = 100;
%let p  = 20; /*Nombre de variables total, mais le modele ne cintient que 11 variables*/
%let MC = 1000;


proc iml;

/* ------------------- Paramètres ------------------- */
DGP_TYPE = &DGP_TYPE;
n  = &n;
p  = &p;
MC = &MC;

chooses  = {"AIC" "SBC" "BIC" "CV" "CP" "ADJRSQ" "AICC"};
nChoose  = ncol(chooses);

/* vecteur de coefficients reels qu'on essa */
beta_11 = j(p, 1, 0);
beta_11[1]  =  0.8;
beta_11[2]  = -0.7;
beta_11[3]  =  0.3;
beta_11[4]  =  1.2;
beta_11[5]  = -3.0;
beta_11[6]  =  0.5;
beta_11[7]  = -1.0;
beta_11[8]  =  0.9;
beta_11[9]  = -0.4;
beta_11[10] =  0.6;
beta_11[11] = -0.8;
/* les 9 dernier =0 */

/*Preparer tous les DGP, pour pouvoir les appeler dans une macro*/

start DGP_H1(X, b);
   return( X*b + normal(j(nrow(X),1,0))*0.5 );
finish;

start DGP_ruptures(X, b);
   n_r = nrow(X);
   t   = floor(n_r/3);
   b_alt = b;
   b_alt[1:3] = {-1, 1, -0.5};
   Y = j(n_r,1,0);
   do i=1 to n_r;
      if i<=t then Y[i] = X[i,]*b_alt;
      else         Y[i] = X[i,]*b;
   end;
   return( Y + normal(j(n_r,1,0))*0.1 );
finish;

start DGP_outliers(X, b);
   n_o = nrow(X);
   eps = normal(j(n_o,1,0))*0.5;
   U   = uniform(j(n_o,1,0));
   idx = loc(U>0.9);
   if ncol(idx)>0 then eps[idx] = normal(j(ncol(idx),1,0))*10;
   return( X*b + eps );
finish;


/*On d�finit des vrais coefficients qui seront une r�f�rence pour evaluer nos modeles*/
beta_true = j(p,1,0);
if DGP_TYPE ^= 1 then beta_true = beta_11;

Res = j(MC*nChoose, 9, .);
row = 1;


/*La boucle de simulation MonteCarlo, qui contient une boucle conditionnelle sur le DGP choisi*/

do iter = 1 to MC;

   /* Données */
   X      = normal(j(n,p,0));
   X_test = normal(j(n,p,0));

   do k=1 to p;
      X[,k]      = (X[,k]-mean(X[,k]))/std(X[,k]);
      X_test[,k] = (X_test[,k]-mean(X_test[,k]))/std(X_test[,k]);
   end;

   /* Choix du DGP */
   if      DGP_TYPE = 1 then do;
      Y = normal(j(n,1,0))*0.5;
      Y_test = normal(j(n,1,0))*0.5;
   end;
   else if DGP_TYPE = 2 then do;
      Y = DGP_H1(X, beta_11);
      Y_test = DGP_H1(X_test, beta_11);
   end;
   else if DGP_TYPE = 3 then do;
      Y = DGP_ruptures(X, beta_11);
      Y_test = DGP_ruptures(X_test, beta_11);
   end;
   else if DGP_TYPE = 4 then do;
      X[,2] = X[,1] + normal(j(n,1,0))*0.01;
      Y     = X*beta_11 + normal(j(n,1,0))*0.1;
      Y_test= X_test*beta_11 + normal(j(n,1,0))*0.1;
   end;
   else if DGP_TYPE = 5 then do;
      Y = DGP_outliers(X, beta_11);
      Y_test = DGP_outliers(X_test, beta_11);
   end;

   TrainMat = X || Y;
   create TRAIN_DS from TrainMat[colname=("X1":"X20" || "Y")];
   append from TrainMat;
   close TRAIN_DS;

   call symputx("method", "&METHOD");

   /*Boucke sur les Algorithmes, et pour chaque algortihme on choisit dedans le critere d'arret, et le critere
   de selection sera choisi par la boucle dans la Macro*/

   do c = 1 to nChoose;

      current_c = chooses[c];
      call symputx("crit", current_c);

      submit;
          

          %if &METHOD = FORWARD %then %do;
              proc glmselect data=TRAIN_DS plots=Coefficients;
                 model Y = X1-X20 / selection=forward(choose=&crit stop=CV) stb;
                 ods output ParameterEstimates = BetaHat;
              run;
          %end;

          %if &METHOD = BACKWARD %then %do;
              proc glmselect data=TRAIN_DS plots=Coefficients;
                 model Y = X1-X20 / selection=backward(choose=&crit stop=CV) stb;
                 ods output ParameterEstimates = BetaHat;
              run;
          %end;

          %if &METHOD = LASSO %then %do;
              proc glmselect data=TRAIN_DS plots=Coefficients;
                 model Y = X1-X20 / selection=lasso(choose=&crit STOP=cv) stb;
                 ods output ParameterEstimates = BetaHat;
              run;
          %end;

          %if &METHOD = LARS %then %do;
              proc glmselect data=TRAIN_DS plots=Coefficients;
                 model Y = X1-X20 / selection=lars(choose=&crit) stb;
                 ods output ParameterEstimates = BetaHat;
              run;
          %end;

          %if &METHOD = ELASTICNET %then %do;
              proc glmselect data=TRAIN_DS plots=Coefficients;
                 model Y = X1-X20 / selection=elasticnet(choose=&crit stop=sbc) stb;
                 ods output ParameterEstimates = BetaHat;
              run;
          %end;

          
      endsubmit;

      /* Lecture des coefficients */
      use BetaHat;
          read all var {Effect}   into EffNames;
          read all var {Estimate} into EstVals;
      close BetaHat;

      beta_hat = j(p,1,0);
      do j = 1 to p;
          idx = loc(EffNames = "X"+strip(char(j)));
          if ncol(idx)>0 then beta_hat[j] = EstVals[idx];
      end;

      active = (beta_hat^=0);

      /*Stocker des indices binaires pour l'overfitting et l'underfitting a chaque it�ration*/
      idx_true = loc(beta_true^=0);
      idx_zero = loc(beta_true=0);

      if ncol(idx_true) > 0 then is_under = (ncol(loc(active[idx_true]=0)) > 0);
      else is_under = 0;

      if ncol(idx_zero) > 0 then is_over = (ncol(loc(active[idx_zero]=1)) > 0);
      else is_over = 0;

      u_f=0; o_f=0; m_f=0; p_f=0;
      if (is_under=0 & is_over=0) then p_f = 1;
      else if (is_under=1 & is_over=0) then u_f = 1;
      else if (is_under=0 & is_over=1) then o_f = 1;
      else if (is_under=1 & is_over=1) then m_f = 1;

      Res[row,] = ssq(beta_hat-beta_true) ||
                  ssq(Y - X*beta_hat)/n ||
                  ssq(Y_test - X_test*beta_hat)/n ||
                  u_f || o_f || m_f || p_f ||
                  c || iter;
      row = row + 1;

   end;

end;


/*Les tables des Outputs*/

create FinalData from Res[colname={
   "Bias2" "MSE_tr" "MSE_te"
   "Under" "Over" "Mixed" "Perfect"
   "ChooseIdx" "Iter"}];
append from Res;
close FinalData;

quit;



proc format;
   value cF
      1 = "AIC"
      2 = "SBC"
      3 = "BIC"
      4 = "CV"
      5 = "CP"
      6 = "ADJRSQ"
      7 = "AICC";
run;


/* La table des moyennes des indices*/
proc summary data=FinalData nway;
   class ChooseIdx;
   var Bias2 Under Over Mixed Perfect MSE_te;
   output out=Stats mean=;
run;

title "Tableau récapitulatif : moyennes par critère CHOOSE";
proc print data=Stats noobs;
   format ChooseIdx cF.;
run;


/* Table pour les graphiques */
proc transpose data=Stats 
               out=LongData(rename=(_NAME_=Metrique COL1=Valeur));
   by ChooseIdx;
   var Under Over Mixed Perfect;
run;


/* Visulaiser le tableau des taux et des probabilites*/
ods graphics on;
ods listing;
ods exclude none;

title "STEPWISE &METHOD – DGP &DGP_TYPE – Comparaison des critères";
proc sgplot data=LongData;
   styleattrs datacolors=(cx6688aa cxcc6666 cx66aaaa cx998866 cxaa88cc);
   vbar ChooseIdx / response=Valeur group=Metrique groupdisplay=cluster datalabel;
   format ChooseIdx cF.;
   yaxis grid label="Probabilité (Somme = 1)";
   xaxis label="Critère CHOOSE";
   keylegend / title="Métrique";
run;



/*MSE Test en barres (moyenne)*/
title "Performance de prédiction : MSE Test (moyenne par critère)";
proc sgplot data=FinalData;
   vbar ChooseIdx / response=MSE_te stat=mean datalabel;
   format ChooseIdx cF.;
   xaxis label="Critère CHOOSE";
   yaxis grid label="MSE Test (moyenne)";
run;

title;
