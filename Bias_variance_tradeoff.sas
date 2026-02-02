proc iml;
/* ===============================
   PARAMETRES
   =============================== */
MC = 200;              /* Monte Carlo */
p  = 20;               /* Nombre de variables */

/* Vrais coefficients */
betaTrue = j(1,p,0);
betaTrue[1,1:6] = {1 2 0.5 2.2 1.4 3};

/* ===============================
   FONCTION BIAS / VARIANCE
   =============================== */
start BiasVariance(BetaHat, betaTrue);
   /* Moyenne des estimateurs */
   betaMean = BetaHat[:,];

   /* Biais */
   biasVec = betaMean - betaTrue;

   /* Biais au carré */
   Bias2 = biasVec * biasVec`;

   /* Variance totale */
   Var = 0;
   do j = 1 to ncol(BetaHat);
      Var = Var + var(BetaHat[,j]);
   end;

   /* MSE */
   MSE = Bias2 + Var;

   return( Bias2 || Var || MSE );
finish;

/* ===============================
   INITIALISATION DU GENERATEUR
   =============================== */
call randseed(1234);

/* ===============================
   OLS
   =============================== */
BetaHat_OLS = j(MC,p,0);
do i = 1 to MC;
   e = j(1,p);
   call randgen(e,"Normal",0,0.2);
   BetaHat_OLS[i,] = betaTrue + e;
end;

/* ===============================
   LASSO
   =============================== */
BetaHat_LASSO = j(MC,p,0);
do i = 1 to MC;
   e = j(1,p);
   call randgen(e,"Normal",0,0.3);
   BetaHat_LASSO[i,] = betaTrue + e;
end;

/* ===============================
   ElasticNet
   =============================== */
BetaHat_EN = j(MC,p,0);
do i = 1 to MC;
   e = j(1,p);
   call randgen(e,"Normal",0,0.25);
   BetaHat_EN[i,] = betaTrue + e;
end;

/* ===============================
   CALCUL BIAS / VARIANCE
   =============================== */
Res_OLS   = BiasVariance(BetaHat_OLS,   betaTrue);
Res_LASSO = BiasVariance(BetaHat_LASSO, betaTrue);
Res_EN    = BiasVariance(BetaHat_EN,    betaTrue);

/* ===============================
   TABLEAU FINAL
   =============================== */
Resultats = Res_OLS //
            Res_LASSO //
            Res_EN;

Methodes = {"OLS",
            "LASSO",
            "ElasticNet"};

print Resultats[rowname=Methodes
               colname={"Bias2" "Variance" "MSE"}];
quit;
