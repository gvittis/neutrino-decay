inline double square(double x)
{
    return x*x;
}
 
  int en2 = 200, flu2=7;
  double miniplus[200][7];

// This function loads MiniBooNE flux

  int loadplus(void)
  {
   FILE *file;
   int q, qi,i,j;
   int state;
   // Read first file "MINIBOONEplus.dat" for MiniBooNE flux
   q=0;  //counter of lines of the file
   file = fopen("./flux/MINIBOONEplus.dat", "r");
   while (q < en2) 
   {
    for(qi=0; qi<flu2;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &miniplus[q][qi]); // scan the column functions
   //  printf("%g \n",miniplus[q][qi]);
     }
    q++;
   }
// for(i=0; i <= 10; i++){ for(j=0; j <= 10; j++){printf("%g \n ",miniplus[j][i]);} printf("\n");} 
 printf("%g \n",miniplus[199][0]); 
   fclose(file);
  
   return 0;
  }

/*inline double square(double x)
{
    return x*x;
} */

/***************************************************************************
 *                     Return the interpolated flux                        *
 ***************************************************************************/

double minip_flux(double E, int f /* 1 a 6 */)
{
  int col = f;
  int n_steps = 200;
  double result;

  int k=0;
  while (k <= n_steps  &&  miniplus[k][0] < E)
    k++;
  if (k <= 0 || k > n_steps)
    return 0.0;
  else
  {
    double E_lo    = miniplus[k-1][0];
    double E_up    = miniplus[k][0];
    double flux_lo = miniplus[k-1][col];
    double flux_up = miniplus[k][col];
    result  = flux_lo + (E - E_lo)*(flux_up - flux_lo)/(E_up - E_lo);
  }
  return result;

}

/***************************************************************************
 *                   Probability functions                                 *
 ***************************************************************************/

double /*complex*/ decay_minip_1(double x, double y, int f, double Erec, double L){
     double /*complex*/ decay=0;
     double En;
     double binwidth = 0.05;
     double /*complex*/ dPdE1=0;
             
for (En = Erec + binwidth; En <  miniplus[199][0]; En= En+binwidth)
{
  dPdE1 = 0.5*2*Erec/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));
  
  decay = decay + dPdE1*minip_flux(En,f)*binwidth;
}
//printf("%g \n",decay); 
     return decay;
}

double /*complex*/ decay_minip_2(double x, double y, int f, double Erec, double L){
     double /*complex*/ decay=0;
     double En;
     double binwidth = 0.05;
     double /*complex*/ dPdE2=0;
     
for (En = Erec + binwidth; En < miniplus[199][0]; En= En+binwidth)
{
  dPdE2 = 0.5*2*(En-Erec)/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE2*minip_flux(En,f)*binwidth;
 }
     return decay;
}

