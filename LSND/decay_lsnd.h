
  int en = 501, flu=7;
  double lsnd[501][7];

// This function loads LSND flux

  int loadflux(void)
  {
   FILE *file;
   int q, qi;
   int state;
   // Read first file "flux.dat" for LSND flux
   q=0;  //counter of lines of the file
   file = fopen("./flux/flux2.dat", "r");
   while (q < en) 
   {
     // state=fscanf(file,"%lf %lf", &a[q], &b[q]); //Scan first two parameters of each column (variables)
    for(qi=0; qi<flu;qi++) // counter of columns of the file ntorres gives the number of columns
     { 
      state=fscanf(file,"%lf", &lsnd[q][qi]); // scan the column functions
     }
    q++;
   }
   fclose(file);
  
   return 0;
  }


/***************************************************************************
 *                     Return the interpolated flux                        *
 ***************************************************************************/

double lsnd_flux(double E, int f /* 1 a 6 */)
{
  int col = f;
  int n_steps = 500;
  double result;

  int k=0;
  while (k <= n_steps  &&  lsnd[k][0] < E)
    k++;
  if (k <= 0 || k > n_steps)
    return 0.0;
  else
  {
    double E_lo    = lsnd[k-1][0];
    double E_up    = lsnd[k][0];
    double flux_lo = lsnd[k-1][col];
    double flux_up = lsnd[k][col];
    result  = flux_lo + (E - E_lo)*(flux_up - flux_lo)/(E_up - E_lo);
  }
  return result;
}

/***************************************************************************
 *                   Probability functions                                 *
 ***************************************************************************/

double complex decay_lsnd_1(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 0.377;
     double complex dPdE1=0;
             
for (En = Erec + binwidth; En < lsnd[500][0]; En= En+binwidth)
{
  dPdE1 = 0.5*2*Erec/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE1*lsnd_flux(En,f)*binwidth;
 }
     return decay;
}

double complex decay_lsnd_2(double x, double y, int f, double Erec, double L){
     double complex decay=0;
     double En;
     double binwidth = 0.377;
     double complex dPdE2=0;
     
for (En = Erec + binwidth; En < lsnd[500][0]; En= En+binwidth)
{
  dPdE2 = 0.5*2*(En-Erec)/(En*En)*x*(1-exp(-5.0677*square(y)*L/(16*PI*En)));

  decay = decay + dPdE2*lsnd_flux(En,f)*binwidth;
 }
     return decay;
}

/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/
double th12;
double th13;
double th23;
double deltacp;
double sdm;
double ldm;

double u;
double mg;

/***************************************************************************
 * Store oscillation parameters in internal data structures.
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{

  th12    = glbGetOscParams(p, GLB_THETA_12);
  th13    = glbGetOscParams(p, GLB_THETA_13);
  th23    = glbGetOscParams(p, GLB_THETA_23);
  deltacp = glbGetOscParams(p, GLB_DELTA_CP);
  sdm     = glbGetOscParams(p, GLB_DM_21);   
  ldm     = glbGetOscParams(p, GLB_DM_31);

  u     = glbGetOscParams(p, GLB_U);
  mg    = glbGetOscParams(p, GLB_MG);

  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
int my_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, th12, GLB_THETA_12);
  glbSetOscParams(p, th13, GLB_THETA_13);
  glbSetOscParams(p, th23, GLB_THETA_23);
  glbSetOscParams(p, deltacp, GLB_DELTA_CP);
  glbSetOscParams(p, sdm, GLB_DM_21);
  glbSetOscParams(p, ldm, GLB_DM_31);

  glbSetOscParams(p, u, GLB_U);
  glbSetOscParams(p, mg, GLB_MG);

  return 0;
}

/***************************************************************************
 * Calculate oscillation probabilities.                                    *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:            The buffer where the probabilities are to be stored     *
 *   cp_sign:      +1 if probalities for neutrinos are requested, -1 for   *
 *                 anti-neutrinos.                                         *
 *   E:            The neutrino energy in GeV                              *
 *   psteps:       Number of constant density layers in the matter profile *
 *   length:       The lengths of these layers in km                       *
 *   density:      The individual densities of these layers in g/cm^3      *
 *   filter_sigma: Width of low-pass filter as given in the AEDL file      *
 ***************************************************************************/
int my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  int i, j;
  double L;
  
 /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;
  
  /* Calculate total baseline */
  L = 0.0;
  for (i=0; i < psteps; i++)
    L += length[i];

 /*anu_e*/    P[1][0] = decay_lsnd_1(u,mg,5,E,L) /*muon*/+ decay_lsnd_2(u,mg,2,E,L) /*pion*/;
     
     P[1][1] = (square(1-u)+square(u)*exp(-5.0677*square(mg)*L/(16*PI*E)))*lsnd_flux(E,2);
    

 return 0;
}

