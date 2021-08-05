/* RATES AND CHI-SQUARED FUNCTIONS FOR LSND */

inline double likelihood(double true_rate, double fit_rate)
{
  double res;
  res = fit_rate - true_rate;
  if (true_rate > 0)
  {
    if (fit_rate <= 0.0)
      res = 1e100;
    else
      res += true_rate * log(true_rate/fit_rate);
  }
    else
    res = fabs(res);

  return 2.0 * res;
}

/* Square of real number */
inline double square(double x)
{
    return x*x;
}

/*inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqrt(square(sqr_sigma));
  else
    return 0.0;
}*/

double chiLSND(int exp, int rule, int np, double *x, double *errors, void* user_data)
 {
  double chi2 = 0.0;
 // double *signal_fit_rates = glbGetSignalFitRatePtr(0, 0);
  double fit_rate, signal_norm, bg_norm;
  int ew_low, ew_high,i,k;
 
  double bg_rates_LSND[]={0.0639,0.463,0.9261,1.5648,1.7245,1.9959,2.5867,2.7783,2.5388,1.7724,0.463};
  double DATA_rates_LSND[]={1.2774,0.016,-0.5269,4.5506,3.4329,6.3869,7.792,8.2391,6.6743,9.6442,3.5128};
//  double sigma_LSND[]={1.3572,0.974,1.6606,2.7624,2.6027,3.1136,3.4169,3.401,3.449,3.6885,2.2034};  

   signal_norm = 1.0 + x[0];
   bg_norm     = 1.0 + x[1];
   
   glbGetEnergyWindowBins(0, 0, &ew_low, &ew_high);

/*****************************************************************
Implementação do erro na forma do espectro de energia
******************************************************************/
   int n_bins = glbGetNumberOfBins(0);
   double signal_fit_rates[n_bins];
   double emin, emax;
   glbGetEminEmax(0, &emin, &emax);
   glbShiftEnergyScale(x[2], glbGetSignalFitRatePtr(0, 0), signal_fit_rates, n_bins, emin, emax);


/******************************************************************/

    for (i=ew_low; i <= ew_high; i++)
     {
    
    fit_rate  = (signal_norm*signal_fit_rates[i] + bg_norm*bg_rates_LSND[i]);
    chi2 += likelihood(DATA_rates_LSND[i], fit_rate/*,fit_rate*//*,sigma_LSND[i]*/);
    
     }
 
   for(k=0; k < np; k++)
    {
    chi2 += square(x[k]/errors[k]);
    }

  return chi2;

}

