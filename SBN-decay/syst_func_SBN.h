/* RATES AND CHI-SQUARED FUNCTIONS FOR LSND */

/* Square of real number */
inline double square(double x)
{
    return x*x;
}

inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}

/*inline double likelihood(double true_rate, double fit_rate, double sigma)
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
}*/



/***************************************************************************
 * Calculate chi^2 using Pedro's function
 ***************************************************************************/

double chiSBN(int exp, int rule, int np, double *x, double *errors, void* user_data)
 {
  int n_bins = glbGetNumberOfBins(0);
  double *true_rates_N = glbGetRuleRatePtr(2, 0);
  double *true_rates_M = glbGetRuleRatePtr(1, 0);
  double *true_rates_F = glbGetRuleRatePtr(0, 0);
  double *signal_fit_rates_N = glbGetSignalFitRatePtr(2, 0);
  double *signal_fit_rates_M = glbGetSignalFitRatePtr(1, 0);
  double *signal_fit_rates_F = glbGetSignalFitRatePtr(0, 0);
  double signal_norm_N, signal_norm_M, signal_norm_F, bg_norm_F, bg_norm_M, bg_norm_N;
  int ew_low, ew_high;
  double fit_rate;
  double chi2 = 0.0;
  int i;

  double bg_N[] = {1432.89, 1747.82, 1839.79, 1767.84, 1590.96, 1496.06, 1662.78, 1422.64, 1521.9, 1074.33, 1919.6};
  double bg_M[] = {108.887, 100.725, 93.7305, 89.7675, 87.669, 80.4405, 91.71, 76.476, 77.3325, 55.9575, 116.58};
  double bg_F[] = {184.895, 180.205, 173.963, 158.396, 153.707, 137.624, 160.664, 127.466, 133.365, 103.075, 204.35};


  glbGetEnergyWindowBins(0, 0, &ew_low, &ew_high);

  /* Loop over all bins in energy window */
  signal_norm_F = 1.0 + x[0];
  signal_norm_M = 1.0 + x[0];
  signal_norm_N = 1.0 + x[0];
  bg_norm_F = 1.0 + x[1];
  bg_norm_M = 1.0 + x[1];
  bg_norm_N = 1.0 + x[1];
  
  if(exp==0){
  for (i=ew_low; i <= ew_high; i++)
  {
    /* Statistical part of chi^2 for far detector */
    fit_rate  = signal_norm_F * signal_fit_rates_F[i] + bg_norm_F * bg_F[i];
    chi2 += likelihood((true_rates_F[i] + bg_F[i]), fit_rate, (true_rates_F[i] + bg_F[i]));

/* Statistical part of chi^2 for near detector */
    fit_rate  = signal_norm_M * signal_fit_rates_M[i] + bg_norm_M * bg_M[i];
    chi2 += likelihood((true_rates_M[i] + bg_M[i]), fit_rate, (true_rates_M[i] + bg_M[i]) );

    /* Statistical part of chi^2 for near detector */
    fit_rate  = signal_norm_N * signal_fit_rates_N[i] + bg_norm_N * bg_N[i];
    chi2 += likelihood((true_rates_N[i] + bg_N[i]), fit_rate, (true_rates_N[i] + bg_N[i]));
  }

  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < np; i++)
    chi2 += square(x[i] / errors[i]);
   }

   else {chi2=0;}

  return chi2;
}
