/* RATES AND CHI-SQUARED FUNCTIONS FOR MINIBOONE */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

struct matrizes_auziliares
{ 

	// essa vai ser a matrix M
	gsl_matrix *Matrix_1;
	int N1;

	// essa vai ser uma matrix auxiliar
	gsl_matrix *Matrix_2;
	int N2;

	// essa vai ser outra matrix auxiliar
	gsl_matrix *Matrix_3;
	int N3;

        // essa vai ser outra matrix auxiliar
	gsl_matrix *Matrix_4;
	int N4;

	//este é uma variável auxiliar para a inversao
	gsl_permutation *p;
};

typedef struct matrizes_auziliares matrizes_auziliares;


void load_matrix_M(matrizes_auziliares *input)
{

	int N_mat=30;
        int N_mat2=19;

	input->Matrix_1=gsl_matrix_alloc(N_mat, N_mat);
	input->Matrix_2=gsl_matrix_alloc(N_mat, N_mat);
	input->Matrix_3=gsl_matrix_alloc(N_mat2, N_mat2);
        input->Matrix_4=gsl_matrix_alloc(N_mat2, N_mat2);

	// Aqui cria um atalho para você acessar os elementos, se precisar, use _Matrix[i][j]
	double (*_Matrix_1)[N_mat] = (double (*)[N_mat]) gsl_matrix_ptr(input->Matrix_1, 0, 0);


	FILE *file;
	file = fopen("./covmat.dat", "r");
	int state;
	double aux;

	for(int i=0; i<30;i++) // counter of columns of the file
	{ 
		for(int j=0; j<30;j++) // counter of columns of the file
		{

			state=fscanf(file,"%lf", &aux); // scan the column functions
			if(i<N_mat && j<N_mat)
			{
				_Matrix_1[i][j]=aux;
			}

		}
	}

	input->N1=N_mat;
	input->N2=N_mat;
	input->N3=N_mat2;
        input->N4=N_mat2;

       input->p=gsl_permutation_alloc(N_mat2);


}

void free_matrizes_auxiliares(matrizes_auziliares *input)
{

	gsl_matrix_free(input->Matrix_1);
	gsl_matrix_free(input->Matrix_2);
	gsl_matrix_free(input->Matrix_3);
	gsl_matrix_free(input->Matrix_4);
        gsl_permutation_free(input->p);
}



double chiMB_anue(int exp, int rule, int np, double *x, double *errors, void* user_data)
{  
	double chi2 = 0.0;
	double *signal_fit_rates = glbGetSignalFitRatePtr(0,0);

	double signal_rate1;
 
	double DATA_rates_MB_anue[]={122.0, 70.0, 65.0, 43.0, 57.0, 39.0, 37.0, 23.0, 22.0, 30.0, 43.0};
     	double bg_rates_MB_anue[]={90.289907, 53.077595, 57.098801, 32.937945, 43.159072, 34.174322, 36.383542, 28.737807, 22.339305, 26.509072, 42.697791}; 
        double pred_muon_a[] = {9998.957451, 13461.301884, 11298.240453, 7604.960449, 4331.886940, 2125.537108, 891.222608, 336.987112};
        double obs_muon_a[] = {9481.0, 13581.0, 11308.0, 7667.0, 4682.0, 2371.0, 985.0, 380.0};

        double vec_P[] = {signal_fit_rates[0], signal_fit_rates[1], signal_fit_rates[2], signal_fit_rates[3], signal_fit_rates[4], signal_fit_rates[5], signal_fit_rates[6], signal_fit_rates[7], signal_fit_rates[8], signal_fit_rates[9], signal_fit_rates[10], bg_rates_MB_anue[0], bg_rates_MB_anue[1], bg_rates_MB_anue[2], bg_rates_MB_anue[3], bg_rates_MB_anue[4], bg_rates_MB_anue[5], bg_rates_MB_anue[6], bg_rates_MB_anue[7], bg_rates_MB_anue[8], bg_rates_MB_anue[9], bg_rates_MB_anue[10], pred_muon_a[0], pred_muon_a[1], pred_muon_a[2], pred_muon_a[3], pred_muon_a[4], pred_muon_a[5], pred_muon_a[6], pred_muon_a[7]};

        double vec_T[] = {DATA_rates_MB_anue[0] - (signal_fit_rates[0] + bg_rates_MB_anue[0]), DATA_rates_MB_anue[1] - (signal_fit_rates[1] + bg_rates_MB_anue[1]), DATA_rates_MB_anue[2] - (signal_fit_rates[2] + bg_rates_MB_anue[2]), DATA_rates_MB_anue[3] - (signal_fit_rates[3] + bg_rates_MB_anue[3]), DATA_rates_MB_anue[4] - (signal_fit_rates[4] + bg_rates_MB_anue[4]), DATA_rates_MB_anue[5] - (signal_fit_rates[5] + bg_rates_MB_anue[5]), DATA_rates_MB_anue[6] - (signal_fit_rates[6] + bg_rates_MB_anue[6]), DATA_rates_MB_anue[7] - (signal_fit_rates[7] + bg_rates_MB_anue[7]), DATA_rates_MB_anue[8] - (signal_fit_rates[8] + bg_rates_MB_anue[8]), DATA_rates_MB_anue[9] - (signal_fit_rates[9] + bg_rates_MB_anue[9]), DATA_rates_MB_anue[10] - (signal_fit_rates[10] + bg_rates_MB_anue[10]),/*start numu osc here*/ pred_muon_a[0]-obs_muon_a[0],pred_muon_a[1]-obs_muon_a[1],pred_muon_a[2]-obs_muon_a[2],pred_muon_a[3]-obs_muon_a[3],pred_muon_a[4]-obs_muon_a[4],pred_muon_a[5]-obs_muon_a[5],pred_muon_a[6]-obs_muon_a[6],pred_muon_a[8]-obs_muon_a[8] /*end numu osc here*/};

	matrizes_auziliares *my_matrix=((matrizes_auziliares *) user_data);


	double (*_Matrix_1)[my_matrix->N1] = (double (*)[my_matrix->N1]) gsl_matrix_ptr(my_matrix->Matrix_1, 0, 0);
	double (*_Matrix_2)[my_matrix->N2] = (double (*)[my_matrix->N2]) gsl_matrix_ptr(my_matrix->Matrix_2, 0, 0);
	double (*_Matrix_3)[my_matrix->N3] = (double (*)[my_matrix->N3]) gsl_matrix_ptr(my_matrix->Matrix_3, 0, 0);
        double (*_Matrix_4)[my_matrix->N4] = (double (*)[my_matrix->N4]) gsl_matrix_ptr(my_matrix->Matrix_4, 0, 0);

	int s;

	for(int i=0; i<my_matrix->N1; i++)
	{
              //  if(i<11){signal_rate1 = (signal_fit_rates[i]);}
              //  else{signal_rate1=0.0;}

		for(int j=0; j<my_matrix->N1; j++)
		{

			//fit_rate2 = (signal_fit_rates[j] + bg_rates_MB_anue[j]);
			if(i==j && i<11){
				_Matrix_2[i][j]=_Matrix_1[i][j]*vec_P[i]*vec_P[j]+signal_fit_rates[i];}
			else{
				_Matrix_2[i][j]=_Matrix_1[i][j]*vec_P[i]*vec_P[j];}
		}

	}

       for(int k=0; k<11; k++) {for(int l=0; l<11; l++) {_Matrix_4[k][l] = _Matrix_2[k][l] + _Matrix_2[k+11][l] + _Matrix_2[k][l+11] + _Matrix_2[k+11][l+11];}}

       for(int k=0; k<11; k++) {for(int l=11; l<19; l++) {_Matrix_4[k][l] = _Matrix_2[k][l+11] + _Matrix_2[k+11][l+11];}}   
     
       for(int k=11; k<19; k++) {for(int l=0; l<11; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l] + _Matrix_2[k+11][l+11];}}
        
       for(int k=11; k<19; k++) {for(int l=11; l<19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11];}}

	// Notice: This function modifies the auxiliar matix Matrix_2. DONT USE MATRIX_1 !!
	gsl_linalg_LU_decomp(my_matrix->Matrix_4, my_matrix->p, &s); 


	double det=gsl_linalg_LU_det(my_matrix->Matrix_4, s);

	gsl_linalg_LU_invert(my_matrix->Matrix_4, my_matrix->p, my_matrix->Matrix_3);

        for(int i=0; i< 11/*my_matrix->N3*/; i++)
	{
		for(int j=0; j< 11 /*my_matrix->N3*/; j++)
		{
			chi2=chi2+(vec_T[i])*(_Matrix_3[i][j])*(vec_T[j]);
		}
	}

  for(int i=0; i< 11/*my_matrix->N3*/; i++)
	{
		for(int j=11; j< 19 /*my_matrix->N3*/; j++)
		{
			chi2=chi2+(vec_T[i])*(_Matrix_3[i][j])*(vec_T[j]);
		}
	}

  for(int i=11; i< 19/*my_matrix->N3*/; i++)
	{
		for(int j=0; j< 11 /*my_matrix->N3*/; j++)
		{
			chi2=chi2+(vec_T[i])*(_Matrix_3[i][j])*(vec_T[j]);
		}
	}

  return chi2;

}

