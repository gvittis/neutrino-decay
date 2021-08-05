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

        // essa vai ser outra matrix auxiliar
	gsl_matrix *Matrix_5;
	int N5;

        // essa vai ser outra matrix auxiliar
	gsl_matrix *Matrix_6;
	int N6;

	//este é uma variável auxiliar para a inversao
	gsl_permutation *p;

        //este é uma variável auxiliar para a inversao
	gsl_permutation *g;
};

typedef struct matrizes_auziliares matrizes_auziliares;


void load_matrix_M(matrizes_auziliares *input)
{

	int N_mat=30;
        int N_mat2=19;
        int N_mat5=8;

	input->Matrix_1=gsl_matrix_alloc(N_mat, N_mat);
	input->Matrix_2=gsl_matrix_alloc(N_mat, N_mat);
	input->Matrix_3=gsl_matrix_alloc(N_mat2, N_mat2);
        input->Matrix_4=gsl_matrix_alloc(N_mat2, N_mat2);
        input->Matrix_5=gsl_matrix_alloc(N_mat5, N_mat5);
        input->Matrix_6=gsl_matrix_alloc(N_mat5, N_mat5);

	// Aqui cria um atalho para você acessar os elementos, se precisar, use _Matrix[i][j]
	double (*_Matrix_1)[N_mat] = (double (*)[N_mat]) gsl_matrix_ptr(input->Matrix_1, 0, 0);


	FILE *file;
	file = fopen("./covanu.dat", "r");
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
        input->N5=N_mat5;
        input->N6=N_mat5;

       input->p=gsl_permutation_alloc(N_mat2);
       input->g=gsl_permutation_alloc(N_mat5);


}

void free_matrizes_auxiliares(matrizes_auziliares *input)
{

	gsl_matrix_free(input->Matrix_1);
	gsl_matrix_free(input->Matrix_2);
	gsl_matrix_free(input->Matrix_3);
	gsl_matrix_free(input->Matrix_4);
        gsl_matrix_free(input->Matrix_5);
	gsl_matrix_free(input->Matrix_6);
        gsl_permutation_free(input->p);
        gsl_permutation_free(input->g);
}

/*I am just putting the antineutrino variable in terms of neutrino ones to make the calculation easier*/

double chiMB_anue(int exp, int rule, int np, double *x, double *errors, void* user_data)
{  
	double chi2 = 0.0;
        double chi2mm = 0.0;
	double *signal_fit_rates = glbGetSignalFitRatePtr(0,0);
 
	double DATA_rates_MB_nue[]={122.0, 70.0, 65.0, 43.0, 57.0, 39.0, 37.0, 23.0, 22.0, 30.0, 43.0};
     	double bg_rates_MB_nue[]={90.289907, 53.077595, 57.098801, 32.937945, 43.159072, 34.174322, 36.383542, 28.737807, 22.339305, 26.509072, 42.697791}; 
        double pred_muon[] = {9998.957451, 13461.301884, 11298.240453, 7604.960449, 4331.886940, 2125.537108, 891.222608, 336.987112};
        double obs_muon[] = {9481.0, 13581.0, 11308.0, 7667.0, 4682.0, 2371.0, 985.0, 380.0};

        double vec_P[] = {signal_fit_rates[0], signal_fit_rates[1], signal_fit_rates[2], signal_fit_rates[3], signal_fit_rates[4], signal_fit_rates[5], signal_fit_rates[6], signal_fit_rates[7], signal_fit_rates[8], signal_fit_rates[9], signal_fit_rates[10], bg_rates_MB_nue[0], bg_rates_MB_nue[1], bg_rates_MB_nue[2], bg_rates_MB_nue[3], bg_rates_MB_nue[4], bg_rates_MB_nue[5], bg_rates_MB_nue[6], bg_rates_MB_nue[7], bg_rates_MB_nue[8], bg_rates_MB_nue[9], bg_rates_MB_nue[10], pred_muon[0], pred_muon[1], pred_muon[2], pred_muon[3], pred_muon[4], pred_muon[5], pred_muon[6], pred_muon[7]};

        double vec_T[] = {DATA_rates_MB_nue[0] - (signal_fit_rates[0] + bg_rates_MB_nue[0]), DATA_rates_MB_nue[1] - (signal_fit_rates[1] + bg_rates_MB_nue[1]), DATA_rates_MB_nue[2] - (signal_fit_rates[2] + bg_rates_MB_nue[2]), DATA_rates_MB_nue[3] - (signal_fit_rates[3] + bg_rates_MB_nue[3]), DATA_rates_MB_nue[4] - (signal_fit_rates[4] + bg_rates_MB_nue[4]), DATA_rates_MB_nue[5] - (signal_fit_rates[5] + bg_rates_MB_nue[5]), DATA_rates_MB_nue[6] - (signal_fit_rates[6] + bg_rates_MB_nue[6]), DATA_rates_MB_nue[7] - (signal_fit_rates[7] + bg_rates_MB_nue[7]), DATA_rates_MB_nue[8] - (signal_fit_rates[8] + bg_rates_MB_nue[8]), DATA_rates_MB_nue[9] - (signal_fit_rates[9] + bg_rates_MB_nue[9]), DATA_rates_MB_nue[10] - (signal_fit_rates[10] + bg_rates_MB_nue[10]),/*start numu osc here*/ obs_muon[0]-pred_muon[0],obs_muon[1]-pred_muon[1],obs_muon[2]-pred_muon[2],obs_muon[3]-pred_muon[3],obs_muon[4]-pred_muon[4],obs_muon[5]-pred_muon[5],obs_muon[6]-pred_muon[6],obs_muon[7]-pred_muon[7] /*end numu osc here*/};

        double dm[]={obs_muon[0]-pred_muon[0],obs_muon[1]-pred_muon[1],obs_muon[2]-pred_muon[2],obs_muon[3]-pred_muon[3],obs_muon[4]-pred_muon[4],obs_muon[5]-pred_muon[5],obs_muon[6]-pred_muon[6],obs_muon[7]-pred_muon[7]};

	matrizes_auziliares *my_matrix=((matrizes_auziliares *) user_data);


	double (*_Matrix_1)[my_matrix->N1] = (double (*)[my_matrix->N1]) gsl_matrix_ptr(my_matrix->Matrix_1, 0, 0);
	double (*_Matrix_2)[my_matrix->N2] = (double (*)[my_matrix->N2]) gsl_matrix_ptr(my_matrix->Matrix_2, 0, 0);
	double (*_Matrix_3)[my_matrix->N3] = (double (*)[my_matrix->N3]) gsl_matrix_ptr(my_matrix->Matrix_3, 0, 0);
        double (*_Matrix_4)[my_matrix->N4] = (double (*)[my_matrix->N4]) gsl_matrix_ptr(my_matrix->Matrix_4, 0, 0);
        double (*_Matrix_5)[my_matrix->N5] = (double (*)[my_matrix->N5]) gsl_matrix_ptr(my_matrix->Matrix_5, 0, 0);
        double (*_Matrix_6)[my_matrix->N6] = (double (*)[my_matrix->N6]) gsl_matrix_ptr(my_matrix->Matrix_6, 0, 0);

	int s;
        int t;

	for(int i=0; i<my_matrix->N1; i++)
	{
              //  if(i<11){signal_rate1 = (signal_fit_rates[i]);}
              //  else{signal_rate1=0.0;}

		for(int j=0; j<my_matrix->N1; j++)
		{

			//fit_rate2 = (signal_fit_rates[j] + bg_rates_MB_nue[j]);
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

       for(int k=0; k<8; k++) {for(int l=0; l<8; l++) {_Matrix_5[k][l] = _Matrix_4[k+11][l+11];}}

        gsl_linalg_LU_decomp(my_matrix->Matrix_5, my_matrix->g, &t); 
	gsl_linalg_LU_invert(my_matrix->Matrix_5, my_matrix->g, my_matrix->Matrix_6);


	// Notice: This function modifies the auxiliar matix Matrix_2. DONT USE MATRIX_1 !!
	//double det=gsl_linalg_LU_det(my_matrix->Matrix_4, s);
        gsl_linalg_LU_decomp(my_matrix->Matrix_4, my_matrix->p, &s); 
	gsl_linalg_LU_invert(my_matrix->Matrix_4, my_matrix->p, my_matrix->Matrix_3);

        for(int i=0; i< 19 /*my_matrix->N3*/; i++)
	{
		for(int j=0; j< 19 /*my_matrix->N3*/; j++)
		{
			chi2=chi2+(vec_T[i])*(_Matrix_3[i][j])*(vec_T[j]);
		}
	}

       for(int i=0; i<8 /*my_matrix->N3*/; i++)
	{
		for(int j=0; j< 8/*my_matrix->N3*/; j++)
		{
			chi2mm=chi2mm+(dm[i])*(_Matrix_6[i][j])*(dm[j]);
		}
	}


  return chi2 - chi2mm;

}

