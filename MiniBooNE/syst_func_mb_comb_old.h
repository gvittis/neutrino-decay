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

	int N_mat=60;
        int N_mat2=38;

	input->Matrix_1=gsl_matrix_alloc(N_mat, N_mat);
	input->Matrix_2=gsl_matrix_alloc(N_mat, N_mat);
	input->Matrix_3=gsl_matrix_alloc(N_mat2, N_mat2);
        input->Matrix_4=gsl_matrix_alloc(N_mat2, N_mat2);

	// Aqui cria um atalho para você acessar os elementos, se precisar, use _Matrix[i][j]
	double (*_Matrix_1)[N_mat] = (double (*)[N_mat]) gsl_matrix_ptr(input->Matrix_1, 0, 0);


	FILE *file;
	file = fopen("./miniboone_full_fractcovmatrix_combined_lowe.txt", "r");
	int state;
	double aux;

	for(int i=0; i<60;i++) // counter of columns of the file
	{ 
		for(int j=0; j<60;j++) // counter of columns of the file
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
	double *signal_fit_rates = glbGetSignalFitRatePtr(1,0);
        double *signal_fit_rates_a = glbGetSignalFitRatePtr(0,0);

	//double signal_rate1;
 
	
        double DATA_rates_MB_nue[]={497.0, 283.0, 313.0, 167.0, 184.0, 163.0, 140.0, 115.0, 97.0, 98.0, 130.0};
     	double bg_rates_MB_nue[]={361.002334, 216.002142, 239.436776, 127.517957, 179.035344, 133.901816, 139.020389, 113.446978, 81.204519, 98.603919, 137.953204}; 
        double pred_muon[] = {38564.217639, 59339.405335, 53069.519495, 37171.337542, 23002.153188, 12423.361945, 6012.845025, 2801.295291};
        double obs_muon[] = {37676.0, 59515.0, 53126.0, 37050.0, 22150.0, 11478.0, 5374.0, 2547.0};

        double DATA_rates_MB_anue[]={122.0, 70.0, 65.0, 43.0, 57.0, 39.0, 37.0, 23.0, 22.0, 30.0, 43.0};
     	double bg_rates_MB_anue[]={90.289907, 53.077595, 57.098801, 32.937945, 43.159072, 34.174322, 36.383542, 28.737807, 22.339305, 26.509072, 42.697791}; 
        double pred_muon_a[] = {9998.957451, 13461.301884, 11298.240453, 7604.960449, 4331.886940, 2125.537108, 891.222608, 336.987112};
        double obs_muon_a[] = {9481.0, 13581.0, 11308.0, 7667.0, 4682.0, 2371.0, 985.0, 380.0};

        double vec_P_a[] = {signal_fit_rates[0], signal_fit_rates[1], signal_fit_rates[2], signal_fit_rates[3], signal_fit_rates[4], signal_fit_rates[5], signal_fit_rates[6], signal_fit_rates[7], signal_fit_rates[8], signal_fit_rates[9], signal_fit_rates[10], bg_rates_MB_nue[0], bg_rates_MB_nue[1], bg_rates_MB_nue[2], bg_rates_MB_nue[3], bg_rates_MB_nue[4], bg_rates_MB_nue[5], bg_rates_MB_nue[6], bg_rates_MB_nue[7], bg_rates_MB_nue[8], bg_rates_MB_nue[9], bg_rates_MB_nue[10], pred_muon[0], pred_muon[1], pred_muon[2], pred_muon[3], pred_muon[4], pred_muon[5], pred_muon[6], pred_muon[7], signal_fit_rates_a[0], signal_fit_rates_a[1], signal_fit_rates_a[2], signal_fit_rates_a[3], signal_fit_rates_a[4], signal_fit_rates_a[5], signal_fit_rates_a[6], signal_fit_rates_a[7], signal_fit_rates_a[8], signal_fit_rates_a[9], signal_fit_rates_a[10], bg_rates_MB_anue[0], bg_rates_MB_anue[1], bg_rates_MB_anue[2], bg_rates_MB_anue[3], bg_rates_MB_anue[4], bg_rates_MB_anue[5], bg_rates_MB_anue[6], bg_rates_MB_anue[7], bg_rates_MB_anue[8], bg_rates_MB_anue[9], bg_rates_MB_anue[10], pred_muon_a[0], pred_muon_a[1], pred_muon_a[2], pred_muon_a[3], pred_muon_a[4], pred_muon_a[5], pred_muon_a[6], pred_muon_a[7]};


        double vec_T_a[] = {DATA_rates_MB_nue[0] - (signal_fit_rates[0] + bg_rates_MB_nue[0]), DATA_rates_MB_nue[1] - (signal_fit_rates[1] + bg_rates_MB_nue[1]), DATA_rates_MB_nue[2] - (signal_fit_rates[2] + bg_rates_MB_nue[2]), DATA_rates_MB_nue[3] - (signal_fit_rates[3] + bg_rates_MB_nue[3]), DATA_rates_MB_nue[4] - (signal_fit_rates[4] + bg_rates_MB_nue[4]), DATA_rates_MB_nue[5] - (signal_fit_rates[5] + bg_rates_MB_nue[5]), DATA_rates_MB_nue[6] - (signal_fit_rates[6] + bg_rates_MB_nue[6]), DATA_rates_MB_nue[7] - (signal_fit_rates[7] + bg_rates_MB_nue[7]), DATA_rates_MB_nue[8] - (signal_fit_rates[8] + bg_rates_MB_nue[8]), DATA_rates_MB_nue[9] - (signal_fit_rates[9] + bg_rates_MB_nue[9]), DATA_rates_MB_nue[10] - (signal_fit_rates[10] + bg_rates_MB_nue[10]),/*start numu osc here*/ pred_muon[0]-obs_muon[0],pred_muon[1]-obs_muon[1],pred_muon[2]-obs_muon[2],pred_muon[3]-obs_muon[3],pred_muon[4]-obs_muon[4],pred_muon[5]-obs_muon[5],pred_muon[6]-obs_muon[6],pred_muon[8]-obs_muon[8]  /*end numu osc here*/, DATA_rates_MB_anue[0] - (signal_fit_rates_a[0] + bg_rates_MB_anue[0]), DATA_rates_MB_anue[1] - (signal_fit_rates_a[1] + bg_rates_MB_anue[1]), DATA_rates_MB_anue[2] - (signal_fit_rates_a[2] + bg_rates_MB_anue[2]), DATA_rates_MB_anue[3] - (signal_fit_rates_a[3] + bg_rates_MB_anue[3]), DATA_rates_MB_anue[4] - (signal_fit_rates_a[4] + bg_rates_MB_anue[4]), DATA_rates_MB_anue[5] - (signal_fit_rates_a[5] + bg_rates_MB_anue[5]), DATA_rates_MB_anue[6] - (signal_fit_rates_a[6] + bg_rates_MB_anue[6]), DATA_rates_MB_anue[7] - (signal_fit_rates_a[7] + bg_rates_MB_anue[7]), DATA_rates_MB_anue[8] - (signal_fit_rates_a[8] + bg_rates_MB_anue[8]), DATA_rates_MB_anue[9] - (signal_fit_rates_a[9] + bg_rates_MB_anue[9]), DATA_rates_MB_anue[10] - (signal_fit_rates_a[10] + bg_rates_MB_anue[10]),/*start numu osc here*/ pred_muon_a[0]-obs_muon_a[0],pred_muon_a[1]-obs_muon_a[1],pred_muon_a[2]-obs_muon_a[2],pred_muon_a[3]-obs_muon_a[3],pred_muon_a[4]-obs_muon_a[4],pred_muon_a[5]-obs_muon_a[5],pred_muon_a[6]-obs_muon_a[6],pred_muon_a[7]-obs_muon_a[7] /*end numu osc here*/};

//pred_muon_a[0]-obs_muon_a[0],pred_muon_a[1]-obs_muon_a[1],pred_muon_a[2]-obs_muon_a[2],pred_muon_a[3]-obs_muon_a[3],pred_muon_a[4]-obs_muon_a[4],pred_muon_a[5]-obs_muon_a[5],pred_muon_a[6]-obs_muon_a[6],pred_muon_a[8]-obs_muon_a[8]

//pred_muon[0]-obs_muon[0],pred_muon[1]-obs_muon[1],pred_muon[2]-obs_muon[2],pred_muon[3]-obs_muon[3],pred_muon[4]-obs_muon[4],pred_muon[5]-obs_muon[5],pred_muon[6]-obs_muon[6],pred_muon[8]-obs_muon[8]

	matrizes_auziliares *my_matrix=((matrizes_auziliares *) user_data);


	double (*_Matrix_1)[my_matrix->N1] = (double (*)[my_matrix->N1]) gsl_matrix_ptr(my_matrix->Matrix_1, 0, 0);
	double (*_Matrix_2)[my_matrix->N2] = (double (*)[my_matrix->N2]) gsl_matrix_ptr(my_matrix->Matrix_2, 0, 0);
	double (*_Matrix_3)[my_matrix->N3] = (double (*)[my_matrix->N3]) gsl_matrix_ptr(my_matrix->Matrix_3, 0, 0);
        double (*_Matrix_4)[my_matrix->N4] = (double (*)[my_matrix->N4]) gsl_matrix_ptr(my_matrix->Matrix_4, 0, 0);

	int s;

	for(int i=0; i< my_matrix->N1; i++)
	{
		for(int j=0; j< my_matrix->N1; j++)
		{
			if(i==j && 0<= i < 11){
                           _Matrix_2[i][j]=_Matrix_1[i][j]*vec_P_a[i]*vec_P_a[j]+signal_fit_rates[i];}
                        if(i==j && 30<= i < 41){
                           _Matrix_2[i][j]=_Matrix_1[i][j]*vec_P_a[i]*vec_P_a[j]+signal_fit_rates_a[i];}
			else{
				_Matrix_2[i][j]=_Matrix_1[i][j]*vec_P_a[i]*vec_P_a[j];}
		}

	}

       /* Bloco Neutrino (k=0-30  l=0-30)*/
       for(int k=0; k<11; k++) {for(int l=0; l<11; l++) {_Matrix_4[k][l] = _Matrix_2[k][l] + _Matrix_2[k+11][l] + _Matrix_2[k][l+11] + _Matrix_2[k+11][l+11];}}
       for(int k=11; k<19; k++) {for(int l=0; l<11; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l] + _Matrix_2[k+11][l+11];}}
       for(int k=0; k<11; k++) {for(int l=11; l<19; l++) {_Matrix_4[k][l] = _Matrix_2[k][l+11] + _Matrix_2[k+11][l+11];}}   
       for(int k=11; k<19; k++) {for(int l=11; l<19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11];}}

       /* Bloco Neutrino-Antineutrino (k=0-30  l=30-60)  */
       for(int k=0; k<11; k++) {for(int l=0+19; l<11+19; l++) {_Matrix_4[k][l] = _Matrix_2[k][l+11] + _Matrix_2[k+11][l+11] + _Matrix_2[k][l+11+11] + _Matrix_2[k+11][l+11+11];}} 
       for(int k=11; k<19; k++) {for(int l=0+19; l<11+19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11] + _Matrix_2[k+11][l+11+11];}} 
       for(int k=0; k<11; k++) {for(int l=11+19; l<19+19; l++) {_Matrix_4[k][l] = _Matrix_2[k][l+11+11] + _Matrix_2[k+11][l+11+11];}}  
       for(int k=11; k<19; k++) {for(int l=11+19; l<19+19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11+11];}}

       /* Bloco Antineutrino-Neutrino (k=30-60  l=0-30)  */
       for(int k=0+19; k<11+19; k++) {for(int l=0; l<11; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l] + _Matrix_2[k+11+11][l] + _Matrix_2[k+11][l+11] + _Matrix_2[k+11+11][l+11];}}
       for(int k=11+19; k<19+19; k++) {for(int l=0; l<11; l++) {_Matrix_4[k][l] = _Matrix_2[k+11+11][l] + _Matrix_2[k+11+11][l+11];}}
       for(int k=0+19; k<11+19; k++) {for(int l=11; l<19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11] + _Matrix_2[k+11+11][l+11];}}   
       for(int k=11+19; k<19+19; k++) {for(int l=11; l<19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11+11][l+11];}}

       /* Bloco Antieutrino (k=30-60  l=30-60)*/
       for(int k=0+19; k<11+19; k++) {for(int l=0+19; l<11+19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11] + _Matrix_2[k+11+11][l+11] + _Matrix_2[k+11][l+11+11] + _Matrix_2[k+11+11][l+11+11];}}
       for(int k=11+19; k<19+19; k++) {for(int l=0+19; l<11+19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11+11][l+11] + _Matrix_2[k+11+11][l+11+11];}}
       for(int k=0+19; k<11+19; k++) {for(int l=11+19; l<19+19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11][l+11+11] + _Matrix_2[k+11+11][l+11+11];}}   
       for(int k=11+19; k<19+19; k++) {for(int l=11+19; l<19+19; l++) {_Matrix_4[k][l] = _Matrix_2[k+11+11][l+11+11];}}

 /*for(int k=0; k<11; k++)
        { 
           for(int l=0; l<11; l++)
           {   
           _Matrix_4[k][l] = _Matrix_2[k][l] + _Matrix_2[k+11][l] + _Matrix_2[k][l+11] + _Matrix_2[k+11][l+11];
           _Matrix_4[k][l+19] = _Matrix_2[k][l+30] + _Matrix_2[k+11][l+30] + _Matrix_2[k][l+11+30] + _Matrix_2[k+11][l+11+30];
           _Matrix_4[k+19][l] = _Matrix_2[k+30][l] + _Matrix_2[k+11+30][l] + _Matrix_2[k+30][l+11] + _Matrix_2[k+11+30][l+11];
           _Matrix_4[k+19][l+19] = _Matrix_2[k+30][l+30] + _Matrix_2[k+11+30][l+30] + _Matrix_2[k+30][l+11+30] + _Matrix_2[k+11+30][l+11+30];
           }
        }

        for(int k=0; k<8; k++)
        { 
           for(int l=0; l<11; l++)
           {   
           _Matrix_4[k+11][l] = _Matrix_2[k+22][l] + _Matrix_2[k+22][l+11];
           _Matrix_4[l][k+11] = _Matrix_2[l][k+22] + _Matrix_2[l+11][k+22];
            
           _Matrix_4[k+11][l+19] = _Matrix_2[k+22][l+30] + _Matrix_2[k+22][l+11+30];
           _Matrix_4[l][k+11+19] = _Matrix_2[l][k+22+30] + _Matrix_2[l+11][k+22+30];
            
           _Matrix_4[k+11+19][l] = _Matrix_2[k+22+30][l] + _Matrix_2[k+22+30][l+11];
           _Matrix_4[l+19][k+11] = _Matrix_2[l+30][k+22] + _Matrix_2[l+11+30][k+22];
          
           _Matrix_4[k+11+19][l+19] = _Matrix_2[k+22+30][l+30] + _Matrix_2[k+22+30][l+11+30];
           _Matrix_4[l+19][k+11+19] = _Matrix_2[l+30][k+22+30] + _Matrix_2[l+11+30][k+22+30];
           }
        }

        for(int ki=0; ki<8; ki++)
        { 
           for(int kj=0; kj<8; kj++)
           {   
           _Matrix_4[ki+11][kj+11] = _Matrix_2[ki+22][kj+22];
           _Matrix_4[ki+11][kj+11+19] = _Matrix_2[ki+22][kj+22+30];
           _Matrix_4[ki+11+19][kj+11] = _Matrix_2[ki+22+30][kj+22];
           _Matrix_4[ki+11+19][kj+11+19] = _Matrix_2[ki+22+30][kj+22+30];
           }
        }*/

	// Notice: This function modifies the auxiliar matix Matrix_2. DONT USE MATRIX_1 !!
	gsl_linalg_LU_decomp(my_matrix->Matrix_4, my_matrix->p, &s); 


	double det=gsl_linalg_LU_det(my_matrix->Matrix_4, s);

	gsl_linalg_LU_invert(my_matrix->Matrix_4, my_matrix->p, my_matrix->Matrix_3);

	for(int i=0; i< 11/*my_matrix->N3*/; i++)
	{
		for(int j=0; j<11 /*my_matrix->N3*/; j++)
		{
		chi2=chi2+(vec_T_a[i])*(_Matrix_3[i][j])*(vec_T_a[j]) + (vec_T_a[i+19])*(_Matrix_3[i+19][j+19])*(vec_T_a[j+19]) + (vec_T_a[i+19])*(_Matrix_3[i+19][j])*(vec_T_a[j]) + (vec_T_a[i])*(_Matrix_3[i][j+19])*(vec_T_a[j+19]) ;
		}
	}

        for(int i=11; i<19 /*my_matrix->N3*/; i++)
	{
		for(int j=0; j< 11/*my_matrix->N3*/; j++)
		{
	        chi2=chi2+(vec_T_a[i])*(_Matrix_3[i][j])*(vec_T_a[j])+ (vec_T_a[i+19])*(_Matrix_3[i+19][j+19])*(vec_T_a[j+19]) + (vec_T_a[i+19])*(_Matrix_3[i+19][j])*(vec_T_a[j]) + (vec_T_a[i])*(_Matrix_3[i][j+19])*(vec_T_a[j+19]);
		}
	}

      for(int i=0; i<11 /*my_matrix->N3*/; i++)
	{
		for(int j=11; j< 19/*my_matrix->N3*/; j++)
		{
		chi2+(vec_T_a[i])*(_Matrix_3[i][j])*(vec_T_a[j])+ (vec_T_a[i+19])*(_Matrix_3[i+19][j+19])*(vec_T_a[j+19]) + (vec_T_a[i+19])*(_Matrix_3[i+19][j])*(vec_T_a[j]) + (vec_T_a[i])*(_Matrix_3[i][j+19])*(vec_T_a[j+19]);
		}
	}

  return chi2;

}

