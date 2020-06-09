#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
#include<omp.h>
using namespace std;

double** metodoGauss(double** A,int m, int n, double* resul);

int getDimension(const char* name,int dim){
    int rows=0;
    int cols=0;
    FILE *archivo;
 	char caracteres[1000];
    char * token;

 	archivo = fopen(name,"r");
 	
 	while (feof(archivo) == 0){
 		fgets(caracteres,100,archivo);
 		token = strtok(caracteres, ",");
        cols=0;
        while( token != NULL ) {
            token = strtok(NULL, ",");
            cols++;
        }
        rows++;
    }

    fclose(archivo);

    if(dim==0){
        return rows;
    }else{
        return cols;
    }
}

double** getMatriz(const char* name){
    int rows=getDimension(name,0);
    int cols=getDimension(name,1);
    int i,j;
    FILE *archivo;
 	double** mat;
 	char caracteres[1000];
    char * token;

 	archivo = fopen(name,"r");
 	
    mat = new double*[rows];
    for(int x=0; x<rows; x++){
        mat[x] = new double[cols];
    }
    i=0;
 	while (feof(archivo) == 0){
 		fgets(caracteres,100,archivo);
 		token = strtok(caracteres, ",");
        j=0;
        while( token != NULL ) {
            mat[i][j] = atof(token);
            //printf("%s ", token );
            token = strtok(NULL, ",");
            j++;
        }    
        i++;        
    }
    fclose(archivo);
    return mat;
}

void mostrarVector(double *V,int m){
		cout<<"\n[ \t";
	for (int i=0; i<m; i++)
    {
		cout<<V[i]<<"\t";
    }	
    cout<<" ]"<<endl;
    
}		

void mostrarMatriz(double** A,int m,int n){
	cout<<endl;
	for (int i=0; i<m; i++){
        cout<<"[\t ";
		for(int j=0;j<n;j++){
			cout<<A[i][j]<<"\t";		
		}
		cout<<" ]\n";
    }
}

int main(){
    double **M = getMatriz("matriz.txt");
    int rowsM1=getDimension("matriz.txt",0);
    int colsM1=getDimension("matriz.txt",1);
    double** resul;
    double* campos;
    campos = new double[rowsM1];
    resul = new double*[rowsM1];
    for (int i = 0; i < rowsM1; i++)
    {
        resul[i] = new double[colsM1];
    }
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        resul = metodoGauss(M,rowsM1,colsM1,campos);
    }
    cout<<"\t.:MATRIZ RESULTANTE:."<<endl;
    mostrarMatriz(resul,rowsM1,colsM1);
    cout<<endl;
    cout<<endl;
    cout<<"\t.:RESULTADOS DEL SISTEMA DE ECUACIONES:."<<endl;
    mostrarVector(campos,rowsM1);

    return 0;
}

double** metodoGauss(double** A,int m, int n,double* resul){
	double **matriz;
    matriz = new double* [m];
    for (int i = 0; i < m; i++) {
        matriz[i] = new double[n];
    }
	matriz=A;
	for (int k = 0; k < m - 1; k++){
        #pragma omp for
        for (int i = k + 1; i < m; i++){
            double factor = matriz[i] [k] / matriz[k] [k];
            
 	        for (int j = k; j <= m; j++){
                matriz[i] [j] = matriz[i] [j]- matriz[k] [j] * factor;
            }
        }
    }
    resul[m-1] = matriz[m-1] [m] / matriz[m-1][m-1];
	
    for (int f = m - 2; f >= 0; f--){
        double sum = 0;
        for (int i = m - 1; i > f; i--){
            sum += matriz[f][i] * resul[i];
        }
        resul[f] = (matriz[f][m] - sum) / matriz[f] [f];
    }
	return matriz;
}