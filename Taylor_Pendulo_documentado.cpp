// Este c�digo implementa el m�todo de Taylor a orden P
// para el problema del p�ndulo  f�sico
// La solucion es mandada al archivo de texto Pendulo.txt con el siguiente formato
// t     W[0]    W[1]    W[2]    W[3]
// donde W[0] y W[1] representan el �ngulo y la velocidad angular y, W[2]  y  W[3]
// son variables auxiliares para hacer de la ecuacion del p�ndulo un sistema de 4 
// ecuaciones de primer orden con flujo polinomial.
// Autores: Jaime Burgos Garc�a, Sim�n Rodr�guez Rodr�guez
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>

using namespace std;

double **Coeff_Taylor; // Puntero a un arreglo que contendr� los coeficientes de la serie de Taylor
double h;              // Tama�o de paso
int ORDER,GDLE;
double Taylor(double t, double *W);

//GDLE=Grados de Libertad Efectivos
//ORDER+1=Orden del m�todo de Taylor

int main()
{
	clock_t t_ini, t_fin; // varibable para almacenar tiempos de inicio y final de ejecuci�n
    double secs;          // tiempo de ejecuci�n en segundos
	
	
	t_ini = clock();  //Toma el tiempo inicial
	
	
	int i,j;
	double t=0;
	double *W;      //puntero a un arreglo que contendr� las soluciones al tiempo t
	
	GDLE=4;
	ORDER=11;
	ORDER=ORDER+1;
	
	W=new double [GDLE]; // Arreglo que contendr� las soluciones a un tiempo ti
	
	//Las siguientes dos lineas asignan la memoria para la tabla de coeficientes de Taylor
	Coeff_Taylor=new double *[GDLE];
	for(i=0;i<GDLE;i++) Coeff_Taylor[i]=new double[ORDER];
	
	
	ofstream datos; //crea el flujo para el archivo donde se guardar� la soluci�n
	datos.open("Pendulo.txt"); //La soluci�n se manda a un archivo de texto Pendulo.txt
	
	datos.setf(ios::fixed);
	datos.setf(ios::showpoint);
	datos.precision(16);
	
	//Tama�o de paso y condiciones iniciales
	
	
	
	h=0.1;   
	W[0]=1;
	W[1]=0;
	W[2]=sin(W[0]);
	W[3]=cos(W[0]);
	
	
	datos<<t<<"    "<<W[0]<<"    "<<W[1]<<"    "<<W[2]<<"    "<<W[3]<<endl;  //salida a archivo
	for(i=0;i<100;i++)
	{
		t=Taylor(t,W); //Llamada a la funci�n Taylor que calcula la soluci�n del siguinte paso de tiempo
		datos<<t<<"    "<<W[0]<<"    "<<W[1]<<"    "<<W[2]<<"    "<<W[3]<<endl; //salida a archivo
	}
	
	datos.close();
	
	t_fin = clock();     //Toma el tiempo final
	
	secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;   //Tiempo de ejecuci�n en segundos
	cout<<"\n\nTiempo de ejecucion: "<<secs<<"secs"<<endl<<endl;
	
	system("PAUSE");
	
}

//Funci�n que calcula la soluci�n para la EDO del pendulo simple
double Taylor(double t, double *W)
{
	char a;
	int i,j;
	
	
	for(i=0;i<GDLE;i++) for(j=0;j<ORDER;j++)  Coeff_Taylor[i][j]=0.0; //Limpia la tabla de coeficientes
	for(i=0;i<GDLE;i++)  Coeff_Taylor[i][0]=W[i];  //Fija los coeficientes de Taylor de orden cero a las condiciones iniciales
	
	// C�lculo de los coeficientes de la expanci�n
	for(i=1;i<ORDER;i++)
	{
	    Coeff_Taylor[0][i]= Coeff_Taylor[1][i-1]/i;

	    Coeff_Taylor[1][i]=-Coeff_Taylor[2][i-1]/i;
	    
		for(j=0;j<i;j++) Coeff_Taylor[2][i]=Coeff_Taylor[2][i]+Coeff_Taylor[1][i-1-j]*Coeff_Taylor[3][j];
		Coeff_Taylor[2][i]=Coeff_Taylor[2][i]/i;
		
		for(j=0;j<i;j++) Coeff_Taylor[3][i]=Coeff_Taylor[3][i]-Coeff_Taylor[1][i-1-j]*Coeff_Taylor[2][j];
		Coeff_Taylor[3][i]=Coeff_Taylor[3][i]/i;
	}
	
	// C�lcula de la soluci�n 
	for(i=0;i<GDLE;i++)
	   for(j=1;j<ORDER;j++)
	   	   W[i]=W[i]+Coeff_Taylor[i][j]*pow(h,j);
	return t+h;  // regresa el siguiente paso de tiempo
}
