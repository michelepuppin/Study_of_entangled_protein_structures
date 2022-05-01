#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

using namespace std;
using namespace Eigen;

struct pos {double x=0.0L; double y=0.0L;	double z=0.0L; double bfact=0.0L;};
struct sotmat {double H1=0.0L; double H2=0.0L; double H3=0.0L;	double H4=0.0L; double H5=0.0L; double H6=0.0L; double H7=0.0L;	double H8=0.0L; double H9=0.0L;};

int main () {
//Leggo il file del tipo (x, y, z, b-factor)
vector <double> v;
double temp=0.0;

fstream file("inputfile.txt");

if (file.is_open())
	{
		while (file>>temp)
		{
		v.push_back(temp);
		}
		file.close();
	}
else
{
	cout << "Impossibile aprire il file" << endl;
}

int size=v.size();
int righe=size/4;
int i=0;

//Creo un vettore di strutture posizione e lo riempio
vector <pos> res(righe);

while (i<righe)
	{
		res[i].x=v.at(4*i);
		res[i].y=v.at(4*i+1);
		res[i].z=v.at(4*i+2);
		res[i].bfact=v.at(4*i+3);
		i++;
	}

//Costruisco la matrice di dimensione NxN con elementi le sottomatrici 3x3 date dalla struttura 
double gamma=0.0L;
double raggioc=0.0L;
int var=0;
double R0=0.0L;
double distance=0.0L;
sotmat M[righe][righe];

//Scelgo il valore del raggio di cutoff 
raggioc=18.0;

for (int i=0; i<righe; i++)
{
	for (int j=0; j<righe; j++)
	{
		R0=(res[i].x-res[j].x)*(res[i].x-res[j].x)+(res[i].y-res[j].y)*(res[i].y-res[j].y)+(res[i].z-res[j].z)*(res[i].z-res[j].z);
		distance=sqrt(R0);

		var=i-j;
		if(distance <= raggioc) //E' stata posta come distanza di cutoff 18A
		{
			if((var==-1)||(var==1)) {gamma=-10.0L;}
        	else {gamma=-1.0L;}
		}
		else {gamma=0.0L;} 

		M[i][j].H1=(gamma/R0)*(res[i].x-res[j].x)*(res[i].x-res[j].x);
		M[i][j].H2=(gamma/R0)*(res[i].x-res[j].x)*(res[i].y-res[j].y);
		M[i][j].H3=(gamma/R0)*(res[i].x-res[j].x)*(res[i].z-res[j].z);
		M[i][j].H4=(gamma/R0)*(res[i].x-res[j].x)*(res[i].y-res[j].y);
		M[i][j].H5=(gamma/R0)*(res[i].y-res[j].y)*(res[i].y-res[j].y);
		M[i][j].H6=(gamma/R0)*(res[i].y-res[j].y)*(res[i].z-res[j].z);
		M[i][j].H7=(gamma/R0)*(res[i].x-res[j].x)*(res[i].z-res[j].z);
		M[i][j].H8=(gamma/R0)*(res[i].y-res[j].y)*(res[i].z-res[j].z);
		M[i][j].H9=(gamma/R0)*(res[i].z-res[j].z)*(res[i].z-res[j].z);
	}
}

//Creo una matrice nellambiente Eigen 
MatrixXd W(3*righe,3*righe);

double diag1=0.0L;
double diag2=0.0L;
double diag3=0.0L;
double diag4=0.0L;
double diag5=0.0L;
double diag6=0.0L;
double diag7=0.0L;
double diag8=0.0L;
double diag9=0.0L;

for (int i=0; i<righe; i++)
	{
	for (int j=0; j<righe; j++)
		{
		if (i!=j)
			{
			W(3*i,3*j)=M[i][j].H1;
			W(3*i,3*j+1)=M[i][j].H2;
			W(3*i,3*j+2)=M[i][j].H3;
			W(3*i+1,3*j)=M[i][j].H4;
			W(3*i+1,3*j+1)=M[i][j].H5;
			W(3*i+1,3*j+2)=M[i][j].H6;
			W(3*i+2,3*j)=M[i][j].H7;
			W(3*i+2,3*j+1)=M[i][j].H8;
			W(3*i+2,3*j+2)=M[i][j].H9;
			}
			else 
			{
			for(int j=0; j<(righe); j++)
				{
				if (i!=j)
					{
					diag1+=M[i][j].H1;
					diag2+=M[i][j].H2;
					diag3+=M[i][j].H3;
					diag4+=M[i][j].H4;
					diag5+=M[i][j].H5;
					diag6+=M[i][j].H6;
					diag7+=M[i][j].H7;
					diag8+=M[i][j].H8;
					diag9+=M[i][j].H9;
					}
				}
			W(3*i,3*j)=(-1)*diag1;
			W(3*i,3*j+1)=(-1)*diag2;
			W(3*i,3*j+2)=(-1)*diag3;
			W(3*i+1,3*j)=(-1)*diag4;
			W(3*i+1,3*j+1)=(-1)*diag5;
			W(3*i+1,3*j+2)=(-1)*diag6;
			W(3*i+2,3*j)=(-1)*diag7;
			W(3*i+2,3*j+1)=(-1)*diag8;
			W(3*i+2,3*j+2)=(-1)*diag9;
			
			diag1=0.0L;
			diag2=0.0L;
			diag3=0.0L;
			diag4=0.0L;
			diag5=0.0L;
			diag6=0.0L;
			diag7=0.0L;
			diag8=0.0L;
			diag9=0.0L;
			}
	}
}

//Stampo la matrice completa
ofstream fileout99;
fileout99.open("matrice.txt");

if (fileout99.is_open())
{
	for (int i=0; i<(3*righe); i++)
	{
		for (int j=0; j<(3*righe); j++)
		{
			fileout99 << W(i,j) << " ";
		}
	fileout99 << endl;
	}
	
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

//Calcolo l'area
double area1=0.0L;

for (int i=0; i<righe; i++)
	{
		area1+=res[i].bfact;
	}

//Scrivo il file con cui fare il grafico del B-factor
ofstream fileout2;
fileout2.open("bfactor.txt");

if (fileout2.is_open())
{
	for (int i=0; i<righe; i++)
	{
		fileout2 << i << " " << res[i].bfact << endl;
	}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

//Calcolo autovalori e autovettori tramite il pacchetto Eigen sfruttando il fatto che la matrice è simmatrica e reale e quindi autoaggiunta
SelfAdjointEigenSolver<MatrixXd> s(W);
s.compute(W);

//Scrivo gli autovalori in un file
ofstream fileout7;
fileout7.open("autovalori.txt");

if (fileout7.is_open())
{
	fileout7 << s.eigenvalues() << endl;
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

//Scrivo gli autovettori su un file
/*ofstream fileout13;
fileout13.open("autovettori.txt");

if (fileout13.is_open())
{
	fileout13 << s.eigenvalues()[6] << " " << s.eigenvalues()[7] << " " << s.eigenvalues()[121] << endl;
	
	for(int i=0; i<(3*righe); i++)
		{
		//for(int k=0 k<(3*righe); k++)
		//	{
			fileout13 << s.eigenvectors().col(6)[i] << " " << s.eigenvectors().col(7)[i] << " " << s.eigenvectors().col(121)[i] << endl;
		//	}
		//fileout13 << endl;
		}		
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}*/

//Calcolo i deltaRi
double deltaRik[righe][3*righe];
MatrixXd P(3*righe,3*righe);
VectorXd U(3);

for(int i=0; i<(righe); i++)
{
	for(int k=0; k<(3*righe); k++)
		{
		U(0)=s.eigenvectors().col(k)[3*i];
		U(1)=s.eigenvectors().col(k)[3*i+1];
		U(2)=s.eigenvectors().col(k)[3*i+2];
		P=U*U.transpose();
		deltaRik[i][k]=(1/(sqrt((s.eigenvalues()[k])*(s.eigenvalues()[k]))))*(P.trace());
		}
}

//N.B.: sopra ho ottenuto una matrice della forma
// deltaR1  (k=1,k=2,...,k=3N)
// deltaR2	(k=1,k=2,...,k=3N)	      
// deltaR3	(k=1,k=2,...,k=3N)            
// ...........................
// deltaRN	(k=1,k=2,...,k=3N)

//Sommo sui k tralasciando i primi sei autovalori
vector<double> deltaRi(righe);

for(int i=0; i<(righe); i++)
{
	for(int k=6; k<(3*righe); k++)
		{
		deltaRi.at(i)+=deltaRik[i][k];
		}
}

//Calcolo l'area
double area2=0.0L;

for (int i=0; i<righe; i++)
	{
		area2+=deltaRi.at(i);
	}

//Scrivo il file con cui fare il grafico del delraR2 in cui "normalizzo" moltiplicando per l'area del bfactor e dividento per l'area di deltaR2
ofstream fileout4;
fileout4.open("deltaRi.txt");

if (fileout4.is_open())
{
	for (int i=0; i<righe; i++)
	{
		fileout4 << i << " " << (area1/area2)*(deltaRi.at(i)) << endl;
	}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

//Scrivo un file per visualizzare i profili di mobilità separatamente per ogni modo
ofstream fileout33;
fileout33.open("deltaRikfiss.txt");

if (fileout33.is_open())
{
	for (int i=0; i<righe; i++)
	{
		fileout33 << i << " " << deltaRik[i][6] << " " << deltaRik[i][7] << " " << deltaRik[i][8] << " " << deltaRik[i][9] << " " << deltaRik[i][10] << " " << deltaRik[i][11] << endl;
	}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

//Calcolo il Kk
//Invece di sommare sui k ora sommo sugli i
double norm=0.0L;
double somma=0.0L;
double alfa=0.0L;
double kappa[3*righe];
double N=0.0L;
N=righe*1.0L;

for(int k=0; k<(3*righe); k++)
{
	for(int i=0; i<(righe); i++)
		{
		norm+=deltaRik[i][k];		 
		}

	alfa=(1/norm);
	
	for(int i=0; i<(righe); i++)
		{
		somma+=alfa*(deltaRik[i][k])*log(alfa*(deltaRik[i][k]));		 
		}	
	
	kappa[k]=(exp((-1)*somma))/N;
	norm=0.0L;
	alfa=0.0L;
	somma=0.0L;
}

//Stampo i K_k
ofstream fileout9;
fileout9.open("kappa.txt");

if (fileout9.is_open())
{
	for (int i=6; i<(3*righe); i++)
	{
		fileout9 << i << " " << kappa[i] << endl;
	}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

return 0;
}