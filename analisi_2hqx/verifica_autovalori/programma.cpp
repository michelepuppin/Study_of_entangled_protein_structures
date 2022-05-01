#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

using namespace std;
using namespace Eigen;

struct pos {double x=0.0L; double y=0.0L;	double z=0.0L;};
struct sotmat {double H1=0.0L; double H2=0.0L; double H3=0.0L;	double H4=0.0L; double H5=0.0L; double H6=0.0L; double H7=0.0L;	double H8=0.0L; double H9=0.0L;};

double linknumb (int righe, int i1, int i2, int j1, int j2, vector <pos> res);
pos centrmass (int righe, vector <pos> res);

int main () {
//Legge un file del tipo (x, y, z, b-factor)
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

//Crea un vettore di strutture pos
vector <pos> res(righe);
vector <double> bfact(righe);

while (i<righe)
	{
		res[i].x=v.at(4*i);
		res[i].y=v.at(4*i+1);
		res[i].z=v.at(4*i+2);
		bfact.at(i)=v.at(4*i+3);
		i++;
	}

double gamma=0.0L;
double raggioc=0.0L;
int var=0;
double R0=0.0L;
double distance=0.0L;
sotmat M[righe][righe];

MatrixXd W(3*righe,3*righe);
double diag1=0.0L, diag2=0.0L, diag3=0.0L, diag4=0.0L, diag5=0.0L, diag6=0.0L, diag7=0.0L, diag8=0.0L, diag9=0.0L;
SelfAdjointEigenSolver<MatrixXd> s(W);
vector <int> counteri; //Vettore per salvare il primo residuo 
vector <int> counterj; //Vettore per salvare il secondo residuo 
raggioc=18.0; 
int bloc1=0; 
int bloc2=0; //Residui tra i quali mettere di volta in volta la molla più forte

for (int i=0; i<righe; i++) //Individuo le coppie di residui non consecutive che si trovano ad una distanza interiore del raggio di cutoff 
{
	for (int j=0; j<righe; j++)
	{
		R0=(res[i].x-res[j].x)*(res[i].x-res[j].x)+(res[i].y-res[j].y)*(res[i].y-res[j].y)+(res[i].z-res[j].z)*(res[i].z-res[j].z);
		distance=sqrt(R0);

		var=i-j;
		if(j>i+1 && distance <= raggioc) //Faccio in modo che le coppie non vengano contate due volte
			{
				counteri.push_back(i);
				counterj.push_back(j);
			}

	}
}

ofstream fileout;
fileout.open("avl.txt");

if (fileout.is_open())
{
for (int p=0; p<(counteri.size()); p++ ) //Mettendo, una alla volta, la molla rafforzata tra ognuna delle coppie di reidui individuata ricalcolo gli autovalori.
{
bloc1=counteri.at(p);
bloc2=counterj.at(p);
for (int i=0; i<righe; i++)
{
	for (int j=0; j<righe; j++)
	{
		R0=(res[i].x-res[j].x)*(res[i].x-res[j].x)+(res[i].y-res[j].y)*(res[i].y-res[j].y)+(res[i].z-res[j].z)*(res[i].z-res[j].z);
		distance=sqrt(R0);

		var=i-j;
		if((i==bloc1&&j==bloc2)||(i==bloc2&&j==bloc1)) //Chiude la coppia di residui scelti con una molla più forte
		{
			gamma=-10.0L;
		}
		else
		{
			if(distance <= raggioc)
			{
				if((var==-1)||(var==1)) {gamma=-10.0L;}
				else {gamma=-1.0L;}
			}
			else {gamma=0.0L;} 
		}

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

s.compute(W);

cout << bloc1 << " " << bloc2 << endl;

fileout << "(" << bloc1 << "-" << bloc2 << ") " << s.eigenvalues()[8] << " " << s.eigenvalues()[9] << endl;
}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}
return 0;
}

