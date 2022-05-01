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
vector <double> bfact(righe);

while (i<righe)
	{
		res[i].x=v.at(4*i);
		res[i].y=v.at(4*i+1);
		res[i].z=v.at(4*i+2);
		bfact.at(i)=v.at(4*i+3);
		i++;
	}

cout << "Ho letto il file di input" << endl;

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
		if(distance <= raggioc)
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

cout << "Ho calcolato la matrice Hessiana con R_C=" << raggioc << " A e l'ho stampata nel file matrice.txt" << endl;

//Calcolo l'area
double area1=0.0L;

for (int i=0; i<righe; i++)
	{
		area1+=bfact.at(i);
	}

//Scrivo il file con cui fare il grafico del B-factor
ofstream fileout2;
fileout2.open("bfactor.txt");

if (fileout2.is_open())
{
	for (int i=0; i<righe; i++)
	{
		fileout2 << i << " " << bfact.at(i) << endl;
	}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

cout << "Ho scritto il file bfact.txt" << endl;

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

cout << "Ho calcolato autovettori e autovalori , quest'ultimi stampati nel file autovalori.txt" << endl;
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

cout << "Ho calcolato il profilo di mobilita' e l'ho stampato nel file deltaRi.txt" << endl;

//Scrivo un file per visualizzare i profili di mobilità separatamente per ogni modo
/*ofstream fileout33;
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
}*/

//Calcolo il Kk
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

cout << "Ho calcolato i gradi di collettività e li ho stampati nel file kappa.txt" << endl;

//Calcolo il linking number
int i1=0;
int i2=0;
int j1a=0;
int j2a=0;
int j1b=0;
int j2b=0;
int j1c=0;
int j2c=0;

i1=10;
i2=39;
j1a=40;
j2a=78;
j1b=42;
j2b=80;
j1c=42;
j2c=75;

cout << "Ho calcolato il linking number per la configurazione base per i tre loop: Ga=" << linknumb(righe, i1, i2, j1a, j2a, res) << ", Gb=" << linknumb(righe, i1, i2, j1b, j2b, res) << ", Gc=" << linknumb(righe, i1, i2, j1c, j2c, res) << endl;

//Genero le configurazioni deformate di prova
vector <pos> def(righe);
double smax=0.0L;
smax=60;
int modo=0;
modo=9;

for (int i=0; i<righe; i++)
	{
		def[i].x=res[i].x+(smax)*(s.eigenvectors().col(modo)[3*i])/(sqrt(s.eigenvalues()[modo]));	
		def[i].y=res[i].y+(smax)*(s.eigenvectors().col(modo)[3*i+1])/(sqrt(s.eigenvalues()[modo]));
		def[i].z=res[i].z+(smax)*(s.eigenvectors().col(modo)[3*i+2])/(sqrt(s.eigenvalues()[modo]));
	}

pos cm;
pos cmk;
cm=centrmass(righe, res);
cmk=centrmass(righe, def);

double dist=0.0L;
double distk=0.0L;

dist=sqrt(pow(cm.x-res[50].x,2)+pow(cm.y-res[50].y,2)+pow(cm.z-res[50].z,2));
distk=sqrt(pow(cmk.x-def[50].x,2)+pow(cmk.y-def[50].y,2)+pow(cmk.z-def[50].z,2));

cout << "Ho calcolato la deformazione di prova con s_max=" << smax << " per il modo k=" << modo << " ottenendo una differenza fra la distanza cm-res50 e cm-res50 deformati del " << 100*(abs(distk-dist))/dist << "% della distanza cm-res50" << endl;

//Faccio variare s nell'intervallo [-smax;smax] e calcolo i linking number per ogni s per i tre loop
vector <double> ss(200);
vector <double> Ga(200);
vector <double> Gb(200);
vector <double> Gc(200);
vector <double> dist1039(200);
vector <double> dist4280(200);
vector <double> dist4275(200);

double variazione=0.0L;
variazione=2*smax/200;

for (int i=0; i<200; i++)
{
ss.at(i)=variazione*i-smax;
}

for (int l=0; l<200; l++)
{
	for (int i=0; i<righe; i++)
	{
		def[i].x=res[i].x+(ss.at(l))*(s.eigenvectors().col(modo)[3*i])/(sqrt(s.eigenvalues()[modo]));	
		def[i].y=res[i].y+(ss.at(l))*(s.eigenvectors().col(modo)[3*i+1])/(sqrt(s.eigenvalues()[modo]));
		def[i].z=res[i].z+(ss.at(l))*(s.eigenvectors().col(modo)[3*i+2])/(sqrt(s.eigenvalues()[modo]));
	}
Ga.at(l)=linknumb(righe, i1, i2, j1a, j2a, def);
Gb.at(l)=linknumb(righe, i1, i2, j1b, j2b, def);
Gc.at(l)=linknumb(righe, i1, i2, j1c, j2c, def);
dist1039.at(l)=sqrt(pow(def[10].x-def[39].x,2)+pow(def[10].y-def[39].y,2)+pow(def[10].z-def[39].z,2));
dist4280.at(l)=sqrt(pow(def[42].x-def[80].x,2)+pow(def[42].y-def[80].y,2)+pow(def[42].z-def[80].z,2));
dist4275.at(l)=sqrt(pow(def[42].x-def[75].x,2)+pow(def[42].y-def[75].y,2)+pow(def[42].z-def[75].z,2));
}

//Stampo il file dei G per il grafico 
ofstream fileout999;
fileout999.open("gvars.txt");

if (fileout999.is_open())
{
	for (int i=0; i<(200); i++)
	{
		fileout999 << ss.at(i) << " " << Ga.at(i) << " " << Gb.at(i) << " " << Gc.at(i) << " " << dist1039.at(i) << " " << dist4280.at(i) << " " << dist4275.at(i)<< endl;
	}
}
else 
{
	cout << "Errore nella scrittura del file." << endl;
}

cout << "Ho calcolato i valori di G al variare di s e li ho stampati ne file gvars.txt" << endl;

return 0;
}





















//========================================================================================================================================
//Funzione che calcola il linking number
double linknumb (int righe, int i1, int i2, int j1, int j2, vector <pos> res) 
{
vector <pos> Ri(righe);
vector <pos> dRi(righe);
vector <pos> Rj(righe);
vector <pos> dRj(righe);

for (int i=i1; i<(i2+1); i++)
{
	Ri[i].x=((res[i].x)+(res[i+1].x))/2;
	Ri[i].y=((res[i].y)+(res[i+1].y))/2;
	Ri[i].z=((res[i].z)+(res[i+1].z))/2;
}

for (int i=i1; i<(i2+1); i++)
{
	dRi[i].x=(res[i+1].x)-(res[i].x);
	dRi[i].y=(res[i+1].y)-(res[i].y);
	dRi[i].z=(res[i+1].z)-(res[i].z);
}

for (int i=j1; i<(j2+1); i++)
{
	Rj[i].x=((res[i].x)+(res[i+1].x))/2;
	Rj[i].y=((res[i].y)+(res[i+1].y))/2;
	Rj[i].z=((res[i].z)+(res[i+1].z))/2;
}

for (int i=j1; i<(j2+1); i++)
{
	dRj[i].x=(res[i+1].x)-(res[i].x);
	dRj[i].y=(res[i+1].y)-(res[i].y);
	dRj[i].z=(res[i+1].z)-(res[i].z);
}

double Gij=0.0L;
double modulo=0.0L;

for (int i=i1; i<i2; i++)
{
	for (int j=j1; j<j2; j++)
	{
		modulo=sqrt(pow((Ri[i].x-Rj[j].x),2)+pow((Ri[i].y-Rj[j].y),2)+pow((Ri[i].z-Rj[j].z),2));
		Gij+=((Ri[i].x-Rj[j].x)*(dRi[i].y*dRj[j].z-dRi[i].z*dRj[j].y) + (Ri[i].y-Rj[j].y)*(dRi[i].z*dRj[j].x-dRi[i].x*dRj[j].z) + (Ri[i].z-Rj[j].z)*(dRi[i].x*dRj[j].y-dRi[i].y*dRj[j].x))/(pow(modulo,3));
	}
}

return Gij/(4*M_PI);
}

//Funzione che calcola le coordinate del centro di massa
pos centrmass (int righe, vector <pos> res)
{
pos R;
R.x=0.0L;
R.y=0.0L;
R.z=0.0L;

for (int i=0; i<righe; i++)
{
	R.x+=res[i].x;
	R.y+=res[i].y;
	R.z+=res[i].z;
}

R.x=R.x/righe;
R.y=R.y/righe;
R.z=R.z/righe;

return R;
}
