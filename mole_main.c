#include "molecule.h"
using namespace std;

/**************************************************************/	
int main()
{
	molecule acet("acetal.txt", 0);
	acet.print_geom();
	double **R_bonds=new double* [acet.natom];
	for(int i=0; i<acet.natom; i++)
	{
		R_bonds[i]=new double [acet.natom];
	}
	for(int i=0; i<acet.natom; i++)
	{
		for(int j=0; j<acet.natom; j++)
		{
			R_bonds[i][j]=acet.bond(i,j);
			cout<<R_bonds[i][j]<<"   ";
		}
		cout<<endl;
	}
	double **ex=new double* [acet.natom];
	double **ey=new double* [acet.natom];
	double **ez=new double* [acet.natom];
	for(int i=0; i< acet.natom;i++)
	{
		ex[i]=new double[acet.natom];
		ey[i]=new double[acet.natom];
		ez[i]=new double[acet.natom];
	}
	for( int i=0; i<acet.natom; i++)
	{
		for(int j=0; j<i ;j++)
		{
			ex[i][j]=ex[j][i]=-(acet.geom[i][0]-acet.geom[j][0])/R_bonds[i][j];
			ey[i][j]=ey[j][i]=-(acet.geom[i][1]-acet.geom[j][1])/R_bonds[i][j];
			ez[i][j]=ez[j][i]=-(acet.geom[i][2]-acet.geom[j][2])/R_bonds[i][j];
		}
	}
	int arlen=acet.natom;
	cout<<endl;
	acet.print(ex,arlen,arlen);

//three dimensional array to store the bond angles
	double ***phi = new double** [acet.natom];
	for(int i=0; i< acet.natom; i++)
	{
		phi[i]=new double* [acet.natom];
		for(int j=0; j<acet.natom; j++)
		{
			phi[i][j]=new double[acet.natom];
		}
	}
	for(int i=0; i< acet.natom; i++)
	{
		for (int j=0; j<i; j++)
		{
			for(int k=0; k<j; k++)
			{
				if (R_bonds[i][j]<4.0 && R_bonds[j][k] < 4.0)
				{
				     phi[i][j][k]=acet.angle(i,j,k);
					cout<<i<<","<<j<<","<<k<<","<<R_bonds[i][j]<<" "<<R_bonds[j][k]<<" "<<phi[i][j][k]*(180.0/acos(-1.0))<<endl;
				}
			}
		}
	}
		
	



//deleting the memory allocated to the bond vector
	for(int i=0; i < acet.natom; i++) 
		delete[] R_bonds[i];
	delete[] R_bonds;
//deleting the memory allocated to the angle vector
	for(int i=0; i<acet.natom; i++)
	{
		delete[] ex[i];
		delete[] ey[i];
		delete[] ez[i];
	}
	delete[] ex;
	delete[] ey;
	delete[] ez;
}
