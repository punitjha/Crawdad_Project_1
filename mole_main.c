#include "molecule.h"
#include "mass.h"
using namespace std;
/**************************************************************/	
int main()
{
	molecule acet("acetal.txt", 0);
	cout<<"We now print the geometry"<<endl;
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
			//cout<<R_bonds[i][j]<<"   ";
		}
	}
	cout<<endl<<"Printing the bond distances"<<endl;
	acet.print(R_bonds, acet.natom,acet.natom);
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
	//cout<<endl<< " We now print the x-components "<<endl;
	//acet.print(ex,arlen,arlen);

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
	cout<<endl<<"Printing the bond angles"<<endl;
	for(int i=0; i< acet.natom; i++)
	{
		for (int j=0; j<i; j++)
		{
			for(int k=0; k<j; k++)
			{
				if (R_bonds[i][j]<4.0 && R_bonds[j][k] < 4.0)
				{
				     phi[i][j][k]=acet.angle(i,j,k);
					printf("%2d %2d %2d %10.4f\n",i,j,k,phi[i][j][k]*(180.0/acos(-1.0)));
				}
			}
		}
}
//now we calculate the out of plane angles
	cout<<endl<<"Out of plane angles"<<endl;
	for (int i=0; i<acet.natom; i++)
	{
		for (int k=0; k<acet.natom; k++)
		{
			for(int j=0; j<acet.natom; j++)
			{
				for(int l=0; l<j; l++)
				{
					if (R_bonds[i][k] < 4.0 && R_bonds[j][k]< 4.0 && R_bonds [k][l]< 4.0 
							&& i!=j && i!=k && i!= l && j!=k && k!=l ) 
					{ 
						printf(" %2d %2d %2d %2d %10.5f\n", i,j,k,l,acet.outplane(i,j,k,l)*(180.0/acos(-1.0)));
					}
				  }
			 }
		}
	}
//now we calculate the torsion/dihedral angles
	cout<<endl<<"The torsion angles"<<endl;
	for (int i=0; i<acet.natom; i++)
	{
		for (int j=0; j<i; j++)
		{
			for(int k=0; k<j; k++)
			{
				for(int l=0; l<k; l++)
				{
					if (R_bonds[i][j] < 4.0 && R_bonds[j][k]< 4.0 && R_bonds [k][l]< 4.0 ) 
					{ 
						printf(" %2d %2d %2d %2d %10.5f\n", i,j,k,l,acet.torsion(i,j,k,l)*(180.0/acos(-1.0)));
					}
				  }
			 }
		}
	}

//center of mass translation
	double x_cm=0.0;
	double y_cm=0.0;
	double z_cm=0.0;
	double sum=0.0;
	for (int i=0; i<acet.natom; i++)
	{
		x_cm=x_cm+mass[acet.zvals[i]]*acet.geom[i][0];
		y_cm=y_cm+mass[acet.zvals[i]]*acet.geom[i][1];
		z_cm=z_cm+mass[acet.zvals[i]]*acet.geom[i][2];
		sum=sum+mass[acet.zvals[i]];
	}
	cout<<endl<<"X_COM"<<' '<<x_cm/sum<<endl;
	cout<<"Y_COM"<<' '<<y_cm/sum<<endl;
	cout<<"Z_COM"<<' '<<z_cm/sum<<endl;

//translating the molecule to the centre of mass
	acet.translate(-x_cm/sum, -y_cm/sum, -z_cm/sum);
	cout<<endl<<"Translating it to the centre of mass"<<endl;
	acet.print_geom();
//Calculating the principle moments of inertia
	arma::mat Inrt= arma::zeros(3,3);
	for (int i=0; i<acet.natom; i++)
	{
		Inrt(0,0)+= mass[acet.zvals[i]]*(acet.geom[i][1]*acet.geom[i][1]+acet.geom[i][2]*acet.geom[i][2]);
		Inrt(1,1)+= mass[acet.zvals[i]]*(acet.geom[i][0]*acet.geom[i][0]+acet.geom[i][2]*acet.geom[i][2]);
		Inrt(2,2)+= mass[acet.zvals[i]]*(acet.geom[i][0]*acet.geom[i][0]+acet.geom[i][1]*acet.geom[i][1]);
		Inrt(0,1)+= mass[acet.zvals[i]]*acet.geom[i][0]*acet.geom[i][1];
		Inrt(0,2)+= mass[acet.zvals[i]]*acet.geom[i][0]*acet.geom[i][2];
		Inrt(1,2)+= mass[acet.zvals[i]]*acet.geom[i][1]*acet.geom[i][2];

	}
	Inrt(1,0)=Inrt(0,1);
	Inrt(2,0)=Inrt(0,2);
	Inrt(2,1)=Inrt(1,2);
	cout<<endl;
	Inrt.print("The moment of inertia matrix/tensors");
	cout<<endl;
	arma::vec eg;
	arma::mat ev;
	arma::eig_sym(eg,ev,Inrt);
	cout<<endl;
	eg.print("These are the eigenvalues in amu bohr^2");
	double conversion=0.529177*0.529177;
	cout<<endl<<"These are the eigenvalues in amu AA"<<endl<<eg*conversion<<endl;
	conversion=1.6605402E-24*0.529177249E-8*0.529177249E-8;
	cout<<endl<<"These are the eigenvalues in gg  cm^2"<<endl<<eg*conversion<<endl;
//classifying the rotor
	if (acet.natom == 2)
		cout<<"The molecule is diatomic"<<endl;
	else if (eg(0)<1e-4)
		cout<<"The molecule is linear"<<endl;
	else if ((fabs(eg(0)-eg(1))< 1e-4) && (fabs(eg(1)-eg(2))< 1e-4))
		cout<<"The molecule is a spherical top"<<endl;
	else if ((fabs(eg(0)-eg(1)) < 1e-4) && (fabs(eg(1)-eg(2))> 1e-4))
		cout<<"The molecule is a oblate symmetric top"<<endl;
	else if ((fabs(eg(0)-eg(1))> 1e-4 )&& (fabs(eg(1)-eg(2))< 1e-4))
		cout<<"The molecule is a prolate symmetric top"<<endl;
	else cout<<"The molecule is a asymmetric top"<<endl;	
//calculating the rotational constants of the molecules
	double plank=6.6260755e-27;
	double pi=acos(-1.0);
	double sp=3.0e10;
	double A=plank/(8.0*pi*pi*sp*eg(0)*1.6605402E-24*0.529177249E-8*0.529177249E-8);
	double B=plank/(8.0*pi*pi*sp*eg(1)*1.6605402E-24*0.529177249E-8*0.529177249E-8);
	double C=plank/(8.0*pi*pi*sp*eg(2)*1.6605402E-24*0.529177249E-8*0.529177249E-8);
	cout<<"\nRotational Constants in cm -1:\n" <<"A="<<A<<"\tB="<<B<<"\tC="<<C<<endl;
	double A_1=1.0e-6*plank/(8.0*pi*pi*eg(0)*1.6605402E-24*0.529177249E-8*0.529177249E-8);
	double B_1=1e-6*plank/(8.0*pi*pi*eg(1)*1.6605402E-24*0.529177249E-8*0.529177249E-8);
	double C_1=1e-6*plank/(8.0*pi*pi*eg(2)*1.6605402E-24*0.529177249E-8*0.529177249E-8);
	cout<<"\nRotational Constants in MHz is:\n" <<"A="<<A_1<<"\tB="<<B_1<<"\tC="<<C_1<<endl;
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
//deleting the memory allocated to the phi vector
	for (int i=0; i<acet.natom; i++)
	{
		for(int j=0; j<acet.natom;j++)
		{
			delete[] phi[i][j];	
		}
	delete[] phi[i];
	}
	delete[] phi;
}
