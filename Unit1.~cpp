//---------------------------------------------------------------------------
#include <vcl\vcl.h>
#pragma hdrstop
#include <math.h>
#include <strstrea.h>
#include <fstream.h>
#include <iostream.h>
#include <fstream.h>
#include <conio.h>
#include "Unit1.h"
#include <complex.h>
//---------------------------------------------------------------------------
#pragma resource "*.dfm"
TForm1 *Form1;
int const
// Number of atoms
Nat = 1000, Nat2 = Nat/2;

double const
pi=3.1415926535897932385,
// Time step
dt=0.005, dt2=dt*dt,
// Discreteness parameter
alpha=-1.0/6.0,
beta=1.0/120.0,
ADB=1.0,
teta=0.257,
omegaDB=2.250419120,
delta=0.3,

// Calculate until time tMAX
tMAX =9900250.0,

// Parameters of Shtormer method
k1=7.0/6.0, k2=5.0/12.0, k3=1.0/3.0, k4=1.0/12.0;

int X1,X2,X3,found;
double e1do,e2do,e3do,e1po,e2po,e3po,TotEn,Kin,Pot,omega,tmax,tmaxP,xmax,xmaxP,
TTpast,DBT,DBomega,DBampl,
P1i,P2i,P1f,P2f,E1i,E2i,E1f,E2f,Ap,Am,Af,Adot,MaxA,Xp,Xm,Xf,Xdot,Rp,Rm,Rf,Rdot,Sp,Sm,Sf,Sdot,Alpha4,XX1,YY1,YY2,V,Amplitudek,MaxKinEner,MinKinEner=100.0;
// Arrays
double *u0,*u1,*u2,*u3,*p,*f0,*f1,*f2,*f3,*analit,*EnerOfAtoms;
double PartParam[5000][3],PartParam0[5000][3];
double KinEner[Nat], PotEner[Nat], TotEner[Nat];

double sqr (double x) {return x*x;}
double cube (double x) {return x*x*x;}
double DET (double a, double b, double c, double d, double e, double f, double g, double h, double i)
           {return a*(e*i-h*f)-d*(b*i-h*c)+g*(b*f-e*c);}
double Ch(double x) {return ( exp(x)+exp(-x) )/2.0;}
double Sh(double x) {return ( exp(x)-exp(-x) )/2.0;}
double Th(double x) {return ( exp(x)-exp(-x) )/( exp(x)+exp(-x) );}
double KinkIM(double x, double t);
double GAMMA(double x);
//double Solution(double x, double t) {return bkRunf(x,t);}
double Solution(double x, double t) {return KinkIM(x,t);}
double PotEnergy (double x[]);
double KinEnergy (double x0[],double x1[],
                  double x2[],double x3[]);
void EnergyOfAtoms (double x0[],double x1[],
                  double x2[],double x3[]);
double EnerFromTo (double x0[],double x1[],
                   double x2[],double x3[],
                   int From, int To);
double MomtFromTo (double x0[],double x1[],
                   double x2[],double x3[],
                   int From, int To);
void ForceCulc (double x[], double y[]);
void PrintOutData();
void ShowMesh();
void ShowAtoms(int, double z[]);
void ShowEnergy(int, double z[]);
void Add_HIST_E(double Al, double Ar);
void Add_HIST_P(double Al, double Ar);
double MR11(double r);

//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner){
u0 = new double [Nat*8]; u1 = new double [Nat*8];
u2 = new double [Nat*8]; u3 = new double [Nat*8];
f0 = new double [Nat*8]; f1 = new double [Nat*8];
f2 = new double [Nat*8]; f3 = new double [Nat*8];
 p = new double [Nat*8]; analit = new double [Nat*8];
 EnerOfAtoms = new double [Nat*8];
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button1Click(TObject *Sender){

ofstream OFSOee("Energies.txt",ios::out);
OFSOee.precision(10);
ofstream OFSOdb("Number_of_DB.txt",ios::out);
OFSOdb.precision(10);
ofstream OFSOdbL("DB_List.txt",ios::out);
OFSOdbL.precision(10);
ofstream OFSOed("Energy_distr.txt",ios::out);
OFSOed.precision(10);

double omegaMAX=2.2360679775;
double omegaMIN=1.0;

// Zero into files
for(int i = 0; i<Nat; i++){
   u0[i]=0.0; u1[i]=0.0; u2[i]=0.0; u3[i]=0.0;
   f0[i]=0.0; f1[i]=0.0; f2[i]=0.0; f3[i]=0.0;
    p[i]=0.0;
}

//**********************************************
//  INITIAL  CONDITIONS
//**********************************************
// Gamma-point mode
double AmplZBM=0.3;
for (int i=0; i<Nat; i++){
   int Rand=random(2000);
   double R=(Rand-1000)/1000.0;
   u0[i]=AmplZBM+R*1.0e-14;
   u1[i]=AmplZBM+R*1.0e-14;
   u2[i]=AmplZBM+R*1.0e-14;
   u3[i]=AmplZBM+R*1.0e-14;
}

//**********************************************
//  MD  CALCULATION
//**********************************************

double t = 0.0;
long nit = 0;
int RandAv=500;
double AverKin=0.0;
int NumAver=0;

while (t<=tMAX+0.0000001){
    //  Stormer procedure
    ForceCulc(u0,f0);
    for(int i=0; i<Nat; i++){
       p[i] = 2.0*u0[i]-u1[i]
            + dt2*(k1*f0[i]-k2*f1[i]+k3*f2[i]-k4*f3[i]);
       }

    // Replace
    for(int i=0; i<Nat; i++){
       u3[i] = u2[i]; u2[i] = u1[i]; u1[i] = u0[i]; u0[i] = p[i];
       f3[i] = f2[i]; f2[i] = f1[i]; f1[i] = f0[i];
    }

    // Show solution
    if (nit==500000){
        Kin=KinEnergy(u0,u1,u2,u3);
        Pot=PotEnergy(u0);
        TotEn = Kin+Pot;
        EnergyOfAtoms(u0,u1,u2,u3);
        ShowMesh();
        ShowAtoms(clRed,u0);
        ShowEnergy(clBlue,EnerOfAtoms);
        // Localization parameter
        double Sum1, Sum2, Localization, MinE, MaxE;
        Sum1=0.0;
        Sum2=0.0;
        Localization=0.0;
        MinE=100000.0;
        MaxE=-100000.0;
        for(int i=0; i<Nat; i++){
           Sum1=Sum1+sqr(EnerOfAtoms[i]);
           Sum2=Sum2+EnerOfAtoms[i];
           if (EnerOfAtoms[i]>MaxE) MaxE=EnerOfAtoms[i];
           if (EnerOfAtoms[i]<MinE) MinE=EnerOfAtoms[i];
        }
        Localization=(Sum1/(1.0*Nat))/(sqr(Sum2)/(1.0*Nat));
        Form1->Canvas->Brush->Color = clBtnFace;
        Form1->Canvas->TextOut(25,5,"Kin   =");
        Form1->Canvas->TextOut(95,5,Kin);
        Form1->Canvas->TextOut(25,20,"Pot   =");
        Form1->Canvas->TextOut(95,20,Pot);
        Form1->Canvas->TextOut(25,35,"TotEn   =");
        Form1->Canvas->TextOut(95,35,TotEn);
        Form1->Canvas->TextOut(25,50,"Loc=                                                       ");
        Form1->Canvas->TextOut(95,50,Localization);
        Form1->Canvas->TextOut(25,65,"DelE=                                                      ");
        Form1->Canvas->TextOut(95,65,MaxE-MinE);
        Form1->Canvas->TextOut(25,80,"t   =");
        Form1->Canvas->TextOut(95,80,t);
        AverKin=AverKin/(1.0*NumAver);
        OFSOee.width(12); OFSOee << t;       OFSOee << " ";  // 1
        OFSOee.width(12); OFSOee << Pot;     OFSOee << " ";  // 2
        OFSOee.width(12); OFSOee << Kin;     OFSOee << " ";  // 3
        OFSOee.width(12); OFSOee << AverKin; OFSOee << " ";  // 4
        OFSOee.width(12); OFSOee << TotEn;   OFSOee << " ";  // 5
        OFSOee.width(12); OFSOee << Localization; OFSOee << " ";  // 6
        OFSOee.width(12); OFSOee << MaxE-MinE; OFSOee << " ";  // 6
        OFSOee << "\n";
        AverKin=0.0;
        NumAver=0;
        nit=0;
        // Count DB
        Kin=KinEnergy(u0,u1,u2,u3);
        Pot=PotEnergy(u0);
        TotEn = (Kin+Pot)/(1.0*Nat);
        EnergyOfAtoms(u0,u1,u2,u3);
        int DB4=0,DB8=0,DB16=0,DB32=0,DB64=0,DB128=0;
        for(int i=1; i<Nat-1; i++){
           if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
              if ((EnerOfAtoms[i]>4.0*TotEn)&&(EnerOfAtoms[i]<8.0*TotEn)) DB4+=1;
           }
           if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
              if ((EnerOfAtoms[i]>8.0*TotEn)&&(EnerOfAtoms[i]<16.0*TotEn)) DB8+=1;
           }
           if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
              if ((EnerOfAtoms[i]>16.0*TotEn)&&(EnerOfAtoms[i]<32.0*TotEn)) DB16+=1;
           }
           if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
              if ((EnerOfAtoms[i]>32.0*TotEn)&&(EnerOfAtoms[i]<64.0*TotEn)) DB32+=1;
           }
           if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
              if ((EnerOfAtoms[i]>64.0*TotEn)&&(EnerOfAtoms[i]<128.0*TotEn)) DB64+=1;
           }
           if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
              if (EnerOfAtoms[i]>128.0*TotEn) DB128+=1;
           }
        }// for(int i=1; i<Nat-1; i++)
        OFSOdb.width(12); OFSOdb << t;              OFSOdb << " ";  // 1
        OFSOdb.width(12); OFSOdb << DB4/(1.0*Nat);  OFSOdb << " ";  // 2
        OFSOdb.width(12); OFSOdb << DB8/(1.0*Nat);  OFSOdb << " ";  // 3
        OFSOdb.width(12); OFSOdb << DB16/(1.0*Nat); OFSOdb << " ";  // 4
        OFSOdb.width(12); OFSOdb << DB32/(1.0*Nat); OFSOdb << " ";  // 5
        OFSOdb.width(12); OFSOdb << DB64/(1.0*Nat); OFSOdb << " ";  // 6
        OFSOdb.width(12); OFSOdb << DB128/(1.0*Nat); OFSOdb << "\n";// 7
    }// show solution

    // Average KinEn
    if (fmod(nit,500+RandAv) == 0){
        Kin=KinEnergy(u0,u1,u2,u3);
        AverKin=AverKin+Kin;
        NumAver=NumAver+1;
        RandAv=random(1000);
    }// for Screen view
                         
    t += dt;
    nit +=1;

}//while (t<=tMAX+0.0000001)

// Count DB
Kin=KinEnergy(u0,u1,u2,u3);
Pot=PotEnergy(u0);
TotEn = (Kin+Pot)/(1.0*Nat);
EnergyOfAtoms(u0,u1,u2,u3);
for(int i=0; i<Nat; i++){
        OFSOed.width(12); OFSOed << i; OFSOed << " ";               // 1
        OFSOed.width(12); OFSOed << EnerOfAtoms[i]; OFSOed << "\n"; // 2
}

int DBN=0;
for(int i=1; i<Nat-1; i++){
        if ((EnerOfAtoms[i]>EnerOfAtoms[i-1])&&(EnerOfAtoms[i]>EnerOfAtoms[i+1])){
                if (EnerOfAtoms[i]>4.0*TotEn){
                        DBN+=1;
                 // count DB energy
                        double DBen=EnerOfAtoms[i];
                        int k=1;
                        while (EnerOfAtoms[i+k-1]>EnerOfAtoms[i+k]){
                                DBen+=EnerOfAtoms[i+k];
                                k+=1;
                        }
                        k=-1;
                        while (EnerOfAtoms[i+k+1]>EnerOfAtoms[i+k]){
                                DBen+=EnerOfAtoms[i+k];
                                k-=1;
                        }
                        OFSOdbL.width(12); OFSOdbL << DBN;  OFSOdbL << " ";  // 1
                        OFSOdbL.width(12); OFSOdbL << i;    OFSOdbL << " ";  // 2
                        OFSOdbL.width(12); OFSOdbL << DBen; OFSOdbL << "\n"; // 3
                }
        }
}// for(int i=1; i<Nat-1; i++)


//Cloze the file
OFSOee.close();
}
//---------------------------------------------------------------------------

//************************************
//************************************
//************************************
//************************************
//      PROCEDURES
//************************************
//************************************
//************************************
//************************************


// Calculation of forces with forth derivative
void ForceCulc (double u[], double f[]){
  for(int i=1; i<Nat-1; i++){
      double u2=u[i]*u[i];
      f[i] = (u[i-1]-2.0*u[i]+u[i+1]) - u[i]*(1.0+alpha*u2+beta*u2*u2);
  }
  f[0] = (u[Nat-1]-2.0*u[0]+u[1]) - u[0]*(1.0+alpha*u[0]*u[0]+beta*u[0]*u[0]*u[0]*u[0]);
  f[Nat-1] = (u[Nat-2]-2.0*u[Nat-1]+u[0]) - u[Nat-1]*(1.0+alpha*u[Nat-1]*u[Nat-1]+beta*u[Nat-1]*u[Nat-1]*u[Nat-1]*u[Nat-1]);
}


// Potential energy of the chain
double PotEnergy (double x[])
{double Pot = 0.0;
for(int i=0; i<Nat-1; i++)
   { Pot += 0.5*sqr(x[i+1]-x[i])+0.5*sqr(x[i])+0.25*alpha*sqr(x[i]*x[i])+(1.0/6.0)*beta*sqr(x[i]*x[i]*x[i]); }
Pot += 0.5*sqr(x[Nat-1]-x[0])+0.5*sqr(x[Nat-1])+0.25*alpha*sqr(x[Nat-1]*x[Nat-1])+(1.0/6.0)*beta*sqr(x[Nat-1]*x[Nat-1]*x[Nat-1]);
return Pot;}

// Kinetic energy of the chain
double KinEnergy (double x0[],double x1[],double x2[],double x3[])
{double Kin = 0.0;
for(int i=0; i<Nat; i++)
  { //Kin += 0.5*sqr((2.0*x0[i]+3.0*x1[i]-6.0*x2[i]+x3[i])/(6.0*dt));
  Kin += 0.5*sqr((11.0*x0[i]-18.0*x1[i]+9.0*x2[i]-2.0*x3[i])/(6.0*dt));
  }
return Kin;}


// Energy of atoms
void EnergyOfAtoms (double x0[],double x1[],double x2[],double x3[]){
  double Kin = 0.0;
  double Pot = 0.0;
  for(int i=1; i<Nat-1; i++){
    Kin = 0.5*sqr((11.0*x0[i]-18.0*x1[i]+9.0*x2[i]-2.0*x3[i])/(6.0*dt));
    Pot = 0.25*sqr(x0[i]-x0[i-1])+0.25*sqr(x0[i+1]-x0[i])+0.5*sqr(x0[i])+0.25*alpha*sqr(x0[i]*x0[i])+(1.0/6.0)*beta*sqr(x0[i]*x0[i]*x0[i]);
    EnerOfAtoms[i]=Kin+Pot;
  }
  EnerOfAtoms[0]=0.5*sqr((11.0*x0[0]-18.0*x1[0]+9.0*x2[0]-2.0*x3[0])/(6.0*dt))
  +0.25*sqr(x0[0]-x0[Nat-1])+0.25*sqr(x0[1]-x0[0])+0.5*sqr(x0[0])+0.25*alpha*sqr(x0[0]*x0[0])+(1.0/6.0)*beta*sqr(x0[0]*x0[0]*x0[0]);
  EnerOfAtoms[Nat-1]=0.5*sqr((11.0*x0[Nat-1]-18.0*x1[Nat-1]+9.0*x2[Nat-1]-2.0*x3[Nat-1])/(6.0*dt))
  +0.25*sqr(x0[Nat-1]-x0[Nat-2])+0.25*sqr(x0[0]-x0[Nat-1])+0.5*sqr(x0[Nat-1])+0.25*alpha*sqr(x0[Nat-1]*x0[Nat-1])+(1.0/6.0)*beta*sqr(x0[Nat-1]*x0[Nat-1]*x0[Nat-1]);
}


/*
// PrintOut the numerical data
void PrintOutData()
{Form1->Canvas->Brush->Color = clBtnFace;
Form1->Canvas->TextOut(25,80,"Time t =");
Form1->Canvas->TextOut(80,80,t+t0);
Form1->Canvas->TextOut(25,95,"Potential energy per atom =");
Form1->Canvas->TextOut(170,95,AverPot);
Form1->Canvas->TextOut(25,110,"Kinetic   energy per atom =");
Form1->Canvas->TextOut(170,110,AverKin);
Form1->Canvas->TextOut(25,125,"Total     energy per atom =");
Form1->Canvas->TextOut(170,125,AverEn);
Form1->Canvas->TextOut(25,140,"Lost      energy per atom =");
Form1->Canvas->TextOut(170,140,AverEn0-AverEn);}
*/

// Show the frame for the picture of atoms
void ShowMesh()
{Form1->Canvas->Brush->Color = clWhite;
Form1->Canvas->Rectangle(1,270,670,630);
Form1->Canvas->MoveTo(40,300);
Form1->Canvas->LineTo(40,600);
Form1->Canvas->LineTo(640,600);
Form1->Canvas->LineTo(640,300);
Form1->Canvas->LineTo(40,300);
Form1->Canvas->MoveTo(40,375);
Form1->Canvas->LineTo(640,375);
Form1->Canvas->MoveTo(40,450);
Form1->Canvas->LineTo(640,450);
Form1->Canvas->MoveTo(40,525);
Form1->Canvas->LineTo(640,525);
Form1->Canvas->MoveTo(90,300);
Form1->Canvas->LineTo(90,600);
Form1->Canvas->MoveTo(140,300);
Form1->Canvas->LineTo(140,600);
Form1->Canvas->MoveTo(190,300);
Form1->Canvas->LineTo(190,600);
Form1->Canvas->MoveTo(240,300);
Form1->Canvas->LineTo(240,600);
Form1->Canvas->MoveTo(290,300);
Form1->Canvas->LineTo(290,600);
Form1->Canvas->MoveTo(340,300);
Form1->Canvas->LineTo(340,600);
Form1->Canvas->MoveTo(390,300);
Form1->Canvas->LineTo(390,600);
Form1->Canvas->MoveTo(440,300);
Form1->Canvas->LineTo(440,600);
Form1->Canvas->MoveTo(490,300);
Form1->Canvas->LineTo(490,600);
Form1->Canvas->MoveTo(540,300);
Form1->Canvas->LineTo(540,600);
Form1->Canvas->MoveTo(590,300);
Form1->Canvas->LineTo(590,600);
Form1->Canvas->TextOut(40,608,"1");
Form1->Canvas->TextOut(630,608,Nat);
Form1->Canvas->TextOut(330,608,Nat2);
Form1->Canvas->TextOut(10,298,"+12.6");
Form1->Canvas->TextOut(10,373,"+6.28");
Form1->Canvas->TextOut(10,448," 0.00");
Form1->Canvas->TextOut(10,523,"-6.28");
Form1->Canvas->TextOut(10,598,"-12.6");
}

// Show displacements of atoms
void ShowAtoms(int Col, double x0[]){
  int x,y;
  for (int i=0; i<Nat; i++){
      x = 40 + floor(i*600.0/Nat);
      y = 450 - floor(50.0*pi*x0[i]*1.0/(4.0*pi));
      Form1->Canvas->Pixels[x][y] = Col;
      Form1->Canvas->Pixels[x+1][y] = Col;
      Form1->Canvas->Pixels[x][y+1] = Col;
      Form1->Canvas->Pixels[x+1][y+1] = Col;
  }
}

// Show energy of atoms
void ShowEnergy(int Col, double x0[]){
  int x,y;
  for (int i=0; i<Nat; i++){
      x = 40 + floor(i*600.0/Nat);
      y = 450 - floor(1.5*pi*x0[i]*150.0/(4.0*pi));
      Form1->Canvas->Pixels[x][y] = Col;
      Form1->Canvas->Pixels[x+1][y] = Col;
      Form1->Canvas->Pixels[x][y+1] = Col;
      Form1->Canvas->Pixels[x+1][y+1] = Col;
  }
}





