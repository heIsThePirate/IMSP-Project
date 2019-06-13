//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <math.h>
#include <strstrea.h>
#include <fstream.h>
#include <iostream.h>
#include <conio.h>
#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

//FILE NAMES
const char
ReadDisplFrom[]    = "Kink_with_IM_h_0_5_Nat_300_E_0_5.txt", // File to read the Initial Conditions
WriteDisplTo[]     = "Kink_b.txt", // File to write initial conditions
WriteEigValuesTo[] = "EigenValues_h_0_20_Nat_400_x0_0_5h.txt"; // File to write initial conditions
ofstream OFSOe(WriteEigValuesTo,ios::out);

int const
// Number of atoms
Nat = 100, Nat2 = Nat/2,
NDF=Nat;

double const
pi=3.1415926535897932385,
// Time step
dt=0.01, dt2=dt*dt,
// Discreteness parameter
AtomMass=1.0,

// Calculate from time tMIN to time tMAX
tMIN = 0.0,
tMAX = 1500.1,
tMAXMD = 50000.0,

// Parameters of Shtormer method
k1=7.0/6.0, k2=5.0/12.0, k3=1.0/3.0, k4=1.0/12.0;

double TotEn,TotEnDo,shift,fric=10.3,h,h2,C,dk,xk,deltak,
Left,Right,t,Kin,Pot;
long nit;
// Arrays
double *u0,*u1,*u2,*u3,*p,*f0,*f1,*f2,*f3,*analit;
double
Are[NDF][NDF],MASS[NDF][NDF],
A[NDF+1][NDF+1],b[NDF+1][NDF+1],
x[NDF+1][NDF+1],d[NDF+1],eigv[NDF+1];

double Bond(double r);      // Inter-particle potential
double DBond(double r);     // Derivative of Inter-particle potential
double DDBond(double r);    // Second derivative of Inter-particle potential
double Subs(double r);      // Substrate potential
double DSubs(double r);     // Derivative of Substrate potential
double DDSubs(double r);    // Second derivative of Substrate potential
void Stiff(double h);
double sqr (double x) {return x*x;}
// FI4
double SpectrumFI4 (double k, double h) {return sqrt(2.0+4.0*sqr(sin(0.5*k))/(h*h));}
double ForceFI4 (double ul, double um, double ur, double h)
       {return (ul-2.0*um+ur)/(h*h)+um-um*um*um;}
double PotEnFI4 (double um, double ur, double h)
       {return 0.5*sqr((ur-um)/h)+0.25*sqr(1.0-um*um);}
double StifLeftFI4 (double ul, double um, double ur, double h) {return -1.0/(h*h);}
double StifMiddleFI4 (double ul, double um, double ur, double h) {return 2.0/(h*h)-1.0+3.0*sqr(um);}
double StifRightFI4 (double ul, double um, double ur, double h) {return -1.0/(h*h);}
// FI4 Speight
double SpectrumFI4S (double k, double h) {return sqrt(2.0+(4.0-2.0*h*h)*sqr(sin(0.5*k))/(h*h));}
double ForceFI4S (double ul, double um, double ur, double h)
       {return (ul-2.0*um+ur)/(h*h)
                 +(1.0/6.0)*((ur+2.0*um)*(1.0-(ur*ur+ur*um+um*um)/3.0)
                           + (ul+2.0*um)*(1.0-(ul*ul+ul*um+um*um)/3.0));}
double PotEnFI4S (double um, double ur, double h)
       {return 0.5*sqr((ur-um)/h)+0.25*sqr(1.0-sqr(0.5*ur+0.5*um));}
double StifLeftFI4S (double ul, double um, double ur, double h)
       {return -1.0/(h*h)-(1.0-sqr(um+ul))/6.0;}
double StifMiddleFI4S (double ul, double um, double ur, double h)
       {return 2.0/(h*h)-(4.0-2.0*sqr(um)-sqr(um+ul)-sqr(um+ur))/6.0;}
double StifRightFI4S (double ul, double um, double ur, double h)
       {return -1.0/(h*h)-(1.0-sqr(um+ur))/6.0;}

// FI4 Kevre 1
double SpectrumFI4K1 (double k, double h) {return sqrt(2.0+4.0*(1.0-h*h)*sqr(sin(0.5*k))/(h*h));}
double ForceFI4K1 (double ul, double um, double ur, double h)
       {return (ul-2.0*um+ur)/(h*h)+(1.0/4.0)*(ur+ul)*(2.0-ur*ur-ul*ul);}
double PotEnFI4K1 (double um, double ur, double h)
       {return 0.5*sqr((ur-um)/h)+0.25*sqr(1.0-sqr(0.5*ur+0.5*um));}
double StifLeftFI4K1 (double ul, double um, double ur, double h)
       {return -1.0/(h*h)-0.25*(2.0-sqr(ul)-sqr(ur))+0.5*ul*(ul+ur);}
double StifMiddleFI4K1 (double ul, double um, double ur, double h)
       {return 2.0/(h*h);}
double StifRightFI4K1 (double ul, double um, double ur, double h)
       {return -1.0/(h*h)-0.25*(2.0-sqr(ul)-sqr(ur))+0.5*ur*(ul+ur);}

// FI4 Kevre 2 (Al)
double SpectrumFI4K2 (double k, double h) {return sqrt(2.0+4.0*sqr(sin(0.5*k))/(h*h));}
double ForceFI4K2 (double ul, double um, double ur, double h)
       {return (ul-2.0*um+ur)/(h*h)+(1.0/2.0)*(ur+ul)*(1.0-um*um);}
double PotEnFI4K2 (double um, double ur, double h)
       {return 0.5*sqr((ur-um)/h)+0.25*sqr(1.0-sqr(0.5*ur+0.5*um));}
double StifLeftFI4K2 (double ul, double um, double ur, double h)
       {return -1.0/(h*h)-0.5*(1.0-sqr(um));}
double StifMiddleFI4K2 (double ul, double um, double ur, double h)
       {return 2.0*(1.0/(h*h)+0.5*(1.0-sqr(um)))-1.0+sqr(um)+um*(ul+ur);}
double StifRightFI4K2 (double ul, double um, double ur, double h)
       {return -1.0/(h*h)-0.5*(1.0-sqr(um));}


// FI4 Kevre 3
double SpectrumFI4K3(double k, double h) {return sqrt(2.0+(4.0-2.0*h*h)*sqr(sin(0.5*k))/(h*h));}
double ForceFI4K3(double ul, double um, double ur, double h)
       {return (ul-2.0*um+ur)/(h*h)+(1.0/8.0)*(ur+ul)*(4.0-ul*ul-2.0*um*um-ur*ur);}
double PotEnFI4K3(double um, double ur, double h)
       {return 0.5*sqr((ur-um)/h)+0.25*sqr(1.0-sqr(0.5*ur+0.5*um));}
double StifLeftFI4K3(double ul, double um, double ur, double h)
       {return -1.0/(h*h)-0.125*(4.0-sqr(ul)-2.0*sqr(um)-sqr(ur))+0.25*ul*(ur+ul);}
double StifMiddleFI4K3(double ul, double um, double ur, double h)
       {return 2.0/(h*h)+0.5*um*(ul+ur);}
double StifRightFI4K3(double ul, double um, double ur, double h)
       {return -1.0/(h*h)-0.125*(4.0-sqr(ul)-2.0*sqr(um)-sqr(ur))+0.25*ur*(ur+ul);}

// FI4 Kevre 4
double SpectrumFI4K4(double k, double h) {return sqrt(2.0+(4.0-2.0*h*h)*sqr(sin(0.5*k))/(h*h));}
double ForceFI4K4(double ul, double um, double ur, double h)
       {return (ul-2.0*um+ur)/(h*h)+um-0.5*um*um*(ul+ur);}
double PotEnFI4K4(double um, double ur, double h)
       {return 0.5*sqr((ur-um)/h)+0.25*sqr(1.0-sqr(0.5*ur+0.5*um));}
double StifLeftFI4K4(double ul, double um, double ur, double h)
       {return -1.0/(h*h)+0.5*um*um;}
double StifMiddleFI4K4(double ul, double um, double ur, double h)
       {return 2.0/(h*h)-um*um-1.0+um*(ul+um+ur);}
double StifRightFI4K4(double ul, double um, double ur, double h)
       {return -1.0/(h*h)+0.5*um*um;}


double Ch(double x) {return (exp(x)+exp(-x))/2.0;}
double Sh(double x) {return (exp(x)-exp(-x))/2.0;}
double Arth(double x) {return 0.5*log((1.0+x)/(1.0-x));}
double Q(double x, double t) {return (x-xk-dk*t)*deltak;}
//double knf(double x, double t) {return -1.0+(4.0/pi)*atan(exp((Arth(0.5*h)*2.0/h)*Q(x,t)));}
double knf(double x, double t) {return -1.0+(4.0/pi)*atan(exp(Q(x,t)));}
double trivial(double x, double t) {return -1.0;}
double Solution(double x, double t) {return knf(x,t);}
//double Solution(double x, double t) {return trivial(x,t);}
double Func(double u[], double t);
double Zero(double u[]);
double PotEnergy (double x[], double h);
double KinEnergy (double x0[],double x1[]);
double KinEnergyStormer (double x0[],double x1[],double x2[],double x3[]);
void ForceCulc (double x[], double y[], double h);
void PrintOutData();
void ShowMesh();
void ShowAtoms(int, double z[]);
void ReadIniCond();
void WriteIniCond();
void PrintAre(int LeftMargin, int TopMargin);
void PrintX(int LeftMargin, int TopMargin);
void PrintEigv(int LeftMargin, int TopMargin);
void Mass();
void MassBig();
void StiffBig();
void Jacobi(int n);
void Ordering();
void SaveEigenValues();

//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner){
u0 = new double [Nat*8]; u1 = new double [Nat*8];
u2 = new double [Nat*8]; u3 = new double [Nat*8];
f0 = new double [Nat*8]; f1 = new double [Nat*8];
f2 = new double [Nat*8]; f3 = new double [Nat*8];
 p = new double [Nat*8]; analit = new double [Nat*8];
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender){

ofstream OFSOkp("KinkPosition.txt",ios::out);
OFSOkp.precision(6);
ofstream OFSOkv("KinkVelocity.txt",ios::out);
OFSOkv.precision(6);

ofstream OFSOd("EnergyDeacay_h_1_05_to_1_25.txt",ios::out);
OFSOd.precision(10);

ofstream OFSOe("EnergyOfLattice.txt",ios::out);
OFSOe.precision(10);

ofstream OFSO("OmegaMAXandMIN.txt",ios::out);
OFSO.precision(10);
OFSOe.precision(10);

ofstream OFSOfp("ForPicture.txt",ios::out);

for (int i=1; i<401; i++){
    double H=0.01*i;
    OFSO.width(12); OFSO << H; OFSO << " ";
    OFSO.width(12); OFSO << SpectrumFI4K4(0.0,H); OFSO << " ";
    OFSO.width(12); OFSO << SpectrumFI4K4(pi,H); OFSO << " ";
    OFSO << "\n";
}

h=1.0;

while (h<1.00001){

OFSOd.width(12); OFSOd << h; OFSOd << "  ";

h2=h*h; C=1.0/(h*h);
// KINK
dk=0.0;
xk=0.5*h;
deltak=1.0/sqrt(1.0-dk*dk);

double Tsave,DTsave,Wmin,Wmax;
Tsave=0.0;
DTsave=1000.0;

Wmin=SpectrumFI4K4(0.0,h);
Wmax=SpectrumFI4K4(pi,h);
Form1->Canvas->TextOut(25,5,"h =");
Form1->Canvas->TextOut(80,5,h);
Form1->Canvas->TextOut(25,20,"Wmin =");
Form1->Canvas->TextOut(80,20,Wmin);
Form1->Canvas->TextOut(25,35,"Wmax =");
Form1->Canvas->TextOut(80,35,Wmax);

// Zero into files
// Zero into files
for(int i = 0; i<Nat; i++){
   u0[i]=-0.0; u1[i]=-0.0; u2[i]=-0.0; u3[i]=-0.0;
   f0[i]=0.0; f1[i]=0.0; f2[i]=0.0; f3[i]=0.0; p[i]=0.0;
}

//**********************************************
//  INITIAL  CONDITIONS
//**********************************************
double t0=tMIN;
Func(u1,t0-1.0*dt);
Func(u0,t0-0.0*dt);

// Initial conditions from file
//ReadIniCond();

// For boundary conditions
Left=-1.0;
Right=1.0;

//**********************************************
//  RELAXATIONAL  MD  CALCULATION
//**********************************************
t = t0;
nit = 0;
//while (waitfinish < 1){
while (t<tMAX){
    //  Stormer procedure
    ForceCulc(u0,f0,h);
    for(int i=1; i<Nat-1; i++)
       p[i]=2.0*u0[i]-u1[i]+f0[i]*dt2-fric*dt*(u0[i]-u1[i]);
       p[0]=Left;
       p[Nat-1]=Right;
//       p[Nat2]=0.0;
    // Replace
    for(int i=0; i<Nat; i++){
       u1[i] = u0[i]; u0[i] = p[i];
    }

    // Screen view
    if (fmod(nit,1000) == 0){
        ShowMesh();
        ShowAtoms(clRed,u0);
        double Kin,Pot;
        Pot = PotEnergy(u1,h);
        Kin = KinEnergy(u0,u1);
        Form1->Canvas->Brush->Color = clBtnFace;
        Form1->Canvas->TextOut(25,50,"t =");
        Form1->Canvas->TextOut(80,50,t);
        Form1->Canvas->TextOut(25,65,"Pot =");
        Form1->Canvas->TextOut(80,65,Pot);
        Form1->Canvas->TextOut(25,80,"Kin =");
        Form1->Canvas->TextOut(80,80,Kin);
        Form1->Canvas->TextOut(25,95,"Tot =");
        Form1->Canvas->TextOut(80,95,Kin+Pot);
    }// for Screen view

    t += dt;
    nit +=1;
}//while ( waitfinish == 0)

ShowMesh();
ShowAtoms(clRed,u1);
Zero(analit);
Func(analit,t);
ShowAtoms(clBlack,analit);


//******************************************
// Dispersion curves culculation ***********
//******************************************

// CALCULATION OF STIFFNESS MATRIX *********

// Formation of STIFF Short (small)
Stiff(h);
//PrintAre(300, 50);
// Mass matrix (small)
Mass();
// Formation of stiffness matrix (big)
StiffBig();
// Mass matrix (big)
MassBig();


// Save stifness matrix re
ofstream OFSO5("STIFFre.txt",ios::out);
OFSO5.precision(8);
for (int i=1; i<=NDF; i++){
   for (int j=1; j<=NDF; j++){
       OFSO5.width(16); OFSO5 << A[i][j];
   }
   OFSO5 << "\n";
}
OFSO5.close();

// Save mass matrix
ofstream OFSO7("MASS.txt",ios::out);
OFSO7.precision(8);
for (int i=1; i<=NDF; i++){
   for (int j=1; j<=NDF; j++){
       OFSO7.width(16); OFSO7 << b[i][j];
   }
   OFSO7 << "\n";
}
OFSO7.close();

fstream InFSOe("e.txt",ios::in);
InFSOe.precision(20);
for (int i=1; i<=NDF; i++){
   InFSOe.width(27); InFSOe >> eigv[i];
}

fstream InFSOB("R.txt",ios::in);
InFSOB.precision(20);
for (int i=1; i<=NDF; i++)
for (int j=1; j<=NDF; j++){
   InFSOB.width(27); InFSOB >> x[i][j];
}

// SHIFTING BY shift
//shift = 25.0;
//for (int i=1; i<=NDF; i++)
//for (int j=1; j<=NDF; j++)
//    A[i][j] += shift*b[i][j];
//Jacobi(NDF);
// BACK SHIFTING BY shift
//for (int i=1; i<=NDF; i++)
//    eigv[i] += -shift;

// Ordering of eigen frequencies
Ordering();

// Transformation of units
//for (int i=1; i<=NDF; i++){
//    if (eigv[i]>=0.0) eigv[i]=sqrt(eigv[i]);
//    else eigv[i]=-sqrt(-eigv[i]);
//}

PrintEigv(250, 50);
//PrintX(250, 50);

// Save Eigen-values
SaveEigenValues();

for (int i=1; i<=NDF; i++){
    double Wmi,Wma;
    if (Wmax>Wmin){ Wma=Wmax; Wmi=Wmin; }
    else          { Wma=Wmin; Wmi=Wmax; }
    if (eigv[i]<Wmi){
       OFSOfp.width(12); OFSOfp << h;       OFSOfp << " ";
       OFSOfp.width(12); OFSOfp << eigv[i]; OFSOfp << "\n";
    }
    if (eigv[i]>Wma){
       OFSOfp.width(12); OFSOfp << h;       OFSOfp << " ";
       OFSOfp.width(12); OFSOfp << eigv[i]; OFSOfp << "\n";
    }
}

//**********************************************
//  INITIAL  CONDITIONS
//**********************************************
// Restore files
for(int i = 0; i<Nat; i++){
   u1[i]=u0[i]; u2[i]=u0[i]; u3[i]=u0[i];
   f0[i]=0.0; f1[i]=0.0; f2[i]=0.0; f3[i]=0.0; p[i]=0.0;
}
t0=tMIN;

// add Goldstone Mode
double AMPL,PHASE,Pot,Kin;
int ModeNumber;
AMPL=0.3;
ModeNumber=1;
// Fing IM amplitude to have IM energy = ENERGYwanted

for (int i=0; i<Nat; i++){
    u0[i]+=AMPL*x[i+1][ModeNumber]*(t0+0.0*dt);
    u1[i]+=AMPL*x[i+1][ModeNumber]*(t0+1.0*dt);
    u2[i]+=AMPL*x[i+1][ModeNumber]*(t0+2.0*dt);
    u3[i]+=AMPL*x[i+1][ModeNumber]*(t0+3.0*dt);
}

// add Internal Mode
PHASE=0.0;
AMPL=0.0;
ModeNumber=2;
for (int i=0; i<Nat; i++){
    u0[i]+=AMPL*x[i+1][ModeNumber]*sin(eigv[ModeNumber]*(t0-0.0*dt+PHASE));
    u1[i]+=AMPL*x[i+1][ModeNumber]*sin(eigv[ModeNumber]*(t0-1.0*dt+PHASE));
    u2[i]+=AMPL*x[i+1][ModeNumber]*sin(eigv[ModeNumber]*(t0-2.0*dt+PHASE));
    u3[i]+=AMPL*x[i+1][ModeNumber]*sin(eigv[ModeNumber]*(t0-3.0*dt+PHASE));
}
Pot = PotEnergy(u1,h);
Kin = KinEnergyStormer(u0,u1,u2,u3);


ForceCulc(u0,f1,h);
ForceCulc(u0,f2,h);
ForceCulc(u0,f3,h);

// Initial conditions from file
//ReadIniCond();

// Write Initial conditions
WriteIniCond();


//**********************************************
//  MD  CALCULATION
//**********************************************
t = tMIN;
nit = 0;
int priznak=0;

OFSOd.width(12); OFSOd << h; OFSOd << " ";
OFSOd.width(12); OFSOd << t; OFSOd << " ";
OFSOd.width(12); OFSOd << Pot+Kin-8.0; OFSOd << " ";
OFSOd << "\n";
Tsave+=DTsave;

//while (waitfinish < 1){
while (t<tMAXMD){
    //  Stormer procedure
    ForceCulc(u0,f0,h);
    for(int i=1; i<Nat-1; i++)
       p[i] = 2.0*u0[i]-u1[i]+dt2*(k1*f0[i]-k2*f1[i]+k3*f2[i]-k4*f3[i]);
       p[0]=Left;
       p[Nat-1]=Right;
    // Replace
    for(int i=0; i<Nat; i++){
       u3[i] = u2[i]; u2[i] = u1[i]; u1[i] = u0[i]; u0[i] = p[i];
       f3[i] = f2[i]; f2[i] = f1[i]; f1[i] = f0[i];
    }

    // Screen view
    if (fmod(nit,100000) == 0){
        ShowMesh();
        ShowAtoms(clRed,u0);
        double Kin,Pot;
        Pot = PotEnergy(u1,h);
        Kin = KinEnergyStormer(u0,u1,u2,u3);
        Form1->Canvas->Brush->Color = clBtnFace;
        Form1->Canvas->TextOut(25,50,"t =");
        Form1->Canvas->TextOut(80,50,t);
        Form1->Canvas->TextOut(25,65,"Pot =");
        Form1->Canvas->TextOut(80,65,Pot);
        Form1->Canvas->TextOut(25,80,"Kin =");
        Form1->Canvas->TextOut(80,80,Kin);
        Form1->Canvas->TextOut(25,95,"Tot =");
        Form1->Canvas->TextOut(80,95,Kin+Pot);
        OFSOe.width(12); OFSOe << t; OFSOe << "  ";
        OFSOe.width(12); OFSOe << Pot; OFSOe << "  ";
        OFSOe.width(12); OFSOe << Kin; OFSOe << "  ";
        OFSOe.width(12); OFSOe << Pot+Kin; OFSOe << "\n";
        if (t>Tsave){
            OFSOd.width(12); OFSOd << Pot+Kin-8.0; OFSOd << "  ";
            Tsave+=DTsave;
        }
    }// for Screen view

    // Find kink position and velocity
    double kinkvel, tL,tR;
    if (fmod(nit,100) == 0){
       int kinkpos=0;
       for (int i=3; i<Nat-3; i++){
           if (u0[i]*u0[i+1]<0.0){
              kinkpos=i;
           }
       }
       if (kinkpos>Nat-20){
          if (priznak==0){
             priznak=1; tR=t;
          }
          else {
             tL=tR; tR=t;
             kinkvel=h*(Nat-40)/(tR-tL);
             OFSOkv.width(12); OFSOkv << t;       OFSOkv << "  ";
             OFSOkv.width(12); OFSOkv << kinkvel; OFSOkv << "\n";
          }
          for(int i=Nat-40; i<Nat; i++){
             u3[i-Nat+40] = u3[i]; u2[i-Nat+40] = u2[i]; u1[i-Nat+40] = u1[i]; u0[i-Nat+40] = u0[i];
             f3[i-Nat+40] = f3[i]; f2[i-Nat+40] = f2[i]; f1[i-Nat+40] = f1[i]; f0[i-Nat+40] = f0[i];
          }
          for(int i=39; i<Nat; i++){
             u3[i] = 1.0; u2[i] = 1.0; u1[i] = 1.0; u0[i] = 1.0;
             f3[i] = 0.0; f2[i] = 0.0; f1[i] = 0.0; f0[i] = 0.0;
          }
       }
       OFSOkp.width(12); OFSOkp << t;         OFSOkp << "  ";
       OFSOkp.width(12); OFSOkp << kinkpos*h; OFSOkp << "\n";
    }// for Find kink position



    t += dt;
    nit +=1;
}//while ( t<tMAX)

OFSOd << "\n";

//ShowMesh();
ShowAtoms(clRed,u1);
Zero(analit);
Func(analit,t);
ShowAtoms(clBlack,analit);

h+=0.05;
}// while (h<1.9)


//Cloze the file
OFSO.close();
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

// Substrate potential
double Subs(double r){
    return 1.0-cos(r);
}

// Derivative of Substrate potential
double DSubs(double r){
    return sin(r);
}

// Second derivative of Substrate potential
double DDSubs(double r){
    return cos(r);
}

// Inter-particle potential
double Bond(double r){
     return sqr(r)/(2.0*h*h);
}

// Derivative of Inter-particle potential
double DBond(double r){
     return r/(h*h);
}

// Second derivative of Inter-particle potential
double DDBond(double r){
     return 1.0/(h*h);
}

void ReadIniCond(){
    fstream InFSO(ReadDisplFrom,ios::in);
    InFSO.precision(20);
    for (int i=0; i<Nat; i++){
       InFSO.width(27); InFSO >> u0[i];
       InFSO.width(27); InFSO >> u1[i];
       InFSO.width(27); InFSO >> u2[i];
       InFSO.width(27); InFSO >> u3[i];
       InFSO.width(27); InFSO >> f1[i];
       InFSO.width(27); InFSO >> f2[i];
       InFSO.width(27); InFSO >> f3[i];
    }
}

void WriteIniCond(){
    ofstream OFSOd(WriteDisplTo,ios::out);
    OFSOd.precision(20);
    for (int i=0; i<Nat; i++){
       OFSOd.width(27); OFSOd << u0[i];
       OFSOd.width(27); OFSOd << u1[i];
       OFSOd.width(27); OFSOd << u2[i];
       OFSOd.width(27); OFSOd << u3[i];
       OFSOd.width(27); OFSOd << f1[i];
       OFSOd.width(27); OFSOd << f2[i];
       OFSOd.width(27); OFSOd << f3[i];
    }
}


double Func(double u[], double t){
   double xx;
   for (int i=0; i<Nat; i++){
       xx=(i-Nat2)*h; u[i]+=Solution(xx,t);
   }
   return *u;
} // Func


double Zero(double u[]){
   for (int i=0; i<Nat; i++) u[i]=0.0;
   return *u;
} // Zero

/*
// Calculation of forces
void ForceCulc (double u[], double f[]){
for(int i=1; i<Nat-1; i++)
        f[i] = 2.0*(sin(0.5*(u[i+1]-u[i]))-sin(0.5*(u[i]-u[i-1])))/h2
             - 0.5*(sin(0.5*(u[i+1]+u[i]))+sin(0.5*(u[i]+u[i-1])));
        f[0] = 2.0*(sin(0.5*(u[1]-u[0]))-sin(0.5*(u[0]-u[Nat-1])))/h2
             - 0.5*(sin(0.5*(u[1]+u[0]))+sin(0.5*(u[0]+u[Nat-1])));
        f[Nat-1] = 2.0*(sin(0.5*(u[0]-u[Nat-1]))-sin(0.5*(u[Nat-1]-u[Nat-2])))/h2
                 - 0.5*(sin(0.5*(u[0]+u[Nat-1]))+sin(0.5*(u[Nat-1]+u[Nat-2])));
//   f[i] = C * (u[i-1]-2.0*u[i]+u[i+1]) - sin(u[i]);
//   f[0] = C * (u[Nat-1]-2.0*u[0]+u[1]) - sin(u[0]);
//   f[Nat-1] = C * (u[Nat-2]-2.0*u[Nat-1]+u[0]) - sin(u[Nat-1]);
}
*/

// Calculation of forces
void ForceCulc (double u[], double f[], double h){
for(int i=1; i<Nat-1; i++)
        f[i] = ForceFI4K4 (u[i-1], u[i], u[i+1], h);
        f[0] = 0.0;
        f[Nat-1] = 0.0;
}


// Potential energy of the chain
double PotEnergy (double x[], double h){
  double Pot = 0.0;
  for (int i=0; i<Nat-1; i++)
      Pot += PotEnFI4S (x[i], x[i+1], h);
  return h*Pot;
}

// Kinetic energy of the chain
double KinEnergy (double x0[],double x1[]){
  double Kin = 0.0;
  for (int i=0; i<Nat; i++) Kin += 0.5*sqr((x0[i]-x1[i])/dt);
  return h*Kin;
}

// Kinetic energy of the chain
double KinEnergyStormer (double x0[],double x1[],double x2[],double x3[]){
  double Kin = 0.0;
  for (int i=0; i<Nat; i++) Kin += 0.5*sqr((2.0*x0[i]+3.0*x1[i]-6.0*x2[i]+x3[i])/(6.0*dt));
  return h*Kin;
}


// Show the frame for the picture of atoms
void ShowMesh(){
Form1->Canvas->Brush->Color = clWhite;
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
Form1->Canvas->TextOut(10,298,"+2.00");
Form1->Canvas->TextOut(10,373,"+1.00");
Form1->Canvas->TextOut(10,448," 0.00");
Form1->Canvas->TextOut(10,523,"-1.00");
Form1->Canvas->TextOut(10,598,"-2.00");
}

// Show displacements of atoms
void ShowAtoms(int Col, double x0[]){
int x,y;
for(int i=0; i<Nat; i++)
{x = 40 + floor(i*600.0/Nat);
y = 450 - floor(1.0*x0[i]*150.0/2.0);
Form1->Canvas->Pixels[x][y] = Col;
Form1->Canvas->Pixels[x+1][y] = Col;
Form1->Canvas->Pixels[x][y+1] = Col;
Form1->Canvas->Pixels[x+1][y+1] = Col;}}



void Mass(){
  for (int i=0; i<NDF; i++)
  for (int j=0; j<NDF; j++)
     MASS[i][j] = 0.0;
  for (int i=0; i<NDF; i++){
     MASS[i][i] = AtomMass;
  }
}

void MassBig(){
  for (int i=1; i<=NDF; i++)
  for (int j=1; j<=NDF; j++)
      b[i][j] = 0.0;
  for (int i=1; i<=NDF; i++){
      b[i][i] = AtomMass;
  }
}

void StiffBig(){
  for (int i=1; i<=NDF; i++)
  for (int j=1; j<=NDF; j++){
      A[i][j] = Are[i-1][j-1];
  }
}

void PrintAre(int LeftMargin, int TopMargin){
    for (int i=0; i<NDF; i++){
       for (int j=0; j<NDF; j++){
           double Print;
           Print=(1.0*ceil(Are[i][j]*10000.0))/10000.0;
           Form1->Canvas->TextOut(LeftMargin+45*j,TopMargin+15*i,Print);
       }
    }
}


void PrintX(int LeftMargin, int TopMargin){
    for (int i=1; i<=NDF; i++){
       for (int j=1; j<=NDF; j++){
           double Print;
           Print=(1.0*ceil(x[i][j]*100.0))/100.0;
           Form1->Canvas->TextOut(LeftMargin+45*j,TopMargin+15*(i-1),Print);
       }
    }
}

void PrintEigv(int LeftMargin, int TopMargin){
    for (int i=1; i<=NDF; i++){
       double Print;
       Print=(1.0*ceil(eigv[i]*10000.0))/10000.0;
       Form1->Canvas->TextOut(LeftMargin,TopMargin+15*(i-1),Print);
    }
}


void SaveEigenValues(){
  for (int i=1; i<=NDF; i++){
      OFSOe << h; OFSOe << " "; OFSOe << eigv[i]; OFSOe << "\n";
  }
}


void Jacobi(int n){
//variables
int k,nsmax,ifpr,error,i,j,jj,nsweep,nr,
    jp1,jn1,kp1,kn1;
double rtol,eps,eptola,eptolb,akk,ajj,ab,cneck,sqch,
       d1,d2,den,ca,cg,aj,ak,bk,bj,xj,xk,tol,dif,
       epsa,epsb,bb;

rtol=0.000000001; nsmax=16; ifpr=1; error=1;

   error=0;
   k=0;
   for (i=1; i<=n; i++)
   {
     if((A[i][i]<=0.0)||(b[i][i]<=0.0)) k=1;
     if(k>0) goto dva;
     d[i]=A[i][i]/b[i][i];
     eigv[i]=d[i];
   dva: }
   if(k==1) goto odin;
   for (i=1; i<=n; i++)
   {
     for (j=1; j<=n; j++) x[i][j]=0.0;
     x[i][i]=1.0;
   }
   if(n==1) goto tri;
   nsweep=0;
   nr=n-1;
   sorok: nsweep=nsweep+1;

   if(ifpr!=0) {cout << " nsweep= "; cout << nsweep;}

   eps=1.0;for (i=1; i<=nsweep; i++) eps=eps*0.01; eps=eps*eps;
   for (j=1; j<=nr; j++) { jj=j+1;

   if (ifpr!=0) {cout << " j= "; cout << j;cout << " nsweep= "; cout << nsweep;}

    for (k=jj; k<=n; k++) {
     eptola=(A[j][k]*A[j][k])/(A[j][j]*A[k][k]);
     eptolb=(b[j][k]*b[j][k])/(b[j][j]*b[k][k]);
     if((eptola<eps)&&(eptolb<eps)) goto dvaodinnol;
     akk=A[k][k]*b[j][k]-b[k][k]*A[j][k];
     ajj=A[j][j]*b[j][k]-A[j][k]*b[j][j];
     ab =A[j][j]*b[k][k]-A[k][k]*b[j][j];
     cneck=(ab*ab+4.0*akk*ajj)/4.0;
     if(cneck<0.0) goto odin;
     sqch=sqrt(cneck);d1=ab/2.0+sqch;d2=ab/2.0-sqch;den=d1;
     if((sqrt(d2*d2))>(sqrt(d1*d1))) den=d2;
     if(den==0) { ca=0.0;cg=-A[j][k]/A[k][k]; }
     else { ca=akk/den;cg=-ajj/den;}
     if(n!=2) { jp1=j+1;jn1=j-1;kp1=k+1;kn1=k-1;
     if(jn1>=1) { for (i=1; i<=jn1; i++) {
      aj=A[i][j];bj=b[i][j];ak=A[i][k];bk=b[i][k];
      A[i][j]=aj+cg*ak;A[i][k]=ak+ca*aj;
      b[i][j]=bj+cg*bk;b[i][k]=bk+ca*bj;}}
     if(kp1<=n) { for (i=kp1; i<=n; i++) {
      aj=A[j][i];bj=b[j][i];ak=A[k][i];bk=b[k][i];
      A[j][i]=aj+cg*ak;b[j][i]=bj+cg*bk;
      A[k][i]=ak+ca*aj;b[k][i]=bk+ca*bj;}}
     if(jp1<=kn1) { for (i=jp1; i<=kn1; i++) {
      aj=A[j][i];bj=b[j][i];ak=A[i][k];bk=b[i][k];
      A[j][i]=aj+cg*ak;b[j][i]=bj+cg*bk;
     A[i][k]=ak+ca*aj;b[i][k]=bk+ca*bj;}}}
     ak=A[k][k];bk=b[k][k];A[k][k]=ak+2.0*ca*A[j][k]+ca*ca*A[j][j];
     b[k][k]=bk+2.0*ca*b[j][k]+ca*ca*b[j][j];
     A[j][j]=A[j][j]+2.0*cg*A[j][k]+cg*cg*ak;
     b[j][j]=b[j][j]+2.0*cg*b[j][k]+cg*cg*bk;
     A[j][k]=0.0;b[j][k]=0.0;
     for (i=1; i<=n; i++) {
      xj=x[i][j];xk=x[i][k];
      x[i][j]=xj+cg*xk;x[i][k]=xk+ca*xj;}
    dvaodinnol: }}
    for (i=1; i<=n; i++) {
     if((A[i][i]<=0.0)||(b[i][i]<=0.0)) goto odin;
     eigv[i]=A[i][i]/b[i][i];}
    j=0;for (i=1; i<=n; i++) {
     tol=rtol*d[i];dif=sqrt((eigv[i]-d[i])*(eigv[i]-d[i]));
     if(dif>tol) j=1;}
    if(j==1)  goto dvavosemnol; eps=rtol*rtol;
    for (j=1; j<=nr; j++){ i=0;
     if(i>0)  goto dvapatnol; jj=j+1; for (k=jj; k<=n; k++) {
      epsa=(A[j][k]*A[j][k])/(A[j][j]*A[k][k]);
      epsb=(b[j][k]*b[j][k])/(b[j][j]*b[k][k]);
      if((epsa>=eps)||(epsb>=eps))  i=1; dvapatnol:}}
      if(i==1) goto dvavosemnol;
    dvapatpat:for (i=1; i<=n; i++) { for (j=1; j<=n; j++) {
     A[j][i]=A[i][j];b[j][i]=b[i][j];}}
    for (j=1; j<=n; j++) { bb=sqrt(b[j][j]);
     for (k=1; k<=n; k++) x[k][j]=x[k][j]/bb; } goto tri;
    dvavosemnol: for (i=1; i<=n; i++) d[i]=eigv[i];
    if (nsweep<nsmax) goto sorok;
    goto dvapatpat;
    odin: cout << " ERROR SOL.STOP: MATR NOT POSITIV DEF";
       error=1;
    tri:
} // Jacobi();

/*
void Ordering(){
for (int j=1; j<=NDF; j++)
for (int i=2; i<=NDF; i++){
    if (eigv[i-1]>eigv[i]){
        double aaa = eigv[i];
        eigv[i] = eigv[i-1];
        eigv[i-1] = aaa;
    }
}
}// for void Ordering()
*/

void Ordering(){
for (int g=2; g<=NDF; g++)
for (int f=2; f<=NDF; f++){
    if (eigv[f]<eigv[f-1]){
       double p=eigv[f];
       eigv[f]=eigv[f-1];
       eigv[f-1]=p;
       for (int i=1; i<=NDF; i++){
          double p=x[i][f];
          x[i][f]=x[i][f-1];
          x[i][f-1]=p;
       }
    }
}
}// for void Ordering()


void Stiff(double h){
int mm;
double Sinus,Cosinus,Rx,Rx2,Ry,Ry2,RLen,RLen2,RLen3,ShortF,ShortFp,A11,A12,A22;
for (int i=0; i<NDF; i++)
for (int j=0; j<NDF; j++){
    Are[i][j] = 0.0;
}
for (int nx=1; nx<Nat-1; nx++){
    Are[nx][nx-1]=StifLeftFI4K4   (u0[nx-1], u0[nx], u0[nx+1], h);
    Are[nx][nx]  =StifMiddleFI4K4 (u0[nx-1], u0[nx], u0[nx+1], h);
    Are[nx][nx+1]=StifRightFI4K4  (u0[nx-1], u0[nx], u0[nx+1], h);
} // for nx
Are[0][0]=Are[1][1];
Are[0][1]=Are[1][2];
Are[Nat-1][Nat-1]=Are[Nat-2][Nat-2];
Are[Nat-1][Nat-2]=Are[Nat-2][Nat-3];
} // void ShortStiff()



