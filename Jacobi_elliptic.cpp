#include<iostream>
#include<algorithm>
#include<cmath>
#include<fstream>
using namespace std;

#define sqr(x) x*x
#define pi 3.14159265358979323846

double sn, cn, dn, JacobiPeriod;

void JacobiElliptic (double z, double m){
     double a[1000], b[1000], c[1000], fi[1000];
     double eps=2.0e-14;
     a[0]=1.0, b[0]=sqrt(1.0-m), c[0]=sqrt(m);
     int N=0;
     while (sqrt(sqr(c[N]))>eps){
        double a1,b1,c1;
        a1=0.5*(a[N]+b[N]);
        b1=sqrt(a[N]*b[N]);
        c1=0.5*(a[N]-b[N]);
        N+=1;
        a[N]=a1;
        b[N]=b1;
        c[N]=c1;
     }// while (sqrt(sqr(c))>eps)
     int N2=1;
     for (int i=0; i<N; i++) N2=N2*2;
     fi[N]=N2*a[N]*z;
     for (int i=0; i<N; i++){
         double arg=(c[N-i]/a[N-i])*sin(fi[N-i]);
         fi[N-i-1]=0.5*(fi[N-i]+atan(arg/sqrt(1.0-arg*arg)));
     }
     sn=sin(fi[0]);
     cn=cos(fi[0]);
     dn=cos(fi[0])/cos(fi[1]-fi[0]);
     JacobiPeriod=2.0*pi/a[N];
}// For void JacobiElliptic (double z, double m)

int main(){
    ofstream fout;
    int cross = 0, counting;
    double time=0.0,dt = 0.005, u_prev = 0.0, u = 0, T2, T1;
    fout.open("u_vs_tt.txt");
    if(fout){
        while (time<10.0){
            double A,AA,P,M,q1,q2,q3,a1,a3,a5;
            A=4.9; AA=A*A;
            a1=3.0;
            a3=1.0/6.0;
            a5=1.0/120.0;
            q1=a1+a3*AA+a5*AA*AA;
            q2=6.0*a1+3.0*a3*AA+ 2.0*a5*AA*AA;
            q3=4.0*a1+3.0*a3*AA+ 2.0*a5*AA*AA;
            P=sqrt(sqrt(q1*q2/6.0))*time;
            M=0.5-0.25*q3*sqrt(3.0/(2*q1*q2));
            JacobiElliptic(P,M);
            u_prev = u;
            u=(A*cn)/sqrt(cn*cn+sqrt(6.0*q1/q2)*sn*sn*dn*dn);
            if(cross<3){
                if(u_prev*u < 0){
                    if(cross == 0){
                        T1 = time;
                    }
                    cross += 1;

                    if(cross == 3){
                        T2 = time;
                    }
                }
            }

            fout<<time<<"\t"<<u<<endl;

            //OFSOes.width(12); OFSOes << time;       OFSOes << " ";  // 1
            //OFSOes.width(12); OFSOes << u;     OFSOes << "\n";  // 2
            time += dt;
        }
    }
    cout<<(T2 - T1)<<endl;
    fout.close();
}
