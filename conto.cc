// Conro.cin is the C++ routine to make the contourplot in the EPS format
// It can be called as conto(o,g,w,v,X,Y,M,N, Lev ,-p,p);
// FILE *o is the output file; it should be opened for writing and the EPS header should be already there.
// (The header can be written with routine ado.cin ).
// double *g is single-dimensional array of length (M+1)*(N+1) used to transfer values of the function to the routine. 
// g[m+(M+1)*n] is interpreted as value of the function at the point of grid with numbers m,n;
// -1<m<M+1; -1<n<N+1
// double *w is working array of length (M+1)*(N+1).
// char *v is working array of length (M+1)*(N+1), 
// v it is used to store the mark of each sell as "visited" to avoid drawing the same line twice.
// double *X is array of length M+1; the abscissas of the grid points should be stored there at the calling of conto.
// double *Y is array of length N+1; the ordinates of the grid points should be stored there at the calling of conto.
// int M is number of cells along abscissas; number of the grid points along x axis is M+1
// int N is number of cells along ordinates; number of the grid points along y axis is N+1
// double Lev, level to be drawn. (At a single call, the only one level is drawn)
// double p and -p, should be something of type double; p and -p determine the interval used for plotting: 
// values smaller than Lev-p or greater than Lev+p are interpreted as "singularities" of the function.
// (only once I used non–symmetric limits, but I still keep this option, some expression may be placed instead of -p)
// Routine 'conto is used in generator Tetre2215.cc to plot the contours of tetration.
// Please let me know if any problem with this routine.
// Copyleft 2008-2011 by Dmitrii Kouznetsov.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define DB double
#include"ado.cc"
#define DB double
#define DO(x,y) for(x=0;x<y;x++)
#define M(x,y) fprintf(o,"%6.4f %6.4f M\n",1.*(x),1.*(y));
#define L(x,y) fprintf(o,"%6.4f %6.4f L\n",1.*(x),1.*(y));
//#define o(x,y) fprintf(o,"%5.3f %5.3f o\n",1.*(x),1.*(y));
#define Mxy M(x,y)
#define Lxy L(x,y)
#define f(m,n) F[(m)*N1+(n)]
#define z(m,n) Z[(m)*N1+(n)] 
#define zmn z(m,n)
#define zMn z(m+1,n)
#define zmN z(m,n+1)
#define zMN z(m+1,n+1)
#define fmn f(m,n)
#define fMn f(m+1,n)
#define fmN f(m,n+1)
#define fMN f(m+1,n+1)
#define Xm X[m] 
#define XM X[m+1] 
#define Yn Y[n]
#define YN Y[n+1]
#define bdpq {b=f(m,n+1);d=f(m+1,n+1);p=f(m,n);q=f(m+1,n);}
#define UPP 1
#define LEF 2
#define DOW 3
#define RIG 4
DB drift(FILE *o,DB *F,char *Z,DB *X,DB *Y,int M,int N,int m,int n,int K)
{int M1=M+1,N1=N+1; DB b,d, p,q, x,y, B,D,P,Q; int mO=m,nO=n;
 //printf("drift: K=%2d, m=%2d n=%2d  \n",K,m,n);
if(K==UPP) goto Up;
if(K==LEF) goto Le;
if(K==DOW) goto Do;
if(K==RIG) goto Ri;

Up:if(zmn=='b'||zmn=='d'||zmn=='p'||zmn=='q'||zmn=='o'){return 0.;}
    if(zmn=='|'||zmn=='+'||zmn=='/'||zmn=='L'){if(m==mO&&n==nO)fprintf(o,"C\n");return 0.;}
bdpq; 
 //printf("Up: m=%2d n=%2d bdpq=%5.2f %5.2f %5.2f %5.2f\n",m,n,b,d,p,q); //getchar();
if(b*d<0){x=Xm+(XM-Xm)*b/(b-d);y=YN;Lxy;if(zmn=='-')zmn='+';else zmn='|';n++;
						    if(n>=N) return 0.; goto Up;} 
if(q*d<0){y=Yn+(YN-Yn)*q/(q-d);x=XM;Lxy;zmn='/';m++;if(m>=M)return 0; goto Ri;}
if(b*p<0){y=Yn+(YN-Yn)*p/(p-b);x=Xm;Lxy;zmn='L';m--;if(m<0)return 0; goto Le;}
 //printf("handle zero, m=%2d n=%2d\n", m,n);
if(d*d==0){ zmn='o';m++;n++;L(X[m],Y[n]); //printf("Go to UR"); 
 		goto UR;}
if(b*b==0){zmn='o';n++; L(X[m],Y[n]);m--; goto UL;}
 //end Up

Le:if(zmn=='b'||zmn=='d'||zmn=='p'||zmn=='q'||zmn=='o'){return 0.;}
    if(zmn=='-'||zmn=='+'||zmn=='/'||zmn=='L'){if(m==mO&&n==nO)fprintf(o,"C\n");return 0.;}
bdpq;
 //printf("Le: m=%2d n=%2d bdpq=%5.2f %5.2f %5.2f %5.2f\n",m,n,b,d,p,q); 
 if(b*p<0){y=Yn+(YN-Yn)*p/(p-b);x=Xm;Lxy;if(zmn=='|') zmn='+';
 				   else zmn='-';m--;if(m<0)return 0; goto Le;}
 if(p*q<0){x=Xm+(XM-Xm)*p/(p-q);y=Yn;Lxy;zmn='/';n--;if(n<0)return 0; goto Do;}
 if(b*d<0){x=Xm+(XM-Xm)*b/(b-d);y=YN;Lxy;zmn='L';n++;if(n>=N)return 0; goto Up;}
 //printf("Le Handles zero\n");
 if(p*p==0){zmn='o';L(Xm,Yn);m--;n--;if(m<0||n<0) return 0;// printf("go to LD\n");
		goto LD;}
 if(b*b==0){zmn='o';n++;L(Xm,Yn);m--;if(m<0||n>=N)return 0;// printf("go to LU\n");
 		goto LU;}
 //end Le
 Do://come to cell m,n from the up and expect to go down.
    if(zmn=='b'||zmn=='d'||zmn=='p'||zmn=='q'||zmn=='o') return 0.;
    if(zmn=='|'||zmn=='+'||zmn=='/'||zmn=='L'){if(m==mO&&n==nO)fprintf(o,"C\n");return 0.;}
 bdpq; //printf("Do: m=%2d n=%2d bdpq=%5.2f %5.2f %5.2f %5.2f\n",m,n,b,d,p,q); 
if(p*q<0){x=Xm+(XM-Xm)*p/(p-q);y=Yn;Lxy;if(zmn=='-')zmn='+';
				else zmn='|';n--;if(n<0)return 0; goto Do;}
if(b*p<0){y=Yn+(YN-Yn)*p/(p-b);x=Xm;Lxy;zmn='/';m--;if(m<0) return 0; goto Le;}
if(q*d<0){y=Yn+(YN-Yn)*q/(q-d);x=XM;Lxy;zmn='L';m++;if(m>=M) return 0; goto Ri;}
//printf("handle zero\n");
if(p*p==0){zmn='o';L(Xm,Yn);m--;n--; goto DL;}
if(q*q==0){zmn='o';m++;L(Xm,Yn);n--; goto DR;}
//end Do
Ri: //expect to go right..
   if(zmn=='b'||zmn=='d'||zmn=='p'||zmn=='q'||zmn=='o')return 0.;
   if(zmn=='-'||zmn=='+'||zmn=='/'||zmn=='L'){if(m==mO&&n==nO)fprintf(o,"C\n");return 0.;}
bdpq;//printf("Ri: m=%2d n=%2d bdpq=%5.2f %5.2f %5.2f %5.2f\n",m,n,b,d,p,q); 
if(d*q<0){y=Yn+(YN-Yn)*q/(q-d);x=XM;Lxy;if(zmn=='-') zmn='+';
					 else zmn='|';m++;if(m>=M)return 0; goto Ri;}
if(b*d<0){x=Xm+(XM-Xm)*b/(b-d);y=YN;Lxy;zmn='/';n++;if(n>=N) return 0; goto Up;}
if(p*q<0){x=Xm+(XM-Xm)*p/(p-q);y=Yn;Lxy;zmn='L';n--;if(n<0)return 0; goto Do;}
//printf("handle zero\n");
if(d*d==0){zmn='o';m++;n++;L(Xm,Yn); goto UR;} 
if(q*q==0){zmn='o';m++;L(Xm,Yn);n--; goto DR;}
//if(n<0 ||m>=M) return 0; goto Ri;} 
return 0;
//end Ri
DL:
LD: //printf("LD m=%2d n=%2d (may be negative)\n",m,n);
//came to the cell (m,n) from upper right corner.
//This cell may exist at the mesh; check this option first.
if(m<0&&n<0) return 0; //corner of the mesh;
if(m<0){m++;bdpq;if(p*q<=0 && zmn==' ') goto Do; return 0;}
if(n<0){n++;bdpq;if(b*p<=0 && zmn==' ') goto Le; return 0;}
bdpq; // pri("inside. bdpq=%5.2f %5.2f %5.2f %5.2f\n", b,d,p,q); 
		if(p*q<=0 && zmn==' ') goto Do;
		if(b*p<=0 && zmn==' ') goto Le;
Q=f(m+2,n); //pri("Q=%5.2f\n",Q);
if(Q*q<=0){m++; 
			    if(zmn==' ') 
					{//printf("go to Do, m=%2d n=%2d\n",m,n);
					goto Do;} return 0;}
B=f(m,n+2); if(B*b<=0){n++; if(zmn==' ') goto Le; return 0;}
return 0;

LU:UL:
//printf("UL: m=%2d n=%2d ( may be out of mesh)\n",m,n); //come from right down.
if(m<0&&n>=N) return 0; //corner of the mesh;
if(m<0){m++;bdpq;if(b*d<=0 && zmn==' ') goto Up; return 0;}
if(n>=N){n--;bdpq;if(b*p<=0 && zmn==' ') goto Le; return 0;}
bdpq; //pri("inside. bdpq=%5.2f %5.2f %5.2f %5.2f\n", b,d,p,q);
                if(b*p<=0 && zmn==' ') goto Le;
                if(b*d<=0 && zmn==' ') goto Up;
D=f(m+2,n+1); //pri("D=%5.2f\n",D);
if(D*d<=0){m++;if(zmn==' ') goto Up;}
	//Q=f(m-1,n+1); //pr("Q=%5.2f\n",Q);
	if(Q*q<=0){m--;if(zmn==' ') goto Le;}
P=f(m,n-1); //pri("P=%5.2f\n",P);
if(P*p<=0){n--;if(zmn==' ') goto Le;}
//if(p*p==0){L(X[m],Y[n]) zmn='-'; Z[m*N1+n-1]='-'; m--;n--; goto LD;}
return 0;

RU:
UR: //printf("UR: m=%2d n=%2d\n",m,n); //come from left down. May be out of mesh.
if(m>=M&&n>=N) return 0; //corner of the mesh;
if(m>=M){m--;bdpq;if(b*d<=0 && zmn==' ') goto Up; return 0;} //more lines
//   if(m>=M){m--;bdpq;if(b*d< 0 && zmn==' ') goto Up; return 0;}  // less lines
if(n>=N){n--;bdpq;if(d*q<=0 && zmn==' ') goto Ri; return 0;} //more lines
//   if(n>=N){n--;bdpq;if(d*q<=0 && zmn==' ') goto Ri; return 0;} //less
bdpq; //pri("inside. bdpq=%5.2f %5.2f %5.2f %5.2f\n", b,d,p,q); 
		if(d*q<=0 && zmn==' ') goto Ri;
		if(b*d<=0 && zmn==' ') goto Up;
B=f(m-1,n+1); //pri("Q=%5.2f\n",Q);
if(B*b<=0){m--; if(zmn==' ') goto Up; return 0;}
Q=f(m+1,n-1); //pri("D=%5.2f\n",D);
if(Q*q<=0){n--; if(zmn==' ') goto Ri; return 0;}
return 0;

DR: RD: //printf("RD: m=%3d n=%2d\n",m,n);
if(m>=M&&n<0) return 0; //corner of the mesh;
if(m>=M){m--;bdpq;if(p*q<=0 && zmn==' ') goto Do; return 0;}
if(n<0 ){n++;bdpq;if(b*p<=0 && zmn==' ') goto Ri; return 0;}
bdpq; //pri("inside. bdpq=%5.2f %5.2f %5.2f %5.2f\n", b,d,p,q); 
		if(p*q<=0 && zmn==' ') goto Do;
		if(d*q<=0 && zmn==' ') goto Ri;
Q=f(m+2,n); //pri("Q=%5.2f\n",Q);
if(Q*q<=0){m++; if(zmn==' ') 
			    {//printf("go to Do, m=%2d n=%2d\n",m,n);
				goto Do;}
  return 0;}
D=f(m+1,n+2); //pri("D=%5.2f\n",D);
	 if(D*d<=0){n++; if(zmn==' ') goto Ri; return 0;}

return 0;}//end drift

DB conto(FILE *o,DB *G,DB *F,char *Z,DB *X, DB *Y,int M,int N,DB L,DB L1,DB L2)
{int m,n; int M1=M+1,N1=N+1;  DB w,	b,d, p,q, x,y;
// printf("conto (copyleft 2008 by Dmitrii Kouznetsov) draws level L=%6.3f\n",L);
//printf("conto draws L=%6.3f\n",L);

//for(n=N;n>=0;n--){DO(m,M1)printf("%5.2f",G[m*N1+n]); printf("\n");}
//getchar();

//DO(m,M1){ M(X[m],Y[0]);L(X[m],Y[N]);}
//DO(n,N1){ M(X[0],Y[n]);L(X[M],Y[n]);}
//fprintf(o,".001 W 0 0 0 RGB S\n");

DO(m,M1)DO(n,N1)z(m,n)=' ';
DO(m,M1)DO(n,N1)
{ w=G[m*N1+n]-L;
  if(L1<w && w<L2) { F[m*N1+n]=w;}
  else{ //o(X[m],Y[n]);
	F[m*N1+n]=0; z(m,n)='p';if(m> 0) z(m-1,n  )='q';
	if(n>0) {z(m,n-1)='b';  if(m> 0 )z(m-1,n-1)='d';}
      }
}

//for(n=N;n>=0;n--){DO(m,M1)pri("%5.2f",F[m*N1+n]); printf("\n");}
//for(n=N;n>=0;n--){DO(m,M1)printf("%2c" , Z[m*N1+n]); printf("\n");}
//getchar();

//printf("Z1 Z2= %c %d   %c %d \n",Z[1],Z[1],Z[2],Z[2]);

//DB t,u,v;//Begin with singularities
// if singularity at the down-left of the cell

DO(m,M-1)
DO(n,N-1)
{ if(zmn=='p')
  {bdpq; //pri("Sp: %2d %2d %5.2f %5.2f %5.2f %5.3f\n",m,n,b,d,p,q);
   if(b*d<0&&z(m,n+1)==' '){x=Xm+(XM-Xm)*b/(b-d);y=YN;Mxy; drift(o,F,Z,X,Y,M,N,m,n+1,UPP);}
   if(q*d<0&&z(m+1,n)==' '){y=Yn+(YN-Yn)*q/(q-d);x=XM;Mxy; drift(o,F,Z,X,Y,M,N,m+1,n,RIG);}
   if(d*d==0&&z(m+1,n+1)==' '){                    M(XM,YN); drift(o,F,Z,X,Y,M,N,m+1,n+1,RIG);}
  }
}

// Check for singularity at down-right
for(m=1;m<M;m++)
DO(n,N-1)
{if(zmn=='q') //how about to go up-left?
 {bdpq; //pri("Sq: %2d %2d %5.2f %5.2f %5.2f %5.3f\n",m,n,b,d,p,q);
  if(b*d<0&&z(m,n+1)==' '){x=Xm+(XM-Xm)*b/(b-d);y=YN;Mxy;drift(o,F,Z,X,Y,M,N,m,n+1,UPP);}
  if(b*p<0&&z(m-1,n)==' '){y=Yn+(YN-Yn)*p/(p-b);x=Xm;Mxy;drift(o,F,Z,X,Y,M,N,m-1,n,LEF);}
  if(b*b==0&&z(m-1,n+1)==' '){                    M(Xm,YN);drift(o,F,Z,X,Y,M,N,m-1,n+1,LEF);}
  }
}

//Check for singularity at the top left corger. How about to go down-right?
DO(m,M-1)
for(n=1;n<N;n++)
{ if(zmn=='b')
  {bdpq; //pri("Sb: %2d %2d %5.2f %5.2f %5.2f %5.3f\n",m,n,b,d,p,q);
   if(q*p<0&&z(m,n-1)==' '){x=Xm+(XM-Xm)*p/(p-q);y=Yn;Mxy;drift(o,F,Z,X,Y,M,N,m,n-1,DOW);}
   if(q*d<0&&z(m+1,n)==' '){y=Yn+(YN-Yn)*q/(q-d);x=XM;Mxy;drift(o,F,Z,X,Y,M,N,m+1,n,RIG);}
   if(q*q==0&&z(m-1,n-1)==' '){                  M(XM,Yn);drift(o,F,Z,X,Y,M,N,m+1,n-1,RIG);}
  }
}

for(m=1;m<M;m++)
for(n=1;n<N;n++)
{ if(zmn=='d') // singularity at the up-right corner of this sell. go down-left?
  {bdpq; //pri("Sd: %2d %2d %5.2f %5.2f %5.2f %5.3f\n",m,n,b,d,p,q);
   if(p*q<0&&z(m,n-1)==' '){x=Xm+(XM-Xm)*p/(p-q);y=Yn;Mxy;drift(o,F,Z,X,Y,M,N,m,n-1,DOW);}
   if(p*b<0&&z(n-1,n)==' '){y=Yn+(YN-Yn)*p/(p-b);x=Xm;Mxy;drift(o,F,Z,X,Y,M,N,m-1,n,LEF);}
   if(p*p==0&&z(m-1,n-1)==' '){                  M(Xm,Yn);drift(o,F,Z,X,Y,M,N,m-1,n-1,LEF);}
} }

 //Trace the margin of the domain

n=0; // printf("n=%3d\n",n);
DO(m,M)
 { if(zmn==' ')	{ bdpq;
		  if(p*q<=0){  if(p>q || q>p)
			    {  x=Xm+(XM-Xm)*p/(p-q);
			       y=Yn; Mxy;
			       drift(o,F,Z,X,Y,M,N,m,n,UPP);}}
		}
 }
n=N-1; //printf("n=%3d\n",n);
DO(m,M)
 {  if(zmn==' '){bdpq;
		  if(b*d<=0){ if(b>d || d>b) 
			    {	x=Xm+(XM-Xm)*b/(b-d);
				y=YN; Mxy;
			  	drift(o,F,Z,X,Y,M,N,m,n,DOW);}}
		}  
 }

m=0; //printf("m=%3d\n",m);
DO(n,N)
 {  if(zmn==' '){bdpq;
		  if(b*p<=0) { if(p>b || p<b) 
			     {	y=Yn+(YN-Yn)*p/(p-b); 
				x=Xm; Mxy;
			  	drift(o,F,Z,X,Y,M,N,m,n,RIG);}}
		}
 }

// This was repaired
m=M-1; //printf("m=%3d\n",m);
DO(n,N)
 {  if(zmn==' '){bdpq;
		  if(d*q<=0) { if( q>d || q<d  )  
			     { 	y=Yn+(YN-Yn)*q/(q-d);
				x=XM; Mxy;
			  	drift(o,F,Z,X,Y,M,N,m,n,LEF);}}
		}
 }


//Check if any loops inside the domain
// for(n=1;n<N-1;n++) 
// for(m=1;m<M-1;m++)
for(n=N-1;n>0;n--) 
for(m=M-1;m>0;m--)
 { 
   if(zmn==' ')	{ bdpq; 
		 if(d*q<0){     y=Yn+(YN-Yn)*q/(q-d); x=XM;Mxy;
				  drift(o,F,Z,X,Y,M,N,m+1,n,RIG);
   			    if(zmn==' '){ Mxy;
				  drift(o,F,Z,X,Y,M,N,m,n,LEF);
					}
			  }
		}
 }
//   if(zmn==' ')	{ if(b*d<0){ 	x=Xm+(XM-Xm)*b/(b-d);	 y=YN; Mxy;
//				  drift(o,F,Z,X,Y,M,N,m,n,UPP);}}}

return 0;}
//end 
#undef Mxy
#undef Lxy
#undef f
#undef z
#undef zmn
#undef zMn
#undef zmN
#undef zMN
#undef fmn
#undef fMn
#undef fmN
#undef fMN
#undef Xm
#undef XM
#undef Yn
#undef YN
#undef bdpq
#undef UPP
#undef LEF
#undef DOW
#undef RIG
/* End of routine */

