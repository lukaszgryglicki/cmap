#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#define DB double
#define DO(x,y) for(x=0;x<y;x++)
using namespace std;
typedef complex<double> z_type;
#define Re(x) x.real()
#define Im(x) x.imag()
#define I z_type(0.,1.)

#include "libfparser.h"
#include "conto.cc"

int eps_cont(const char* fn, int M, int N, 
		DB sx, DB ex, DB sy, DB ey, 
		int ngx, int ngy, 
		DB mir, DB mii, DB mim, DB mar, DB mai, DB mam, 
		DB lev_rf, DB lev_rt, DB lev_ri, 
		DB lev_if, DB lev_it, DB lev_ii, 
		DB lev_mf, DB lev_mt, DB lev_mi, 
		int dr, int di, int dm, 
		complex<double> (*fp)(complex<double>))
{
    printf("fn=%s M=%d N=%d sx=%f ex=%f sy=%f ey=%f ngx=%d ngy=%d mir=%f mii=%f mim=%f mar=%f mai=%f mam=%f lrf=%f lrt=%f lri=%f lif=%f lit=%f lii=%f lmf=%f lmt=%f lmi=%f dr=%d di=%d dm=%d f=%p\n", 
	    fn, M, N, sx, ex, sy, ey, ngx, ngy, mir, mii, mim, mar, mai, mam, lev_rf, lev_rt, lev_ri, lev_if, lev_it, lev_ii, lev_mf, lev_mt, lev_mi, dr, di, dm, fp);

    FILE *o;
    int m, M1=M+1;
    int n, N1=N+1;
    DB X[M1], Y[N1];
    DB *g1, *g2, *g3, *w1, *w2, *w3;
    DB x, y, p, q;
    z_type z1, z2;
    char *v1, *v2, *v3;

    g1 = new DB[M1*N1];
    g2 = new DB[M1*N1];
    g3 = new DB[M1*N1];
    w1 = new DB[M1*N1];
    w2 = new DB[M1*N1];
    w3 = new DB[M1*N1];
    v1 = new char[M1*N1];
    v2 = new char[M1*N1];
    v3 = new char[M1*N1];
    
    DB scaX = DB(M)/(ex-sx);
    DB scaY = DB(N)/(ey-sy);
    DB cx   = (ex+sx)/2.;
    DB cy   = (ey+sy)/2.;
    DB traX = DB(M)/2.-cx*scaX;
    DB traY = DB(N)/2.-cy*scaY;
    
    o = fopen(fn,"w");
    ado(o,M,N);

    fprintf(o,"%f %f translate\n%f %f scale\n", traX, traY, scaX, scaY);

    DO(m,M1) X[m] = sx + ((ex-sx)*m)/M;
    DO(n,N1) Y[n] = sy + ((ey-sy)*n)/N;

    DO(m,M1) DO(n,N1) g1[m*N1+n] = g2[m*N1+n] = g3[m*N1+n] = 0.;

    DO(m,M1)
    {
	x=X[m];
	DO(n,N1)
	{
	    y = Y[n]; 
	    z1 = z_type(x,y);
	    z2 = fp(z1);
	    p = Re(z2);
	    q = Re(z1);
	    if (p > mir && p < mar)  g1[m*N1+n] = p;
	    else g1[m*N1+n] = q;
	    p = Im(z2);
	    q = Im(z1);
	    if (p > mii && p < mai)  g2[m*N1+n] = p;
	    else g2[m*N1+n] = q;
	    p = abs(z2);
	    q = abs(z1);
	    if (p > mim && p < mam)  g3[m*N1+n] = p;
	    else g3[m*N1+n] = q;
	}
    }

    if (ngx > 0) for (x=sx;x<=ex;x+=(ex-sx)/(DB)ngx) 
    {
	M(x,sy)L(x,ey)
    }
    if (ngy > 0) for (y=sy;y<=ey;y+=(ey-sy)/(DB)ngy) 
    {
	M(sx,y)L(ex,y)
    }
    fprintf(o,".001 W 0 0 0 RGB S\n");

    fprintf(o,"1 setlinejoin 1 setlinecap\n");

    if (dm)
    {
	p = (mam-mim)/2.;
	for (q=lev_mf;q<=lev_mt;q+=lev_mi) conto(o,g3,w3,v3,X,Y,M,N,q,-p,p);
	fprintf(o,".001 W 0 0.5 0 RGB S\n");
    }
    
    if (di)
    {
	p = (mai-mii)/2.;
	for (q=lev_if;q<=lev_it;q+=lev_ii) conto(o,g2,w2,v2,X,Y,M,N,q,-p,p);
	fprintf(o,".001 W 0 0 0.8 RGB S\n");
    }
    
    if (dr)
    {
	p = (mar-mir)/2.;
	for (q=lev_rf;q<=lev_rt;q+=lev_ri) conto(o,g1,w1,v1,X,Y,M,N,q,-p,p);
	fprintf(o,".001 W 0.7 0 0 RGB S\n");
    }
    
    fprintf(o,"showpage\n%cTrailer",'%');
    fclose(o);

    delete [] g1;
    delete [] g2;
    delete [] g3;
    delete [] w1;
    delete [] w2;
    delete [] w3;
    delete [] v1;
    delete [] v2;
    delete [] v3;

    return 0;
}

int main(int lb, char** par)
{
    if (lb < 29)
    {
	printf("Args: %d < 29\nArgs given and in [] example values\n", lb);
	printf("file_out[out.eps] width[600] height[600]\n");
	printf("start_z_real[-2] end_z_real[2] start_z_imag[-2] end_z_imag[2]\n");
	printf("grid_lines_real[10] grid_lines_imag[10]\n");
	printf("min_real_val[-100] min_imag_val[-100] min_modulo_val[0]\n");
	printf("max_real_val[100] max_imag_val[100] max_modulo_val[200]\n");
	printf("level_real_from[-5] level_real_to[5] level_real_inc[1]\n");
	printf("level_imag_from[-5] level_imag_to[5] level_imag_inc[1]\n");
	printf("level_mod_from[0] level_mod_to[10] level_mod_inc[1]\n");
	printf("draw_real[1] draw_imag[1] draw_modulo[1] func_def['cos(x)+i']\n");
	return 1;
    }
		
    fpar_function(par[28]);
    if (!fpar_ok())
    {
	printf("\"%s\" - function definition parse error\n", par[28]);
	return 2;
    }

    /*
    complex<double> z, fz;
    for (double r=-2.;r<=2.;r+=1.)
    for (double i=-2.;i<=2.;i+=1.)
    {
	z = complex<double>(r, i);
	fz = fpar_f(z);
	printf("f(%lf,%lfi) = %lf,%lfi\n", Re(z), Im(z), Re(fz), Im(fz));
    }
    return 0;
    */

    eps_cont(
	    	par[1], atoi(par[2]), atoi(par[3]), 			// fn, pixels_w, pixels_h
		atof(par[4]), atof(par[5]), atof(par[6]), atof(par[7]),	// sx, ex, sy, ey
		atoi(par[8]), atoi(par[9]),				// grid_x, grid_y 
		atof(par[10]), atof(par[11]), atof(par[12]),		// min_r, min_i, min_m
		atof(par[13]), atof(par[14]), atof(par[15]),		// max_r, max_i, max_m
		atof(par[16]), atof(par[17]), atof(par[18]),		// lev_r_f, lev_r_t, lev_r_d 
		atof(par[19]), atof(par[20]), atof(par[21]),		// lev_i_f, lev_i_t, lev_i_d 
		atof(par[22]), atof(par[23]), atof(par[24]),		// lev_m_f, lev_m_t, lev_m_d 
		atoi(par[25]), atoi(par[26]), atoi(par[27]), 		// draw_r, draw_i, draw_m
		fpar_f);
    return 0;
}


