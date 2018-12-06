#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void diffr( double * x, double * xr, int m, int n );
void diffs( double * x, double * xr, int m, int n );
void exportGrid(double * x, double * y, int m, int n);
void exportSolu(double * x, double * y, double * u, int m, int n, int step);

int main()
{
	int i, j, m, n, ind, st, maxsteps;
	double *x, *y, *xr, *yr, *xs, *ys, *u, *res, *det, *ur, *us, *ux, *uy;
	double r, s, alpha, pi, dt;
	m = 30;
	n = 100;
	alpha = 1;
	maxsteps = 150;
	pi = 3.1415926535;
	dt = 0.001;

	/* Allocate memory */
	x = (double *)malloc( m*n*sizeof(double));
	y = (double *)malloc( m*n*sizeof(double));
	xr = (double *)malloc( m*n*sizeof(double));
	yr = (double *)malloc( m*n*sizeof(double));
	xs = (double *)malloc( m*n*sizeof(double));
	ys = (double *)malloc( m*n*sizeof(double));
	u = (double *)malloc( m*n*sizeof(double));
	res = (double *)malloc( m*n*sizeof(double));
	det = (double *)malloc( m*n*sizeof(double));
	ur = (double *)malloc( m*n*sizeof(double));
	us = (double *)malloc( m*n*sizeof(double));
	ux = (double *)malloc( m*n*sizeof(double));
	uy = (double *)malloc( m*n*sizeof(double));

	/* Generate grid */
	for( j= 0 ; j<n ; j++ )
	for( i= 0 ; i<m ; i++ )
	{
		ind = i + m*j;
		r = (i)/(m-1.0);
		s = (j)/(n-4.0);
		x[ind] = (1+2*r)*cos(2*pi*s);
		y[ind] = (1+2*r)*sin(2*pi*s);
	}

	/* Compute metric */
	diffr( x, xr, m, n );
	diffr( y, yr, m, n );
	diffs( x, xs, m, n );
	diffs( y, ys, m, n );

	for( ind = 0 ; ind < m*n ; ind++ )
	det[ind] = xr[ind]*ys[ind]-xs[ind]*yr[ind];

	/* Get initial data */
	for( ind = 0 ; ind < m*n ; ind++ )
	u[ind] = 1 + exp(-alpha*( (x[ind]-2)*(x[ind]-2) + y[ind]*y[ind] ) );

	/* Export grid */
	exportGrid(x, y, m, n);

	/* Main loop for time stepping */
	for( st = 0 ; st < maxsteps ; st++ )
	{
		/* Compute the second derivatives */
		diffr( u, ur, m, n );
		diffs( u, us, m, n );

		for( ind = 0 ; ind < m*n ; ind++ )
		{
			ux[ind] =( ur[ind]*ys[ind]-us[ind]*yr[ind])/det[ind];
			uy[ind] =(-ur[ind]*xs[ind]+us[ind]*xr[ind])/det[ind];
		}

		diffr( ux, ur, m, n );
		diffs( ux, us, m, n );

		for( ind = 0 ; ind < m*n ; ind++ )
		res[ind] =( ur[ind]*ys[ind]-us[ind]*yr[ind])/det[ind];

		diffr( uy, ur, m, n );
		diffs( uy, us, m, n );

		for( ind = 0 ; ind < m*n ; ind++ )
		u[ind] = u[ind] + dt*( res[ind] +
			( -ur[ind]*xs[ind]+us[ind]*xr[ind])/det[ind] );

		/* Inner and outer boundaries */
		for( j=0 ; j < n ; j++ )
		{
			u[m*j] = 2;
			u[m-1+m*j] = 1;
		}

		/* Periodic boundary */
		for( i=0 ; i < m ; i++ )
		{
			u[ i ] = u[i + m*(n-4) ];
			u[ i +m ] = u[i + m*(n-3) ];
			u[ i +m*(n-2) ] = u[i + m*2];
			u[ i +m*(n-1) ] = u[i + m*3];
		}

		/* Export solution */
		if( st % 10 == 0)
			exportSolu(x, y, u, m, n, st);
	}

	/* Write out solution and grid */
	//for( ind=0 ; ind < n*m ; ind++ )
	//printf("%g %g %g \n" , x[ind], y[ind], u[ind] );

	/* Should give back memory here */
}

/* ---------- Procedure to compute r-derivative --------- */
void diffr( double * x, double * xr, int m, int n )
{
	int i, j, ind;
	for( j= 0 ; j<n ; j++ )
	{
		i = 0;
		ind = i + m*j;
		xr[ind] = x[ind+1]-x[ind];

		for( i= 1 ; i<m-1 ; i++ )
		{
			ind = i + m*j;
			xr[ind] = 0.5*(x[ind+1]-x[ind-1]);
		}

		i = m-1;
		ind = i + m*j;
		xr[ind] = x[ind]-x[ind-1];
	}
}

/* ---------- Procedure to compute s-derivative --------- */
void diffs( double * x, double * xs, int m, int n )
{
	int i, j, ind;
	for( i= 0 ; i<m ; i++ )
	{
		j = 0;
		ind = i + m*j;
		xs[ind] = x[ind+m]-x[ind];

		for( j= 1 ; j<n-1 ; j++ )
		{
			ind = i + m*j;
			xs[ind] = 0.5*(x[ind+m]-x[ind-m]);
		}

		ind = n-1;
		ind = i + m*j;
		xs[ind] = x[ind]-x[ind-m];
	}
}

/* ---------- Procedure to export grid ------------------ */
void exportGrid(double * x, double * y, int m, int n)
{
	FILE *pF;	
	char str[20] = "grid.vtk";
	if ((pF = fopen(str, "w")) == NULL) 
	{
    		printf("Cannot open file \"%s\" ", str);
    		exit(-1);
	}

	fprintf(pF,"# vtk DataFile Version 2.0 \n");
  	fprintf(pF,"LES slice \n");
  	fprintf(pF,"ASCII \n");
  	fprintf(pF,"DATASET STRUCTURED_GRID \n");

  	fprintf(pF,"DIMENSIONS %d %d %d \n", m, n, 1 );
	fprintf(pF,"POINTS %d float \n", m * n);

	for( int ind=0 ; ind < n*m ; ind++ )
	fprintf(pF,"%f %f %f \n" , x[ind], y[ind], 0.0 );

	fclose( pF );
}

void exportSolu(double * x, double * y, double * u, int m, int n, int step)
{
	FILE *pF;
	char str[20];	
	sprintf(str, "%7.7d", step);
  	char str2[20] = ".vtk";
  	sprintf(str2,strcat(str,str2));
  	char str3[20] = "solution_";
  	sprintf(str2,strcat(str3,str2));
  	if ((pF = fopen(str2, "w")) == NULL) 
	{
    		printf("Cannot open file \"%s\" ", str2);
    		exit(-1);
	}

	fprintf(pF,"# vtk DataFile Version 2.0 \n");
  	fprintf(pF,"LES slice \n");
  	fprintf(pF,"ASCII \n");
  	fprintf(pF,"DATASET STRUCTURED_GRID \n");

  	fprintf(pF,"DIMENSIONS %d %d %d \n", m, n, 1 );
	fprintf(pF,"POINTS %d float \n", m * n);

	for( int ind=0 ; ind < n*m ; ind++ )
	fprintf(pF,"%f %f %f \n" , x[ind], y[ind], 0.0 );
	
	fprintf(pF,"POINT_DATA %d \n", m * n );
  	fprintf(pF,"SCALARS u float \n");
	fprintf(pF,"LOOKUP_TABLE default \n");
	for( int ind=0 ; ind < n*m ; ind++ )
	fprintf(pF,"%f \n" ,u[ind] );

	fclose( pF );
}
