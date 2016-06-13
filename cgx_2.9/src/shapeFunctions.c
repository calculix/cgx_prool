/* bug somewere */
/* dirty hag in shape10tet_(): without sprintf(buffer,"dirtyhag\n") the code does not calc xs[] (???) valgrind did not show anything */

#include <extUtil.h>

void shapeHe8(double pcoords[3],double shp[20])
{
  double  xi,et,ze,omg,omh,omr,opg,oph,opr;

  xi=pcoords[0];
  et=pcoords[1];
  ze=pcoords[2];
  //
  //     shape functions and their glocal derivatives
  //
  omg=1.0-xi;
  omh=1.0-et;
  omr=1.0-ze;
  opg=1.0+xi;
  oph=1.0+et;
  opr=1.0+ze;
  //
  //     shape functions
  //
  shp[0]=omg*omh*omr/8.0;
  shp[1]=opg*omh*omr/8.0;
  shp[2]=opg*oph*omr/8.0;
  shp[3]=omg*oph*omr/8.0;
  shp[4]=omg*omh*opr/8.0;
  shp[5]=opg*omh*opr/8.0;
  shp[6]=opg*oph*opr/8.0;
  shp[7]=omg*oph*opr/8.0;
}

void shapeW6(double pcoords[3],double shp[15])
{
  double  xi,et,ze,a;

  xi=pcoords[0];
  et=pcoords[1];
  ze=pcoords[2];
  a=1.-xi-et;
  //
  //
  //     shape functions
  //
  shp[0]=0.5*a *(1.-ze); 
  shp[1]=0.5*xi*(1.-ze); 
  shp[2]=0.5*et*(1.-ze); 
  shp[3]=0.5*a *(1.+ze); 
  shp[4]=0.5*xi*(1.+ze); 
  shp[5]=0.5*et*(1.+ze);   
}

void shapeW15(double pcoords[3],double shp[15])
{
  double  xi,et,ze,a;

  xi=pcoords[0];
  et=pcoords[1];
  ze=pcoords[2];
  a=1.-xi-et;
  //
  //
  //     shape functions
  //
  shp[0]=-0.5*a*(1.0-ze)*(2.0*xi+2.0*et+ze);
  shp[1]=0.5*xi*(1.0-ze)*(2.0*xi-2.0-ze);
  shp[2]=0.5*et*(1.0-ze)*(2.0*et-2.0-ze);
  shp[3]=-0.5*a*(1.0+ze)*(2.0*xi+2.0*et-ze);
  shp[4]=0.5*xi*(1.0+ze)*(2.0*xi-2.0+ze);
  shp[5]=0.5*et*(1.0+ze)*(2.0*et-2.0+ze);
  shp[6]=2.0*xi*a*(1.0-ze);
  shp[7]=2.0*xi*et*(1.0-ze);
  shp[8]=2.0*et*a*(1.0-ze);
  shp[9]=2.0*xi*a*(1.0+ze);
  shp[10]=2.0*xi*et*(1.0+ze); 
  shp[11]=2.0*et*a*(1.0+ze);
  shp[12]= a*(1.0-ze*ze);
  shp[13]=xi*(1.0-ze*ze);
  shp[14]=et*(1.0-ze*ze);
}

//    shape functions for a 10-node quadratic
//    isoparametric tetrahedral element. 0<=xi,et,ze<=1,xi+et+ze<=1.
void shapeTet10(double pcoords[3],double shp[20])
{
  double a;
  double xi,et,ze;

  xi=pcoords[0];
  et=pcoords[1];
  ze=pcoords[2];
  //
  //     shape functions
  //
  a=1.0-xi-et-ze;
  shp[0]=(2.0*a-1.0)*a;
  shp[1]=xi*(2.0*xi-1.0);
  shp[2]=et*(2.0*et-1.0);
  shp[3]=ze*(2.0*ze-1.0);
  shp[4]=4.0*xi*a;
  shp[5]=4.0*xi*et;
  shp[6]=4.0*et*a;
  shp[7]=4.0*ze*a;
  shp[8]=4.0*xi*ze;
  shp[9]=4.0*et*ze;
}



void shapeHe20(double pcoords[3],double weights[20])
{
  //VTK needs parametric coordinates to be between "pcoords" (0,1). Isoparametric
  //shape functions are formulated between (-1,1). Here we do a
  //coordinate system conversion from (0,1) to (-1,1).
  double r = 2.0*(pcoords[0]-0.5);
  double s = 2.0*(pcoords[1]-0.5);
  double t = 2.0*(pcoords[2]-0.5);

  double rm = 1.0 - r;
  double rp = 1.0 + r;
  double sm = 1.0 - s;
  double sp = 1.0 + s;
  double tm = 1.0 - t;
  double tp = 1.0 + t;
  double r2 = 1.0 - r*r;
  double s2 = 1.0 - s*s;
  double t2 = 1.0 - t*t;

  //The eight corner points
  weights[0] = 0.125 * rm * sm * tm * (-r - s - t - 2.0);
  weights[1] = 0.125 * rp * sm * tm * ( r - s - t - 2.0);
  weights[2] = 0.125 * rp * sp * tm * ( r + s - t - 2.0);
  weights[3] = 0.125 * rm * sp * tm * (-r + s - t - 2.0);
  weights[4] = 0.125 * rm * sm * tp * (-r - s + t - 2.0);
  weights[5] = 0.125 * rp * sm * tp * ( r - s + t - 2.0);
  weights[6] = 0.125 * rp * sp * tp * ( r + s + t - 2.0);
  weights[7] = 0.125 * rm * sp * tp * (-r + s + t - 2.0);

  //The mid-edge nodes
  weights[8] =  0.25 * r2 * sm * tm;
  weights[9] =  0.25 * s2 * rp * tm;
  weights[10] = 0.25 * r2 * sp * tm;
  weights[11] = 0.25 * s2 * rm * tm;
  weights[12] = 0.25 * r2 * sm * tp;
  weights[13] = 0.25 * s2 * rp * tp;
  weights[14] = 0.25 * r2 * sp * tp;
  weights[15] = 0.25 * s2 * rm * tp;
  weights[16] = 0.25 * t2 * rm * sm;
  weights[17] = 0.25 * t2 * rp * sm;
  weights[18] = 0.25 * t2 * rp * sp;
  weights[19] = 0.25 * t2 * rm * sp;
}




void shape4q(double xi, double et, double *xl, double *xsj)
{
  int i,j,k;
  double xs[3][2], shp[4][8];

      shp[0][0]=-(1.-et)/4.;
      shp[0][1]=(1.-et)/4.;
      shp[0][2]=(1.+et)/4.;
      shp[0][3]=-(1.+et)/4.;

      shp[1][0]=-(1.-xi)/4.;
      shp[1][1]=-(1.+xi)/4.;
      shp[1][2]=(1.+xi)/4.;
      shp[1][3]=(1.-xi)/4.;
/*
!     computation of the local derivative of the global coordinates
!     (xs)
*/
      for (i=0; i<3; i++)
      {
        for (j=0; j<2; j++)
        {
          xs[i][j]=0.;
          for (k=0; k<4; k++)  {
	    //  printf(" i:%d j:%d k:%d xs:%f xl:%f shp:%f \n",i,j,k, xs[i][j], xl[], shp[j][k] );
	    xs[i][j]=xs[i][j]+xl[k*3+i]*shp[j][k];
          }
	}
      }
/*
!     computation of the jacobian vector
*/
      xsj[0]=xs[1][0]*xs[2][1]-xs[2][0]*xs[1][1];
      xsj[1]=xs[0][1]*xs[2][0]-xs[2][1]*xs[0][0];
      xsj[2]=xs[0][0]*xs[1][1]-xs[1][0]*xs[0][1];
}



void shape6tri(double xi, double et, double *xl, double *xsj)
{
  int i,j,k;
  double xs[3][2], shp[4][8];

      shp[0][0]=4.*(xi+et)-3.;
      shp[0][1]=4.*xi-1.;
      shp[0][2]=0.;
      shp[0][3]=4.*(1.-2.*xi-et);
      shp[0][4]=4.*et;
      shp[0][5]=-4.*et;

      shp[1][0]=4.*(xi+et)-3.;
      shp[1][1]=0.;
      shp[1][2]=4.*et-1.;
      shp[1][3]=-4.*xi;
      shp[1][4]=4.*xi;
      shp[1][5]=4.*(1.-xi-2.*et);
/*
!     computation of the local derivative of the global coordinates
!     (xs)
*/
      for (i=0; i<3; i++)
      {
        for (j=0; j<2; j++)
        {
          xs[i][j]=0.;
          for (k=0; k<6; k++)  {
	    //  printf(" i:%d j:%d k:%d xs:%f xl:%f shp:%f \n",i,j,k, xs[i][j], xl[], shp[j][k] );
	    xs[i][j]=xs[i][j]+xl[k*3+i]*shp[j][k];
          }
	}
      }
/*
!     computation of the jacobian vector
*/
      xsj[0]=xs[1][0]*xs[2][1]-xs[2][0]*xs[1][1];
      xsj[1]=xs[0][1]*xs[2][0]-xs[2][1]*xs[0][0];
      xsj[2]=xs[0][0]*xs[1][1]-xs[1][0]*xs[0][1];
}



void shape8q(double xi, double et, double *xl, double *xsj)
{
  int i,j,k;
  double xs[3][2], shp[4][8];

/*
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
!     local derivatives of the shape functions: xi-derivative
      xl local coordinates of the element nodes [nr][xyz]
*/
      shp[0][0]=(1.-et)*(2.*xi+et)/4.;
      shp[0][1]=(1.-et)*(2.*xi-et)/4.;
      shp[0][2]=(1.+et)*(2.*xi+et)/4.;
      shp[0][3]=(1.+et)*(2.*xi-et)/4.;
      shp[0][4]=-xi*(1.-et);
      shp[0][5]=(1.-et*et)/2.;
      shp[0][6]=-xi*(1.+et);
      shp[0][7]=-(1.-et*et)/2.;
/*
!     local derivatives of the shape functions: eta-derivative
*/
      shp[1][0]=(1.-xi)*(2.*et+xi)/4.;
      shp[1][1]=(1.+xi)*(2.*et-xi)/4.;
      shp[1][2]=(1.+xi)*(2.*et+xi)/4.;
      shp[1][3]=(1.-xi)*(2.*et-xi)/4.;
      shp[1][4]=-(1.-xi*xi)/2.;
      shp[1][5]=-et*(1.+xi);
      shp[1][6]=(1.-xi*xi)/2.;
      shp[1][7]=-et*(1.-xi);
/*
!     computation of the local derivative of the global coordinates
!     (xs)
*/
      for (i=0; i<3; i++)
      {
        for (j=0; j<2; j++)
        {
          xs[i][j]=0.;
          for (k=0; k<8; k++)  {
	    //  printf(" i:%d j:%d k:%d xs:%f xl:%f shp:%f \n",i,j,k, xs[i][j], xl[], shp[j][k] );
	    xs[i][j]=xs[i][j]+xl[k*3+i]*shp[j][k];
          }
	}
      }
/*
!     computation of the jacobian vector
*/
      xsj[0]=xs[1][0]*xs[2][1]-xs[2][0]*xs[1][1];
      xsj[1]=xs[0][1]*xs[2][0]-xs[2][1]*xs[0][0];
      xsj[2]=xs[0][0]*xs[1][1]-xs[1][0]*xs[0][1];
}


/* shape8h.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/* Subroutine */ 
int shape8h_(xi, et, ze, xl, xsj)
double *xi, *et, *ze, *xl, *xsj;
{
    static int i__, j, k;
    static double xs[9]	/* was [3][3] */, shp[60]	/* was [3][20]
	     */;


/*     shape functions and derivatives for a 8-node linear isoparametric */
/*     solid element */





/*     shape functions and their glocal derivatives */

/*     local derivatives of the shape functions: xi-derivative */

    /* Parameter adjustments */
    xl -= 4;

    /* Function Body */
    shp[0] = -(1. - *et) * (*ze + 1.) / 8.;
    shp[3] = (1. - *et) * (*ze + 1.) / 8.;
    shp[6] = (1. - *et) * (1. - *ze) / 8.;
    shp[9] = -(1. - *et) * (1. - *ze) / 8.;
    shp[12] = -(*et + 1.) * (*ze + 1.) / 8.;
    shp[15] = (*et + 1.) * (*ze + 1.) / 8.;
    shp[18] = (*et + 1.) * (1. - *ze) / 8.;
    shp[21] = -(*et + 1.) * (1. - *ze) / 8.;

/*     local derivatives of the shape functions: eta-derivative */

    shp[1] = -(1. - *xi) * (*ze + 1.) / 8.;
    shp[4] = -(*xi + 1.) * (*ze + 1.) / 8.;
    shp[7] = -(*xi + 1.) * (1. - *ze) / 8.;
    shp[10] = -(1. - *xi) * (1. - *ze) / 8.;
    shp[13] = (1. - *xi) * (*ze + 1.) / 8.;
    shp[16] = (*xi + 1.) * (*ze + 1.) / 8.;
    shp[19] = (*xi + 1.) * (1. - *ze) / 8.;
    shp[22] = (1. - *xi) * (1. - *ze) / 8.;

/*     local derivatives of the shape functions: zeta-derivative */

    shp[2] = (1. - *xi) * (1. - *et) / 8.;
    shp[5] = (*xi + 1.) * (1. - *et) / 8.;
    shp[8] = -(*xi + 1.) * (1. - *et) / 8.;
    shp[11] = -(1. - *xi) * (1. - *et) / 8.;
    shp[14] = (1. - *xi) * (*et + 1.) / 8.;
    shp[17] = (*xi + 1.) * (*et + 1.) / 8.;
    shp[20] = -(*xi + 1.) * (*et + 1.) / 8.;
    shp[23] = -(1. - *xi) * (*et + 1.) / 8.;

/*     computation of the local derivative of the global coordinates */
/*     (xs) */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    xs[i__ + j * 3 - 4] = 0.;
	    for (k = 1; k <= 8; ++k) {
		xs[i__ + j * 3 - 4] += xl[i__ + k * 3] * shp[j + k * 3 - 4];
	    }
	}
    }

/*     computation of the jacobian determinant */

    *xsj = xs[0] * (xs[4] * xs[8] - xs[7] * xs[5]) - xs[3] * (xs[1] * xs[8] - 
	    xs[7] * xs[2]) + xs[6] * (xs[1] * xs[5] - xs[4] * xs[2]);

    return(0);
} /* shape8h_ */




/* shape20h.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/* Subroutine */ 
int shape20h_(xi, et, ze, xl, xsj)
double *xi, *et, *ze, *xl, *xsj;
{
    static int i__, j, k;
    static double xs[9]	/* was [3][3] */, shp[60]	/* was [3][20]
	     */;


/*     shape functions and derivatives for a 20-node quadratic */
/*     isoparametric brick element. Watch out: 0<=xi<=2, -1<=et,ze<=1 !!! */





/*     shape functions and their glocal derivatives */

/*     local derivatives of the shape functions: xi-derivative */

    /* Parameter adjustments */
    xl -= 4;

    /* Function Body */
    shp[0] = *et * (*et * ((double)-.125 - *ze * (double).125) + (double).25 + *
	    ze * (*ze * (double).125 + (double).375)) - (double).125 + *ze * ((
	    double)-.25 - *ze * (double).125) + *xi * (*et * ((double)-.25 - *ze 
	    * (double).25) + (double).25 + *ze * (double).25);
    shp[3] = *et * (*et * (*ze * (double).125 + (double).125) + (double).25 + *
	    ze * ((double).125 - *ze * (double).125)) - (double).375 + *ze * (*
	    ze * (double).125 - (double).25) + *xi * (*et * ((double)-.25 - *ze *
	     (double).25) + (double).25 + *ze * (double).25);
    shp[6] = *et * (*et * ((double).125 - *ze * (double).125) + (double).25 + *
	    ze * ((double)-.125 - *ze * (double).125)) - (double).375 + *ze * (*
	    ze * (double).125 + (double).25) + *xi * (*et * (*ze * (double).25 - 
	    (double).25) + (double).25 - *ze * (double).25);
    shp[9] = *et * (*et * (*ze * (double).125 - (double).125) + (double).25 + *
	    ze * (*ze * (double).125 - (double).375)) - (double).125 + *ze * ((
	    double).25 - *ze * (double).125) + *xi * (*et * (*ze * (double).25 - 
	    (double).25) + (double).25 - *ze * (double).25);
    shp[12] = *et * (*et * ((double)-.125 - *ze * (double).125) - (double).25 + *
	    ze * ((double)-.375 - *ze * (double).125)) - (double).125 + *ze * ((
	    double)-.25 - *ze * (double).125) + *xi * (*et * (*ze * (double).25 
	    + (double).25) + (double).25 + *ze * (double).25);
    shp[15] = *et * (*et * (*ze * (double).125 + (double).125) - (double).25 + *
	    ze * (*ze * (double).125 - (double).125)) - (double).375 + *ze * (*
	    ze * (double).125 - (double).25) + *xi * (*et * (*ze * (double).25 + 
	    (double).25) + (double).25 + *ze * (double).25);
    shp[18] = *et * (*et * ((double).125 - *ze * (double).125) - (double).25 + *
	    ze * (*ze * (double).125 + (double).125)) - (double).375 + *ze * (*
	    ze * (double).125 + (double).25) + *xi * (*et * ((double).25 - *ze * 
	    (double).25) + (double).25 - *ze * (double).25);
    shp[21] = *et * (*et * (*ze * (double).125 - (double).125) - (double).25 + *
	    ze * ((double).375 - *ze * (double).125)) - (double).125 + *ze * ((
	    double).25 - *ze * (double).125) + *xi * (*et * ((double).25 - *ze * 
	    (double).25) + (double).25 - *ze * (double).25);
    shp[24] = *et * ((double)-.5 - *ze * (double).5) + (double).5 + *ze * (double)
	    .5 + *xi * (*et * (*ze * (double).5 + (double).5) - (double).5 - *ze 
	    * (double).5);
    shp[27] = *et * (*ze * (*ze * (double).25) - (double).25) + (double).25 + *
	    ze * (-(*ze) * (double).25);
    shp[30] = *et * (*ze * (double).5 - (double).5) + (double).5 - *ze * (double)
	    .5 + *xi * (*et * ((double).5 - *ze * (double).5) - (double).5 + *ze 
	    * (double).5);
    shp[33] = *et * (*ze * (-(*ze) * (double).25) + (double).25) - (double).25 + 
	    *ze * (*ze * (double).25);
    shp[36] = *et * (*ze * (double).5 + (double).5) + (double).5 + *ze * (double)
	    .5 + *xi * (*et * ((double)-.5 - *ze * (double).5) - (double).5 - *
	    ze * (double).5);
    shp[39] = *et * (*ze * (-(*ze) * (double).25) + (double).25) + (double).25 + 
	    *ze * (-(*ze) * (double).25);
    shp[42] = *et * ((double).5 - *ze * (double).5) + (double).5 - *ze * (double)
	    .5 + *xi * (*et * (*ze * (double).5 - (double).5) - (double).5 + *ze 
	    * (double).5);
    shp[45] = *et * (*ze * (*ze * (double).25) - (double).25) - (double).25 + *
	    ze * (*ze * (double).25);
    shp[48] = *et * (*et * (*ze * (double).25 + (double).25)) - (double).25 - *
	    ze * (double).25;
    shp[51] = *et * (*et * ((double)-.25 - *ze * (double).25)) + (double).25 + *
	    ze * (double).25;
    shp[54] = *et * (*et * (*ze * (double).25 - (double).25)) + (double).25 - *
	    ze * (double).25;
    shp[57] = *et * (*et * ((double).25 - *ze * (double).25)) - (double).25 + *
	    ze * (double).25;

/*     local derivatives of the shape functions: eta-derivative */

    shp[1] = *et * (*ze * (double).5 + (double).5) + *ze * ((double)-.25 - *ze * 
	    (double).25) + *xi * (*et * ((double)-.25 - *ze * (double).25) + (
	    double).25 + *ze * (*ze * (double).125 + (double).375) + *xi * ((
	    double)-.125 - *ze * (double).125));
    shp[4] = *xi * (*et * (*ze * (double).25 + (double).25) + (double).25 + *ze *
	     ((double).125 - *ze * (double).125) + *xi * ((double)-.125 - *ze * (
	    double).125));
    shp[7] = *xi * (*et * ((double).25 - *ze * (double).25) + (double).25 + *ze *
	     ((double)-.125 - *ze * (double).125) + *xi * (*ze * (double).125 - (
	    double).125));
    shp[10] = *et * ((double).5 - *ze * (double).5) + *ze * ((double).25 - *ze * 
	    (double).25) + *xi * (*et * (*ze * (double).25 - (double).25) + (
	    double).25 + *ze * (*ze * (double).125 - (double).375) + *xi * (*ze *
	     (double).125 - (double).125));
    shp[13] = *et * (*ze * (double).5 + (double).5) + *ze * (*ze * (double).25 + 
	    (double).25) + *xi * (*et * ((double)-.25 - *ze * (double).25) - (
	    double).25 + *ze * ((double)-.375 - *ze * (double).125) + *xi * (*ze 
	    * (double).125 + (double).125));
    shp[16] = *xi * (*et * (*ze * (double).25 + (double).25) - (double).25 + *ze 
	    * (*ze * (double).125 - (double).125) + *xi * (*ze * (double).125 + (
	    double).125));
    shp[19] = *xi * (*et * ((double).25 - *ze * (double).25) - (double).25 + *ze 
	    * (*ze * (double).125 + (double).125) + *xi * ((double).125 - *ze * (
	    double).125));
    shp[22] = *et * ((double).5 - *ze * (double).5) + *ze * (*ze * (double).25 - 
	    (double).25) + *xi * (*et * (*ze * (double).25 - (double).25) - (
	    double).25 + *ze * ((double).375 - *ze * (double).125) + *xi * ((
	    double).125 - *ze * (double).125));
    shp[25] = *xi * ((double)-.5 - *ze * (double).5 + *xi * (*ze * (double).25 + 
	    (double).25));
    shp[28] = *xi * (*ze * (*ze * (double).25) - (double).25);
    shp[31] = *xi * (*ze * (double).5 - (double).5 + *xi * ((double).25 - *ze * (
	    double).25));
    shp[34] = *ze * (*ze * (double).5) - (double).5 + *xi * (*ze * (-(*ze) * (
	    double).25) + (double).25);
    shp[37] = *xi * (*ze * (double).5 + (double).5 + *xi * ((double)-.25 - *ze * 
	    (double).25));
    shp[40] = *xi * (*ze * (-(*ze) * (double).25) + (double).25);
    shp[43] = *xi * ((double).5 - *ze * (double).5 + *xi * (*ze * (double).25 - (
	    double).25));
    shp[46] = *ze * (-(*ze) * (double).5) + (double).5 + *xi * (*ze * (*ze * (
	    double).25) - (double).25);
    shp[49] = *et * (*ze * (double)-1. - (double)1.) + *xi * (*et * (*ze * (
	    double).5 + (double).5));
    shp[52] = *xi * (*et * ((double)-.5 - *ze * (double).5));
    shp[55] = *xi * (*et * (*ze * (double).5 - (double).5));
    shp[58] = *et * (*ze * (double)1. - (double)1.) + *xi * (*et * ((double).5 - 
	    *ze * (double).5));

/*     local derivatives of the shape functions: zeta-derivative */

    shp[2] = *et * (*et * (double).25 - (double).25 - *ze * (double).5) + *ze * (
	    double).5 + *xi * (*et * ((double).375 - *et * (double).125 + *ze * (
	    double).25) - (double).25 - *ze * (double).25 + *xi * (*et * (double)
	    -.125 + (double).125));
    shp[5] = *xi * (*et * (*et * (double).125 + (double).125 - *ze * (double).25)
	     - (double).25 + *ze * (double).25 + *xi * (*et * (double)-.125 + (
	    double).125));
    shp[8] = *xi * (*et * ((double)-.125 - *et * (double).125 - *ze * (double)
	    .25) + (double).25 + *ze * (double).25 + *xi * (*et * (double).125 - 
	    (double).125));
    shp[11] = *et * ((double).25 - *et * (double).25 - *ze * (double).5) + *ze * 
	    (double).5 + *xi * (*et * (*et * (double).125 - (double).375 + *ze * 
	    (double).25) + (double).25 - *ze * (double).25 + *xi * (*et * (double)
	    .125 - (double).125));
    shp[14] = *et * (*et * (double).25 + (double).25 + *ze * (double).5) + *ze * 
	    (double).5 + *xi * (*et * ((double)-.375 - *et * (double).125 - *ze *
	     (double).25) - (double).25 - *ze * (double).25 + *xi * (*et * (
	    double).125 + (double).125));
    shp[17] = *xi * (*et * (*et * (double).125 - (double).125 + *ze * (double)
	    .25) - (double).25 + *ze * (double).25 + *xi * (*et * (double).125 + 
	    (double).125));
    shp[20] = *xi * (*et * ((double).125 - *et * (double).125 + *ze * (double)
	    .25) + (double).25 + *ze * (double).25 + *xi * (*et * (double)-.125 
	    - (double).125));
    shp[23] = *et * ((double)-.25 - *et * (double).25 + *ze * (double).5) + *ze *
	     (double).5 + *xi * (*et * (*et * (double).125 + (double).375 - *ze *
	     (double).25) + (double).25 - *ze * (double).25 + *xi * (*et * (
	    double)-.125 - (double).125));
    shp[26] = *xi * (*et * (double)-.5 + (double).5 + *xi * (*et * (double).25 - 
	    (double).25));
    shp[29] = *xi * (*et * (*ze * (double).5) - *ze * (double).5);
    shp[32] = *xi * (*et * (double).5 - (double).5 + *xi * (*et * (double)-.25 + 
	    (double).25));
    shp[35] = *et * (*ze * (double)1.) + *ze * (double)-1. + *xi * (*et * (*ze *
	     (double)-.5) + *ze * (double).5);
    shp[38] = *xi * (*et * (double).5 + (double).5 + *xi * (*et * (double)-.25 - 
	    (double).25));
    shp[41] = *xi * (*et * (*ze * (double)-.5) - *ze * (double).5);
    shp[44] = *xi * (*et * (double)-.5 - (double).5 + *xi * (*et * (double).25 + 
	    (double).25));
    shp[47] = *et * (*ze * (double)-1.) + *ze * (double)-1. + *xi * (*et * (*ze 
	    * (double).5) + *ze * (double).5);
    shp[50] = *et * (*et * (double)-.5) + (double).5 + *xi * (*et * (*et * (
	    double).25) - (double).25);
    shp[53] = *xi * (*et * (*et * (double)-.25) + (double).25);
    shp[56] = *xi * (*et * (*et * (double).25) - (double).25);
    shp[59] = *et * (*et * (double).5) - (double).5 + *xi * (*et * (*et * (
	    double)-.25) + (double).25);

/*     computation of the local derivative of the global coordinates */
/*     (xs) */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    xs[i__ + j * 3 - 4] = 0.;
	    for (k = 1; k <= 20; ++k) {
		xs[i__ + j * 3 - 4] += xl[i__ + k * 3] * shp[j + k * 3 - 4];
	    }
	}
    }

/*     computation of the jacobian determinant */
    *xsj = xs[0] * (xs[4] * xs[8] - xs[7] * xs[5]) - xs[3] * (xs[1] * xs[8] - 
	    xs[7] * xs[2]) + xs[6] * (xs[1] * xs[5] - xs[4] * xs[2]);

    return 0;
} /* shape20h_ */




/* shape10tet.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/* Subroutine */ 
int shape10tet_(xi, et, ze, xl, xsj)
double *xi, *et, *ze, *xl, *xsj;
{
    static int i__, j, k;
    static double xs[9]	/* was [3][3] */, shp[30]	/* was [3][10]
	     */;

    char buffer[10];

/*     shape functions and derivatives for a 10-node quadratic */
/*     isoparametric tetrahedral element. 0<=xi,et,ze<=1,xi+et+ze<=1. */





/*     shape functions and their glocal derivatives */

/*     local derivatives of the shape functions: xi-derivative */

    /* Parameter adjustments */
    xl -= 4;

    /* Function Body */
    shp[0] = 1. - (1. - *xi - *et - *ze) * 4.;
    shp[3] = *xi * 4. - 1.;
    shp[6] = 0.;
    shp[9] = 0.;
    shp[12] = (1. - *xi * 2. - *et - *ze) * 4.;
    shp[15] = *et * 4.;
    shp[18] = *et * -4.;
    shp[21] = *ze * -4.;
    shp[24] = *ze * 4.;
    shp[27] = 0.;

/*     local derivatives of the shape functions: eta-derivative */

    shp[1] = 1. - (1. - *xi - *et - *ze) * 4.;
    shp[4] = 0.;
    shp[7] = *et * 4. - 1.;
    shp[10] = 0.;
    shp[13] = *xi * -4.;
    shp[16] = *xi * 4.;
    shp[19] = (1. - *xi - *et * 2. - *ze) * 4.;
    shp[22] = *ze * -4.;
    shp[25] = 0.;
    shp[28] = *ze * 4.;

/*     local derivatives of the shape functions: zeta-derivative */

    shp[2] = 1. - (1. - *xi - *et - *ze) * 4.;
    shp[5] = 0.;
    shp[8] = 0.;
    shp[11] = *ze * 4. - 1.;
    shp[14] = *xi * -4.;
    shp[17] = 0.;
    shp[20] = *et * -4.;
    shp[23] = (1. - *xi - *et - *ze * 2.) * 4.;
    shp[26] = *xi * 4.;
    shp[29] = *et * 4.;
	    


/*     computation of the local derivative of the global coordinates */
/*     (xs) */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	  /* dirty hag: without sprintf() the code does not calc xs[] (???) valgrind did not show anything */
	  sprintf(buffer,"dirtyhag\n");
	    xs[i__ + j * 3 - 4] = 0.;
	    for (k = 1; k <= 10; k++) {
		xs[i__ + j * 3 - 4] += xl[i__ + k * 3] * shp[j + k * 3 - 4];
		//printf("xs[%d]=%f %d\n", i__ + j * 3 - 4, xs[i__ + j * 3 - 4], i__ + k * 3);
	    }
	}
    }

/*     computation of the jacobian determinant */
    //printf("xsj:%f %f\n", xs[0],xs[8]);
    *xsj = xs[0] * (xs[4] * xs[8] - xs[7] * xs[5]) - xs[3] * (xs[1] * xs[8] - 
	    xs[7] * xs[2]) + xs[6] * (xs[1] * xs[5] - xs[4] * xs[2]);

    return 0;
} /* shape10tet_ */



/* Subroutine */ int shape15h_(xi, et, ze, xl, xsj)
double *xi, *et, *ze, *xl, *xsj;
{
    static double a;
    static int i__, j, k;
    static double sh[3], xs[9], pshp[60]	/* was [3][3] */, xsi[9]	/* 
	    was [3][3] */;
    double *shp;

    shp=&pshp[0];

/*     shape functions and derivatives for a 15-node quadratic */
/*     isoparametric wedge element. 0<=xi,et<=1,-1<=ze<=1,xi+et<=1. */

/*     iflag=1: calculate only the value of the shape functions */
/*     iflag=2: calculate the value of the shape functions and */
/*              the Jacobian determinant */
/*     iflag=3: calculate the value of the shape functions, the */
/*              value of their derivatives w.r.t. the global */
/*              coordinates and the Jacobian determinant */


/*     Copyright (c) 2003 WB */

/*     Written February 2003 on the basis of the Guido's shape function files */



/*     shape functions and their glocal derivatives */

    /* Parameter adjustments */
    shp -= 5;
    xl -= 4;

    /* Function Body */
    a = 1. - *xi - *et;

/*     shape functions */

    shp[8] = a * (float)-.5 * ((float)1. - *ze) * (*xi * (float)2. + *et * (
	    float)2. + *ze);
    shp[12] = *xi * (float).5 * ((float)1. - *ze) * (*xi * (float)2. - (float)
	    2. - *ze);
    shp[16] = *et * (float).5 * ((float)1. - *ze) * (*et * (float)2. - (float)
	    2. - *ze);
    shp[20] = a * (float)-.5 * (*ze + (float)1.) * (*xi * (float)2. + *et * (
	    float)2. - *ze);
    shp[24] = *xi * (float).5 * (*ze + (float)1.) * (*xi * (float)2. - (float)
	    2. + *ze);
    shp[28] = *et * (float).5 * (*ze + (float)1.) * (*et * (float)2. - (float)
	    2. + *ze);
    shp[32] = *xi * (float)2. * a * ((float)1. - *ze);
    shp[36] = *xi * (float)2. * *et * ((float)1. - *ze);
    shp[40] = *et * (float)2. * a * ((float)1. - *ze);
    shp[44] = *xi * (float)2. * a * (*ze + (float)1.);
    shp[48] = *xi * (float)2. * *et * (*ze + (float)1.);
    shp[52] = *et * (float)2. * a * (*ze + (float)1.);
    shp[56] = a * ((float)1. - *ze * *ze);
    shp[60] = *xi * ((float)1. - *ze * *ze);
    shp[64] = *et * ((float)1. - *ze * *ze);


/*     local derivatives of the shape functions: xi-derivative */

    shp[5] = ((float)1. - *ze) * (float).5 * (*xi * (float)4. + *et * (float)
	    4. + *ze - (float)2.);
    shp[9] = ((float)1. - *ze) * (float).5 * (*xi * (float)4. - *ze - (float)
	    2.);
    shp[13] = 0.;
    shp[17] = (*ze + (float)1.) * (float).5 * (*xi * (float)4. + *et * (float)
	    4. - *ze - (float)2.);
    shp[21] = (*ze + (float)1.) * (float).5 * (*xi * (float)4. + *ze - (float)
	    2.);
    shp[25] = 0.;
    shp[29] = ((float)1. - *ze) * (float)2. * ((float)1. - *xi * (float)2. - *
	    et);
    shp[33] = *et * (float)2. * ((float)1. - *ze);
    shp[37] = *et * (float)-2. * ((float)1. - *ze);
    shp[41] = (*ze + (float)1.) * (float)2. * ((float)1. - *xi * (float)2. - *
	    et);
    shp[45] = *et * (float)2. * (*ze + (float)1.);
    shp[49] = *et * (float)-2. * (*ze + (float)1.);
    shp[53] = -((float)1. - *ze * *ze);
    shp[57] = (float)1. - *ze * *ze;
    shp[61] = 0.;

/*     local derivatives of the shape functions: eta-derivative */

    shp[6] = ((float)1. - *ze) * (float).5 * (*xi * (float)4. + *et * (float)
	    4. + *ze - (float)2.);
    shp[10] = 0.;
    shp[14] = ((float)1. - *ze) * (float).5 * (*et * (float)4. - *ze - (float)
	    2.);
    shp[18] = (*ze + (float)1.) * (float).5 * (*xi * (float)4. + *et * (float)
	    4. - *ze - (float)2.);
    shp[22] = 0.;
    shp[26] = (*ze + (float)1.) * (float).5 * (*et * (float)4. + *ze - (float)
	    2.);
    shp[30] = *xi * (float)-2. * ((float)1. - *ze);
    shp[34] = *xi * (float)2. * ((float)1. - *ze);
    shp[38] = ((float)1. - *ze) * (float)2. * ((float)1. - *xi - *et * (float)
	    2.);
    shp[42] = *xi * (float)-2. * (*ze + (float)1.);
    shp[46] = *xi * (float)2. * (*ze + (float)1.);
    shp[50] = (*ze + (float)1.) * (float)2. * ((float)1. - *xi - *et * (float)
	    2.);
    shp[54] = -((float)1. - *ze * *ze);
    shp[58] = 0.;
    shp[62] = (float)1. - *ze * *ze;

/*     local derivatives of the shape functions: zeta-derivative */

    shp[7] = a * (*xi + *et + *ze - (float).5);
    shp[11] = *xi * (-(*xi) + *ze + (float).5);
    shp[15] = *et * (-(*et) + *ze + (float).5);
    shp[19] = a * (-(*xi) - *et + *ze + (float).5);
    shp[23] = *xi * (*xi + *ze - (float).5);
    shp[27] = *et * (*et + *ze - (float).5);
    shp[31] = *xi * -2 * a;
    shp[35] = *xi * -2 * *et;
    shp[39] = *et * -2 * a;
    shp[43] = *xi * 2 * a;
    shp[47] = *xi * 2 * *et;
    shp[51] = *et * 2 * a;
    shp[55] = a * -2 * *ze;
    shp[59] = *xi * -2 * *ze;
    shp[63] = *et * -2 * *ze;

/*     computation of the local derivative of the global coordinates */
/*     (xs) */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    xs[i__ + j * 3 - 4] = 0.;
	    for (k = 1; k <= 15; ++k) {
		xs[i__ + j * 3 - 4] += xl[i__ + k * 3] * shp[j + (k << 2)];
	    }
	}
    }

/*     computation of the jacobian determinant */

    *xsj = xs[0] * (xs[4] * xs[8] - xs[7] * xs[5]) - xs[3] * (xs[1] * xs[8] - 
	    xs[7] * xs[2]) + xs[6] * (xs[1] * xs[5] - xs[4] * xs[2]);

/*     computation of the global derivative of the local coordinates */
/*     (xsi) (inversion of xs) */

    xsi[0] = (xs[4] * xs[8] - xs[5] * xs[7]) / *xsj;
    xsi[3] = (xs[6] * xs[5] - xs[3] * xs[8]) / *xsj;
    xsi[6] = (xs[3] * xs[7] - xs[4] * xs[6]) / *xsj;
    xsi[1] = (xs[7] * xs[2] - xs[1] * xs[8]) / *xsj;
    xsi[4] = (xs[0] * xs[8] - xs[2] * xs[6]) / *xsj;
    xsi[7] = (xs[6] * xs[1] - xs[0] * xs[7]) / *xsj;
    xsi[2] = (xs[1] * xs[5] - xs[2] * xs[4]) / *xsj;
    xsi[5] = (xs[3] * xs[2] - xs[0] * xs[5]) / *xsj;
    xsi[8] = (xs[0] * xs[4] - xs[1] * xs[3]) / *xsj;

/*     computation of the global derivatives of the shape functions */

    for (k = 1; k <= 15; ++k) {
	for (j = 1; j <= 3; ++j) {
	    sh[j - 1] = shp[(k << 2) + 1] * xsi[j * 3 - 3] + shp[(k << 2) + 2]
		     * xsi[j * 3 - 2] + shp[(k << 2) + 3] * xsi[j * 3 - 1];
	}
	for (j = 1; j <= 3; ++j) {
	    shp[j + (k << 2)] = sh[j - 1];
	}
    }

    return 0;
} /* shape15w_ */


/* Subroutine */ int shape6h_(xi, et, ze, xl, xsj)
double *xi, *et, *ze, *xl, *xsj;
{
    static double a;
    static int i__, j, k;
    static double sh[3], xs[9], pshp[60]	/* was [3][3] */, xsi[9]	/* 
	    was [3][3] */;
    double *shp;

    shp=&pshp[0];


/*     shape functions and derivatives for a 6-node linear */
/*     isoparametric wedge element. 0<=xi,et<=1,xi+et<=1,-1<=ze<=1. */

/*     iflag=1: calculate only the value of the shape functions */
/*     iflag=2: calculate the value of the shape functions and */
/*              the Jacobian determinant */
/*     iflag=3: calculate the value of the shape functions, the */
/*              value of their derivatives w.r.t. the global */
/*              coordinates and the Jacobian determinant */


/*    Copyright (c) 2003 WB */

/*    Written January 2003 on the basis of the Guido's shape function files */





/*     shape functions and their glocal derivatives */

    /* Parameter adjustments */
    shp -= 5;
    xl -= 4;

    /* Function Body */
    a = 1. - *xi - *et;

/*     shape functions */

    shp[8] = a * .5 * (1. - *ze);
    shp[12] = *xi * .5 * (1. - *ze);
    shp[16] = *et * .5 * (1. - *ze);
    shp[20] = a * .5 * (*ze + 1.);
    shp[24] = *xi * .5 * (*ze + 1.);
    shp[28] = *et * .5 * (*ze + 1.);

/*     local derivatives of the shape functions: xi-derivative */

    shp[5] = (1. - *ze) * -.5;
    shp[9] = (1. - *ze) * .5;
    shp[13] = 0.;
    shp[17] = (*ze + 1.) * -.5;
    shp[21] = (*ze + 1.) * .5;
    shp[25] = 0.;

/*     local derivatives of the shape functions: eta-derivative */

    shp[6] = (1. - *ze) * -.5;
    shp[10] = 0.;
    shp[14] = (1. - *ze) * .5;
    shp[18] = (*ze + 1.) * -.5;
    shp[22] = 0.;
    shp[26] = (*ze + 1.) * .5;

/*     local derivatives of the shape functions: zeta-derivative */

    shp[7] = a * -.5;
    shp[11] = *xi * -.5;
    shp[15] = *et * -.5;
    shp[19] = a * .5;
    shp[23] = *xi * .5;
    shp[27] = *et * .5;


/*     computation of the local derivative of the global coordinates */
/*     (xs) */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    xs[i__ + j * 3 - 4] = 0.;
	    for (k = 1; k <= 6; ++k) {
		xs[i__ + j * 3 - 4] += xl[i__ + k * 3] * shp[j + (k << 2)];
	    }
	}
    }

/*     computation of the jacobian determinant */

    *xsj = xs[0] * (xs[4] * xs[8] - xs[7] * xs[5]) - xs[3] * (xs[1] * xs[8] - 
	    xs[7] * xs[2]) + xs[6] * (xs[1] * xs[5] - xs[4] * xs[2]);

/*     computation of the global derivative of the local coordinates */
/*     (xsi) (inversion of xs) */

    xsi[0] = (xs[4] * xs[8] - xs[5] * xs[7]) / *xsj;
    xsi[3] = (xs[6] * xs[5] - xs[3] * xs[8]) / *xsj;
    xsi[6] = (xs[3] * xs[7] - xs[4] * xs[6]) / *xsj;
    xsi[1] = (xs[7] * xs[2] - xs[1] * xs[8]) / *xsj;
    xsi[4] = (xs[0] * xs[8] - xs[2] * xs[6]) / *xsj;
    xsi[7] = (xs[6] * xs[1] - xs[0] * xs[7]) / *xsj;
    xsi[2] = (xs[1] * xs[5] - xs[2] * xs[4]) / *xsj;
    xsi[5] = (xs[3] * xs[2] - xs[0] * xs[5]) / *xsj;
    xsi[8] = (xs[0] * xs[4] - xs[1] * xs[3]) / *xsj;

/*     computation of the global derivatives of the shape functions */

    for (k = 1; k <= 6; ++k) {
	for (j = 1; j <= 3; ++j) {
	    sh[j - 1] = shp[(k << 2) + 1] * xsi[j * 3 - 3] + shp[(k << 2) + 2]
		     * xsi[j * 3 - 2] + shp[(k << 2) + 3] * xsi[j * 3 - 1];
	}
	for (j = 1; j <= 3; ++j) {
	    shp[j + (k << 2)] = sh[j - 1];
	}
    }

    return 0;
} /* shape6w_ */


/* e_c3d.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
 e_c3d8.f gauss.f in ccx_1.8/src
*/
/* Subroutine */ 
int e_c3d_volu_(double *xl,char *elty, double *volu, double *cg )
{
  int j;
  double weight,weights[20], gauscg[3],gausvolu, pcoords[3];

    /* xi,et,ze */
    static double gauss3d1[3]	/* was [3][1] */ = { 0.,0.,0. };
    static double gauss3d2[24]	/* was [3][8] */ = { 
	    -.577350269189626,-.577350269189626,.577350269189626,
	    .577350269189626,-.577350269189626,.577350269189626,
	    -.577350269189626,-.577350269189626,-.577350269189626,
	    .577350269189626,-.577350269189626,-.577350269189626,
	    -.577350269189626,.577350269189626,.577350269189626,
	    .577350269189626,.577350269189626,.577350269189626,
	    -.577350269189626,.577350269189626,-.577350269189626,
	    .577350269189626,.577350269189626,-.577350269189626 };
    static double gauss3d3[81]	/* was [3][27] */ = { 
	    -.774596669241483,-.774596669241483,.774596669241483,0.,
	    -.774596669241483,.774596669241483,.774596669241483,
	    -.774596669241483,.774596669241483,-.774596669241483,
	    -.774596669241483,0.,0.,-.774596669241483,0.,.774596669241483,
	    -.774596669241483,0.,-.774596669241483,-.774596669241483,
	    -.774596669241483,0.,-.774596669241483,-.774596669241483,
	    .774596669241483,-.774596669241483,-.774596669241483,
	    -.774596669241483,0.,.774596669241483,0.,0.,.774596669241483,
	    .774596669241483,0.,.774596669241483,-.774596669241483,0.,0.,0.,
	    0.,0.,.774596669241483,0.,0.,-.774596669241483,0.,
	    -.774596669241483,0.,0.,-.774596669241483,.774596669241483,0.,
	    -.774596669241483,-.774596669241483,.774596669241483,
	    .774596669241483,0.,.774596669241483,.774596669241483,
	    .774596669241483,.774596669241483,.774596669241483,
	    -.774596669241483,.774596669241483,0.,0.,.774596669241483,0.,
	    .774596669241483,.774596669241483,0.,-.774596669241483,
	    .774596669241483,-.774596669241483,0.,.774596669241483,
	    -.774596669241483,.774596669241483,.774596669241483,
	    -.774596669241483 };
    static double gauss3d4[12]	/* was [3][4] */ = { .138196601125011,
	    .138196601125011,.138196601125011,.585410196624968,
	    .138196601125011,.138196601125011,.138196601125011,
	    .585410196624968,.138196601125011,.138196601125011,
	    .138196601125011,.585410196624968 };
    static double gauss3d7[]= {
            0.333333333333333,0.333333333333333,-0.577350269189626,
            0.333333333333333,0.333333333333333,0.577350269189626 };
    static double gauss3d8[]= {
      0.166666666666667,0.166666666666667,-0.774596669241483,
      0.666666666666667,0.166666666666667,-0.774596669241483,
      0.166666666666667,0.666666666666667,-0.774596669241483,
      0.166666666666667,0.166666666666667,0.,
      0.666666666666667,0.166666666666667,0.,
      0.166666666666667,0.666666666666667,0.,
      0.166666666666667,0.166666666666667,0.774596669241483,
      0.666666666666667,0.166666666666667,0.774596669241483,
      0.166666666666667,0.666666666666667,0.774596669241483 };

    static double weight3d1[]= { 8. };

    static double weight3d2[]={1.,1.,1.,1.,1.,1.,1.,1.};

    static double weight3d3[]={
      0.171467764060357,0.274348422496571,0.171467764060357,
      0.274348422496571,0.438957475994513,0.274348422496571,
      0.171467764060357,0.274348422496571,0.171467764060357,
      0.274348422496571,0.438957475994513,0.274348422496571,
      0.438957475994513,0.702331961591221,0.438957475994513,
      0.274348422496571,0.438957475994513,0.274348422496571,
      0.171467764060357,0.274348422496571,0.171467764060357,
      0.274348422496571,0.438957475994513,0.274348422496571,
      0.171467764060357,0.274348422496571,0.171467764060357};

    //static double weight3d4[]={0.166666666666667};

    static double weight3d5[]={
      0.041666666666667,0.041666666666667,0.041666666666667,
      0.041666666666667 };
    static double weight3d7[]={ 0.5,0.5};
    static double weight3d8[]={
       0.092592592592593,0.092592592592593,0.092592592592593,
       0.148148148148148,0.148148148148148,0.148148148148148,
       0.092592592592593,0.092592592592593,0.092592592592593 };

    /* System generated locals */
    int i__1;

    /* Builtin functions */
    int s_cmp();

    /* Local variables */
    static int nope;
    //extern /* Subroutine */ int shape10tet_();
    static int mint3d, kk;
    static double et, ze, xi;
    extern /* Subroutine */ int shape8h_(), shape20h_();
    static double xsj;

/*     contains Gauss point information */

/*     gauss3d1: hex, 1-point integration (1 integration point) */
/*     gauss3d2: hex, 2-point integration (8 integration points) */
/*     gauss3d3: hex, 3-point integration (27 integration points) */
/*     gauss3d4: tet, 4 integration points */
/*     gauss3d7: wedge, 2 integration points C3D6 */
/*     gauss3d8: wedge, 9 integration points C3D15 */

    /* Parameter adjustments */
    xl -= 4;

    /* Function Body */

    /* Initialized data */
    weight=*volu=0.;
    for(j=0; j<3; j++) cg[j] = 0.;

/*     number of nodes per element */

      if((compare(elty,"C3D20R",6)==6)||(compare(elty,"C3D20",5)==5)) 
        nope=20;
      else if((compare(elty,"C3D8R",5)==5)||(compare(elty,"C3D8",4)==4))
        nope=8;
      else if (compare(elty,"C3D10",5)==5)
         nope=10;
      else if (compare(elty,"C3D15",5)==5)
         nope=15;
      else if (compare(elty,"C3D6",4)==4)
         nope=6;

/*     number of integration points per element */

    if (compare(elty,"C3D8R",5)==5)
      mint3d = 1;
    else if ((compare(elty,"C3D8",4)==4)||(compare(elty,"C3D20R",6)==6))
      mint3d = 8;
    else if (compare(elty,"C3D20",5)==5)
      mint3d = 27;
    else if (compare(elty,"C3D10",5)==5)
      mint3d = 4;
    else if (compare(elty,"C3D15",5)==5)
      mint3d = 9;
    else if (compare(elty,"C3D6",4)==4)
      mint3d = 2;
    else
      mint3d = 0;

/*     computation of the matrix: loop over the Gauss points */
    i__1 = mint3d;
    for (kk = 1; kk <= i__1; ++kk) {
	if (compare(elty,"C3D8R",5)==5)
	{
	    xi = gauss3d1[kk * 3 - 3];
	    et = gauss3d1[kk * 3 - 2];
	    ze = gauss3d1[kk * 3 - 1];
            weight=weight3d1[kk-1];
	}
        else if ((compare(elty,"C3D8",4)==4)|| (compare(elty,"C3D20R",6)==6))
        {
	    xi = gauss3d2[kk * 3 - 3];
	    if (nope == 20)
            {
		xi += 1.;
	    }
	    et = gauss3d2[kk * 3 - 2];
	    ze = gauss3d2[kk * 3 - 1];
            weight=weight3d2[kk-1];
	}
        else if (compare(elty,"C3D20",5)==5)
        {
	    xi = gauss3d3[kk * 3 - 3] + 1.;
	    et = gauss3d3[kk * 3 - 2];
	    ze = gauss3d3[kk * 3 - 1];
            weight=weight3d3[kk-1];
	}
        else if (compare(elty,"C3D10",5)==5)
        {
	    xi = gauss3d4[kk * 3 - 3];
	    et = gauss3d4[kk * 3 - 2];
	    ze = gauss3d4[kk * 3 - 1];
            weight=weight3d5[kk-1];
	}
        else if (compare(elty,"C3D15",5)==5)
        {
            xi=gauss3d8[kk * 3 - 3];
            et=gauss3d8[kk * 3 - 2];
            ze=gauss3d8[kk * 3 - 1];
            weight=weight3d8[kk-1];
	}
        else if (compare(elty,"C3D6",4)==4)
        {
            xi=gauss3d7[kk * 3 - 3];
            et=gauss3d7[kk * 3 - 2];
            ze=gauss3d7[kk * 3 - 1];
            weight=weight3d7[kk-1];
	}

	//printf("%s %f %f %f w: %f\n", elty, xi, et, ze, weight );

/*           calculation of the shape functions and their derivatives */
/*           in the gauss point */

/*           xl contains the coordinates of the nodes belonging to */
/*           the element */

	/* calc volume */
	if (nope == 20)
        {
	  //printf(" %f %f %f xl: %f %f %f   %f\n", xi, et, ze, xl[0+4],xl[1+4],xl[2+4],xsj );
	    shape20h_(&xi, &et, &ze, &xl[4], &xsj);
	    //exit(0);
	} 
        else if (nope == 8)
        {
	    shape8h_(&xi, &et, &ze, &xl[4], &xsj);
	}
        else if (nope == 10)
        {
	    shape10tet_(&xi, &et, &ze, &xl[4], &xsj);
	}
        else if(nope == 15)
        {
  	    // call shape15w(xi,et,ze,xl,xsj,shp,iflag)
	    shape15h_(&xi, &et, &ze, &xl[4], &xsj);
	}
        else
        {
	  //call shape6w(xi,et,ze,xl,xsj,shp,iflag)
	    shape6h_(&xi, &et, &ze, &xl[4], &xsj);
	}

        gausvolu=weight*xsj;
        *volu=*volu+gausvolu;
	//printf(" gausvolu %f volu %f\n", gausvolu, *volu);

        /* calc CG */
        pcoords[0]=xi*.5;
        pcoords[1]=(et+1)*.5;
        pcoords[2]=(ze+1)*.5;
	if (nope == 20) {
	  //printf(" %f %f %f  %f %f\n",pcoords[0],pcoords[1],pcoords[2], weight, xsj);
          shapeHe20(&pcoords[0], &weights[0]);
	} else if (nope == 8) {
          shapeHe8(&pcoords[0], &weights[0]);
	} else if (nope == 10) {
          shapeTet10(&pcoords[0], &weights[0]);
	} else if (nope == 15) {
          pcoords[0]=xi;
          pcoords[1]=et;
          pcoords[2]=ze;
          shapeW15(&pcoords[0], &weights[0]);
	} else {
          pcoords[0]=xi;
          pcoords[1]=et;
          pcoords[2]=ze;
          shapeW6(&pcoords[0], &weights[0]);
	}
        for(j=0; j<3; j++) gauscg[j] = 0.;
        for(j=0; j<nope; j++)
	{
          gauscg[0]+=xl[4+j*3+0] * weights[j];
          gauscg[1]+=xl[4+j*3+1] * weights[j];
          gauscg[2]+=xl[4+j*3+2] * weights[j];
	  //printf("%d  xyz: %f %f %f weights:%f\n",j, xl[4+j*3],xl[4+j*3+1],xl[4+j*3+2],weights[j]);
	}
	//printf(" gauscg %f %f %f\n",gauscg[0],gauscg[1],gauscg[2]);
        cg[0]+=gauscg[0]*gausvolu;
        cg[1]+=gauscg[1]*gausvolu;
        cg[2]+=gauscg[2]*gausvolu;

    }
    cg[0]/=*volu;
    cg[1]/=*volu;
    cg[2]/=*volu;
    //printf(" cg %f %f %f volu:%f\n",cg[0],cg[1],cg[2],*volu);

    return(1);
} /* e_c3d_volu_ */



/* e_c3d.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/* Subroutine */ 
int e_c3d__(double *xl,char *elty)
{
    /* Initialized data */

    static double gauss3d1[3]	/* was [3][1] */ = { 0.,0.,0. };
    static double gauss3d2[24]	/* was [3][8] */ = { 
	    -.577350269189626,-.577350269189626,.577350269189626,
	    .577350269189626,-.577350269189626,.577350269189626,
	    -.577350269189626,-.577350269189626,-.577350269189626,
	    .577350269189626,-.577350269189626,-.577350269189626,
	    -.577350269189626,.577350269189626,.577350269189626,
	    .577350269189626,.577350269189626,.577350269189626,
	    -.577350269189626,.577350269189626,-.577350269189626,
	    .577350269189626,.577350269189626,-.577350269189626 };
    static double gauss3d3[81]	/* was [3][27] */ = { 
	    -.774596669241483,-.774596669241483,.774596669241483,0.,
	    -.774596669241483,.774596669241483,.774596669241483,
	    -.774596669241483,.774596669241483,-.774596669241483,
	    -.774596669241483,0.,0.,-.774596669241483,0.,.774596669241483,
	    -.774596669241483,0.,-.774596669241483,-.774596669241483,
	    -.774596669241483,0.,-.774596669241483,-.774596669241483,
	    .774596669241483,-.774596669241483,-.774596669241483,
	    -.774596669241483,0.,.774596669241483,0.,0.,.774596669241483,
	    .774596669241483,0.,.774596669241483,-.774596669241483,0.,0.,0.,
	    0.,0.,.774596669241483,0.,0.,-.774596669241483,0.,
	    -.774596669241483,0.,0.,-.774596669241483,.774596669241483,0.,
	    -.774596669241483,-.774596669241483,.774596669241483,
	    .774596669241483,0.,.774596669241483,.774596669241483,
	    .774596669241483,.774596669241483,.774596669241483,
	    -.774596669241483,.774596669241483,0.,0.,.774596669241483,0.,
	    .774596669241483,.774596669241483,0.,-.774596669241483,
	    .774596669241483,-.774596669241483,0.,.774596669241483,
	    -.774596669241483,.774596669241483,.774596669241483,
	    -.774596669241483 };
    static double gauss3d4[12]	/* was [3][4] */ = { .138196601125011,
	    .138196601125011,.138196601125011,.585410196624968,
	    .138196601125011,.138196601125011,.138196601125011,
	    .585410196624968,.138196601125011,.138196601125011,
	    .138196601125011,.585410196624968 };

    /* System generated locals */
    int i__1;

    /* Builtin functions */
    int s_cmp();

    /* Local variables */
    static int nope;
    //extern /* Subroutine */ int shape10tet_();
    static int mint3d, kk;
    static double et, ze, xi;
    extern /* Subroutine */ int shape8h_(), shape20h_();
    static double xsj;

/*     contains Gauss point information */

/*     gauss3d1: hex, 1-point integration (1 integration point) */
/*     gauss3d2: hex, 2-point integration (8 integration points) */
/*     gauss3d3: hex, 3-point integration (27 integration points) */
/*     gauss3d4: tet, 4 integration points */

    /* Parameter adjustments */
    xl -= 4;

    /* Function Body */

/*     number of nodes per element */


      if((compare(elty,"C3D20R",6)==6)||(compare(elty,"C3D20",5)==5)) 
        nope=20;
      else if((compare(elty,"C3D8R",5)==5)||(compare(elty,"C3D8",4)==4))
        nope=8;
      else
         nope=10;

/*     number of integration points per element */

    if (compare(elty,"C3D8R",5)==5)
      mint3d = 1;
    else if ((compare(elty,"C3D8",4)==4)||(compare(elty,"C3D20R",6)==6))
      mint3d = 8;
    else if (compare(elty,"C3D20",5)==5)
      mint3d = 27;
    else
      mint3d = 4;

/*     computation of the matrix: loop over the Gauss points */
    i__1 = mint3d;
    for (kk = 1; kk <= i__1; ++kk) {
	if (compare(elty,"C3D8R",5)==5)
	{
	    xi = gauss3d1[kk * 3 - 3];
	    et = gauss3d1[kk * 3 - 2];
	    ze = gauss3d1[kk * 3 - 1];
	}
        else if ((compare(elty,"C3D8",4)==4)|| (compare(elty,"C3D20R",6)==6))
        {
	    xi = gauss3d2[kk * 3 - 3];
	    if (nope == 20)
            {
		xi += 1.;
	    }
	    et = gauss3d2[kk * 3 - 2];
	    ze = gauss3d2[kk * 3 - 1];
	}
        else if (compare(elty,"C3D20",5)==5)
        {
	    xi = gauss3d3[kk * 3 - 3] + 1.;
	    et = gauss3d3[kk * 3 - 2];
	    ze = gauss3d3[kk * 3 - 1];
	}
        else if (compare(elty,"C3D10",5)==5)
        {
	    xi = gauss3d4[kk * 3 - 3];
	    et = gauss3d4[kk * 3 - 2];
	    ze = gauss3d4[kk * 3 - 1];
	}
	//printf("elem coords %f %f %f %s\n", xi, et, ze,elty );

/*           calculation of the shape functions and their derivatives */
/*           in the gauss point */

/*           xl contains the coordinates of the nodes belonging to */
/*           the element */

	if (nope == 20) {
	    shape20h_(&xi, &et, &ze, &xl[4], &xsj);
	} else if (nope == 8) {
	    shape8h_(&xi, &et, &ze, &xl[4], &xsj);
	} else {
	    shape10tet_(&xi, &et, &ze, &xl[4], &xsj);
	}

/*           check the jacobian determinant */
	// printf("%s %d %d xi:%e et:%e xsj:%e\n", elty, mint3d, kk, xi,et,xsj);

	if (xsj <= 0.) {
	    return(0);
	}
    }

    return(1);
} /* e_c3d__ */


/* e_c3d.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/* Subroutine */ 
int e_c3d_nodes_(double *xl,char *elty, int *elem, double *eqal)
{
    /* Initialized data */

    static double gauss3d1[3]	/* was [3][1] */ = { 0.,0.,0. };
    static double gauss3d2[24]	/* was [3][8] */ = {
        -1.,-1.,-1.,
         1.,-1.,-1.,
        -1., 1.,-1.,
         1., 1.,-1.,
        -1.,-1., 1.,
         1.,-1., 1.,
        -1., 1., 1.,
        1., 1., 1.};
    static double gauss3d3[81]	/* was [3][27] */ = {
        -1.,-1.,-1.,
         0.,-1.,-1.,
         1.,-1.,-1.,
        -1., 0.,-1.,
         0., 0.,-1.,
         1., 0.,-1.,
        -1., 1.,-1.,
         0., 1.,-1.,
         1., 1.,-1.,
        -1.,-1., 0.,
         0.,-1., 0.,
         1.,-1., 0.,
        -1., 0., 0.,
         0., 0., 0.,
         1., 0., 0.,
        -1., 1., 0.,
         0., 1., 0.,
         1., 1., 0.,
        -1.,-1., 1.,
         0.,-1., 1.,
         1.,-1., 1.,
        -1., 0., 1.,
         0., 0., 1.,
         1., 0., 1.,
        -1., 1., 1.,
         0., 1., 1.,
         1., 1., 1.  };

    /* midside nodes not done  */
    static double gauss3d4[12]	/* was [3][4] */ = {
      0.,0.,0.,
      0.,0.,1.,
      0.,1.,0.,
      1.,0.,0. };

    /* System generated locals */
    int i__1;

    /* Builtin functions */
    int s_cmp();

    /* Local variables */
    static int nope;
    //extern /* Subroutine */ int shape10tet_();
    static int mint3d, kk;
    static double et, ze, xi;
    extern /* Subroutine */ int shape8h_(), shape20h_();
    static double xsj;
    double jacobi[27], sum_jacobi, max,min;

/*     contains Gauss point information */

/*     gauss3d1: hex, 1-point integration (1 integration point) */
/*     gauss3d2: hex, 2-point integration (8 integration points) */
/*     gauss3d3: hex, 3-point integration (27 integration points) */
/*     gauss3d4: tet, 4 integration points */

    /* Parameter adjustments */
    xl -= 4;

    /* Function Body */

/*     number of nodes per element */

      if((compare(elty,"C3D20R",6)==6)||(compare(elty,"C3D20",6)==6)) 
        nope=20;
      else if((compare(elty,"C3D8R",5)==5)||(compare(elty,"C3D8",4)==4))
        nope=8;
      else
         nope=10;

/*     number of integration points per element */

    if (compare(elty,"C3D8R",5)==5)
      mint3d = 1;
    else if ((compare(elty,"C3D8",4)==4)||(compare(elty,"C3D20R",6)==6))
      mint3d = 8;
    else if (compare(elty,"C3D20",5)==5)
      mint3d = 27;
    else
      mint3d = 4;

/*     computation of the matrix: loop over the Gauss points */

    i__1 = mint3d;
    sum_jacobi=0;
    for (kk = 1; kk <= i__1; ++kk) {
	if (compare(elty,"C3D8R",5)==5)
	{
	    xi = gauss3d1[kk * 3 - 3];
	    et = gauss3d1[kk * 3 - 2];
	    ze = gauss3d1[kk * 3 - 1];
	}
        else if ((compare(elty,"C3D8",4)==4)|| (compare(elty,"C3D20R",6)==6))
        {
	    xi = gauss3d2[kk * 3 - 3];
	    if (nope == 20)
            {
		xi += 1.;
	    }
	    et = gauss3d2[kk * 3 - 2];
	    ze = gauss3d2[kk * 3 - 1];
	}
        else if (compare(elty,"C3D20",5)==5)
        {
	    xi = gauss3d3[kk * 3 - 3] + 1.;
	    et = gauss3d3[kk * 3 - 2];
	    ze = gauss3d3[kk * 3 - 1];
	}
        else if (compare(elty,"C3D10",5)==5)
        {
	    xi = gauss3d4[kk * 3 - 3];
	    et = gauss3d4[kk * 3 - 2];
	    ze = gauss3d4[kk * 3 - 1];
	}

/*           calculation of the shape functions and their derivatives */
/*           in the gauss point */

/*           xl contains the coordinates of the nodes belonging to */
/*           the element */

	if (nope == 20) {
	    shape20h_(&xi, &et, &ze, &xl[4], &xsj);
	} else if (nope == 8) {
	    shape8h_(&xi, &et, &ze, &xl[4], &xsj);
	} else {
	    shape10tet_(&xi, &et, &ze, &xl[4], &xsj);
	}

/*           check the jacobian determinant */
        /* printf("%s %d %d xi:%e et:%e xsj:%e\n", elty, mint3d, kk, xi,et,xsj); */
	jacobi[kk-1]=xsj;
        sum_jacobi+=xsj;
    }
    /* check the element */
    /* xsj = Vworldcoordinates / Velementcoordinates */
    /* - he8/he20 == 8   */
    /* - tet10    == 1/6 */
    max=-MAX_INTEGER;
    min= MAX_INTEGER;
    for(kk=0; kk<i__1; kk++)
    {
        if(jacobi[kk]>max) max=jacobi[kk];
        if(jacobi[kk]<min) min=jacobi[kk];
        /* printf("%s %d %d jac:%e jac/sumjac:%e\n", elty, mint3d, kk, jacobi[kk], jacobi[kk]/sum_jacobi); */
    }
    if(min<=0.) return(0);
    if(max/min>*eqal)
    {
      printf("elem:%d max:%lf min:%lf max/min:%lf\n",*elem, max,min,max/min);
      return(0);
    }

    return(1);
} /* e_c3d__ */



/* neue attach funktionen (ccx_1.5) */
/* attach.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



/*     CalculiX - A 3-dimensional finite element program */
/*              Copyright (C) 1998-2005 Guido Dhondt */

/*     This program is free software; you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation(version 2); */


/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* Subroutine */ int attach_new(pneigh, pnode, nterms, ratio, dist, elemcoords)
double *pneigh, *pnode;
int *nterms;
double *ratio, *dist;
double *elemcoords;
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    static double aold[9]	/* was [3][3] */;
    static int imin, jmin;
    static double a[9]	/* was [3][3] */;
    static int i__, j;
    static double p[3], etold[9]	/* was [3][3] */, xiold[9]	/* 
	    was [3][3] */, d1, d2, d3, d4, etopt, xiopt;
    static int ii;
    static double et[9]	/* was [3][3] */, xi[9]	/* was [3][3] */;
    extern /* Subroutine */ int distattach_new();
    static double distmin;


/*     ataches node with coordinates in "pnode" to the face containing */
/*     "nterms" nodes with coordinates in field "pneigh" (nterms < 9). */




    /* Parameter adjustments */
    --ratio;
    --pnode;
    pneigh -= 4;

    /* Function Body */
    d1 = .25;
    d2 = .003125;
    d3 = .000039063;
    d4 = .0000001;

/*     initialisation */

    for (i__ = -1; i__ <= 1; ++i__) {
	for (j = -1; j <= 1; ++j) {
	    xi[i__ + j * 3 + 4] = i__ * d1;
	    et[i__ + j * 3 + 4] = j * d1;
	    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &pneigh[4]
		    , &pnode[1], &a[i__ + j * 3 + 4], p, &ratio[1], nterms);
	}
    }

/*     minimizing the distance from the face to the node */

    for (ii = 1; ii <= 1000000; ++ii) {
	distmin = a[4];
	imin = 0;
	jmin = 0;
	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (a[i__ + j * 3 + 4] < distmin) {
		    distmin = a[i__ + j * 3 + 4];
		    imin = i__;
		    jmin = j;
		}
	    }
	}

/*       exit if minimum found */

	if (imin == 0 && jmin == 0) {
	    goto L100;
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		aold[i__ + j * 3 + 4] = a[i__ + j * 3 + 4];
		xiold[i__ + j * 3 + 4] = xi[i__ + j * 3 + 4];
		etold[i__ + j * 3 + 4] = et[i__ + j * 3 + 4];
	    }
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (i__ + imin >= -1 && i__ + imin <= 1 && j + jmin >= -1 && 
			j + jmin <= 1) {
		    a[i__ + j * 3 + 4] = aold[i__ + imin + (j + jmin) * 3 + 4]
			    ;
		    xi[i__ + j * 3 + 4] = xiold[i__ + imin + (j + jmin) * 3 + 
			    4];
		    et[i__ + j * 3 + 4] = etold[i__ + imin + (j + jmin) * 3 + 
			    4];
		} else {
		    xi[i__ + j * 3 + 4] += imin * d1;
		    et[i__ + j * 3 + 4] += jmin * d1;

/* Computing MIN */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = max(d__1,-1.);

		    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &
			    pneigh[4], &pnode[1], &a[i__ + j * 3 + 4], p, &
			    ratio[1], nterms);
/*              write(*,*) a(i,j) */
		}
	    }
	}
    }
L100:

/*     2nd run */
/*     initialisation */

    xiopt = xi[4];
    etopt = et[4];
    for (i__ = -1; i__ <= 1; ++i__) {
	for (j = -1; j <= 1; ++j) {
	    xi[i__ + j * 3 + 4] = xiopt + i__ * d2;
	    et[i__ + j * 3 + 4] = etopt + j * d2;
/* Computing MIN */
	    d__1 = xi[i__ + j * 3 + 4];
	    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
	    d__1 = xi[i__ + j * 3 + 4];
	    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
	    d__1 = et[i__ + j * 3 + 4];
	    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
	    d__1 = et[i__ + j * 3 + 4];
	    et[i__ + j * 3 + 4] = max(d__1,-1.);
	    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &pneigh[4]
		    , &pnode[1], &a[i__ + j * 3 + 4], p, &ratio[1], nterms);
	}
    }

/*     minimizing the distance from the face to the node */

    for (ii = 1; ii <= 1000000; ++ii) {
	distmin = a[4];
	imin = 0;
	jmin = 0;
	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (a[i__ + j * 3 + 4] < distmin) {
		    distmin = a[i__ + j * 3 + 4];
		    imin = i__;
		    jmin = j;
		}
	    }
	}

/*       exit if minimum found */

	if (imin == 0 && jmin == 0) {
	    goto L200;
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		aold[i__ + j * 3 + 4] = a[i__ + j * 3 + 4];
		xiold[i__ + j * 3 + 4] = xi[i__ + j * 3 + 4];
		etold[i__ + j * 3 + 4] = et[i__ + j * 3 + 4];
	    }
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (i__ + imin >= -1 && i__ + imin <= 1 && j + jmin >= -1 && 
			j + jmin <= 1) {
		    a[i__ + j * 3 + 4] = aold[i__ + imin + (j + jmin) * 3 + 4]
			    ;
		    xi[i__ + j * 3 + 4] = xiold[i__ + imin + (j + jmin) * 3 + 
			    4];
		    et[i__ + j * 3 + 4] = etold[i__ + imin + (j + jmin) * 3 + 
			    4];
		} else {
		    xi[i__ + j * 3 + 4] += imin * d2;
		    et[i__ + j * 3 + 4] += jmin * d2;

/* Computing MIN */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = max(d__1,-1.);

		    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &
			    pneigh[4], &pnode[1], &a[i__ + j * 3 + 4], p, &
			    ratio[1], nterms);
/*              write(*,*) a(i,j) */
		}
	    }
	}
    }
L200:

/*     3rd run */
/*     initialisation */

    xiopt = xi[4];
    etopt = et[4];
    for (i__ = -1; i__ <= 1; ++i__) {
	for (j = -1; j <= 1; ++j) {
	    xi[i__ + j * 3 + 4] = xiopt + i__ * d3;
	    et[i__ + j * 3 + 4] = etopt + j * d3;
/* Computing MIN */
	    d__1 = xi[i__ + j * 3 + 4];
	    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
	    d__1 = xi[i__ + j * 3 + 4];
	    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
	    d__1 = et[i__ + j * 3 + 4];
	    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
	    d__1 = et[i__ + j * 3 + 4];
	    et[i__ + j * 3 + 4] = max(d__1,-1.);
	    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &pneigh[4]
		    , &pnode[1], &a[i__ + j * 3 + 4], p, &ratio[1], nterms);
	}
    }

/*     minimizing the distance from the face to the node */

    for (ii = 1; ii <= 1000000; ++ii) {
	distmin = a[4];
	imin = 0;
	jmin = 0;
	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (a[i__ + j * 3 + 4] < distmin) {
		    distmin = a[i__ + j * 3 + 4];
		    imin = i__;
		    jmin = j;
		}
	    }
	}

/*       exit if minimum found */

	if (imin == 0 && jmin == 0) {
	    goto L400;
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		aold[i__ + j * 3 + 4] = a[i__ + j * 3 + 4];
		xiold[i__ + j * 3 + 4] = xi[i__ + j * 3 + 4];
		etold[i__ + j * 3 + 4] = et[i__ + j * 3 + 4];
	    }
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (i__ + imin >= -1 && i__ + imin <= 1 && j + jmin >= -1 && 
			j + jmin <= 1) {
		    a[i__ + j * 3 + 4] = aold[i__ + imin + (j + jmin) * 3 + 4]
			    ;
		    xi[i__ + j * 3 + 4] = xiold[i__ + imin + (j + jmin) * 3 + 
			    4];
		    et[i__ + j * 3 + 4] = etold[i__ + imin + (j + jmin) * 3 + 
			    4];
		} else {
		    xi[i__ + j * 3 + 4] += imin * d3;
		    et[i__ + j * 3 + 4] += jmin * d3;

/* Computing MIN */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = max(d__1,-1.);

		    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &
			    pneigh[4], &pnode[1], &a[i__ + j * 3 + 4], p, &
			    ratio[1], nterms);
/*              write(*,*) a(i,j) */
		}
	    }
	}
    }
L400:

/*     4th run */
/*     initialisation */

    xiopt = xi[4];
    etopt = et[4];
    for (i__ = -1; i__ <= 1; ++i__) {
	for (j = -1; j <= 1; ++j) {
	    xi[i__ + j * 3 + 4] = xiopt + i__ * d4;
	    et[i__ + j * 3 + 4] = etopt + j * d4;
/* Computing MIN */
	    d__1 = xi[i__ + j * 3 + 4];
	    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
	    d__1 = xi[i__ + j * 3 + 4];
	    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
	    d__1 = et[i__ + j * 3 + 4];
	    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
	    d__1 = et[i__ + j * 3 + 4];
	    et[i__ + j * 3 + 4] = max(d__1,-1.);
	    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &pneigh[4]
		    , &pnode[1], &a[i__ + j * 3 + 4], p, &ratio[1], nterms);
	}
    }

/*     minimizing the distance from the face to the node */

    for (ii = 1; ii <= 1000000; ++ii) {
	distmin = a[4];
	imin = 0;
	jmin = 0;
	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (a[i__ + j * 3 + 4] < distmin) {
		    distmin = a[i__ + j * 3 + 4];
		    imin = i__;
		    jmin = j;
		}
	    }
	}

/*       exit if minimum found */

	if (imin == 0 && jmin == 0) {
	    goto L300;
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		aold[i__ + j * 3 + 4] = a[i__ + j * 3 + 4];
		xiold[i__ + j * 3 + 4] = xi[i__ + j * 3 + 4];
		etold[i__ + j * 3 + 4] = et[i__ + j * 3 + 4];
	    }
	}

	for (i__ = -1; i__ <= 1; ++i__) {
	    for (j = -1; j <= 1; ++j) {
		if (i__ + imin >= -1 && i__ + imin <= 1 && j + jmin >= -1 && 
			j + jmin <= 1) {
		    a[i__ + j * 3 + 4] = aold[i__ + imin + (j + jmin) * 3 + 4]
			    ;
		    xi[i__ + j * 3 + 4] = xiold[i__ + imin + (j + jmin) * 3 + 
			    4];
		    et[i__ + j * 3 + 4] = etold[i__ + imin + (j + jmin) * 3 + 
			    4];
		} else {
		    xi[i__ + j * 3 + 4] += imin * d4;
		    et[i__ + j * 3 + 4] += jmin * d4;

/* Computing MIN */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = xi[i__ + j * 3 + 4];
		    xi[i__ + j * 3 + 4] = max(d__1,-1.);
/* Computing MIN */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = min(d__1,1.);
/* Computing MAX */
		    d__1 = et[i__ + j * 3 + 4];
		    et[i__ + j * 3 + 4] = max(d__1,-1.);

		    distattach_new(&xi[i__ + j * 3 + 4], &et[i__ + j * 3 + 4], &
			    pneigh[4], &pnode[1], &a[i__ + j * 3 + 4], p, &
			    ratio[1], nterms);
/*              write(*,*) a(i,j) */
		}
	    }
	}
    }
L300:

    distattach_new(&xi[4], &et[4], &pneigh[4], &pnode[1], &a[4], p, &ratio[1], 
	    nterms);

    for (i__ = 1; i__ <= 3; ++i__) {
	pnode[i__] = p[i__ - 1];
    }

    *dist = a[4];

    elemcoords[0] = xi[4];
    elemcoords[1] = et[4];
    return 0;
} /* attach_ */




/* distattach.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/




/*     CalculiX - A 3-dimensional finite element program */
/*              Copyright (C) 1998-2005 Guido Dhondt */

/*     This program is free software; you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation(version 2); */


/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* Subroutine */ int distattach_new(xig, etg, pneigh, pnode, a, p, ratio, nterms)
double *xig, *etg, *pneigh, *pnode, *a, *p, *ratio;
int *nterms;
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;

    /* Builtin functions */
    int s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static int i__, j;
    static double dummy, et, xi;



/*     calculates the distance between the node with coordinates */
/*     in "pnode" and the node with local coordinates xig and etg */
/*     in a face described by "nterms" nodes with coordinates */
/*     in pneigh */




    /* Parameter adjustments */
    --ratio;
    --p;
    --pnode;
    pneigh -= 4;

    /* Function Body */
    if (*nterms == 3) {
	xi = (*xig + 1.) / 2.;
	et = (*etg + 1.) / 2.;
	if (xi + et > 1.) {
	    dummy = xi;
	    xi = 1. - et;
	    et = 1. - dummy;
	}
	ratio[1] = 1. - xi - et;
	ratio[2] = xi;
	ratio[3] = et;
    } else if (*nterms == 4) {
	xi = *xig;
	et = *etg;
	ratio[1] = (1. - xi) * (1. - et) / 4.;
	ratio[2] = (xi + 1.) * (1. - et) / 4.;
	ratio[3] = (xi + 1.) * (et + 1.) / 4.;
	ratio[4] = (1. - xi) * (et + 1.) / 4.;
    } else if (*nterms == 6) {
	xi = (*xig + 1.) / 2.;
	et = (*etg + 1.) / 2.;
	if (xi + et > 1.) {
	    dummy = xi;
	    xi = 1. - et;
	    et = 1. - dummy;
	}
	ratio[1] = (.5 - xi - et) * 2. * (1. - xi - et);
	ratio[2] = xi * (xi * 2. - 1.);
	ratio[3] = et * (et * 2. - 1.);
	ratio[4] = xi * 4. * (1. - xi - et);
	ratio[5] = xi * 4. * et;
	ratio[6] = et * 4. * (1. - xi - et);
    } else if (*nterms == 8) {
	xi = *xig;
	et = *etg;
	ratio[1] = (1. - xi) * (1. - et) * (-xi - et - 1.) / 4.;
	ratio[2] = (xi + 1.) * (1. - et) * (xi - et - 1.) / 4.;
	ratio[3] = (xi + 1.) * (et + 1.) * (xi + et - 1.) / 4.;
	ratio[4] = (1. - xi) * (et + 1.) * (-xi + et - 1.) / 4.;
	ratio[5] = (1. - xi * xi) * (1. - et) / 2.;
	ratio[6] = (xi + 1.) * (1. - et * et) / 2.;
	ratio[7] = (1. - xi * xi) * (et + 1.) / 2.;
	ratio[8] = (1. - xi) * (1. - et * et) / 2.;
    } else {
	printf ("ERROR in distattach\n");
    }

/*     calculating the position in the face */

    for (i__ = 1; i__ <= 3; ++i__) {
	p[i__] = 0.;
	i__1 = *nterms;
	for (j = 1; j <= i__1; ++j) {
	    p[i__] += ratio[j] * pneigh[i__ + j * 3];
	}
    }

/*     calculating the distance */

/* Computing 2nd power */
    d__1 = pnode[1] - p[1];
/* Computing 2nd power */
    d__2 = pnode[2] - p[2];
/* Computing 2nd power */
    d__3 = pnode[3] - p[3];
    *a = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;

    return 0;
} /* distattach_new */

