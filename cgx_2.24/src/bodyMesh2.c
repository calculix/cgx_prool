/* --------------------------------------------------------------------  */
/*                          CALCULIX                                     */
/*                   - GRAPHICAL INTERFACE -                             */
/*                                                                       */
/*     A 3-dimensional pre- and post-processor for finite elements       */
/*              Copyright (C) 1996 Klaus Wittig                          */
/*                                                                       */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; version 2 of           */
/*     the License.                                                      */
/*                                                                       */
/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */
/*                                                                       */
/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */
/* --------------------------------------------------------------------  */


#include <cgx.h>


int bodyMesh2(int *wmax, int *vmax, int *umax, double *x, double *y, double * z)
/* mesht Bodys indem es Scheibchenweise mit surfMesh arbeitet (von body surf1 nach surf2 (in u dir.)) */
{
  int u,v,w,nn;
  double *xn, *yn, *zn;

  //printf("bodyMesh2\n");

  nn=*umax**vmax**wmax;
  if((xn=(double *)malloc(nn*sizeof(double)) )==NULL)
  { printf(" ERROR: malloc failure\n\n"); return(0); }
  if((yn=(double *)malloc(nn*sizeof(double)) )==NULL)
  { printf(" ERROR: malloc failure\n\n"); return(0); }
  if((zn=(double *)malloc(nn*sizeof(double)) )==NULL)
  { printf(" ERROR: malloc failure\n\n"); return(0); }
  
  for (u=1; u<*umax-1; u++)
  {
    /*
    printf(" in:\n");
    */
    for (v=0; v<*vmax; v++)
    {
      for (w=0; w<*wmax; w++)
      {
        xn[v* *wmax +w]=x[u* *vmax* *wmax +v* *wmax +w];
        yn[v* *wmax +w]=y[u* *vmax* *wmax +v* *wmax +w];
        zn[v* *wmax +w]=z[u* *vmax* *wmax +v* *wmax +w];
	//printf("pnt r%d:%d:%d %lf %lf %lf \n", u,v,w, xn[v* *wmax +w],yn[v* *wmax +w],zn[v* *wmax +w]);
      }
    }
    //printf("call surfMesh\n");
    surfMesh( wmax, vmax, xn, yn, zn);
    /*
    printf(" out:\n");
    */
    for (v=1; v<*vmax-1; v++)
    {
      for (w=1; w<*wmax-1; w++)
      {
        x[u* *vmax* *wmax +v* *wmax +w]=xn[v* *wmax +w];
        y[u* *vmax* *wmax +v* *wmax +w]=yn[v* *wmax +w];
        z[u* *vmax* *wmax +v* *wmax +w]=zn[v* *wmax +w];
	//printf("pnt s%d:%d:%d %lf %lf %lf \n", u,v,w, xn[v* *wmax +w],yn[v* *wmax +w],zn[v* *wmax +w]);
      }
    }
    // test
    /*
    printf("seto w%d\n",w);
    for (v=0; v<*vmax; v++)
    {
      for (w=0; w<*wmax; w++)
      {
        printf("pnt b%d:%d:%d %lf %lf %lf \n", u,v,w, xn[v* *wmax +w],yn[v* *wmax +w],zn[v* *wmax +w]);
      }
    }
    printf("setc w%d\n",w);
    */
  }
  free(xn);
  free(yn);
  free(zn);
  return(1);
}
