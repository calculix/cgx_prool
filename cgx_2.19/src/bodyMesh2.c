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
/* mesht Bodys indem es Scheibchenweise mit surfMesh arbeitet */
{
  int u,v,w,nn;
  double *xn, *yn, *zn;

  nn=*umax**vmax**wmax;
  if((xn=(double *)malloc(nn*sizeof(double)) )==NULL)
  { printf(" ERROR: malloc failure\n\n"); return(0); }
  if((yn=(double *)malloc(nn*sizeof(double)) )==NULL)
  { printf(" ERROR: malloc failure\n\n"); return(0); }
  if((zn=(double *)malloc(nn*sizeof(double)) )==NULL)
  { printf(" ERROR: malloc failure\n\n"); return(0); }

    for (w=1; w<*wmax-1; w++)
    {
      /*
	printf(" in:\n");
      */
      for (u=0; u<*umax; u++)
      {
        for (v=0; v<*vmax; v++)
        {
        xn[u* *vmax +v]=x[u* *vmax* *wmax +v* *wmax +w];
        yn[u* *vmax +v]=y[u* *vmax* *wmax +v* *wmax +w];
        zn[u* *vmax +v]=z[u* *vmax* *wmax +v* *wmax +w];
      /*
	printf("n[%d,%d,%d] %lf %lf %lf \n", u,v,w, xn[u* *vmax +v],yn[u* *vmax +v],zn[u* *vmax +v]);
      */
        }
      }
      surfMesh( vmax, umax, xn, yn, zn);
      /*
	printf(" out:\n");
       */
     for (u=1; u<*umax-1; u++)
      {
        for (v=1; v<*vmax-1; v++)
        {
        x[u* *vmax* *wmax +v* *wmax +w]=xn[u* *vmax +v];
        y[u* *vmax* *wmax +v* *wmax +w]=yn[u* *vmax +v];
        z[u* *vmax* *wmax +v* *wmax +w]=zn[u* *vmax +v];
      /*
	printf("n[%d,%d,%d] %lf %lf %lf \n", u,v,w, xn[u* *vmax +v],yn[u* *vmax +v],zn[u* *vmax +v]);
      */
        }
      }
    }
    free(xn);
    free(yn);
    free(zn);
  return(1);
}
