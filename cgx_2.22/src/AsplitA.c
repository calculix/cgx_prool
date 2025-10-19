#include <extUtil.h>


/******************************************************************/
/*   area (a) cut area (b) return -1 if failed                    */
/*  in:                                                           */
/*  pa1-3 points define area a                                    */
/*  pb1-3 points define area b                                    */
/*  out:                                                          */
/*  ps1-2 points derived from vectors betw. split points of area b*/
/******************************************************************/
int AsplitA( double *pa1, double *pa2, double *pa3, double *pb1, double *pb2, double *pb3, double *ps1, double *ps2)
{
  int i,n=0, mini;
  double g,l,gl,ps[3][3], pa[3], pb[3], eu[3], ev[3], eg[3], en[3], pba[3], pbps[3];
  double psps[3][3], minval;
  static double tol1=-0.000001;
  static double tol2=1.000001;


  /* determine the vectors (eu,ev) defining area a */
  v_result( pa1, pa2, pa );
  v_norm( pa, eu );
  v_result( pa1, pa3, pa );
  v_norm( pa, ev );

  /* determine split points ps1, ps2 on the vectors between points defining area b */

  /* split line pb12 */
  v_result( pb1, pa1, pba );
  v_result( pb1, pb2, pb );
  l=v_norm( pb, eg );
  v_prod(eu,ev,en);
  g = AsplitL( pba, eu, ev, eg, en );
  gl=g/l;
  if((gl>=tol1)&&(gl<=tol2))
  {
    /* area a splits line pb12 */
    v_scal( &g, eg, pbps );
    v_add( pb1, pbps, ps[n] );
    n++;
  }

  /* split line pb23 */
  v_result( pb2, pa1, pba );
  v_result( pb2, pb3, pb );
  l=v_norm( pb, eg );
  v_prod(eu,ev,en);
  g = AsplitL( pba, eu, ev, eg, en );
  gl=g/l;
  if((gl>=tol1)&&(gl<=tol2))
  {
    /* area a splits line pb23 */
    v_scal( &g, eg, pbps );
    v_add( pb2, pbps, ps[n] );
    n++;
  }

  /* split line pb31 */
  v_result( pb3, pa1, pba );
  v_result( pb3, pb1, pb );
  l=v_norm( pb, eg );
  v_prod(eu,ev,en);
  g = AsplitL( pba, eu, ev, eg, en );
  gl=g/l;
  if((gl>=tol1)&&(gl<=tol2))
  {
    /* area a splits line pb31 */
    v_scal( &g, eg, pbps );
    v_add( pb3, pbps, ps[n] );
    n++;
  }

  if(n<2) return(-1); /* area a does not intersect area b */
  else     for (i=0; i<3; i++) ps1[i]=ps[0][i];

  if(n==2)
  {
    for (i=0; i<3; i++)
    {
      ps2[i]=ps[1][i];
    }
    return(1);
  }
  else
  {
    //for(i=0; i<3; i++) printf("AsplitA,n=3 ps %f %f %f\n",ps[i][0],ps[i][1],ps[i][2]);
    /* two points have to be very close and have to be averaged */
    v_result( &ps[0][0], &ps[1][0], &psps[0][0] );
    v_result( &ps[1][0], &ps[2][0], &psps[1][0] );
    v_result( &ps[2][0], &ps[0][0], &psps[2][0] );
    minval=MAX_FLOAT;
    for(i=0; i<3; i++)
    {
      ev[i]=v_betrag( &psps[i][0] );
      if(ev[i]<minval) { minval=ev[i]; mini=i; }
    }
    if(mini==0) { v_add(&ps[0][0], &ps[1][0],pa); mini=2; }
    else if(mini==1) { v_add(&ps[1][0], &ps[2][0],pa); mini=0; }
    else { v_add(&ps[2][0], &ps[0][0],pa); mini=1; }
    for (i=0; i<3; i++) ps1[i]=ps[mini][i];
    l=0.5;
    v_scal(&l,pa,ps2);
    //printf(" ps1 %f %f %f\n",ps1[0],ps1[1],ps1[2]);
    //printf(" ps2av %f %f %f\n",ps2[0],ps2[1],ps2[2]);
    return(2);
  }
  return(-1);
}


