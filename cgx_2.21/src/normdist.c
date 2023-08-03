#include <extUtil.h>


double normdist( double *p0, double *p1, double *p2, double *p3 )
/******************************************************************/
/*   return distance of p0 to plane p1p2p3                        */
/******************************************************************/
{
  double v0[3], v1[3], v2[3], v3[3], en[3], fi;
  
  // normalized normal vector
  v_result( p1,p2, v1);
  v_result( p1,p3, v2);
  v_prod( v1, v2, v3 );
  v_norm( v3, en );

  // angle between p1p0 and en
  v_result( p1, p0, v0);
  fi=v_angle(v0, en );

  /*
  v_print(v0);
  v_print(v1);
  v_print(v2);
  v_print(v3);
  v_print(en);
  printf("fi:%f cos(fi):%f v0:%f\n",fi*180./PI, cos(fi), v_betrag(v0));
  */
  
  // normal distance of p0 based on trigonometry G=H*cos(PI-fi)
  return((v_betrag(v0)*cos(fi)));
}
