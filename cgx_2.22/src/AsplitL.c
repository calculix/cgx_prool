#include <extUtil.h>
#define TOL_ACOS  1e-6  // try and error

double AsplitL( double *b, double *eu, double *ev, double *eg, double *en )
/******************************************************************/
/*   Grade (eg) schneidet Ebene (eu,ev) return g                  */
/*   determinante dritter Ordnung                                 */
/*   b= Abstand zwischen den Aufpunkten der Linie und Ebene       */
/*   b=eu*u + ev*v + eg *g  (e: Einheitsvektoren )                */
/* return MAX_FLOAT if eg parallel to euXev                       */
/******************************************************************/
{
  double g, D, Dg, a, c;

  if(abs(v_sprod(eg,en))<TOL_ACOS) return(MAX_FLOAT);

  a = eu[0]*ev[1]*eg[2]+ ev[0]*eg[1]*eu[2]+ eg[0]*eu[1]*ev[2];
  c = eg[0]*ev[1]*eu[2]+ eu[0]*eg[1]*ev[2]+ ev[0]*eu[1]*eg[2];
  D = a-c;

  a = eu[0]*ev[1]* b[2]+ ev[0]* b[1]*eu[2]+  b[0]*eu[1]*ev[2];
  c =  b[0]*ev[1]*eu[2]+ eu[0]* b[1]*ev[2]+ ev[0]*eu[1]* b[2];
  Dg= a-c;
  g = Dg / D;
  return (g);
}


