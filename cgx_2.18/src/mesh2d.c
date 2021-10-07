/* 
   Static Encoding of Unconstrainted Delaunay Triangulation 
               Algorithm by Dr. B. Kaan Karamete 

  Mesh Generator:mesh2d.c
   Must be compiled like: cc mesh2d.c -lm -O -o mesh2d

  Input : mesh2d.dat
  Input Generator: model2d.c with modler2d.h
  Input Format: 
          alpha beta <- local and global adaptation parameters.
          NB         <- number of boundaries
          bindx[1].n1 bindx[1].n2 <- Each boundary start and end nodes
            ..        ..
          bindx[NB].n1 bindx[NB].n2
          NBP        <-  Total number of nodes
          pc[1].x[1]  pc[1].x[2]    pc[1].x[0] <- Each node's coords and
                                                  point spacing value. 
            ..          ..             .. 
            ..          ..             ..
          pc[NBP].x[1] pc[NBP].x[2]  pc[NBP].x[0]  

  Output: mesh2d.out

*/


#include <extUtil.h>


#define TEST 0
#define MAX_LOOPS 1e5
/*
   if maximum number of nodes is exceeded or you have more memory 
   you can increase both MAXPSIZE and/or MAXESIZE 
*/

/*
 replaced by malloc, wittig
#define MAXPSIZE 50001
#define MAXESIZE 100001
*/

#define SUFFERED 50
#define BOUNDLOOPS 100
#define eq       ==
#define or       ||
#define neq      !=
#define and      &&
#define true     1
#define false    0
#define sqr(x) (x)*(x)

/* Data Structure Definitions */

typedef struct point_list {double x[3];int bnd,next;int m;} point;

typedef struct triangle_list 
      {
       int n[4];
       int t[4],nt;
       point xc;
       double R;
      }
triangle;

typedef struct boundary_index {int n1,n2;} boundary;

typedef struct integer_point {long x[3];} ipoint;

typedef struct {
int NP,NE,NBP,NT,NSE,NewE,NE1,NB,MNP,NBND;
point    *pc;
triangle *tc;
int      suf[SUFFERED+1];
triangle ntc[3*SUFFERED+1];
boundary bindx[BOUNDLOOPS];
int near[BOUNDLOOPS],NNEARS;
double alpha, beta;
int nadapt; /*default number of adaptation cycles*/
int bloking;/*default node numbering style */
int smoothing; /*default flag for coordinate smoothing*/
point max,min,a,b;
int sufdel[SUFFERED];
int EXTRA; /* extra nodes to be put in case of missing edges*/ 
} GlobalVariables;

/* Function Prototyping */

void Load_Triangle(GlobalVariables *,int,int,int,int,int);
point Vector(GlobalVariables *,point,point);
double Cross(GlobalVariables *,point,point);
double Dot(GlobalVariables *,point,point);
int Node_In(GlobalVariables *,int,int);
short InCircleTest(GlobalVariables *,point,triangle);
triangle Compute_Circumcircle(GlobalVariables *,triangle);
void Initialize(GlobalVariables *);
void Suffered(GlobalVariables *,int,int);
void Construct_New_Triangles(GlobalVariables *,int);
int Newplaces(GlobalVariables *,int);
int Readfile(GlobalVariables *);
void Write2tecfile(GlobalVariables *);
int Engine(GlobalVariables *,int, int *ppivot);
void Construct_New_Neighbors1(GlobalVariables *);
void Construct_New_Neighbors2(GlobalVariables *);
int Find_Missing_Edges(GlobalVariables *);
void Set_Next_Fields(GlobalVariables *);
int Node_Renumber(GlobalVariables *);

int Insert_Nodes(GlobalVariables *,int *ppivot);
void Compute_Centroids(GlobalVariables *,int tno,double *xc,double *yc);
void Intro(GlobalVariables *,int c,char **v);
void Smooth(GlobalVariables *);

/* Global Variables */
//char *LastUpdate="  Last Update: 4 Dec 1997";
char *LastUpdate="  Last Update: 15. Feb 2019 Wi";

/* Implementation of Functions */


void Load_Triangle(GlobalVariables *gv, int id,int n1,int n2,int n3,int nt)
{ 
  gv->tc[id].n[1]=n1;
  gv->tc[id].n[2]=n2;
  gv->tc[id].n[3]=n3;
  gv->tc[id].nt=nt;
}
 


point Vector(GlobalVariables *gv,point p,point q)
{ int i;
  point pq;
  for(i=1;i<=2;++i) pq.x[i]=q.x[i]-p.x[i];
  return(pq);
}



ipoint Vectori(GlobalVariables *gv,point p,point q)
{ int i;
  ipoint pp,qq; 
  ipoint pq;
  int j;
  for(j=1;j<=2;++j)
  {
    pp.x[j]=(long)(gv->a.x[j]*p.x[j]+gv->b.x[j]);
    qq.x[j]=(long)(gv->a.x[j]*q.x[j]+gv->b.x[j]);
  } 
 
  for(i=1;i<=2;++i) pq.x[i]=qq.x[i]-pp.x[i];
  return(pq);
}



double Cross(GlobalVariables *gv,point p1,point p2)
{
  return((p1.x[1]*p2.x[2]-p1.x[2]*p2.x[1]));
} 



long Crossi(GlobalVariables *gv,ipoint p1,ipoint p2)
{ 
  return((p1.x[1]*p2.x[2]-p1.x[2]*p2.x[1]));
} 



double Dot(GlobalVariables *gv,point p1,point p2)
{
  return((p1.x[1]*p2.x[1]+p1.x[2]*p2.x[2]));
}



int Node_In(GlobalVariables *gv,int node,int prev)
{
  double epsnodein;
  /*  long epsnodein;*/
  int i,j,tt,k,search;
  int n1,n2,n3,n4;
  short notin,ok;
  
  int count=0;

#if TEST
  FILE *t;
#endif

  epsnodein=-1e-8;
  /* epsnodein=0.; */ 

  search=prev;
  do
  { 
    //if(count>MAX_LOOPS) return(-1); else count++;
    if(count>MAX_LOOPS)
    //if(node>62)
    {

#if TEST
  printf(" failed mesh written to mesh2d.out\n"); 
  t=fopen("mesh2d.out","w"); 
  fprintf(t,"*NODE\n");
  for(i=1;i<=gv->NP;++i) 
    fprintf(t,"%d, %lf, %lf, %lf\n", i, gv->pc[i].x[1],gv->pc[i].x[2],0.);
  fprintf(t,"*ELEMENT, TYPE=S3R, ELSET=Ernum \n");
  for(i=1;i<=gv->NE;++i)
    fprintf(t,"%d, %d, %d, %d\n",i,gv->tc[i].n[1], gv->tc[i].n[2],gv->tc[i].n[3]);
  fclose(t);
#endif

      return(-1);
    }
    else count++;
    notin=false;

    for(i=1;i<=3;++i)
    {
      if(Cross(gv,Vector(gv,gv->pc[gv->tc[search].n[i]],gv->pc[gv->tc[search].n[(i%3)+1]]),
           Vector(gv,gv->pc[gv->tc[search].n[i]],gv->pc[node])) < epsnodein)
      { notin=true;
        n1=gv->tc[search].n[i];n2=gv->tc[search].n[(i%3)+1];
        ok=false; 
        for(j=1;j<=gv->tc[search].nt;++j)
        { 
         tt=gv->tc[search].t[j];
         for(k=1;k<=3;++k)
         { 
           n3=gv->tc[tt].n[k]; n4=gv->tc[tt].n[(k%3)+1];          
           if(((n1 eq n3) and (n2 eq n4)) or  
              ((n1 eq n4) and (n2 eq n3)))    
           {
             search=tt;
             ok=true;
             break;
           }
         }
         if(ok) break; 
       } 
       if(ok) break;
      }
    }
  } while(notin); 
  return(search);
} 


  
short InCircleTest(GlobalVariables *gv,point node,triangle tri)
{ double epsradius=-1e-5;
  point point2R ;
  point2R=Vector(gv,node,tri.xc);
  if(Dot(gv,point2R,point2R)<(1.0+epsradius)*tri.R) return(true);
  else return(false);
}



triangle Compute_Circumcircle(GlobalVariables *gv,triangle tri)
{
  int i;
  point pq,pr,pxc,xc;
  double Area, f[4];
  triangle triupdate;

  pq=Vector(gv,gv->pc[tri.n[1]],gv->pc[tri.n[2]]);
  pr=Vector(gv,gv->pc[tri.n[1]],gv->pc[tri.n[3]]);
  Area=Cross(gv,pq,pr);
  for(i=1;i<=3;++i) f[i]=Dot(gv,gv->pc[tri.n[i]],gv->pc[tri.n[i]]);
  xc.x[1]=((f[2]-f[1])*pr.x[2]-(f[3]-f[1])*pq.x[2])/(2.0*Area);
  xc.x[2]=((f[3]-f[1])*pq.x[1]-(f[2]-f[1])*pr.x[1])/(2.0*Area);
  triupdate=tri;
  triupdate.xc=xc;
  pxc=Vector(gv,gv->pc[tri.n[1]],xc);
  triupdate.R=Dot(gv,pxc,pxc);
  return(triupdate);
} 



void Test_Cavity(GlobalVariables *gv,int node)
{
  int i,k,n1,n2,j,l;
  short olmadi;
  long Area;
  for(i=1;i<=gv->NSE;++i) gv->sufdel[i]=1;
  for(i=1;i<=gv->NSE;++i)
  for(k=1;k<=3;++k)
   {
     n1=gv->tc[gv->suf[i]].n[k];
     n2=gv->tc[gv->suf[i]].n[(k%3)+1];
     olmadi=false;
    for(j=1;j<=gv->NSE;++j)
    if(i neq j)
    {
      for(l=1;l<=3;++l)
      if(( (gv->tc[gv->suf[j]].n[l] eq n1) and
           (gv->tc[gv->suf[j]].n[(l%3)+1] eq n2)
          ) or
         ( (gv->tc[gv->suf[j]].n[l] eq n2) and
           (gv->tc[gv->suf[j]].n[(l%3)+1] eq n1) ))
       {olmadi=true;break;}
 
     }
      if(!olmadi)
      { gv->ntc[1].n[1]=n1;
        gv->ntc[1].n[2]=n2;
        gv->ntc[1].n[3]=node;
      
        Area=Crossi(gv,Vectori(gv, gv->pc[gv->ntc[1].n[1]],gv->pc[gv->ntc[1].n[2]]),
                    Vectori(gv,gv->pc[gv->ntc[1].n[1]],gv->pc[gv->ntc[1].n[3]]));
        if(Area<0) gv->sufdel[i]=-1;
      }
   }
}
   


void Suffered(GlobalVariables *gv,int node,int pivot)
{
  int i;
  short olmadi;
  int ind,j;
  /*  int update[SUFFERED],nnse;  */

  gv->NSE=1;
  gv->suf[gv->NSE]=pivot;
  ind=0;
  while(ind<gv->NSE and gv->NSE<SUFFERED)
  {
    ++ind;
    pivot=gv->suf[ind];
    for(i=1;i<=gv->tc[pivot].nt;++i)
    {
      olmadi=false; 
      for(j=1;j<=gv->NSE;++j)
      if(gv->tc[pivot].t[i] eq gv->suf[j]) {olmadi=true;break;}
      if(!olmadi)
        if(InCircleTest(gv,gv->pc[node],gv->tc[gv->tc[pivot].t[i]])) 
          if (gv->tc[gv->tc[pivot].t[i]].t[0] neq 1)
          {
            ++gv->NSE;
            gv->suf[gv->NSE]=gv->tc[pivot].t[i];
          }
    }
  }
  /*
  Test_Cavity(node);
  nnse=0;
  for(i=1;i<=gv->NSE;++i)
   if(gv->sufdel[i] neq -1) { ++nnse; update[nnse]=gv->suf[i];}
  gv->NSE=nnse;
  for(i=1;i<=gv->NSE;++i) gv->suf[i]=update[i];
  */
}

void Construct_New_Triangles(GlobalVariables *gv,int Newnode)
{
  int i,j,k,l;
  short olmadi;
  int n1,n2;

  for(i=1;i<=3*SUFFERED-1;++i) gv->ntc[i].nt=0;
  gv->NewE=0;

  for(i=1;i<=gv->NSE;++i)
    for(k=1;k<=3;++k)
    {
     n1=gv->tc[gv->suf[i]].n[k];
     n2=gv->tc[gv->suf[i]].n[(k%3)+1];
     olmadi=false;
    for(j=1;j<=gv->NSE;++j)
    if(i neq j)
    { 
      for(l=1;l<=3;++l)
      if(( (gv->tc[gv->suf[j]].n[l] eq n1) and
           (gv->tc[gv->suf[j]].n[(l%3)+1] eq n2)
          ) or
         ( (gv->tc[gv->suf[j]].n[l] eq n2) and
           (gv->tc[gv->suf[j]].n[(l%3)+1] eq n1) ))
       {olmadi=true;break;}
 
     }    
      if(!olmadi) { ++gv->NewE;
                    gv->ntc[gv->NewE].n[1]=n1; 
                    gv->ntc[gv->NewE].n[2]=n2; 
                    gv->ntc[gv->NewE].n[3]=Newnode; 
                    gv->ntc[gv->NewE].n[0]=gv->suf[i];
                    gv->ntc[gv->NewE]=Compute_Circumcircle(gv,gv->ntc[gv->NewE]);
                    
                    } 
   }
}  



int Newplaces(GlobalVariables *gv,int newel)
{
  if(newel<=gv->NSE) return(gv->suf[newel]);
  else
  {
    return(gv->NE+(newel-gv->NSE));
  }
}
/* construct neighborhood info for new triangles
  for new triangles 
  if any two nodes of the new triangle=the other new triangles
  or
  the neighborhood triangles of its suffered triangle
  stored at gv->ntc[i].n[0].
*/



void Construct_New_Neighbors1(GlobalVariables *gv)
{
  int i,j,k,l,n1,n2,kk;
  short oldu=0;
  int neybor=0;

  /* check for neighborhood triangles of its suffered triangle*/

  for(i=1;i<=gv->NewE;++i)
    for(k=1;k<=3;++k)
    {
      n1=gv->ntc[i].n[k];
      n2=gv->ntc[i].n[(k%3)+1];
      for (j=1;j<=gv->tc[gv->ntc[i].n[0]].nt;++j)
      { oldu=false;
        neybor=gv->tc[gv->ntc[i].n[0]].t[j];
        for(l=1;l<=3;++l)
        if( ((gv->tc[neybor].n[l] eq n1) and (gv->tc[neybor].n[(l%3)+1] eq n2))
            or
          ((gv->tc[neybor].n[l] eq n2) and (gv->tc[neybor].n[(l%3)+1] eq n1)))
        { oldu=true; break;}
        if(oldu) break;
      }
      if(oldu)
      {
        ++gv->ntc[i].nt;
        gv->ntc[i].t[gv->ntc[i].nt]=neybor;
        for(kk=1;kk<=gv->tc[neybor].nt;++kk)
        if(gv->tc[neybor].t[kk] eq gv->ntc[i].n[0]) break;
        gv->tc[neybor].t[kk]=Newplaces(gv,i);
      
      }   
    }         
}



void Construct_New_Neighbors2(GlobalVariables *gv)
{
  int i,j,k,l,n1,n2;
  short oldu;
  oldu=false;

  for(i=1;i<=gv->NewE;++i)
  for(k=1;k<=3;++k)
  {
    n1=gv->ntc[i].n[k];
    n2=gv->ntc[i].n[(k%3)+1];
    for(j=1;j<=gv->NewE;++j)
    if(i neq j)
    { oldu=false;
      for(l=1;l<=3;++l)
      if( ((gv->ntc[j].n[l] eq n1) and (gv->ntc[j].n[(l%3)+1] eq n2))
          or
          ((gv->ntc[j].n[l] eq n2) and (gv->ntc[j].n[(l%3)+1] eq n1)))
      { oldu=true; break;}
     if (oldu) break;
    }
    if(oldu)
    {
    ++gv->ntc[i].nt;
    gv->ntc[i].t[gv->ntc[i].nt]=Newplaces(gv,j);
    }
  }
  if(gv->NE+gv->NewE-gv->NSE>gv->NE)
  {
#if TEST
    printf("realloc gv->tc to %d elems\n", gv->NE+gv->NewE-gv->NSE+1);
#endif
    if( (gv->tc = (triangle *)realloc((triangle *)gv->tc, (gv->NE+gv->NewE-gv->NSE+1)*sizeof(triangle) )) == NULL )
      { printf(" ERROR: realloc failure in mesh2d() exiting\n\n");  exit(0); }
  }
  for(i=1;i<=gv->NewE;++i)
    gv->tc[Newplaces(gv,i)]=gv->ntc[i];
  gv->NE=gv->NE+gv->NewE-gv->NSE;
}



int Engine(GlobalVariables *gv,int Newnode, int *ppivot)
{
  if( Newnode==-1) { *ppivot=1;  return(0); }

  *ppivot=Node_In(gv,Newnode,*ppivot);
  if (*ppivot==-1) return(-1);
  Suffered(gv,Newnode,*ppivot);
  Construct_New_Triangles(gv,Newnode);
  Construct_New_Neighbors1(gv);
  Construct_New_Neighbors2(gv);
  return(1);
}



void Set_Next_Fields(GlobalVariables *gv)
{
  int i,j;
  for(i=1;i<=gv->NBP;++i) {gv->pc[i].next=0;gv->pc[i].bnd=0;}
  for(i=1;i<=gv->NB;++i)
  {
    for(j=gv->bindx[i].n1;j<=gv->bindx[i].n2;++j)
    {
      gv->pc[j].next=j+1;
      gv->pc[j].bnd=i;
    }
    gv->pc[gv->bindx[i].n2].next=gv->bindx[i].n1;
  }
  for(i=1;i<=gv->NBP;++i) gv->pc[i].m=gv->pc[i].next;
}



void Initialize(GlobalVariables *gv)
{
  int i;
  point boxmax,boxmin;
  double lr;


  /* create initial convex hull */
   /*
  for(i=1;i<=2;++i)
  {
   gv->pc[1].x[i]=gv->min.x[i]-fabs((gv->max.x[i]-gv->min.x[i])*0.2); 
   gv->pc[3].x[i]=gv->max.x[i]+fabs((gv->max.x[i]-gv->min.x[i])*0.2);
  }
   */

  /* length ratio between dir 1 and 2 (wittig, to gain better start-elements for stretched domains) */
  lr=(gv->max.x[1]-gv->min.x[1])/(gv->max.x[2]-gv->min.x[2]);
  i=1;
   gv->pc[1].x[i]=gv->min.x[i]-fabs((gv->max.x[i]-gv->min.x[i])*0.2); 
   gv->pc[3].x[i]=gv->max.x[i]+fabs((gv->max.x[i]-gv->min.x[i])*0.2);
  i=2;
   gv->pc[1].x[i]=gv->min.x[i]-fabs((gv->max.x[i]-gv->min.x[i])*0.2*lr); 
   gv->pc[3].x[i]=gv->max.x[i]+fabs((gv->max.x[i]-gv->min.x[i])*0.2*lr);

  gv->pc[2].x[1]=gv->pc[3].x[1];
  gv->pc[2].x[2]=gv->pc[1].x[2];
  gv->pc[4].x[1]=gv->pc[1].x[1];
  gv->pc[4].x[2]=gv->pc[3].x[2];
 
  Load_Triangle(gv,1,1,2,3,1);
  gv->tc[1].t[1]=2;
  Load_Triangle(gv,2,1,3,4,1);
  gv->tc[2].t[1]=1;
  gv->NP=4; 
  gv->tc[1]=Compute_Circumcircle(gv,gv->tc[1]);
  gv->tc[2]=Compute_Circumcircle(gv,gv->tc[2]);
  gv->NE=2;
  /*
  Points 1-4 and Triangles 1-2 are used for creating  the initial template 
  */
  for(i=1;i<=2;++i)
  { 
   boxmin.x[i]=1.0;
   boxmax.x[i]=10000.0;
   gv->a.x[i]=(boxmax.x[i]-boxmin.x[i])/(gv->pc[3].x[i]-gv->pc[1].x[i]);
   gv->b.x[i]=boxmax.x[i]-gv->a.x[i]*gv->pc[3].x[i];
  }
#if TEST
  puts("**********************************************************");
  printf("\nGenerator in progress starting with %d nodes and\n",gv->NBP-4);
  printf(" with the parameters n1=%d n2=%d n3=%d\n\n",
                   gv->nadapt,gv->bloking,gv->smoothing);
  puts("**********************************************************");
#endif
}



void Remove_Unwanted_Triangles(GlobalVariables *gv)
{
  int i,j;
  int pivot;

#if TEST
  puts("Unwanted triangles are being swept off...");
#endif
  for(i=1;i<=gv->NE;++i)
    for(j=1;j<=3;++j)
      if(gv->tc[i].n[j] <= 4) {gv->tc[i].t[0]=1;++gv->NE1;break;}

  for(i=1;i<=gv->NE;++i)
  if (
     gv->pc[gv->tc[i].n[2]].bnd eq gv->pc[gv->tc[i].n[1]].bnd and
     gv->pc[gv->tc[i].n[3]].bnd eq gv->pc[gv->tc[i].n[1]].bnd and 
     gv->pc[gv->tc[i].n[1]].bnd neq 0)
  { 
    j=1; 
    pivot=gv->tc[i].n[1];
    do
    { 
      if(gv->pc[pivot].next eq gv->tc[i].n[(j%3)+1]) break;
      else if(gv->pc[pivot].next eq gv->tc[i].n[(((j%3)+1)%3)+1])
                        {gv->tc[i].t[0]=1;++gv->NE1;break;}
      else if(j<3) pivot=gv->tc[i].n[++j];
      else pivot=gv->pc[pivot].next;
    }while(true);
  }
} 



int Find_Missing_Edges(GlobalVariables *gv)
{
  int i,j,k;
  // changed to realloc Wittig
  int *missing=NULL;
  gv->MNP=0;

  for(i=1;i<=gv->NE;++i)
    for(j=1;j<=3;++j)
      for(k=1;k<=3;++k)
        if(gv->pc[gv->tc[i].n[j]].m eq gv->tc[i].n[k])
          gv->pc[gv->tc[i].n[j]].m=0;
 
  for(i=1;i<=gv->NP;++i)
  {
    if(gv->pc[i].m neq 0)
    {
      if( (missing = (int *)realloc((int *)missing, (gv->MNP+2)*sizeof(int) )) == NULL )
      { printf(" ERROR: realloc failure in mesh2d()\n\n");      return(0); }
      missing[++gv->MNP]=i;
    } 
  }

  for(i=1;i<=gv->MNP;++i)
  {
    ++gv->NBP;
    ++gv->NP;
    if( (gv->pc = (point *)realloc((point *)gv->pc, (gv->NP+1)*sizeof(point) )) == NULL )
    { printf(" ERROR: realloc failure in mesh2d()\n\n");      return(0); }
    gv->pc[gv->NP].next=gv->pc[missing[i]].next;
    gv->pc[gv->NP].m=gv->pc[missing[i]].next;
    gv->pc[gv->NP].bnd=gv->pc[missing[i]].bnd;
    gv->pc[gv->NP].x[1]=(gv->pc[missing[i]].x[1]+gv->pc[gv->pc[missing[i]].next].x[1])/2.0;
    gv->pc[gv->NP].x[2]=(gv->pc[missing[i]].x[2]+gv->pc[gv->pc[missing[i]].next].x[2])/2.0;
    gv->pc[gv->NP].x[0]=gv->pc[missing[i]].x[0];
#if TEST
    printf("edge %lf %lf\n",gv->pc[gv->NP].x[1],gv->pc[gv->NP].x[2]);
#endif
    if((gv->pc[gv->NP].x[1]==gv->pc[gv->NP-1].x[1])&&(gv->pc[gv->NP].x[2]==gv->pc[gv->NP-1].x[2]))
    {
      /* found the same point again, break */
      gv->NP--;
      gv->NBP--;
      gv->MNP--;
    }
    gv->pc[missing[i]].next=gv->NP;
    gv->pc[missing[i]].m=gv->NP;
    gv->pc[missing[i]].bnd=gv->pc[missing[i]].bnd;
    gv->NBND=gv->bindx[gv->NB].n2+1;
  }
#if TEST
  printf("Detecting and Curing %d Missing Edges...\n",gv->MNP);
#endif
  free(missing);
  gv->EXTRA=gv->EXTRA+gv->MNP;
  return(gv->MNP);
}



void Near_Nodes(GlobalVariables *gv,int tno,int *nn)
{ int i,k,l;
  int pivot[5];
  pivot[1]=tno;
  for(i=1;i<=gv->tc[tno].nt;++i) pivot[i+1]=gv->tc[tno].t[i]; 
  k=0;
  for(l=1;l<=gv->tc[tno].nt+1;++l)
    for(i=1;i<=gv->tc[pivot[l]].nt;++i)
      if(gv->tc[gv->tc[pivot[l]].t[i]].n[0] neq 0)
      {++k;gv->near[k]=gv->tc[gv->tc[pivot[l]].t[i]].n[0];}
  *nn=k;
} 



int Insert_Nodes(GlobalVariables *gv,int *ppivot)
{
  int  i,j;
  point p[5];
  double dm[4],sj;
  short  reject;
  /* short EXCEEDED=false; */
  int nold=0;
  int cnt=0;

  while(cnt<gv->nadapt and (gv->NP-nold)>(0.1*nold))
  { ++cnt;
    nold=gv->NP;
    for(i=1;i<=gv->NE;++i) gv->tc[i].n[0]=0;
    for(i=1;i<=gv->NE;++i)
    if(gv->tc[i].t[0] eq 0)
    {
      reject=0;
      for (j=1;j<=3;++j) p[j]=gv->pc[gv->tc[i].n[j]];
  p[4].x[0]=(p[1].x[0]+p[2].x[0]+p[3].x[0])/3.0;
  p[4].x[1]=(p[1].x[1]+p[2].x[1]+p[3].x[1])/3.0;
  p[4].x[2]=(p[1].x[2]+p[2].x[2]+p[3].x[2])/3.0;
  for (j=1;j<=3;++j)
  dm[j]=sqrt((p[j].x[1]-p[4].x[1])*(p[j].x[1]-p[4].x[1])
           +(p[j].x[2]-p[4].x[2])*(p[j].x[2]-p[4].x[2]));
  for (j=1;j<=3;++j)  if (dm[j]<=(gv->alpha*p[4].x[0])) { reject=1; break;}
  if(!reject)
  {
    Near_Nodes(gv,i,&gv->NNEARS);
    for(j=1;j<=gv->NNEARS;++j)
    {
      sj=sqrt(sqr(gv->pc[gv->near[j]].x[1]-p[4].x[1])+
        sqr(gv->pc[gv->near[j]].x[2]-p[4].x[2]));
      if (sj<=(gv->beta*p[4].x[0])) { reject=1;break;}
    }  
   }
   if(!reject)
   { 
    ++gv->NP;
    if( (gv->pc = (point *)realloc((point *)gv->pc, (gv->NP+1)*sizeof(point) )) == NULL )
    { printf(" ERROR: realloc failure in mesh2d()\n\n");      return(0); }
    gv->pc[gv->NP]=p[4];
    gv->pc[gv->NP].next=0;
    gv->pc[gv->NP].bnd=0;
    gv->tc[i].n[0]=gv->NP;
    // due to malloc deactivated by wittig
    //if(gv->NP eq MAXPSIZE-1) { EXCEEDED=true;break;}
      }
    }
#if TEST
    printf("%d Adaptation: %d nodes",cnt,gv->NP-4);
#endif
    for(i=nold+1;i<=gv->NP;++i) if(Engine(gv,i,ppivot)<0) return(-1);
#if TEST
    printf(" %d triangles\n",gv->NE-gv->NE1);
#endif
    // due to malloc deactivated by wittig
    //if(EXCEEDED) { printf("MAXPSIZE in mesh2d.c exceeded! Stopping... Increase MAXPSIZE if possible\n");break;}
  }
  return(1);
}



void Compute_Centroids(GlobalVariables *gv,int tno,double *xc,double *yc)
{
  int i;
  double sumx,sumy;
  sumx=0.0;sumy=0.0;
  for(i=1;i<=3;++i)
  {
    sumx=sumx+gv->pc[gv->tc[tno].n[i]].x[1];
    sumy=sumy+gv->pc[gv->tc[tno].n[i]].x[2];
  }
  *xc=sumx/3.0;
  *yc=sumy/3.0;
}


 
int Node_Renumber(GlobalVariables *gv)
{
  // replaced by wittig
  //int tri[MAXESIZE],node[MAXPSIZE],nnode[MAXPSIZE],ntri[MAXESIZE];
  int *tri,*node,*nnode,*ntri;
  int i,j,pivot=0,NNE,NNP,cpivot=0,trino,cnode=0,ind,ok;
  double xc,yc,xsc,ysc,minimum,diff;
  int k,m,z,maxz,band[4],maxband,neks=0;
  int piv=0,sakla,stack[500],nstack,maxinter;
  point orta;
#if TEST
  FILE *t;
  puts("Node and Triangle Renumbering...");      
#endif
  NNP=0;
  NNE=1;
  if( (tri = (int *)malloc((gv->NE+1)*sizeof(int) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");
    return(0); }
  if( (ntri = (int *)malloc((gv->NE+2)*sizeof(int) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");
    return(0); }
  if( (node = (int *)malloc((gv->NP+1)*sizeof(int) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");
    return(0); }
  if( (nnode = (int *)malloc((gv->NP+2)*sizeof(int) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");
    return(0); }
  
  for(i=1;i<=gv->NE;++i) tri[i]=0;
  for(i=1;i<=gv->NP;++i) node[i]=0;
  if(gv->bloking eq 3){   
    orta.x[1]=gv->min.x[1]+(gv->max.x[1]-gv->min.x[1])/2.0;
     orta.x[2]=gv->min.x[2]+(gv->max.x[2]-gv->min.x[2])/2.0;}
  
   minimum=sqr(gv->max.x[1]);  
   for(i=1;i<=gv->NE;++i)
     if(gv->tc[i].t[0] neq 1)
     {
        Compute_Centroids(gv,i,&xc,&yc);
         switch (gv->bloking)
          { case 1: diff=sqr(xc-gv->max.x[1]);
                    if(diff<minimum) {minimum=diff;pivot=i;}
                    break;
            case 2: diff=sqr(yc-gv->max.x[2]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break;
            case 3:
                  diff=sqr(xc-orta.x[1])+sqr(yc-orta.x[2]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break; 
            default:
                  diff=sqr(xc-orta.x[1])+sqr(yc-orta.x[2]);
                  if(diff<minimum) {minimum=diff;pivot=i;}
                  break; 
        }
     }  
  tri[pivot]=NNE;
  for(i=1;i<=3;++i) {++NNP;node[gv->tc[pivot].n[i]]=NNP;}
  Compute_Centroids(gv,pivot,&xc,&yc);
  xsc=xc;ysc=yc;
  nstack=0;
  stack[++nstack]=pivot;
  maxinter=0;
  maxband=0;
  do
  { ind=0;
    if(nstack>maxinter) maxinter=nstack;
    for(j=1;j<=nstack;++j)
    { pivot=stack[j];
      sakla=ind;
    for(i=1;i<=gv->tc[pivot].nt;++i)
    { trino=gv->tc[pivot].t[i];
      if((tri[trino] eq 0) and (gv->tc[trino].t[0] neq 1))
      {
        Compute_Centroids(gv,trino,&xc,&yc);
        ++ind;
        switch (gv->bloking)
        { case 1: diff=sqr(xc-xsc);break;
          case 2: diff=sqr(yc-ysc);break;
          case 3: diff=sqr(xc-xsc)+sqr(yc-ysc);break;
          case 4: z=0;
                  for(k=1;k<=3;++k)
                  for(m=1;m<=3;++m)
                    if(gv->tc[pivot].n[k]==gv->tc[trino].n[m])
                      band[++z]=gv->tc[pivot].n[k];
                  for(m=1;m<=3;++m)
                  { ok=0;
                   for(k=1;k<=z;++k)
                   if(band[k] eq gv->tc[trino].n[m]) {ok=1; break;}
                    if(!ok) {neks=gv->tc[trino].n[m];break;}
                   }
               if(node[neks] eq 0) neks=NNP+1; else neks=node[neks];
                maxz=0;
                for(k=1;k<=z;++k) 
                  if(abs(node[band[k]]-neks)>maxz)
                      maxz=abs(node[band[k]]-neks);
                 diff=(sqr(xc-xsc)+sqr(yc-ysc))*(maxz);
                 break;
        default: diff=sqr(xc-xsc)+sqr(yc-ysc);break;
      }
      if(ind eq 1) { minimum=diff; cpivot=trino;piv=pivot;}
      else if (minimum>diff) {minimum=diff; cpivot=trino;piv=pivot;}
    }
  }
  if(sakla eq ind) stack[j]=0;
  }   
  for(i=1;i<=nstack;++i)
    if(stack[i] eq 0)
    { 
      for(j=i+1;j<=nstack;++j)
       stack[j-1]=stack[j];
      --nstack;
     }
    
  if(ind>0)
  {
    stack[++nstack]=cpivot;
    Compute_Centroids(gv,cpivot,&xc,&yc);
    xsc=(xsc*NNE+xc)/(NNE+1);
    ysc=(ysc*NNE+yc)/(NNE+1);
    ++NNE;
    tri[cpivot]=NNE;
    for(i=1;i<=3;++i)
    {  
      ok=0;
      for(j=1;j<=3;++j)
      if(gv->tc[cpivot].n[i] eq gv->tc[piv].n[j])
      { ok=1;break;}
     if(!ok)
     {
      cnode=gv->tc[cpivot].n[i];
      break;
     }
    }
    if(node[cnode] eq 0) {++NNP;node[cnode]=NNP;}

     for(k=1;k<=3;++k) 
      if(abs(node[gv->tc[cpivot].n[k]]-node[gv->tc[cpivot].n[(k%3)+1]])>maxband)
        maxband=abs(node[gv->tc[cpivot].n[k]]-node[gv->tc[cpivot].n[(k%3)+1]]);
   }    
  else
  break;
  }while(NNE<(gv->NE-gv->NE1));
#if TEST
  printf("Maximum interface elements= %d\n",maxinter);
  printf("Maximum band= %d\n",maxband);
#endif
  if(gv->NP-4 eq NNP)
  {
      ; 
#if TEST
   printf("ok.\n"); 
#endif
  }
  else {printf("ERROR: Renumbering could not be implemented properly\n"); 
         return(0);}
  
  for(i=1;i<=gv->NP;++i) 
    nnode[node[i]]=i;
  for(i=1;i<=gv->NE;++i)
    ntri[tri[i]]=i;
  
#if TEST
  t=fopen("mesh2d.out","w"); 
  //fprintf(t,"VARIABLES=x y dpi bnd\n");
  //fprintf(t,"ZONE N=%d,E=%d,F=FEPOINT,ET=TRIANGLE\n",NNP,NNE);
  fprintf(t,"*NODE\n");
  for(i=1;i<=NNP;++i) 
    //fprintf(t,"%lf %lf %lf %d\n", gv->pc[nnode[i]].x[1],gv->pc[nnode[i]].x[2],gv->pc[nnode[i]].x[0],gv->pc[nnode[i]].bnd);
    fprintf(t,"%d, %lf, %lf, %lf\n", nnode[i], gv->pc[nnode[i]].x[1],gv->pc[nnode[i]].x[2],0.);
  
  fprintf(t,"*ELEMENT, TYPE=S3R, ELSET=Ernum \n");
  for(i=1;i<=NNE;++i)
    //fprintf(t,"%d %d %d\n",node[gv->tc[ntri[i]].n[1]], node[gv->tc[ntri[i]].n[2]],node[gv->tc[ntri[i]].n[3]]);
    fprintf(t,"%d, %d, %d, %d\n",i,node[gv->tc[ntri[i]].n[1]], node[gv->tc[ntri[i]].n[2]],node[gv->tc[ntri[i]].n[3]]);
  fclose(t);
  printf("Output file mesh2d.out is generated!\n");
  
  t=fopen("mesh2d.dx","w");
  fprintf(t,"object 1 class array type double rank 1 shape 3 items %d data follows\n",NNP-1);
  for(i=0;i<=NNP-1;++i) 
    fprintf(t,"%lf %lf %lf\n", gv->pc[nnode[i+1]].x[1],gv->pc[nnode[i+1]].x[2],gv->pc[nnode[i+1]].x[0]);
 
  fprintf(t,"object 2 class array type int rank 1 shape 3 items %d data follows\n",NNE-1);
  for(i=0;i<=NNE-1;++i)
    fprintf(t,"%d %d %d\n",node[gv->tc[ntri[i+1]].n[1]], node[gv->tc[ntri[i+1]].n[2]],node[gv->tc[ntri[i+1]].n[3]]);
  fprintf(t,"attribute \"element type\" string \"triangles\"\n");
  fprintf(t,"attribute \"ref\" string \"positions\"\n");
  fprintf(t,"object \"mesh2d\" class field\n");
  fprintf(t,"component \"positions\" value 1\n");
  fprintf(t,"component \"connections\" value 2\n");
  fprintf(t,"end\n");
  fclose(t);
  printf("Output file mesh2d.dx is generated!\n");
 
  
  t=fopen("mesh2d.vtk","w");
  fprintf(t,"# vtk DataFile Version 1.0\n");
  fprintf(t,"Mesh2d Output\n");
  fprintf(t,"ASCII\n\n");
  fprintf(t,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(t,"POINTS %d double\n",NNP);
  for(i=1;i<=NNP;++i) 
    fprintf(t,"%lf %lf 0.0\n", gv->pc[nnode[i]].x[1],gv->pc[nnode[i]].x[2]);
  fprintf(t,"\n"); 
  fprintf(t,"CELLS %d %d\n",NNE,4*NNE);
  for(i=1;i<=NNE;++i)
    fprintf(t,"3 %d %d %d\n",node[gv->tc[ntri[i]].n[1]]-1, node[gv->tc[ntri[i]].n[2]]-1,node[gv->tc[ntri[i]].n[3]]-1);
  fprintf(t,"\n");
  fprintf(t,"CELL_TYPES %d\n",NNE);
  for(i=1;i<=NNE;++i)
    fprintf(t,"5\n\n");
  fclose(t);
  printf("Output file mesh2d.vtk is generated!\n");
  
  t=fopen("ptc.dat","w");
  fprintf(t,"%d\n",NNE);
  for(i=1;i<=NNE;++i)
    { fprintf(t,"%d ",i);
       for(j=1;j<=gv->tc[ntri[i]].nt;++j)
       { trino=gv->tc[ntri[i]].t[j];
         if(gv->tc[trino].t[0] neq 1)
          fprintf(t,"%d ",tri[trino]);
         else fprintf(t,"-1 ");
        }
      fprintf(t,"\n");
    }
  fclose(t);
  printf("Triangles info file ptc.dat is generated!\n");
#endif
  free(tri);
  free(node);
  free(nnode);
  free(ntri);
  return(1);
}



void Intro(GlobalVariables *gv,int c,char *v[])
{
  puts("\nA Delaunay 2D Constrained Unstructured Mesh");
  puts(" Generator developed by B. Kaan Karamete");
  puts(LastUpdate);
  puts(" Usage: mesh2d <n1> <n2> <n3>");
  puts("   n1: number of adaptation cycles<1..10>");
  puts("   n2: centroidal node renumbering style<1..3>");
  puts("   n3: Flag for coordinate Smoothing<0..1>");
  puts(" defaults: n1=8; n2=3; n3=1");    
  if(c>1) gv->nadapt=atoi(v[1]);  if(gv->nadapt<0) gv->nadapt=8;
  if(c>2) gv->bloking=atoi(v[2]); if(!(gv->bloking>=0 and gv->bloking<5)) gv->bloking=3;
  if(c>3) gv->smoothing=atoi(v[3]); if(gv->smoothing neq 0 and gv->smoothing neq 1) 
  gv->smoothing=1;
}



void Smooth(GlobalVariables *gv)
{
  // replaced by wittig
  //int ptc[MAXPSIZE][20];
  int **ptc;
  int i,j,tno,k;
  double xc,yc,xsc,ysc;

#if TEST
  printf("Smoothing three times... ");
#endif
  //if(gv->NP>=MAXPSIZE) { printf(" FATAL ERROR in Smooth gv->NP:%d > MAXPSIZE:%d\n", gv->NP,MAXPSIZE); exit(0); }
  //for(i=1;i<=gv->NP;++i)  ptc[i][0]=0;

  if( (ptc = (int **)malloc((gv->NP+1)*sizeof(int *) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");
    return; }
  for(i=0; i<=gv->NP; i++)
    if( (ptc[i] = (int *)calloc((20),sizeof(int) )) == NULL )
    { printf(" ERROR: realloc failure in mesh2d()\n\n");
      return; }

  for(i=1;i<=gv->NE;++i)
    if(gv->tc[i].t[0] neq 1)
    {
      for(j=1;j<=3;++j)
      {
        // corrected by wittig
        //if(gv->tc[i].n[j]>gv->bindx[gv->NB].n2+4+gv->EXTRA)
        if(gv->tc[i].n[j]>gv->bindx[gv->NB].n2+gv->EXTRA)
        {
          ++ptc[gv->tc[i].n[j]][0];
          tno=ptc[gv->tc[i].n[j]][0];
          if( (ptc[gv->tc[i].n[j]] = (int *)realloc((int *)ptc[gv->tc[i].n[j]], (tno+20)*sizeof(int) )) == NULL )
          { printf(" ERROR: realloc failure in mesh2d()\n\n"); return; }
          ptc[gv->tc[i].n[j]][tno]=i;
        }
      }
    }
  //else printf("gv->tc:%d skipped\n", i);

  for(k=1;k<=3;++k)
    for(i=1;i<=gv->NP;++i)
        // corrected by wittig
      //if(i>gv->bindx[gv->NB].n2+4+gv->EXTRA)
      if(i>gv->bindx[gv->NB].n2+gv->EXTRA)
      {
        xsc=0.0;
        ysc=0.0;
        for(j=1;j<=ptc[i][0];++j)
        {
          Compute_Centroids(gv,ptc[i][j],&xc,&yc);
          xsc=xsc+xc;
          ysc=ysc+yc;
        }
        if(ptc[i][0]>0)
        {
          gv->pc[i].x[1]=xsc/ptc[i][0];
          gv->pc[i].x[2]=ysc/ptc[i][0];
        }
        //else printf("gv->pc:%d skipped, ptc[i][0]:%d\n", i, ptc[i][0]);

      }
  //else printf("gv->pc:%d skipped, gv->bindx[gv->NB].n2:%d gv->EXTRA:%d\n", i, gv->bindx[gv->NB].n2, gv->EXTRA);
  for(i=1;i<=gv->NP;++i)  free(ptc[i]); free(ptc);

}




void Write2tecfile(GlobalVariables *gv)
{
  FILE *t,*f;
  int i,j,trino;

  t=fopen("mesh2d.out","w");
  //fprintf(t,"VARIABLES=x y dpi bnd\n");
  //fprintf(t,"ZONE N=%d,E=%d,F=FEPOINT,ET=TRIANGLE\n",gv->NP,gv->NE-gv->NE1);
  fprintf(t,"*NODE\n");
  for(i=1;i<=gv->NP;++i) 
    //fprintf(t,"%lf %lf %lf %d\n",  gv->pc[i].x[1],gv->pc[i].x[2],gv->pc[i].x[0],gv->pc[i].bnd);
    fprintf(t,"%d, %lf, %lf, %lf\n", i, gv->pc[i].x[1],gv->pc[i].x[2],0.);
  fprintf(t,"*ELEMENT, TYPE=S3R, ELSET=Etec \n");
  for(i=1;i<=gv->NE;++i)
  if(gv->tc[i].t[0] neq 1)
    // fprintf(t,"%d %d %d\n",gv->tc[i].n[1],gv->tc[i].n[2],gv->tc[i].n[3]);
    fprintf(t,"%d, %d, %d, %d\n",i,gv->tc[i].n[1],gv->tc[i].n[2],gv->tc[i].n[3]);
  fclose(t);
  printf("Output file mesh2d.out is generated!\n");

  f=fopen("ptc.dat","w");
  fprintf(f,"%d\n",gv->NE);
  for(i=1;i<=gv->NE;++i)
  { fprintf(f,"%d ",i);
     for(j=1;j<=gv->tc[i].nt;++j)
     { trino=gv->tc[i].t[j];
       if(gv->tc[trino].t[0] neq 1)
        fprintf(f,"%d ",trino);
       else fprintf(f,"-1 ");
      }
    fprintf(f,"\n");
  }
  fclose(f);
  printf("Triangles info file ptc.dat is generated!\n");
}  



int Readfile(GlobalVariables *gv)
{ 
  FILE *in;
  int i,j;
  
  in=fopen("mesh2d.dat","r");
  if(in eq NULL) return(0);
    
  fscanf(in,"%lf %lf\n",&gv->alpha,&gv->beta);
  fscanf(in,"%d\n",&gv->NB);
  for(i=1;i<=gv->NB;++i)
  {
   fscanf(in,"%d %d\n",&gv->bindx[i].n1,&gv->bindx[i].n2);
   gv->bindx[i].n1=gv->bindx[i].n1+4;
   gv->bindx[i].n2=gv->bindx[i].n2+4;
  }
  fscanf(in,"%d\n",&gv->NT); /* gv->NP olacak */
  gv->NBP=gv->NT+4;
  fscanf(in,"%lf %lf %lf\n",&gv->pc[5].x[1],&gv->pc[5].x[2],&gv->pc[5].x[0]);
   gv->min=gv->pc[5];
   gv->max=gv->pc[5];
  for(i=6;i<=gv->NBP;++i)
  {
   fscanf(in,"%lf %lf %lf\n",&gv->pc[i].x[1],&gv->pc[i].x[2],&gv->pc[i].x[0]);
    for(j=1;j<=2;++j)
    {
     if(gv->pc[i].x[j]<gv->min.x[j]) gv->min.x[j]=gv->pc[i].x[j]; 
     if(gv->pc[i].x[j]>gv->max.x[j]) gv->max.x[j]=gv->pc[i].x[j]; 
    }
   } 
  fclose(in);
  gv->NBND=gv->bindx[gv->NB].n2;
  Set_Next_Fields(gv);
  return(1);
}




/*
IN:
_nt: no of uv-points
_nb: no of trimming loops
npc[]: no of uv-points in each loop
pnt_u: u (of all curves)
pnt_v: v

OUT:
pnt_u: u (of all nodes)
pnt_v: v
apre: no of nodes and elements
epre: elements (tr3)
return: number of tri3 or 0 if in error, -1 if the distance between points is 0.
 */



int mesh2d(int *_nt, int _nb, int *npc, double **_pnt_u, double **_pnt_v, int **_pnt_flag, int **_tri3, double _alpha, double _beta, int _nadapt)
{
#if TEST
  FILE *handle;
#endif
  int i,j,n,nexact;
  double *pnt_u, *pnt_v;
  int *pnt_flag;
  int *tri3;
  double *dl=NULL, p0[3], p1[3], p0p1[3];
  int sum_sp, sum_ep;

  double scale_umax;
  double scale_umin;
  double scale_vmax;
  double scale_vmin;
  double scale_u;
  double scale_v;
  double scale_w;

  GlobalVariables *gv;
  int pivot=1;

  /* initialize all global variables */
  if( (gv = (GlobalVariables *)calloc( (1),sizeof(GlobalVariables) )) == NULL )
  { printf(" ERROR: malloc failure in mesh2d()\n\n");      return(0); }
  gv->NP=gv->NE=gv->NBP=gv->NT=gv->NSE=gv->NewE=gv->NE1=gv->NB=gv->MNP=gv->NBND=gv->EXTRA=gv->NNEARS=0;
  if( (gv->pc = (point *)calloc( (6),sizeof(point) )) == NULL )
  { printf(" ERROR: malloc failure in mesh2d()\n\n");      return(0); }
  if( (gv->tc = (triangle *)calloc( (6),sizeof(triangle) )) == NULL )
  { printf(" ERROR: malloc failure in mesh2d()\n\n");      return(0); }
  Engine(gv,-1,&pivot);

  gv->alpha=_alpha;
  gv->beta=_beta;
  gv->nadapt=_nadapt;    /* number of adaptation cycles*/
  gv->bloking=0;   /* node numbering style */
  gv->smoothing=1; /* flag for coordinate smoothing*/
#if TEST
  printf("gv->alpha:%lf gv->beta:%lf gv->nadapt:%d\n", gv->alpha,gv->beta,gv->nadapt);
#endif
  pnt_u=*_pnt_u;
  pnt_v=*_pnt_v;
  pnt_flag=*_pnt_flag;
  tri3=*_tri3;

  /* store in global mesh2d variables */
  gv->NT=*_nt;
  gv->NB=_nb;
  if(gv->NB>=BOUNDLOOPS) { printf(" ERROR: Increase BOUNDLOOPS in mesh2d to more than:%d\n",gv->NB); return(0); }

  /* scale coords to higher values, this code can not deal with small values */
  scale_umax=-MAX_FLOAT;
  scale_umin=MAX_FLOAT;
  scale_vmax=-MAX_FLOAT;
  scale_vmin=MAX_FLOAT;

  for(i=0;i<gv->NT;i++)
  { 
    if(pnt_u[i]>scale_umax) scale_umax=pnt_u[i];
    if(pnt_u[i]<scale_umin) scale_umin=pnt_u[i];
    if(pnt_v[i]>scale_vmax) scale_vmax=pnt_u[i];
    if(pnt_v[i]<scale_vmin) scale_vmin=pnt_u[i];
  }
  scale_u=(scale_umax+scale_umin)/2.;
  scale_v=(scale_vmax+scale_vmin)/2.;
  scale_umax=scale_umax-scale_u;
  scale_vmax=scale_vmax-scale_v;
  scale_umin=scale_umin-scale_u;
  scale_vmin=scale_vmin-scale_v;
  if (scale_umax < (-scale_umin)) scale_umax=(-scale_umin);
  if (scale_vmax < (-scale_vmin)) scale_vmax=(-scale_vmin);
  scale_w=scale_umax;
  if (scale_w < scale_vmax){ scale_w=scale_vmax;}

  scale_w/=0.4; /* nochmal scaliert */
  if (scale_w<=0.) scale_w=1.;

  for(i=0;i<gv->NT;i++)
  { 
    pnt_u[i]=(pnt_u[i]-scale_u)/scale_w;
    pnt_v[i]=(pnt_v[i]-scale_v)/scale_w;
#if TEST
    printf("uv %e %e\n", pnt_u[i],pnt_v[i]);
#endif
  }

  /* startno of the curve */
  sum_sp=1+4;
  gv->bindx[1].n1=5;
  for(i=2;i<=gv->NB;++i)
  {
    sum_sp+=npc[i-2];
    gv->bindx[i].n1=sum_sp;
  }

  /* endno of the curves */
  sum_ep=4;
  for(i=1;i<=gv->NB;++i)
  {
    sum_ep+=npc[i-1];
    gv->bindx[i].n2=sum_ep;
  }

  gv->NBP=gv->NT+4;

  /* calculate the spacing between the uv-points */
  n=0; p0[2]=p1[2]=0.;
  /* alloc dl */
  if( (dl = (double *)malloc((sum_ep)*sizeof(double) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");
    return(0); }
  for(j=1; j<=gv->NB; j++)
  {
    for(i=gv->bindx[j].n1-5;i<gv->bindx[j].n2-5; i++)
    {
      if(i<0) break;
      p0[0]=pnt_u[i];
      p0[1]=pnt_v[i];
      p1[0]=pnt_u[i+1];
      p1[1]=pnt_v[i+1];
      v_result(p0, p1, p0p1);
      dl[i]=v_betrag(p0p1);
      if(dl[i]==0.)
      {
        printf(" ERROR: dist between nodes in uv-coords is %f for %d 1st=%d last:%d\n",dl[i], i, gv->bindx[j].n1-5, gv->bindx[j].n2-5);
        printf(" p1:%e %e\n", p0[0],p0[1]);
        printf(" p2:%e %e\n", p1[0],p1[1]);
        printf(" p1p2:%e %e\n", p0p1[0],p0p1[1]);
        free(dl); return(-1);
      }
    }
    if(i<1) { printf(" ERROR: Problems with dist between nodes in uv-coords, i:%d\n",i); free(dl); return(-1); }
    else
    {
      dl[i]=dl[i-1];
      if(dl[i]==0.) { printf(" ERROR: dist between nodes in uv-coords is %f for %d\n",dl[i],i); free(dl); return(-1); }
    }
  }

  /* store the points in the mesh2d structure */
  gv->pc[5].x[1]=pnt_u[0];
  gv->pc[5].x[2]=pnt_v[0];
  gv->pc[5].x[0]=dl[0];
  gv->min=gv->pc[5];
  gv->max=gv->pc[5];
  if( (gv->pc = (point *)realloc((point *)gv->pc, (gv->NBP+1)*sizeof(point) )) == NULL )
  { printf(" ERROR: realloc failure in mesh2d()\n\n");      return(0); }
  for(i=6;i<=gv->NBP;++i)
  {
    gv->pc[i].x[1]=pnt_u[i-5];
    gv->pc[i].x[2]=pnt_v[i-5];
    gv->pc[i].x[0]=dl[i-5];
    for(j=1;j<=2;++j)
    {
      if(gv->pc[i].x[j]<gv->min.x[j]) gv->min.x[j]=gv->pc[i].x[j]; 
      if(gv->pc[i].x[j]>gv->max.x[j]) gv->max.x[j]=gv->pc[i].x[j]; 
    }
  } 
  free(dl);

#if TEST
  handle = fopen ("mesh2d1.dat", "w");
  fprintf(handle, "curves:%d points:%d\n", gv->NB, gv->NT);
  fprintf(handle, "min1:%lf max1:%lf\n",gv->min.x[1], gv->max.x[1]);
  fprintf(handle, "min2:%lf max2:%lf\n",gv->min.x[1], gv->max.x[1]);
  for(i=1;i<=gv->NB;++i) fprintf(handle, "%d %d\n",gv->bindx[i].n1,gv->bindx[i].n2);
  for(i=5;i<=gv->NBP;++i) fprintf(handle, "%lf %lf %lf\n", gv->pc[i].x[1],gv->pc[i].x[2],gv->pc[i].x[0]);
  fclose(handle);
#endif

  gv->NBND=gv->bindx[gv->NB].n2;
  Set_Next_Fields(gv);

  /* copy from original main() routine */
  Initialize(gv);

  do { ++gv->NP; if(Engine(gv,gv->NP,&pivot)<0) { printf(" Engine1 failed\n"); return(0);}  } while(gv->NP<gv->NBP);

  j=gv->NBP;  
  while((gv->MNP=Find_Missing_Edges(gv)) neq 0)
  {
    if((j-gv->MNP)<=0) { printf(" Meshing failed, please repeat with changed line divisions\n"); return(0); }
    j=gv->MNP;
    for(i=gv->NP-gv->MNP+1;i<=gv->NP;++i) if(Engine(gv,i,&pivot)<0)
    { printf(" Engine2 failed\n"); return(0); }
  }

  for(i=1;i<=gv->NE;++i) gv->tc[i].t[0]=0;
  gv->NE1=0;
  Remove_Unwanted_Triangles(gv);
  if(Insert_Nodes(gv,&pivot)<0) { printf(" Insert_Nodes failed\n"); return(0); }
  if(gv->smoothing) Smooth(gv);
#if TEST
  if(gv->bloking!=0) 
  {
    if(!Node_Renumber(gv)) Write2tecfile(gv);
  } else Write2tecfile(gv);
#endif  

  /* store the results in the interface variables */
  if( (pnt_u = (double *)realloc((double *)pnt_u, (gv->NP+1)*sizeof(double) )) == NULL )
  { printf(" ERROR: realloc failure11 in mesh2d()\n\n");
    return(0); }
  if( (pnt_v = (double *)realloc((double *)pnt_v, (gv->NP+1)*sizeof(double) )) == NULL )
  { printf(" ERROR: realloc failure11 in mesh2d()\n\n");
    return(0); }
  if( (pnt_flag = (int *)realloc((int *)pnt_flag, (gv->NP+1)*sizeof(int) )) == NULL )
  { printf(" ERROR: realloc failure11 in mesh2d()\n\n");
    return(0); }
  if( (tri3 = (int *)realloc((int *)tri3, ((gv->NE+1)*3)*sizeof(int) )) == NULL )
  { printf(" ERROR: realloc failure11 in mesh2d()\n\n");
    return(0); }
  for(i=0;i<4;i++)
  {
    pnt_u[i]=0.;
    pnt_v[i]=0.;
    pnt_flag[i]=0;
  }
  for(i=5;i<=gv->NP;++i)
  { 
    pnt_u[i-1]=gv->pc[i].x[1];
    pnt_v[i-1]=gv->pc[i].x[2];
    pnt_flag[i-1]=gv->pc[i].bnd;
    /* gv->pc[i].x[0],gv->pc[i].bnd */
  }

  /* back to original scale */
  for(i=4;i<gv->NP;++i)
  { 
    pnt_u[i]=(pnt_u[i]*scale_w)+scale_u;
    pnt_v[i]=(pnt_v[i]*scale_w)+scale_v;
  }

  n=0;
  for(i=1;i<=gv->NE;++i)
  { 
    if(gv->tc[i].t[0] neq 1)
    {
      tri3[n++]=gv->tc[i].n[1];
      tri3[n++]=gv->tc[i].n[2];
      tri3[n++]=gv->tc[i].n[3];
    }
  }
  *_nt=gv->NP;
  *_pnt_u=pnt_u;
  *_pnt_v=pnt_v;
  *_tri3=tri3;
  *_pnt_flag=pnt_flag;
  free(gv->pc);
  free(gv->tc);

  /* check the result */
  /* printf("With %d nodes %d triangles are generated!\n",gv->NP-4,gv->NE-gv->NE1); */
  nexact=2*((gv->NP-4)-(gv->NBND-4))+(gv->NBND-4)-2+2*(gv->NB-1); 
  if(nexact neq (gv->NE-gv->NE1)) 
  { 
    errMsg("ERROR: Something wrong in mesh2d()! triangles are generated:%d but should be=%d\n",gv->NE-gv->NE1, nexact);
    return(0);
  }
  return(n/3);
}
