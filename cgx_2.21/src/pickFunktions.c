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


/*
TO DO:
- replace line_ with line_i in qsplitLine() and delete line_() in defineEntity and rename line_i to line_ in all functions.
- extend hitUndo()

- not finalized:
      qmshOperator=0;
*/


#define     TEST            0    /* debugging */

#include <cgx.h>

#define DPICK_BUFFER 10000
#define MAX_BUF      100
#define N1NM_NODES   62          // max number of nodes related to a second node in pre_cut


extern int   width_ini, height_ini;            /* Grafik-Fensterbreite/hoehe */
extern int   width_menu, height_menu;
extern int   w0, w1, activWindow;
extern int   width_w0, height_w0;
extern int   width_w1, height_w1;
extern int   MouseMode;                                   /* status maustasten */
extern double dtx, dty, drx, dry, drz, ds;                 /* Verschiebungen */
extern GLdouble R[4][4];                                   /* Rotationsmatrix */
extern char  surfFlag;                /* zeichne nur Oberflaechenelemente (1), sonst (0)*/
extern char  modelEdgeFlag;           /* zeichne mit Modell-Ecken (1), sonst (0)*/
extern char  drawMode;                /* protokoliert drawFunktion (drawLoad,Light,Animate,set)*/
extern char  frameFlag;               /* mit (1) oder ohne Rahmen um das Grafikfenster */
extern char  pickFlag;
extern char  graphFlag;
extern char  sequenceFlag;                  /* 1: play a sequence of LC */
extern char  vectorFlag;
extern double dx ,dy;                  /* Mauskoordinaten */
extern double centerPnt[3];            /* Rotationszentrum */
extern int   centerNode;            /* Nr of center Node, 0:no centernode */
extern double     gtol;                                    /* geometry tolerance for merging */
extern int     ddiv;
extern int   animList;
extern int   lcase_animList;
extern double   aspectRatio_w1;         /* width_w1/height_w1 */

extern Scale     scale[1];
extern Summen    anz[1];
extern Nodes     *node;
extern Elements  *e_enqire;

extern Datasets *lcase;
extern Faces     *face;
extern Texts     *ntext;

extern Alias     *alias;
extern Sets      *set;                
extern Shapes    *shape;
extern Psets     *pset;                
extern Points    *point;
extern Lines     *line;
extern Lcmb      *lcmb;
extern Gsur      *surf;
extern Gbod      *body;
extern Nurbl     *nurbl;
extern Nurbs     *nurbs;
extern SumGeo    anzGeo[1];
extern SumAsci   sumAsci[1];
extern double     *colNr;

extern char  datin[MAX_LINE_LENGTH];                          /* Input-data-file */

/* global variables for picking */
extern char pickfunc[MAX_LINE_LENGTH];  /* pickfunc either "qenq" "qadd" "qrem" ..  */
char hidefunc[MAX_LINE_LENGTH];         /* hidefunc stores the original commando if a qfunction was started inside another qfunction  */
extern char mode[2];                    /* pickmode */
char hidemode[2];                       /* hidemode stores the original pickmode */
extern GLdouble dx_cur, dy_cur;         /* pick-cursor Area */
extern double backgrndcol_rgb[4];
extern char    buffer[MAX_LINE_LENGTH];        /* common string-buffer */
extern char  printFlag;                     /* printf 1:on 0:off */
extern int  cur_entity;                                       /* aktive entity (component) */
extern char  v_dim;                         /* 1: scalar plot, 2: a 2D vector plot, 3: a 3D vectorplot, 4: a 3D vectorplot with signed vals */
extern int   entity_v[6];                                         /* components of a vector-entity */
extern int  cur_lc;
extern char  addDispFlag;                    /* 0: original node-coordinates, 1: node-coordinates+displacements */
extern char  illumResultFlag;                /* results with illumination */

extern double v[4];                                        /* drehkorrekturen fuer centerPkt */
extern GLdouble dR[4][4];                                  /* dR= R-Rmem fuer center */
extern GLdouble Rmem[4][4];
extern double dtx, dty, dtz;

int anz_threads=1;
int qmshOperator=0;                     /* operator for qdiv in qmsh 0|'*'|'/' */
int qbiaFlag=0;                         /* for qmsh, 1:switch from qlin to qbia */
int qdisFlag=0;                         /* for qdis */
int intersectFlag=0;                    /* for qint, qfill */
double qaddTol=-1;                      /* for qadd, angle tolerance for set completition of adjacent faces */
double filletRadius=-1;                 /* for qfill */
int    shp_pindx=0;                     /* for qshp */


static GLint hits;                             /* number of picked items, must be global because of glutKeyboardFunc(defineDiv) */
GLuint *selectBuf=NULL;                 /* buffer which holds the picked items, must be global because of glutKeyboardFunc(defineDiv) */
extern SpecialSet specialset[1];
extern int       setall;                /* setNr of the default set "all" */

int     pick_zmin;                      /* kleinster z-wert der gepickten items */
int     setNrbuf;                       /* buffer for pick() */
int     setNrdiv=-1;                    /* subset of buffer for qmsh(), stores only directly selected lines for defineDiv() */
int     qmshbuf;                         /* additional entity buffer for qmsh in pick() */
int    *pickdata;                       /* Zwischenspeicher fuer pickdaten */
double    pickbuf;                      /* kurzzeitbuffer */
int     pick_buffer;                    /* groesse von pickdata */
int     qnorCounter=0;                  /* counts number of selected points for qnor() */
int     qcutCounter=0;                  /* counts number of selected nodes for qcut() */
int     qaliCounter=0;                  /* counts number of selected nodes for qali() */
int     xbuf, ybuf;                     /* mauskoordinaten */
int     *selem=NULL;                    /* buffer for qflp, points to surfs */
int     entitybuf=-1;                   /* saves the index of the last created entity for undo */
char    keybuf=0;                       /* saves the type of the last selection */
int     set_highl=-1;                   /* set which stores actual entity */
int     setNr=-1;                       /* selected set */
int     pntNr=-1;                       /* selected point */
int     lineNr=-1;                      /* selected line */
int     surfNr=-1;                      /* selected surf, surface as a target for the point-projection */
int     nurbsNr=-1;                     /* selected Nurbs, used for all touched surfs (s,g) */
int     shapeNr=-1;                     /* selected Shape, used for all touched surfs (s,g) */
int     bodyNr=-1;                      /* selected body */

/* values for movePoint() */
int movedp=-1, pntnr=-1, pnt_undo=-1, moveFlag=0;

#define MAX_LINES 100
int qspl_i, qspl_indx[MAX_LINES]; /* buffer for already splitted lines (qspl) */
GLuint GLubuf[2];

/* qcut, cut */
#define CUTS 20          // max expected cuts per element
typedef struct {
  int anz_n;             // nr of splitpnt's
  double npre[CUTS][8];  // stores splitpnt coordinates (0-2), position (3)-> saved in qcut_nod[].val, element nodes at the splitted edge (4,5),
                         // nodes used for dataset interpolation (6,7, later saved in qcut_nod[].n1 and ~.n2)
                         // where qcut_nod[] is used to interpolate the datasets to the new node .nnew, later saved in qcut_nod[].nr.
  int ipre[CUTS];        // stores the index of the sorted npre (i<anz_n)
  int nnew[CUTS];        // stores the final node at the cut
  int nref[CUTS];        // stores the referenced base element node if the new node is 'close' to it
  int skip;              // will be set to 1 if a second cutParam exists which references at least the same nodes (to prevent double elements in one spot)
  double cgx,cgy,cgz;    // center of cutParam (cg)
} CutParam;
extern Qcut_nodes *qcut_nod;

typedef struct {
  int *n2;
  int *ns;
  int nn;
} Checknode;
Checknode *checknode=NULL;

// reference the midside nodes, used by pre_cut
typedef struct {
  int sum, *n2, *nm;
}N1nm;
N1nm *n1nm=NULL;

/* the copied node-sets which have to be filled with values from new loaded Datasets */
extern CopiedNodeSets copiedNodeSets[1];

/* shape buffer */
Shapes shapebuf;

/* threading */
int            glob_cutElement=0;
sem_t   sem_cute;
extern sem_t   sem_n;


int hitAction( GLuint *name, char *type, int x, int y )
{
  /* name[0]== entity type of picked name   */
  /* name[1]== picked name(Nr)  */
  /* type[0] == requested entity type */
  /* type[1] == often gkey in pick() */

  char buf[MAX_LINE_LENGTH], printbuf[MAX_LINE_LENGTH];
  int i,j,k,p,l;
  int bias_fbd, index;
  static int p_qnor[3], cpFlag=0, shFlag=0;
  static double v_qdis[4][3], v_val;
  double v0[3], v1[3], da,dax,day,daz, r1,r2,rm,lcir, dist;
  int matchFlag, lin[2]={0,0};

  /* highlight the entity if its not qdel or qrem ,qtxt */
  if((compare(pickfunc, "qdel", 4) != 4)&&(compare(pickfunc, "qrem", 4) != 4)&&(compare(pickfunc, "qtxt", 4) != 4))
  {
    buf[0]=name[0];
    buf[1]='\0';
    seta(set_highl, buf, name[1] );
  }

  //printf("type %c %c\n", type[0],type[1]);

  entitybuf=name[1];
  keybuf=name[0];

  if (name[0] == 'e')
  {
    printf ("%d t:%d n: ", name[1], e_enqire[name[1]].type);
    if(e_enqire[name[1]].type==1) k=8;
    else if(e_enqire[name[1]].type==2) k=6;
    else if(e_enqire[name[1]].type==3) k=4;
    else if(e_enqire[name[1]].type==4) k=20;
    else if(e_enqire[name[1]].type==5) k=15;
    else if(e_enqire[name[1]].type==6) k=10;
    else if(e_enqire[name[1]].type==7) k=3;
    else if(e_enqire[name[1]].type==8) k=6;
    else if(e_enqire[name[1]].type==9) k=4;
    else if(e_enqire[name[1]].type==10) k=8;
    else if(e_enqire[name[1]].type==11) k=2;
    else if(e_enqire[name[1]].type==12) k=3;
    else k=0;
    for (j=0; j<k; j++) printf("%d ",e_enqire[name[1]].nod[j]);
    printf("\n"); 
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].elem,set[i].anz_e,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qmsh", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qflp", 4) == 4) flip( type, name[1]);
  }
  else if (name[0] == 'f')
  {
    printf ("%d e:%d s:%d n= ", name[1], face[name[1]].elem_nr, face[name[1]].nr+1 );
      if(face[name[1]].type==7) k=3;
      else if(face[name[1]].type==8) k=6;
      else if(face[name[1]].type==9) k=4;
      else if(face[name[1]].type==10) k=8;
      else if(face[name[1]].type==11) k=2;
      else if(face[name[1]].type==12) k=3;
      else k=0;
      for (j=0; j<k; j++) printf("%d ",face[name[1]].nod[j]);
      printf("\n"); 
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].face,set[i].anz_f,name[1]) >-1))
          {
            sprintf(&printbuf[strlen(printbuf)],"%s(%d)",set[i].name,set[i].index);

            /* is a contact surface referenced ? */
            if(set[i].anz_se)
            {
              sprintf(&printbuf[strlen(printbuf)],"->[ ");
              for(j=0; j<set[i].anz_se; j++) 
                sprintf(&printbuf[strlen(printbuf)],"%s(%d) ",set[set[i].set[j]].name,set[set[i].set[j]].index);
              sprintf(&printbuf[strlen(printbuf)],"]");
            }
            sprintf(&printbuf[strlen(printbuf)],",");
          }
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qdis", 4) == 4)
    {
      delSet(specialset->impc);
      setNrbuf=pre_seta( specialset->impc, "i", 0);
      if(setNrbuf<0) { errMsg(" ERROR: could not create set for qdis\n"); return(0); } 
      seta( setNrbuf, type, name[1]);
      
      /* complete the set by all connected faces which do not violate the tolerance */
      qaddTol=QADDTOL;
      completeFacesByTolerance(set_highl, setNrbuf, qaddTol);
    }
  }
  else if (name[0] == 'h')
  {
    printf ("%s p1=%s p2=%s p3=%s\n", shape[name[1]].name, point[shape[name[1]].p[0]].name, point[shape[name[1]].p[1]].name, point[shape[name[1]].p[2]].name );
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].shp,set[i].anz_sh,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qsur", 4) == 4)
    {
      /* if a a default surf is known apply the shape */
      if(surfNr>-1)
      {
        surf[surfNr].sh=name[1];
        shape_refSurf(surf[surfNr].sh, surfNr);

	/* if it was prepared for illumination, change the interior def of the selected surf */
        if(surf[surfNr].pgn!=NULL)
	{
          free(surf[surfNr].pgn); surf[surfNr].pgn=NULL; surf[surfNr].npgn=0;

          /* create the interior */
          repSurf(surfNr,1);
	}
        printf(" Interiour of surf:%s changed to shape:%s\n", surf[surfNr].name, shape[name[1]].name);
      }
      else
      {
	printf(" Active Shape:%s\n", shape[name[1]].name);
        shapeNr=name[1];
      }
    }
    else if (compare(pickfunc, "qshp", 4) == 4)
    {
	printf(" Active Shape:%s\n", shape[name[1]].name);
        shapeNr=name[1];
        nurbsNr=-1;
    }
  }
  else if (name[0] == 't')
  {
    if(type[1]=='d') ntext[name[1]].node_nr=0;
    if(type[1]=='p') moveText(name[1], x, y);
    if(type[1]=='b') { ntext[name[1]].tx= -1; }
    if(type[1]=='n') { ntext[name[1]].nFlag=!ntext[name[1]].nFlag; }
    if(type[1]=='v') { ntext[name[1]].vFlag=!ntext[name[1]].vFlag; }
    if(type[1]=='t') { ntext[name[1]].tFlag=!ntext[name[1]].tFlag; }
    if(type[1]=='s') { ntext[name[1]].pFlag=!ntext[name[1]].pFlag; }
    if(type[1]=='f')
    {
      ntext[name[1]].fFlag++;
      if(ntext[name[1]].fFlag>2) ntext[name[1]].fFlag=0;
    }
  }
  else if ((name[0] == 'n')&& (node[name[1]].pflag!=1))
  { 
    if(!anz->l) printf ("%d xyz= %lf %lf %lf\n", name[1],
      (node[name[1]].nx* scale->w+scale->x),
      (node[name[1]].ny* scale->w+scale->y),
      (node[name[1]].nz* scale->w+scale->z) );
    else
    {
      if(sequenceFlag) rm=lcase[lcase_animList].dat[animList][name[1]];
      else rm=lcase[cur_lc].dat[cur_entity][name[1]];

      printf ("%d v= %e ", name[1], rm);
      /* search linked values (iexist=2) */
      if(sequenceFlag)
      {
        for(i=0; i<lcase[lcase_animList].ncomps; i++)
          if((i!=animList)&&(lcase[lcase_animList].iexist[i]==2)) printf ("%s:%5.4e ",lcase[lcase_animList].compName[i],lcase[lcase_animList].dat[i][name[1]]);
      }
      else
      {
        for(i=0; i<lcase[cur_lc].ncomps; i++)
          if((i!=cur_entity)&&(lcase[cur_lc].iexist[i]==2)) printf ("%s:%5.4e ",lcase[cur_lc].compName[i],lcase[cur_lc].dat[i][name[1]]);
      }

      printf (" xyz= %lf %lf %lf ",
      (node[name[1]].nx* scale->w+scale->x),
      (node[name[1]].ny* scale->w+scale->y),
      (node[name[1]].nz* scale->w+scale->z) );
    }

    v0[0]=0.;
    v0[1]=node[name[1]].ny* scale->w+scale->y;
    v0[2]=node[name[1]].nz* scale->w+scale->z;
    v1[0]=0.; v1[1]=1.; v1[2]=0.;
    dax=v_angle(v0,v1);
    v0[0]=node[name[1]].nx* scale->w+scale->x;
    v0[1]=0.;
    v0[2]=node[name[1]].nz* scale->w+scale->z;
    v1[0]=0.; v1[1]=0.; v1[2]=1.;
    day=v_angle(v0,v1);
    v0[0]=node[name[1]].nx* scale->w+scale->x;
    v0[1]=node[name[1]].ny* scale->w+scale->y;
    v0[2]=0.;
    v1[0]=1.; v1[1]=0.; v1[2]=0.;
    daz=v_angle(v0,v1);
    v0[0]=node[name[1]].nx* scale->w+scale->x;
    v0[1]=node[name[1]].ny* scale->w+scale->y;
    v0[2]=node[name[1]].nz* scale->w+scale->z;
    printf("axyz= %lf %lf %lf rxyz= %lf %lf %lf\n",
      dax*180./PI,day*180./PI,daz*180./PI,
      sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
      sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
      sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].node,set[i].anz_n,name[1]) >-1))
          {
            sprintf(&printbuf[strlen(printbuf)],"%s(%d)",set[i].name,set[i].index);  

            /* is a contact surface referenced ? */
            if(set[i].anz_se)
            {
              sprintf(&printbuf[strlen(printbuf)],"->[ ");
              for(j=0; j<set[i].anz_se; j++) 
                sprintf(&printbuf[strlen(printbuf)],"%s(%d) ",set[set[i].set[j]].name,set[set[i].set[j]].index);
              sprintf(&printbuf[strlen(printbuf)],"]");
            }
            sprintf(&printbuf[strlen(printbuf)],",");
          }
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdis", 4) == 4)
    {
      if( type[1] == 'h') /* plane */
      {
        cpFlag++;
        v_qdis[cpFlag][0]=node[name[1]].nx* scale->w+scale->x;
        v_qdis[cpFlag][1]=node[name[1]].ny* scale->w+scale->y;
        v_qdis[cpFlag][2]=node[name[1]].nz* scale->w+scale->z;
        if(cpFlag==3) { type[1]=0; cpFlag=0; shFlag++; }
      }
      
      else if( type[1] == 'c') /* set centerpnt */
      {
        type[1]=0;
        cpFlag=1;
        v_qdis[2][0]=node[name[1]].nx* scale->w+scale->x;
        v_qdis[2][1]=node[name[1]].ny* scale->w+scale->y;
        v_qdis[2][2]=node[name[1]].nz* scale->w+scale->z;
      }
      
      else if(qdisFlag>0)
      {
        qdisFlag=0;
	shFlag=0;
        v_qdis[1][0]=node[name[1]].nx* scale->w+scale->x;
        v_qdis[1][1]=node[name[1]].ny* scale->w+scale->y;
        v_qdis[1][2]=node[name[1]].nz* scale->w+scale->z;
        if(cpFlag==1)
        {
          cpFlag=0;  
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]-v_qdis[2][i]; v1[i]=v_qdis[1][i]-v_qdis[2][i]; }
          da=v_angle(v0,v1);
          r1=v_betrag(v0);
          r2=v_betrag(v1);
          rm=(r1+r2)*.5;
          lcir=2.*rm*da;
          printf(" lcir:%lf da:%lf dr:%lf r1:%lf r2:%lf\n", lcir, da*180./PI, r2-r1, r1, r2 );
        }
        else
        {
          da=v_angle(v_qdis[0],v_qdis[1]);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          v0[0]=v1[0]=0.; 
          dax=v_angle(v0,v1);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          v0[1]=v1[1]=0.; 
          day=v_angle(v0,v1);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          v0[2]=v1[2]=0.; 
          daz=v_angle(v0,v1);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          if(!anz->l) 
          {
	    dist=sqrt((v_qdis[1][0]-v_qdis[0][0])*(v_qdis[1][0]-v_qdis[0][0])+
            (v_qdis[1][1]-v_qdis[0][1])*(v_qdis[1][1]-v_qdis[0][1])+
		      (v_qdis[1][2]-v_qdis[0][2])*(v_qdis[1][2]-v_qdis[0][2]));
	    if(dist>=1e-4)
	    {
            printf("dist:%lf dxyz: %lf %lf %lf da:%lf daxyz: %lf %lf %lf dr:%lf drxyz: %lf %lf %lf\n",
		   dist,v_qdis[1][0]-v_qdis[0][0],v_qdis[1][1]-v_qdis[0][1],v_qdis[1][2]-v_qdis[0][2],
            da*180./PI,dax*180./PI,day*180./PI,daz*180./PI,
            sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[1]*v1[1])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
	    }
	    else
	    {
              printf("dist:%e dxyz: %e %e %e da:%lf daxyz: %lf %lf %lf dr:%e drxyz: %e %e %e\n",
		     dist,v_qdis[1][0]-v_qdis[0][0],v_qdis[1][1]-v_qdis[0][1],v_qdis[1][2]-v_qdis[0][2],
            da*180./PI,dax*180./PI,day*180./PI,daz*180./PI,
            sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[1]*v1[1])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
	    }
          }
  	  else
          {
            dist=sqrt((v_qdis[1][0]-v_qdis[0][0])*(v_qdis[1][0]-v_qdis[0][0])+
            (v_qdis[1][1]-v_qdis[0][1])*(v_qdis[1][1]-v_qdis[0][1])+
		 (v_qdis[1][2]-v_qdis[0][2])*(v_qdis[1][2]-v_qdis[0][2])); 
	    if(dist>=1e-4)
	    {
              printf("dist:%lf dv:%lf dxyz: %lf %lf %lf da:%lf daxyz: %lf %lf %lf dr:%lf drxyz: %lf %lf %lf\n",
		   dist,lcase[cur_lc].dat[cur_entity][name[1]]-v_val,
                 v_qdis[1][0]-v_qdis[0][0],v_qdis[1][1]-v_qdis[0][1],v_qdis[1][2]-v_qdis[0][2],
            da*180./PI,dax*180./PI,day*180./PI,daz*180./PI,
            sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[1]*v1[1])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
	    }
	    else
	    {
              printf("dist:%e dv:%e dxyz: %e %e %e da:%lf daxyz: %lf %lf %lf dr:%e drxyz: %e %e %e\n",
		   dist,lcase[cur_lc].dat[cur_entity][name[1]]-v_val,
                 v_qdis[1][0]-v_qdis[0][0],v_qdis[1][1]-v_qdis[0][1],v_qdis[1][2]-v_qdis[0][2],
            da*180./PI,dax*180./PI,day*180./PI,daz*180./PI,
            sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
            sqrt(v1[0]*v1[0]+v1[1]*v1[1])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
	    }
          }
        }
      }
      
      else
      {
	if(qaddTol>-1) // distance of that node to a set of surfaces or faces
	{
	  qaddTol=-1;
	  sprintf(buffer,"%d",name[1]);
	  delSet(specialset->mpc);
          pre_seta( specialset->mpc, "n", buffer);
	  sprintf(buffer,"dist %s %s nor",specialset->mpc,specialset->impc);
	  pre_proj(buffer);
	  //delSet(specialset->mpc);
	  //delSet(specialset->impc);
	}
	else
	{
          qdisFlag=2; shFlag++;
          if(anz->l) v_val=lcase[cur_lc].dat[cur_entity][name[1]];
          v_qdis[0][0]=node[name[1]].nx* scale->w+scale->x;
          v_qdis[0][1]=node[name[1]].ny* scale->w+scale->y;
          v_qdis[0][2]=node[name[1]].nz* scale->w+scale->z;
	}
      }
      
      if(shFlag==2) // calc the normal dist of v_qdis[0] to the shape v_qdis[1] v_qdis[2] v_qdis[3]
      {
        qdisFlag=0;
	shFlag=0;
	printf("normal dist:%lf\n",normdist(&v_qdis[0][0],&v_qdis[1][0],&v_qdis[2][0],&v_qdis[3][0]));
      }
    }
    else if (compare(pickfunc, "qcnt", 4) == 4)
    {
      centerNode=name[1];
      center( node[name[1]].nx, node[name[1]].ny, node[name[1]].nz);
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      glutSetWindow( w1);
      /* Keyboardfunktion wiederherstellen (von pick())  */
      pickFlag=0;
      free(selectBuf); selectBuf=NULL;
      glutKeyboardFunc ( Keyboard );
      glutSetWindow( activWindow );
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) { if(set[setNrbuf].type) seqr( setNrbuf, type, name[1]); else setr( setNrbuf, type, name[1]); }
    else if (compare(pickfunc, "qmsh", 4) == 4) entitybuf=createElem(name[1]);
    else if (compare(pickfunc, "qmov", 4) == 4) moveSet(name, x, y);
    else if (compare(pickfunc, "qnod", 4) == 4) moveNode(name[1], x, y);
    else if (compare(pickfunc, "qpnt", 4) == 4) movePoint(name, x, y);
    else if (compare(pickfunc, "qali", 4) == 4) pre_align( name[1], 0 );
    else if (compare(pickfunc, "qcut", 4) == 4) pre_cut( name[1], type[1] );
    else if (compare(pickfunc, "qtxt", 4) == 4) moveText(createText(name[1], -1, 0 ), x, y );
    else if (compare(pickfunc, "qshp", 4) == 4)
    {
      /* create point on node */
      p= getNewName( buf, "p" );
      printf(" create point:%s %lf %lf %lf\n", buf, node[name[1]].nx, node[name[1]].ny, node[name[1]].nz );
      p=pnt( buf, node[name[1]].nx, node[name[1]].ny, node[name[1]].nz, 0 );
      shapebuf.p[shp_pindx++]=p;
      printf(" %d pnt picked\n",shp_pindx);
      if(shp_pindx>6)
      {
        shp_pindx=0;
      }
    }
  }
  else if (name[0] == 'p')
  {
    printf ("%s xyz= %lf %lf %lf ", point[name[1]].name,
      (point[name[1]].px* scale->w+scale->x),
      (point[name[1]].py* scale->w+scale->y),
      (point[name[1]].pz* scale->w+scale->z) );

      v0[0]=0.;
      v0[1]=point[name[1]].py* scale->w+scale->y;
      v0[2]=point[name[1]].pz* scale->w+scale->z;
      v1[0]=0.; v1[1]=1.; v1[2]=0.;
      dax=v_angle(v0,v1);
      v0[0]=point[name[1]].px* scale->w+scale->x;
      v0[1]=0.;
      v0[2]=point[name[1]].pz* scale->w+scale->z;
      v1[0]=0.; v1[1]=0.; v1[2]=1.;
      day=v_angle(v0,v1);
      v0[0]=point[name[1]].px* scale->w+scale->x;
      v0[1]=point[name[1]].py* scale->w+scale->y;
      v0[2]=0.;
      v1[0]=1.; v1[1]=0.; v1[2]=0.;
      daz=v_angle(v0,v1);
      v0[0]=point[name[1]].px* scale->w+scale->x;
      v0[1]=point[name[1]].py* scale->w+scale->y;
      v0[2]=point[name[1]].pz* scale->w+scale->z;
      printf("axyz= %lf %lf %lf rxyz= %lf %lf %lf\n",
        dax*180./PI,day*180./PI,daz*180./PI,
        sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
        sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
        sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].pnt,set[i].anz_p,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),", set[i].name,set[i].index); 
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qnor", 4) == 4)
    {
      if(qnorCounter==2)
      {
        p_qnor[0]=name[1];
        qnorCounter=0;
        getNewName(buffer,"l");
        if((i=normalLine(buffer, p_qnor, pickbuf))<0)
        { printf("qnor: could not create new line\n"); }
        else { printf(" create line:%s with new end-point %s\n", buffer, point[i].name ); }
      }
      else p_qnor[++qnorCounter]=name[1];
    }
    else if (compare(pickfunc, "qdis", 4) == 4)
    {
      if( type[1] == 'h') /* plane */
      {
        cpFlag++;
        v_qdis[cpFlag][0]=point[name[1]].px* scale->w+scale->x;
        v_qdis[cpFlag][1]=point[name[1]].py* scale->w+scale->y;
        v_qdis[cpFlag][2]=point[name[1]].pz* scale->w+scale->z;
        if(cpFlag==3) { type[1]=0; cpFlag=0; shFlag++; }
      }
      
      else if( type[1] == 'c') /* set centerpnt */
      {
        type[1] =0;
        cpFlag=1;
        v_qdis[2][0]=point[name[1]].px* scale->w+scale->x;
        v_qdis[2][1]=point[name[1]].py* scale->w+scale->y;
        v_qdis[2][2]=point[name[1]].pz* scale->w+scale->z;
      }
      
      else if(qdisFlag>0)
      {
        qdisFlag=0;
	shFlag=0;
        v_qdis[1][0]=point[name[1]].px* scale->w+scale->x;
        v_qdis[1][1]=point[name[1]].py* scale->w+scale->y;
        v_qdis[1][2]=point[name[1]].pz* scale->w+scale->z;
        if(cpFlag==1)
        {
          cpFlag=0;  
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]-v_qdis[2][i]; v1[i]=v_qdis[1][i]-v_qdis[2][i]; }
          da=v_angle(v0,v1);
          r1=v_betrag(v0);
          r2=v_betrag(v1);
          rm=(r1+r2)*.5;
          lcir=2.*rm*da;
          printf(" lcir:%lf da:%lf dr:%lf r1:%lf r2:%lf\n", lcir, da*180./PI, r2-r1, r1, r2 );
        }
        else
        {
          da=v_angle(v_qdis[0],v_qdis[1]);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          v0[0]=v1[0]=0.; 
          dax=v_angle(v0,v1);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          v0[1]=v1[1]=0.; 
          day=v_angle(v0,v1);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          v0[2]=v1[2]=0.; 
          daz=v_angle(v0,v1);
          for(i=0; i<3; i++) { v0[i]=v_qdis[0][i]; v1[i]=v_qdis[1][i]; }
          printf("dist= %lf dxyz= %lf %lf %lf da= %lf daxyz= %lf %lf %lf dr=%lf drxyz= %lf %lf %lf\n",
          sqrt((v_qdis[1][0]-v_qdis[0][0])*(v_qdis[1][0]-v_qdis[0][0])+
          (v_qdis[1][1]-v_qdis[0][1])*(v_qdis[1][1]-v_qdis[0][1])+
          (v_qdis[1][2]-v_qdis[0][2])*(v_qdis[1][2]-v_qdis[0][2])),
          v_qdis[1][0]-v_qdis[0][0],v_qdis[1][1]-v_qdis[0][1],v_qdis[1][2]-v_qdis[0][2],
          da*180./PI, dax*180./PI, day*180./PI, daz*180./PI,
          sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]),
          sqrt(v1[1]*v1[1]+v1[2]*v1[2])-sqrt(v0[1]*v0[1]+v0[2]*v0[2]),
          sqrt(v1[0]*v1[0]+v1[2]*v1[2])-sqrt(v0[0]*v0[0]+v0[2]*v0[2]),
          sqrt(v1[0]*v1[0]+v1[1]*v1[1])-sqrt(v0[0]*v0[0]+v0[1]*v0[1]) );
        }
      }
      
      else
      {
	if(qaddTol>-1) // distance of that node to a set of surfaces or faces
	{
	  qaddTol=-1;
	  delSet(specialset->mpc);
          pre_seta( specialset->mpc, "p", point[name[1]].name);
	  sprintf(buffer,"dist %s %s nor",specialset->mpc,specialset->impc);
	  pre_proj(buffer);
	  //delSet(specialset->mpc);
	  //delSet(specialset->impc);
	}
	else
	{
          qdisFlag=1; shFlag++;
          v_qdis[0][0]=point[name[1]].px* scale->w+scale->x;
          v_qdis[0][1]=point[name[1]].py* scale->w+scale->y;
          v_qdis[0][2]=point[name[1]].pz* scale->w+scale->z;
        }
      }

      if(shFlag==2) // calc the normal dist of v_qdis[0] to the shape v_qdis[1] v_qdis[2] v_qdis[3]
      {
        qdisFlag=0;
	shFlag=0;
	printf("normal dist:%lf\n",normdist(&v_qdis[0][0],&v_qdis[1][0],&v_qdis[2][0],&v_qdis[3][0]));
      }      
    }
    else if (compare(pickfunc, "qlin", 4) == 4)
    {
      if( type[1] == 'b') entitybuf=createLine( point[name[1]].name, 0 ); /* start line */
      if( type[1] == 'g') entitybuf=createLine( point[name[1]].name, 1 ); /* end  line */
      if( type[1] == 'c') entitybuf=createLine( point[name[1]].name, 2 ); /* add centerpnt */
      if( type[1] == 'm') entitybuf=createLine( point[name[1]].name, 3 ); /* midpnt, create centerpnt */
      if( type[1] == 't') entitybuf=createLine( point[name[1]].name, 4 ); /* seq, create set */
      if( type[1] == 'p') entitybuf=createLine( point[name[1]].name, 5 ); /* change length of line (displace) */
    }
    else if (compare(pickfunc, "qcnt", 4) == 4)
    {
      center( point[name[1]].px, point[name[1]].py, point[name[1]].pz);
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      glutSetWindow( w1);
      /* Keyboardfunktion wiederherstellen (von pick())  */
      pickFlag=0;
      free(selectBuf); selectBuf=NULL;
      glutKeyboardFunc ( Keyboard );
      glutSetWindow( activWindow );
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) { if(set[setNrbuf].type) seqr( setNrbuf, type, name[1]); else setr( setNrbuf, type, name[1]); }
    else if (compare(pickfunc, "qmov", 4) == 4) moveSet(name, x, y);
    else if (compare(pickfunc, "qpnt", 4) == 4)
    {
      if((surfNr>-1)||(nurbsNr>-1)||(shapeNr>-1))  /* move points to a surface */
        seta( setNrbuf, type, name[1]);
      else movePoint(name, x, y);
    }
    else if (compare(pickfunc, "qali", 4) == 4) pre_align( name[1], 1 );
    else if (compare(pickfunc, "qcut", 4) == 4) pre_cut( name[1], type[1] );
    else if (compare(pickfunc, "qshp", 4) == 4)
    {
      shapebuf.p[shp_pindx++]=name[1];
      printf(" %d pnt picked\n",shp_pindx);
      if(shp_pindx>6)
      {
        shp_pindx=0;
      }
    }
  }
  else if (name[0] == 'l')
  {
    if (line[name[1]].name == (char *)NULL) return(0);

    if (compare(pickfunc, "qenq", 4) == 4)
    {
      if(type[1]=='c')
      {
        /* lcmb must be identified and added */
        for(i=0; i<anzGeo->c; i++)
        {
          if( lcmb[i].name != (char *)NULL )
          {
            for(j=0; j<lcmb[i].nl; j++) if(lcmb[i].l[j]==name[1])
            {
              printf (" lcmb:%s ", lcmb[i].name );
              for (k=0; k<lcmb[i].nl; k++)
	      {
                printf (" %1c %s", lcmb[i].o[k], line[lcmb[i].l[k]].name );
	      }
              printf (" \n");
            }
          }
        }
        return(1);
      }
    }
    else if (compare(pickfunc, "qadd", 4) == 4)
    {
      if(type[1]=='c')
      {
        /* lcmb must be identified and added */
        for(i=0; i<anzGeo->c; i++)
        {
          if( lcmb[i].name != (char *)NULL )
          {
            for(j=0; j<lcmb[i].nl; j++) if(lcmb[i].l[j]==name[1])
            {
              printf (" lcmb:%s ", lcmb[i].name );
              for (k=0; k<lcmb[i].nl; k++)
	      {
                printf (" %1c %s", lcmb[i].o[k], line[lcmb[i].l[k]].name );
                seta( setNrbuf, "l",lcmb[i].l[k] );
	      }
              printf (" \n");
              seta( setNrbuf, "c", i);
            }
          }
        }
        return(1);
      }
      else  seta( setNrbuf, type, name[1]);
    }
    else if (compare(pickfunc, "qdel", 4) == 4)
    {
      if(type[1]=='c')
      {
        /* lcmb must be identified and added */
        for(i=0; i<anzGeo->c; i++)
        {
          if( lcmb[i].name != (char *)NULL )
          {
            for(j=0; j<lcmb[i].nl; j++) if(lcmb[i].l[j]==name[1])
            {
              printf (" lcmb:%s ", lcmb[i].name );
              seta( setNrbuf, "c", i);
            }
          }
        }
        return(1);
      }
      else  seta( setNrbuf, type, name[1]);
    }
    else if (compare(pickfunc, "qrem", 4) == 4)
    {
      if(type[1]=='c')
      {
        /* lcmb must be identified and added */
        for(i=0; i<anzGeo->c; i++)
        {
          if( lcmb[i].name != (char *)NULL )
          {
            for(j=0; j<lcmb[i].nl; j++) if(lcmb[i].l[j]==name[1])
            {
              printf (" lcmb:%s ", lcmb[i].name );
              for (k=0; k<lcmb[i].nl; k++)
	      {
                printf (" %1c %s", lcmb[i].o[k], line[lcmb[i].l[k]].name );
                setr( setNrbuf, "l",lcmb[i].l[k] );
	      }
              printf (" \n");
              setr( setNrbuf, "c", i);
            }
          }
        }
        return(1);
      }
      else setr( setNrbuf, type, name[1]); 
    }
    else if ((compare(pickfunc, "qbia", 4) == 4)&&(line[name[1]].div>1)) 
    {
      if(type[3]=='c')
        line[name[1]].bias=1./line[name[1]].bias;
      else
      {
        if(line[name[1]].bias<1.)
          line[name[1]].bias = 1./pow((double)pickbuf, (1./((double)line[name[1]].div-1.)));
        else
          line[name[1]].bias = pow((double)pickbuf, (1./((double)line[name[1]].div-1.)));
      }
      repLine(name[1]);
    }
    else if (compare(pickfunc, "qdiv", 4) == 4)
    {
      if((line[name[1]].bias!=1.)&&(pickbuf>1))
      {
        bias_fbd=getBias_fbd(name[1], line);
        // printf(" bias_fbd:%d\n", bias_fbd);

        // code from splitBiasDiv()
        bias_fbd*=10;

        /* in the fbd-format-definition bias is defined as bias(a)=(elem_length(last)/elem_length(first)) */
        /* therefore a re-calculation is necessary because in cgx the definition is: */ 
        /* bias(b) = (elem_length(n+1)/elem_length(n)) => bias(a)= bias(b)**(div-1) */
        if(pickbuf>1)
        {
          if(bias_fbd<0) line[name[1]].bias= 1./pow(((double)bias_fbd*-.1), (1./((double)pickbuf-1.)));
          else           line[name[1]].bias= pow(((double)bias_fbd*.1), (1./((double)pickbuf-1.)));
        }
        else line[name[1]].bias=1.;
      }
      else if(pickbuf==1) line[name[1]].bias=1.;
      if(pickbuf>0) line[name[1]].div=pickbuf;
      repLine(name[1]);
    }

    bias_fbd=getBias_fbd(name[1], line);

    if (line[name[1]].typ=='a')
      printf ("%s p1:%s p2:%s pc:%s div:%d bias:%d bias_el:%lf"
      , line[name[1]].name, point[line[name[1]].p1].name, point[line[name[1]].p2].name
      , point[line[name[1]].trk].name, line[name[1]].div, bias_fbd, line[name[1]].bias  );
    else if (line[name[1]].typ=='s')
      printf ("%s p1:%s p2:%s set:%s div:%d bias:%d bias_el:%lf"
      , line[name[1]].name, point[line[name[1]].p1].name, point[line[name[1]].p2].name
      , set[line[name[1]].trk].name, line[name[1]].div, bias_fbd, line[name[1]].bias  );
    else if (line[name[1]].typ==' ')
      printf ("%s p1:%s p2:%s div:%d bias:%d bias_el:%lf"
      , line[name[1]].name, point[line[name[1]].p1].name, point[line[name[1]].p2].name
      , line[name[1]].div, bias_fbd, line[name[1]].bias );
    printf(" etyp:%d attr:%d lock:%d\n", line[name[1]].etyp, line[name[1]].eattr, line[name[1]].elock);
    
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].line,set[i].anz_l,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qmsh", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qint", 4) == 4) 
    {
      if(intersect( name[1])==0) for (j=0; j<set[setall].anz_l; j++) repLine(set[setall].line[j]);
    }
    else if (compare(pickfunc, "qfil", 4) == 4)
    {
      if(createFillet( name[1], filletRadius)==0) for (j=0; j<set[setall].anz_l; j++) repLine(set[setall].line[j]);
    }
    else if (compare(pickfunc, "qseq", 4) == 4) seta( setNrbuf, type, name[1]);
    else if ((compare(pickfunc, "qspl", 4) == 4)||((compare(pickfunc, "qlin", 4) == 4)&&(type[1]=='s') ))
    {
      for(j=0; j<qspl_i; j++) if(qspl_indx[j]==name[1]) return(0);
      if (qspl_i<MAX_LINES) qspl_indx[qspl_i]=name[1];
      else { printf(" ERROR in hitAction(), to much lines picked. Increase MAX_LINES\n"); return(0); }
      qspl_i++;
      k= qsplitLine( name[1], x, y );
      if(k>-1) seta( setNrbuf, "p", k);
      if(set[setNrbuf].anz_p>1)
      {
        sprintf( buffer,"p %s %lf", set[setNrbuf].name, gtol);
        pre_merge( buffer);
      }
    }
    else if (compare(pickfunc, "qsur", 4) == 4)
    {
      if(type[1]=='l') pickstack(name[1]);
      else strcpy(type, line[name[1]].name);
    }
    else if ((compare(pickfunc, "qlin", 4) == 4)&&(type[1]=='e') )
    {
      if(entitybuf==name[1]) entitybuf=lineNr;
      if(entitybuf<0) return(0);
      if(entitybuf==name[1]) return(0);

      printf ("try to replace line %s by %s\n", line[name[1]].name, line[entitybuf].name);
      if(name[1] !=entitybuf)
      {
        /* versuche eine linie direkt auszutauschten */
        if((line[name[1]].p1==line[entitybuf].p1)&&(line[name[1]].p2==line[entitybuf].p2))
	{
          line[entitybuf].div=line[name[1]].div;
          /* check lcmbs for the line[i]  */
          for (j=0; j<anzGeo->c; j++) if( lcmb[j].name != (char *)NULL )
          {
            for (k=0; k<lcmb[j].nl; k++)
            {
              if (lcmb[j].l[k] == name[1])
	      {
                lcmb[j].l[k]=entitybuf;
                printf ("replace %s in lcmb %s by %s\n", line[name[1]].name, lcmb[j].name, line[entitybuf].name);
	      }
            }
          }
          /* check surfs for the line[i]  */
          for (j=0; j<anzGeo->s; j++) if( surf[j].name != (char *)NULL )
          {
            for (k=0; k<surf[j].nl; k++)
            {
              if ((surf[j].l[k] == name[1])&&(surf[j].typ[k]=='l'))
	      {
                surf[j].l[k]=entitybuf;
                printf ("replace %s in surf %s by %s\n", line[name[1]].name, surf[j].name, line[entitybuf].name);
	      }
            }
          }
          pre_seta( "-delete", "l", line[name[1]].name);
          return(1); 
        }
        else if((line[name[1]].p1==line[entitybuf].p2)&&(line[name[1]].p2==line[entitybuf].p1))
	{
          line[entitybuf].div=line[name[1]].div;
          /* check lcmbs for the line[i]  */
          for (j=0; j<anzGeo->c; j++) if( lcmb[j].name != (char *)NULL )
          {
            for (k=0; k<lcmb[j].nl; k++)
            {
              if (lcmb[j].l[k] == name[1])
              {
                lcmb[j].l[k]=entitybuf;
                printf ("replace %s in lcmb %s by %s\n", line[name[1]].name, lcmb[j].name, line[entitybuf].name);
                if(lcmb[j].o[k]=='+') lcmb[j].o[k]='-';
                else                 lcmb[j].o[k]='+';
              }
            }
          }
          /* check surfs for the line[i]  */
          for (j=0; j<anzGeo->s; j++) if( surf[j].name != (char *)NULL )
          {
            for (k=0; k<surf[j].nl; k++)
            {
              if ((surf[j].l[k] == name[1])&&(surf[j].typ[k]=='l'))
              {
                printf ("replace %s in surf %s by %s\n", line[name[1]].name, surf[j].name, line[entitybuf].name);
                surf[j].l[k]=entitybuf;
                if(surf[j].o[k]=='+') surf[j].o[k]='-';
                else                 surf[j].o[k]='+';
              }
            }
          }
          pre_seta( "-delete", "l", line[name[1]].name);
          return(1); 
	}

        /* untersuche alle lcmbs ob line ein Mitglied ist */
        for (i=0; i<anzGeo->c; i++) if( lcmb[i].name != (char *)NULL )
        {
          for (j=0; j<lcmb[i].nl; j++)
          {
            if( name[1] == lcmb[i].l[j] )
            {
              /* compare the end-points to find out if the lcmb has the same range */
              matchFlag=0;
	      if(lcmb[i].o[0]=='+')
              {  if((line[lcmb[i].l[0]].p1==line[entitybuf].p1)||(line[lcmb[i].l[0]].p1==line[entitybuf].p2)) matchFlag=1; }
	      else
              {  if((line[lcmb[i].l[0]].p2==line[entitybuf].p1)||(line[lcmb[i].l[0]].p2==line[entitybuf].p2)) matchFlag=1; }


	      if(lcmb[i].o[lcmb[i].nl-1]=='-')
	      {if((line[lcmb[i].l[lcmb[i].nl-1]].p1==line[entitybuf].p1)||(line[lcmb[i].l[lcmb[i].nl-1]].p1==line[entitybuf].p2)) matchFlag++; }
	      else
	      { if((line[lcmb[i].l[lcmb[i].nl-1]].p2==line[entitybuf].p1)||(line[lcmb[i].l[lcmb[i].nl-1]].p2==line[entitybuf].p2)) matchFlag++; }


              /*
	      if(matchFlag) printf("direct match o1:%c lcmblp:%s %s lp:%s %s\n", lcmb[i].o[0], point[line[lcmb[i].l[0]].p1].name,point[line[lcmb[i].l[0]].p2].name,point[line[entitybuf].p1].name,point[line[entitybuf].p2].name);
	      if(matchFlag==2) printf("direct match o2:%c lcmblp:%s %s lp:%s %s\n", lcmb[i].o[lcmb[i].nl-1], point[line[lcmb[i].l[lcmb[i].nl-1]].p1].name,point[line[lcmb[i].l[lcmb[i].nl-1]].p2].name,point[line[entitybuf].p1].name,point[line[entitybuf].p2].name);
	      */


       	      if(matchFlag==2)
	      {
                line[entitybuf].div=0;
                for(l=0; l<lcmb[i].nl; l++)
                {
                  line[entitybuf].div+=line[lcmb[i].l[l]].div;
                  pre_seta( "-delete", "l", line[lcmb[i].l[l]].name);
  	        }
                repLine(entitybuf);
                printf("complete edge (lcmb:%s) will be deleted and replaced by the new line:%s \n", lcmb[i].name, line[entitybuf].name );
                if ((lcmb[i].o = (char *)realloc( (char *)lcmb[i].o, (1)*sizeof(char)) ) == NULL )
                { printf("\n\n ERROR: realloc failure in qspl, lcmb.o:%s not changed\n\n",lcmb[i].name ); return(0); }
                if ((lcmb[i].l = (int *)realloc( (int *)lcmb[i].l, (1)*sizeof(int)) ) == NULL )
                { printf("\n\n ERROR: realloc failure in qspl, lcmb.l:%s not changed\n\n", lcmb[i].name); return(0); }
      
                lcmb[i].o[0]='+';
                lcmb[i].l[0]=entitybuf;
                lcmb[i].nl=1;
                printf ("replace lines in lcmb %s by %s\n", lcmb[i].name, line[entitybuf].name);
                return(1); 
              }
            }
          }
        }

        /* no exactly matching lcmb was found, go again over all lcmb and replace a line-sequence in a matching lcmb */
        for (i=0; i<anzGeo->c; i++) if( lcmb[i].name != (char *)NULL )
        {
          for (j=0; j<lcmb[i].nl; j++)
          {
            if( name[1] == lcmb[i].l[j] )
            {
              /* search the first and the last line which match with the new line */
              matchFlag=0;
              for (k=0; k<lcmb[i].nl; k++)
              {
                if(!matchFlag)
		{
	          if(lcmb[i].o[k]=='+')
                  {  if((line[lcmb[i].l[k]].p1==line[entitybuf].p1)||(line[lcmb[i].l[k]].p1==line[entitybuf].p2)) matchFlag++; lin[0]=k;  }
	          else
                  {  if((line[lcmb[i].l[k]].p2==line[entitybuf].p1)||(line[lcmb[i].l[k]].p2==line[entitybuf].p2)) matchFlag++; lin[0]=k;  }
		}
                else
		{
	          if(lcmb[i].o[k]=='-')
                  { if((line[lcmb[i].l[k]].p1==line[entitybuf].p1)||(line[lcmb[i].l[k]].p1==line[entitybuf].p2)) { matchFlag++; lin[1]=k; } }
	          else
		  { if((line[lcmb[i].l[k]].p2==line[entitybuf].p1)||(line[lcmb[i].l[k]].p2==line[entitybuf].p2)) { matchFlag++; lin[1]=k;  } }
		}
	      }

	      /*
	      if(matchFlag)
                printf("o1:%c l:%s p:%s %s lp:%s %s\n", lcmb[i].o[lin[0]], line[lcmb[i].l[lin[0]]].name, point[line[lcmb[i].l[lin[0]]].p1].name,point[line[lcmb[i].l[lin[0]]].p2].name,point[line[entitybuf].p1].name,point[line[entitybuf].p2].name);
	      if(matchFlag==2)
                printf("o2:%c l:%s p:%s %s lp:%s %s\n", lcmb[i].o[lin[1]], line[lcmb[i].l[lin[1]]].name, point[line[lcmb[i].l[lin[1]]].p1].name,point[line[lcmb[i].l[lin[1]]].p2].name,point[line[entitybuf].p1].name,point[line[entitybuf].p2].name);
	      */

       	      if(matchFlag==2) 
	      {
                line[entitybuf].div=0;
                for(l=lin[0]; l<=lin[1]; l++)
                {
                  line[entitybuf].div+=line[lcmb[i].l[l]].div;
                  pre_seta( "-delete", "l", line[lcmb[i].l[l]].name);
  	        }
                repLine(entitybuf);
  
                lcmb[i].l[lin[0]]=entitybuf;
                lcmb[i].o[lin[0]]='+';
                k=lin[0];
                for(l=lin[1]+1; l<lcmb[i].nl; l++)
  	        {
                  k++;
                  lcmb[i].o[k]=lcmb[i].o[l];
                  lcmb[i].l[k]=lcmb[i].l[l];
                }
                lcmb[i].nl=k+1; 
                printf ("replace lines in lcmb %s by %s\n", lcmb[i].name, line[entitybuf].name);
                return(1); 
              }
            }
          }
        }
      }
    }
    else if ((compare(pickfunc, "qlin", 4) == 4)&&(type[1]=='x') )
    {
      printf("Line redefined to be straight\n");
      line[name[1]].typ=' ';
      line[name[1]].trk=-1;
      repLine(name[1]);
    }
    else if (compare(pickfunc, "qlin", 4) == 4) lineNr=name[1];
  }
  /*
  // lcmbs do not exist as graphical entities and can not be selected so far (TBD)
  else if (name[0] == 'c')
  {
      printf("hallo ch\n");
    if (lcmb[name[1]].name == (char *)NULL) return(0);
    else if (compare(pickfunc, "qseq", 4) == 4) convertLCMB( name[1] );
  }
  */
  else if (name[0] == 's')
  {
    if (surf[name[1]].sh<0) printf ("%s %1c BLEND ", surf[name[1]].name, surf[name[1]].ori );
    else printf ("%s %1c %s ", surf[name[1]].name, surf[name[1]].ori, shape[surf[name[1]].sh].name );
    for (j=0; j<surf[name[1]].nl; j++)
    {
      if (surf[name[1]].typ[j]=='l')
      {
        printf (" %1c %s", surf[name[1]].o[j], line[surf[name[1]].l[j]].name );
      }
      else
      {
        printf (" %1c %s", surf[name[1]].o[j], lcmb[surf[name[1]].l[j]].name );
        for (k=0; k<lcmb[surf[name[1]].l[j]].nl; k++ )
        {
          printf (" (%c) (%s)", lcmb[surf[name[1]].l[j]].o[k], line[lcmb[surf[name[1]].l[j]].l[k]].name );
        }
      }
    }
    if( surf[name[1]].eparm != (char *)NULL )
      printf (" etyp:%d attr:%d mshp:%s lock:%d\n", surf[name[1]].etyp, surf[name[1]].eattr, surf[name[1]].eparm, surf[name[1]].elock);
    else
      printf (" etyp:%d attr:%d lock:%d\n", surf[name[1]].etyp, surf[name[1]].eattr, surf[name[1]].elock);
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].surf,set[i].anz_s,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qmsh", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qdis", 4) == 4)
    {
      delSet(specialset->impc);
      setNrbuf=pre_seta( specialset->impc, "i", 0);
      if(setNrbuf<0) { errMsg(" ERROR: could not create set for qdis\n"); return(0); } 
      seta( setNrbuf, type, name[1]);
      qaddTol=0;  // marker
    }
    else if (compare(pickfunc, "qsur", 4) == 4)
    {
      /* change the interior def to BLEND */ 
      if(type[1]=='b')
      {
        j=surf[name[1]].sh;
        if(j>-1) shape[j].ns=iremove(&shape[j].s, shape[j].ns, name[1]);
        if(surf[name[1]].pgn!=NULL)
	{
          surf[name[1]].sh=-1; free(surf[name[1]].pgn); surf[name[1]].pgn=NULL; surf[name[1]].npgn=0;
          repSurf(name[1],1);
	}
        else
        {
          surf[name[1]].sh=-1; free(surf[name[1]].pgn); surf[name[1]].pgn=NULL; surf[name[1]].npgn=0;
	}
        printf (" interior changed to BLEND\n");
      }
      else if(type[1]=='h')
      {
        shapeNr=surf[name[1]].sh;
	if( shapeNr>-1)
	{
          //if(shape[surf[name[1]].sh].type==4) nurbsNr=shape[surf[name[1]].sh].p[0];
          printf (" default shape set to %s\n",shape[shapeNr].name);
	}
        else printf(" surf has no related shape\n");
      }
      else if(type[1]=='S')
      {
        shapeNr=surf[name[1]].sh;
	if( shapeNr>-1)
	{
          if(shape[surf[name[1]].sh].type==4) nurbsNr=shape[surf[name[1]].sh].p[0];
          printf (" default nurbs set to %s (referenced by shape:%s)\n",nurbs[nurbsNr].name, shape[nurbsNr].name);
	}
        else printf(" surf has no related nurbs\n");
      }
      else if(type[1]=='s')
      {
	/* if a Nurbs or a shape was chosen change the interior def of the selected surf */
        if(nurbsNr>-1)
	{
          surf[name[1]].sh=getShapeNr(nurbs[nurbsNr].name);
          shape_refSurf(surf[name[1]].sh, name[1]);
          if(surf[name[1]].sh>-1)
	  {
            if(surf[name[1]].pgn!=NULL)
	    {
              free(surf[name[1]].pgn); surf[name[1]].pgn=NULL; surf[name[1]].npgn=0;

              /* create the interior */
              repSurf(name[1],1);
	    }
            printf (" interior changed to nurbs: %s (referenced by shape:%s)\n",nurbs[nurbsNr].name, shape[surf[name[1]].sh].name );
	  }
	}
        else if(shapeNr>-1)
	{
          surf[name[1]].sh=shapeNr;
          shape_refSurf(surf[name[1]].sh, name[1]);
          if(surf[name[1]].pgn!=NULL)
	  {
            free(surf[name[1]].pgn); surf[name[1]].pgn=NULL; surf[name[1]].npgn=0;

            /* create the interior */
            repSurf(name[1],1);
	  }
          printf (" interior changed to shape: %s\n", shape[surf[name[1]].sh].name );
	}
        else
        {
	  /* define this surf-name as the default name */
          printf(" Active surf:%s\n", surf[name[1]].name);
          surfNr=name[1];
        }
      }  
    }
    else if (compare(pickfunc, "qbod", 4) == 4) pickstack(name[1]);
    else if (compare(pickfunc, "qshp", 4) == 4)
    {
      /* if a a default shape or nurbs is known apply to surf */
      if((nurbsNr>-1)||(shapeNr>-1))
      {
        if(nurbsNr>-1) { surf[name[1]].sh=getShapeNr(nurbs[nurbsNr].name); shape_refSurf(surf[name[1]].sh, name[1]); }
        if(shapeNr>-1) { surf[name[1]].sh=shapeNr; shape_refSurf(surf[name[1]].sh, name[1]); }
        

	/* if it was prepared for illumination, change the interior def of the selected surf */
        if(surf[name[1]].pgn!=NULL)
	{
          free(surf[name[1]].pgn); surf[name[1]].pgn=NULL; surf[name[1]].npgn=0;

          /* create the interior */
          repSurf(name[1],1);
	}
        printf (" interior changed to Shape: %s\n", shape[surf[name[1]].sh].name );
      }
      else
      {
        printf (" ERROR: No active shape\n");
      }
    }
    else if (compare(pickfunc, "qflp", 4) == 4) flip( type, name[1]);
    else if (compare(pickfunc, "qpnt", 4) == 4) 
    {
      /* mark a surface as a target for the point-projection */
      if(surf[name[1]].sh>-1)
      {
        surfNr=name[1];
      }
      else
      {
        errMsg("WARNING: surf:%s is not related to a nurbs, surf can not be used for projection\n", surf[name[1]].name);
      }
    }
  }
  else if (name[0] == 'b')
  {
    printf ("%s %1c", body[name[1]].name, body[name[1]].ori );
    for (j=0; j<body[name[1]].ns; j++)
      printf (" %1c %s", body[name[1]].o[j], surf[body[name[1]].s[j]].name );
    if( body[name[1]].eparm != (char *)NULL )
      printf (" etyp:%d attr:%d mshp:%s lock:%d\n", body[name[1]].etyp, body[name[1]].eattr, body[name[1]].eparm, body[name[1]].elock);
    else
      printf (" etyp:%d attr:%d lock:%d\n", body[name[1]].etyp, body[name[1]].eattr, body[name[1]].elock);
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].body,set[i].anz_b,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qbod", 4) == 4) bodyNr=name[1];
    else if (compare(pickfunc, "qflp", 4) == 4) flip( type, name[1]);
  }
  else if (name[0] == 'L')
  {
    printf ("%s \n", nurbl[name[1]].name);
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].nurl,set[i].anz_nurl,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
  }
  else if (name[0] == 'S')
  {
    printf ("%s \n", nurbs[name[1]].name);
    if (compare(pickfunc, "qenq", 4) == 4)
    {
      sprintf(printbuf, "    in set=");
      index=2;
      for(i=1; i<anz->sets; i++)
        if(!set[i].type)
        {
          set[i].index=index++;
          if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].nurs,set[i].anz_nurs,name[1]) >-1))
            sprintf(&printbuf[strlen(printbuf)],"%s(%d),",set[i].name,set[i].index);  
        }
      sprintf(&printbuf[strlen(printbuf)],"\n"); if(strlen(printbuf)>12) printf("%s", printbuf);
    }
    else if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qpnt", 4) == 4) 
    {
      /* mark a surface as a target for the point-projection */
      nurbsNr=name[1];
    }
    else if (compare(pickfunc, "qshp", 4) == 4)
    {
	printf(" Active Nurbs:%s\n", nurbs[name[1]].name);
        nurbsNr=name[1];
        shapeNr=-1;
    }
  }
  else printf ("\n");
  return(1);
}



int hitUndo( GLuint *name, char *type, int x, int y )
{
  /* name[0]== entity type of picked name   */
  /* name[1]== picked name(Nr)  */
  /* type[0] == requested entity type */
  /* type[1] == often gkey in pick() */

  if (name[0] == 'e')
  {
    printf (" element: %d \n", name[1]);
    /* if (compare(pickfunc, "qdel", 4) == 4) */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else if (name[0] == 'f')
  {
    printf (" face: %d \n", name[1]);
    /* if (compare(pickfunc, "qdel", 4) == 4) */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else if (name[0] == 'h')
  {
    /* if (compare(pickfunc, "qdel", 4) == 4) */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else if (name[0] == 'n')
  {
    /* if (compare(pickfunc, "qdel", 4) == 4) */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else if (name[0] == 'p')
  {
    /*
    if (compare(pickfunc, "qdel", 4) == 4)
    */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else if (name[0] == 'l')
  {
    /*
    if (compare(pickfunc, "qdiv", 4) == 4) { line[name[1]].div=pickbuf; repLine(name[1]); }
    if (compare(pickfunc, "qspl", 4) == 4) qsplitLine( name[1], x, y, type[1]);
    if (compare(pickfunc, "qdel", 4) == 4) setr( setNrbuf, type, name[1]);
    */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else if (name[0] == 's')
  {
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
    /*
    if (compare(pickfunc, "qdel", 4) == 4) 
    */
  }
  else if (name[0] == 'b')
  {
    /* if (compare(pickfunc, "qdel", 4) == 4) */
    if (compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, name[1]);
    else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, name[1]);
  }
  else printf ("\n");
  return(1);
}



int pickstack( GLuint name )
{
  printf (" pickstack: name:%d added\n", name );
  pickdata[0]++;
  if (pick_buffer<=pickdata[0])
  {
    pick_buffer=DPICK_BUFFER+pickdata[0]; 
    if((pickdata=(int *)realloc((int *)pickdata,(pick_buffer)*sizeof(int)))==NULL)
    { errMsg("\n\nERROR: realloc failure in pickstack\n");
    return(-1); }
  }    
  pickdata[pickdata[0]]=name;
  return(pickdata[0]);
}



int processHits( GLint hits, GLuint *buffer, char *type, char *mode, int x, int y )
{
  int i;
  GLuint   *ptr, key;
  GLuint   name[2];
  static GLuint ubuffer;
  float value;

  pick_zmin=(unsigned int)(pow(2,32)-1);
  name[0]=-1;
  name[1]=-1;

  key=type[0];
  ptr=buffer;
  /* ptr pro hit: anz_names, pick_zmin, max_z, name_1, name_n, ... */
  if (hits)
  {
    
    // printf (" hits= %d mode %c type %c %c\n", hits,  mode[0], type[0], type[1]);
    
    if ( mode[0] == 'i')
    {
      for (i=0; i<hits; i++) /* suche hit mit pick_zmin */
      {
        
        //printf ("hitrecord: %u %u %u type:%c entity:%u %d key:%c\n", *(ptr), *(ptr+1), *(ptr+2), *(ptr+3), *(ptr+4), *(ptr+4), key);
	
        ptr+=2;
        if(((unsigned int)ptr[0] <= (unsigned int)pick_zmin) &&(ptr[1] == key) )
        {
          pick_zmin=*ptr; name[0]=ptr[1]; name[1]=*(ptr+2);
	  
          //printf ("pick_zmin:%u typ:%c name:%u key:%c\n", pick_zmin, name[0], name[1], (unsigned int)type[0] );
          
        }
        ptr+=1+*(ptr-2);
      }
      if (name[0]!=(GLuint)-1)
      {
        if (type[3]=='u')  hitUndo( name, type, x, y );
        else               hitAction( name, type, x, y);
      }
    }
    if ( mode[0] == 'a')
    {
      for (i=0; i<hits; i++) /* collect all hits */
      {
      
        //printf ("hitrecord: %u %u %u %c %u\n", *(ptr), *(ptr+1), *(ptr+2), *(ptr+3), *(ptr+4) );
      
        ptr+=2;
        if ((unsigned int)ptr[1] == (unsigned int)key)
        {
          pick_zmin=*ptr; name[0]=ptr[1]; name[1]=*(ptr+2);
           
          //printf ("pick_zmin:%u typ:%1c name:%u key:%1c\n", pick_zmin, name[0], name[1], type[0] );
          
          if (type[3]=='u')  hitUndo( name, type, x, y );
          else               hitAction( name, type, x, y );
        }
        if (key=='a')
        {
          if((unsigned int)ptr[0] <= (unsigned int)pick_zmin)
          {
            pick_zmin=*ptr; name[0]=ptr[1]; name[1]=*(ptr+2);
	  
            //printf ("pick_zmin:%u typ:%c name:%u key:%c\n", pick_zmin, name[0], name[1], (unsigned int)type[0] );
          
          }
        }
        ptr+=1+*(ptr-2);
      }
    }
    if ( mode[0] == 'm') /* search max/min values in given range of nodes */
    {
      if (type[3]=='u')
      {
        ntext[anz->t-1].node_nr=0;
        if (compare(pickfunc, "qdel", 4) == 4) setr( setNrbuf, type, ubuffer);
        else if(compare(pickfunc, "qadd", 4) == 4) setr( setNrbuf, type, ubuffer);
        else if (compare(pickfunc, "qrem", 4) == 4) seta( setNrbuf, type, ubuffer);
        return(1);
      }
      if (type[1]=='h') value=-MAX_FLOAT; else value=MAX_FLOAT;
      name[1]=0;
      for (i=0; i<hits; i++) /* collect all hits */
      {
        ptr+=2;
        if ((unsigned int)ptr[1] == (unsigned int)key)
        {
          name[0]=*(ptr+2); /* buffer, will be redefined as 'n' later */
          if (type[1]=='h')
	  {
	    if(sequenceFlag)
	    {
              if(value < lcase[lcase_animList].dat[animList][name[0]])
	      {
                name[1]=name[0];
                value=lcase[lcase_animList].dat[animList][name[1]];
	      }
	    }
            else
	    {
              if(value < lcase[cur_lc].dat[cur_entity][name[0]])
	      {
                name[1]=name[0];
                value=lcase[cur_lc].dat[cur_entity][name[1]];
	      }
	    }
 	  }
          else
	  {
	    if(sequenceFlag)
	    {
              if(value > lcase[lcase_animList].dat[animList][name[0]])
	      {
                name[1]=name[0];
                value=lcase[lcase_animList].dat[animList][name[1]];
	      }
	    }
            else
	    {
              if(value > lcase[cur_lc].dat[cur_entity][name[0]])
	      {
                name[1]=name[0];
                value=lcase[cur_lc].dat[cur_entity][name[1]];
	      }
	    }
	  }
        }
        ptr+=1+*(ptr-2);
      }

      if(name[1])
      {
        name[0]='n';
        ubuffer=name[1];

        printf ("%d v= %e ", name[1], value);
        /* search linked values (iexist=2) */
        if(sequenceFlag)
        {
          for(i=0; i<lcase[lcase_animList].ncomps; i++)
            if((i!=animList)&&(lcase[lcase_animList].iexist[i]==2)) printf ("%s:%5.4e ",lcase[lcase_animList].compName[i],lcase[lcase_animList].dat[i][name[1]]);
        }
        else
        {
          for(i=0; i<lcase[cur_lc].ncomps; i++)
            if((i!=cur_entity)&&(lcase[cur_lc].iexist[i]==2)) printf ("%s:%5.4e ",lcase[cur_lc].compName[i],lcase[cur_lc].dat[i][name[1]]);
        }
  
        printf (" xyz= %lf %lf %lf \n",
        (node[name[1]].nx* scale->w+scale->x),
        (node[name[1]].ny* scale->w+scale->y),
        (node[name[1]].nz* scale->w+scale->z) );
        createText(name[1], -1, 0 );
        if (compare(pickfunc, "qdel", 4) == 4) seta( setNrbuf, type, name[1]);
        else if(compare(pickfunc, "qadd", 4) == 4) seta( setNrbuf, type, name[1]);
        else if (compare(pickfunc, "qrem", 4) == 4) setr( setNrbuf, type, name[1]);
      }
    }
  }
  else errMsg ("\n");
  return(1);
}



void defineDiv( unsigned char gkey, int x, int y )
{
  int i,j;
  static int ii=0;
  static char type[MAX_LINE_LENGTH], buffer[MAX_LINE_LENGTH];

  if (ii)
  {
    printf("%c",gkey);
    ii--;
    type[ii]=gkey;
  }
  else
  {
    /* 3 digits number? */
    if(gkey==' ')
    {
      printf(" define 3 digits Div: \n");
      ii=3;
    }

    /* 2 digits number? */
    else
    {
      ii=2;
    }
    type[ii]='\0';

    /* 2 digits number */
    if (ii==2)
    {
      printf("%c",gkey);
      ii--;
      type[ii]=gkey;
    }
  }

  if (!ii)
  {
    printf("\n");
    j=0;
    for(i=strlen(type)-1; i>=0; i--) buffer[j++]=type[i];
    //buffer[j]='\0';
    buffer[j]=0;
    i=0;
    pickbuf=atof(buffer);
    glutKeyboardFunc ( pick );
    type[0]='l';
    type[1]=' ';
    type[2]='\0';
    if (compare(pickfunc, "qmsh", 4) == 4)
    {
      if(setNrdiv==-1) return;
      if((!set[setNrbuf].anz_l)&&(!set[setNrbuf].anz_s)) { printf(" WARNING: Nothing selected\n"); return; }
      /* delete the mesh */
      sprintf(buffer,"me %s", set[setNrbuf].name);
      pre_del(buffer);
      for(i=0; i<set[setNrbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[setNrbuf].line[i]].name);
        line[set[setNrbuf].line[i]].nn=0;
        line[set[setNrbuf].line[i]].ne=0;
      }
      for(i=0; i<set[setNrbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[setNrbuf].surf[i]].name);
        surf[set[setNrbuf].surf[i]].nn=0;
        surf[set[setNrbuf].surf[i]].ne=0;
      }
      if(!qbiaFlag)
      {
        /* change division */
        strcpy(pickfunc, "qdiv");
        goPicking(xbuf,ybuf,type);
        for(i=0; i<set[setNrdiv].anz_l; i++)
        { line[set[setNrdiv].line[i]].div=line[entitybuf].div; repLine(set[setNrdiv].line[i]); }
      }
      else
      {
        strcpy(pickfunc, "qbia");
        goPicking(xbuf,ybuf,type);
        for(i=0; i<set[setNrdiv].anz_l; i++)
        {
          line[set[setNrdiv].line[i]].bias=line[entitybuf].bias;
          repLine(set[setNrdiv].line[i]);
        }
      }
      qmshOperator=0;
      qbiaFlag=0;
      strcpy(pickfunc, "qmsh");
      /* mesh again */
      if(set[setNrbuf].anz_s)
      {
        for(i=0; i<set[setNrbuf].anz_s; i++)
        {
          printf("mesh %s\n",surf[set[setNrbuf].surf[i]].name);
        }
        completeSet(set[setNrbuf].name, "do");
        pre_mesh( set[setNrbuf].name );
        sprintf(buffer,"n %s",set[setNrbuf].name);
        pre_merge(buffer);
      }
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
    }
    else goPicking(xbuf,ybuf,type);
  }
}

void defineValue( unsigned char gkey, int x, int y )
{
  static int i=0;
  static char type[MAX_LINE_LENGTH];
  type[i]=gkey;
  i++;
  if (gkey==( char )0xd)
  {
    i=0;
    printf("\n");
    pickbuf=atof(type);
    glutKeyboardFunc ( pick );
    type[0]='p';
    type[1]='p';
    type[2]='\0';
    goPicking(xbuf,ybuf,type);
  }
  else
  {
    printf("%c",gkey);
    fflush(stdout);
  }
}


/* fragt Keyboard ab, was im pick-mode zu tun ist  */
/* wenn "q" dann exit  */
void pick( unsigned char gkey, int x, int y )
{
  int  i,j,n,k;
  char flag;
  static int ii={0}, xcur[2], ycur[2], wtogle={0};
  static double yw[2], xw[2], dist;
  static int anz_lpc[EDGES_PER_SURF], clines[EDGES_PER_SURF][LINES_PER_EDGE], edge[EDGES_PER_SURF], newc[EDGES_PER_SURF];
  char ctyp[EDGES_PER_SURF];
  static char type[MAX_LINE_LENGTH];
  int  icor, nr, shbuf, linbuf;
  static int anz_c=0, anz_newc=0;
  static char name[MAX_LINE_LENGTH], ori[2], blend[MAX_LINE_LENGTH], cori[EDGES_PER_SURF][2];
  static char *lori=NULL;
  static char face[SURFS_PER_BODY][MAX_LINE_LENGTH];
  static char *alreadyChecked=NULL;
  GLuint size_selectBuf;
  double dum;

  dum=(double)anz->n+anz->e+anz->f+anz->g+anz->t+anzGeo->p+anzGeo->l+anzGeo->s+anzGeo->b+anzGeo->sh+anzGeo->nurl+anzGeo->nurs+1000;
  dum*=10;
  if(dum > MAX_INTEGER) size_selectBuf=MAX_INTEGER;
  else size_selectBuf=(GLuint)dum;
  do{
    if( ( selectBuf= (GLuint *)realloc((GLuint *)selectBuf, size_selectBuf * sizeof(GLuint))) == NULL )
    {
      printf ("WARNING: in Pick() is size_selectBuf: %d to large and is reduced\n", size_selectBuf);
      size_selectBuf/=2;
    }
    if(size_selectBuf<100)
    {
      errMsg("\n\n ERROR: realloc Failure in pick()\n\n") ;
      return;
    }
  }while(!selectBuf);
  glSelectBuffer (size_selectBuf, selectBuf);


  /* small UNDO capability */
  if ( gkey == 'u') /* undo last picking-action */
  {

    if (compare(pickfunc, "qcut", 4) == 4) uncut(1);

    else if (compare(pickfunc, "qmov", 4) == 4)
    { 
      GLubuf[1]=-2;
      moveSet(GLubuf, x, y);
    }
    else if (compare(pickfunc, "qnod", 4) == 4)
    { 
      moveNode(-2, x, y);
    }
    else if (compare(pickfunc, "qmsh", 4) == 4)
    { 
      if(entitybuf>=0) { delElem( 1, &entitybuf); entitybuf=-1; }
    }
    else if (compare(pickfunc, "qpnt", 4) == 4)
    { 
      if(entitybuf>=0) { delPnt( 1, &entitybuf); entitybuf=-1; }
      else
      {
        GLubuf[1]=-2;
        movePoint(GLubuf, x, y);
      }
    }
    else if (compare(pickfunc, "qlin", 4) == 4)
    {
      if(entitybuf>=0) { delLine( 1, &entitybuf); entitybuf=-1; }
    }
    else if (compare(pickfunc, "qsur", 4) == 4)
    {
      if(entitybuf>=0) { delSurf( 1, &entitybuf); entitybuf=-1; }
      /* reset counters */
      anz_c=0;
      for (i=0; i<EDGES_PER_SURF; i++) anz_lpc[i]=0;
    }
    else if (compare(pickfunc, "qbod", 4) == 4)
    {
      if(entitybuf>=0) { delBody( 1, &entitybuf); entitybuf=-1; }
    }
    else if (compare(pickfunc, "qenq", 4) == 4)
    { 
      if(mode[0]=='m') ntext[anz->t-1].node_nr=0;
    }
    else if ((compare(pickfunc, "qadd", 4) == 4)||(compare(pickfunc, "qrem", 4) == 4))
    {
      type[3]='u';
      processHits( hits, selectBuf, type, mode, x, y );
      if (type[0]=='e') updateDispLists();
    }
    return;
  }

  type[0]=gkey;
  type[2]='\0';
  type[3]='\0';

  /* QUIT: delete and re-open the temporary set which stores highlited entities, or quit */
  if (( type[0] == 'q')||( type[0] == 'Q')) /* Exit Pick() */
  {
    /* if a q-func was started inside another q-func, go back to 1st q-func */
    if(hidefunc[0]!=0)
    {
      printf(" change back to %s\n", hidefunc);
      strcpy(pickfunc, hidefunc);
      hidefunc[0]=0;
      mode[0]=hidemode[0];
      mode[1]=hidemode[1];
      return;
    }
    delSet(specialset->highl);
    dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
    glutSetWindow( w1);
    /* glutPassiveMotionFunc ( NULL ); */   /* schaltet mitgehenden Fangrahmen aus */
    free(selectBuf); selectBuf=NULL;
    glutKeyboardFunc ( Keyboard );    /* Keyboardfunkt wiederherstellen (von pick())  */
    glutSetWindow( activWindow );
    pickFlag=0;
    qbiaFlag=0;
    qdisFlag=0;
    intersectFlag=0; 
    anz_c=0;
    entitybuf=-1;
    pntNr=-1;  
    lineNr=-1;   
    surfNr=-1;
    nurbsNr=-1;
    shapeNr=-1;
    bodyNr=-1;  
    qaddTol=-1;
    if (compare(pickfunc, "qsur", 4) == 4) 
    {
      anz_c=0;
      for (i=0; i<EDGES_PER_SURF; i++) anz_lpc[i]=0;
    }
    else if (compare(pickfunc, "qflp", 4) == 4) { free(selem); selem=NULL; }
    else if (compare(pickfunc, "qmsh", 4) == 4)
    {
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      /* switch to original nr of thread */
      anz->threads=anz_threads;
    }

    if(getSetNr(specialset->tmp)>-1) delSet(specialset->tmp);

    pickfunc[0]='\0';
    printf (" done\n");
    return;
  }

  /* create a pset for the selected entity and display it in magenta */
  delSet(specialset->highl);

  if ((pset = (Psets *)realloc( (Psets *)pset, (anzGeo->psets+1)*sizeof(Psets)) ) == NULL )
  {
    printf(" ERROR: realloc failure in plot, pset not installed\n\n");
    return;
  }
  set_highl= pre_seta( specialset->highl, "i", 0 );
  pset[anzGeo->psets].nr= set_highl;
  pset[anzGeo->psets].type[1]='\0';
  pset[anzGeo->psets].type[2]='\0';
  pset[anzGeo->psets].width=0;

  if (compare(pickfunc, "qsur", 4) == 4) { if(type[0]=='s') pset[anzGeo->psets].type[0]='s'; else pset[anzGeo->psets].type[0]= 'l'; }
  else if (compare(pickfunc, "qlin", 4) == 4) pset[anzGeo->psets].type[0]= 'p';
  else pset[anzGeo->psets].type[0]= tolower(gkey);

  if (pset[anzGeo->psets].type[0]=='p'||pset[anzGeo->psets].type[0]=='l'||pset[anzGeo->psets].type[0]=='s'||pset[anzGeo->psets].type[0]=='b') pset[anzGeo->psets].type[1]= 'a';

  pset[anzGeo->psets].col=7;
  anzGeo->psets++;
  glutPostRedisplay();

  /* Action: Look what should be done */
  if ( type[0] == 'r') /*  new pick-area, print x,y */
  {
    xcur[ii]=x; ycur[ii]=y;
    if (ii)
    {
      ii=0;
      dx_cur= xcur[1]-xcur[0];
      dx_cur= sqrt( dx_cur*dx_cur );
      dy_cur= ycur[1]-ycur[0];
      dy_cur= sqrt( dy_cur*dy_cur );
      printf (" cursorArea set to dx_cur=%lf dy_cur=%lf\n", dx_cur, dy_cur);
    }
    else
      ii=1;
    return;
  }
  else if ( type[0] == 'w') /* show window koordinates */
  {
    xw[wtogle]=(double)x/width_w1*scale->w*2.*ds * aspectRatio_w1;
    yw[wtogle]=(double)y/width_w1*scale->w*2.*ds * aspectRatio_w1;
    if (wtogle)
    {
      dist=sqrt( (xw[1]-xw[0])*(xw[1]-xw[0]) + (yw[1]-yw[0])*(yw[1]-yw[0]) );
      printf (" P2 xw:%lf yw:%lf \n", xw[wtogle], yw[wtogle] );
      printf (" distance: %e \n", dist);
      wtogle=0;
    }
    else
    {
      printf (" P1 xw:%lf yw:%lf \n", xw[wtogle], yw[wtogle]);
      wtogle=1;
    }
    return;
  }


  else if ((compare(pickfunc, "qenq", 4) == 4)||(compare(pickfunc, "qadd", 4) == 4)||
           (compare(pickfunc, "qdis", 4) == 4)||
           (compare(pickfunc, "qrem", 4) == 4)||(compare(pickfunc, "qdel", 4) == 4))
  {
    if ( type[0] == 't') /*  change to "qtxt" */
    {
      printf ("mode:%c\n", 't');
      strcpy(hidefunc, pickfunc);
      strcpy(pickfunc, "qtxt");
      printf(" change to %s\n", pickfunc );
      hidemode[0]=mode[0];
      hidemode[1]=mode[1];
      mode[0]='i';
      moveText(anz->t-1, x, y);
      return;
    }
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'c')
    {
      /*  center-point for qdis */
      if(compare(pickfunc, "qdis", 4) == 4)
      {
        type[1]=type[0];
        if(qdisFlag!=2) type[0]='p';
        else type[0]='n';
      }
      /*  lcmb in qadd, qdel and qrem */
      if((compare(pickfunc, "qenq", 4) == 4)||(compare(pickfunc, "qadd", 4) == 4)||(compare(pickfunc, "qrem", 4) == 4)||(compare(pickfunc, "qdel", 4) == 4))
      {
        type[1]=type[0];
        type[0]='l';
      }      
    }
    if ( type[0] == 'g')
    {
      /*  pick points to define a plane for qdis */
      if(compare(pickfunc, "qdis", 4) == 4)
      {
        type[1]='h';
        type[0]='p';
      }
    }
    if ( type[0] == 'm')
    {
      /*  pick nodes to define a plane for qdis */
      if(compare(pickfunc, "qdis", 4) == 4)
      {
        type[1]='h';
        type[0]='n';
      }
      else
      {
        /*  search max/min-value at nodes in qenq */
        mode[0]='m';
        printf ("mode:%s\n", mode);
        minus("nt all    \n");
        plus("nt all k   \n");
        return;
      }
    }

    if ( mode[0] == 'm')
    {
      type[1]=type[0];
      type[0]='n';
    }      

    if (compare(pickfunc, "qdel", 4) == 4) {
      if((type[0]=='n')||(type[0]=='e')||(type[0]=='p')||(type[0]=='h')||(type[0]=='l')||(type[0]=='s')||(type[0]=='b')
        ||(type[0]=='L')||(type[0]=='S'))
      {
        if(getSetNr(specialset->zap)>-1) delSet( specialset->zap);
        setNrbuf=pre_seta(specialset->zap, "i", 0 );
        if(setNrbuf<0) { errMsg(" ERROR: could not create set for qdel\n"); return; } 
        goPicking(x,y,type);
        zap(set[setNrbuf].name);
      }
    }
    else goPicking(x,y,type);

    if ((compare(pickfunc, "qadd", 4) == 4)&&(qaddTol>-1))
    {
      /* complete the set by all connected faces which do not violate the tolerance */
      completeFacesByTolerance(set_highl, setNrbuf, qaddTol);
    }

    if ((compare(pickfunc, "qadd", 4) == 4)||(compare(pickfunc, "qrem", 4) == 4)||(compare(pickfunc, "qdel", 4) == 4))
      if (type[0]=='e') updateDispLists();
    type[1]=type[0]=0;
  }

  else if ((compare(pickfunc, "qdiv", 4) == 4)||(compare(pickfunc, "qbia", 4) == 4))
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    else if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }

    else if( type[0] == 'c')
    {
      if(compare(pickfunc, "qbia", 4) == 4) /* change the bias-direction */
      {
        type[3]='c';
        type[0]='l';
      }
    }
    else if( type[0] == ' ')
    {
        xbuf=x; ybuf=y;
        printf(" define 2 digits Div: \n");
        glutKeyboardFunc(defineDiv);
        return;
    }
    else
    {
      /* look if we have a valid integer > 0, then pick a line */
      type[1] = ' ';
      icor=atoi(type);
      if (( icor <1 )||( icor >9 ))
      {
        printf(" Key not known\n");
        return;
      }
      type[0]='l';
      pickbuf=icor;
    }
    goPicking(x,y,type);
  }

  else if (compare(pickfunc, "qcnt", 4) == 4)
  {
    if((type[0]=='n')||(type[0]=='p'))
      goPicking(x,y,type);
  }
  else if (compare(pickfunc, "qnor", 4) == 4)
  {
    if(type[0]=='p')
    {
      if(qnorCounter>1)
      {
        xbuf=x; ybuf=y;
        printf(" Displacement: \n");
        glutKeyboardFunc(defineValue);
        return;
      }
      goPicking(x,y,type);
    }
  }

  else if (compare(pickfunc, "qsur", 4) == 4)
  {
    /* create a gsur */
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }

    if(type[0]=='b')
    {
      type[0]='s'; type[1]='b';
      goPicking(x,y,type);
    }
    else if(type[0]=='c') { printf (" default shape, nurbs and surface reseted\n"); nurbsNr=shapeNr=surfNr=-1; }
    else if(type[0]=='l') { type[0]='l'; type[1]='l'; goPicking(x,y,type); }
    else if(type[0]=='s') { type[0]='s'; type[1]='s'; goPicking(x,y,type); }
    else if(type[0]=='h') { type[0]='s'; type[1]='h'; goPicking(x,y,type); }
    else if(type[0]=='S') { type[0]='s'; type[1]='S'; goPicking(x,y,type); }
    else if(type[0]=='g')
    {
      strcpy(ori, "+");

      /* if several lines where picked at once, generate a surf for unstructured meshing */
      if(pickdata[0])
      {
        printf ("generate surf ");

        if(pickdata[0]>=EDGES_PER_SURF)
	{
          errMsg(" ERROR: more(%d) than %d edges defined, start again\n", pickdata[0], EDGES_PER_SURF);
          /* reset counters */
          anz_c=0;
          pickdata[0]=0;
          return;
	}

        /* to avoid multiple lines store them in a temporary set */
        linbuf=pre_seta("-qsurl","i",0);
        for (i=0; i<pickdata[0]; i++) seta (linbuf,"l",pickdata[i+1]);
        anz_c=set[linbuf].anz_l;
        for(i=0; i<anz_c; i++)
        {
          strcpy( cori[i], "+");
          edge[i]=set[linbuf].line[i];
          ctyp[i]='l';
        }
        delSet("-qsurl");

        /* try to generate a new surface in any case */
        getNewName( name, "s" );
        entitybuf= surface_i( name, ori[0], shapeNr, anz_c, &cori[0][0], edge, ctyp );

        if (entitybuf>-1)
        {
          /* if a surface was pre-selected to be updated delete the new one and continue */
          if(surfNr>-1)
          {
            delSurf( 1, &entitybuf );
            /* surfNr will be used as a name for the next generated surface and the original surface will be destroyed. The buffer is then reseted */
            strcpy(name,surf[surfNr].name); ori[0]=surf[surfNr].ori; shbuf=surf[surfNr].sh; surfNr=-1;
            entitybuf= surface_i( name, ori[0], shbuf, anz_c, &cori[0][0], edge, ctyp );
            if (entitybuf>-1) printf (" name: %s elty:%d\n", name, surf[entitybuf].etyp );
            else printf(" failed\n");
          }
          else
	  {
            printf (" name: %s elty: tr6u\n", name );

            /* set the element type to tr6u */
            surf[entitybuf].etyp=8;
            surf[entitybuf].eattr=-1;
            if(surfToShape(entitybuf)<0)
            {
              errMsg("\n WARNING: shape could not be defined. Therefore elty tr6u will not be possible if no shape or nurbs will be assigned manually. Check if all points exist in one common plane.\n\n");
	    }
	  }
	}
        else printf(" failed\n");

        /* reset counters */
        anz_c=0;
        pickdata[0]=0;

        return;
      }

      /* look if we have the right amount of edges */
      if(anz_c<3)
      {
          errMsg(" ERROR: less(%d) than %d edges defined, start again\n", anz_c, 3);
          /* reset counters */
          anz_c=0;
          pickdata[0]=0;
        return;
      }
      if( anz_c>=EDGES_PER_SURF-1)
      {
        errMsg("ERROR: more(%d) than %d edges defined, start again\n", anz_c, EDGES_PER_SURF-1);
        anz_c=0;
        for (i=0; i<EDGES_PER_SURF; i++)
          anz_lpc[i]=0;
        return;
      }
      /* look if all edges are propperly defined */
      for (i=0; i<anz_c; i++)
      {
        if (anz_lpc[i]<1)
	{
          printf("ERROR: define edge Nr:%d\n", i+1);
          return;
        }
      }
      printf ("generate surf ");
      anz_newc=0; 
      for (i=0; i<anz_c; i++)
      {
        strcpy( cori[i], "+");
        if (anz_lpc[i]==1)                 /* we have a single line */
	{
          nr=clines[i][0];
          edge[i]=nr;
          ctyp[i]='l';
        }
        else                               /* we need a lcmb */
	{
          /* do we already have a suitable lcmb? */
          flag=0;
          for (j=0; j<anzGeo->c; j++ )
	  {
            if((lcmb[j].name != (char *)NULL)&&(lcmb[j].nl==anz_lpc[i]))   /* same amount of lines */
	    {
              if (printFlag) printf ("check lcmb:%s \n", lcmb[j].name);
	     
              if((alreadyChecked=(char *)realloc((char *)alreadyChecked,(anz_lpc[i])*sizeof(char)))==NULL)
              { printf(" ERROR: realloc failure in pick()\n\n"); return; }
	      for (k=0; k<anz_lpc[i]; k++) alreadyChecked[k]=0; /* reset */
              for (n=0; n<lcmb[j].nl; n++)
	      {
                flag=0;
                for (k=0; k<anz_lpc[i]; k++)
		{
		  /*
                  printf ("c:%s == l:%s\n", line[lcmb[j].l[n]].name, line[clines[i][k]].name);
		  */
                  if ((lcmb[j].l[n]==clines[i][k])&&(alreadyChecked[k]==0))
                  {
                    flag=1;               /* share that line */
                    alreadyChecked[k]=1;
                  }
		}
                if (!flag) goto not_equal_lcmb;
	      }
	     
              if (printFlag) printf ("equal\n");
	     
              break;
	    }
            else flag=0;
            not_equal_lcmb:;
	  }

          if (!flag)  /* no lcmb was found, so create one */
	  {
            if((lori=(char *)realloc((char *)lori, (anz_lpc[i])*sizeof(char)) ) == NULL )
            { printf(" ERROR: realloc failure in pick()\n\n"); return; }
            for (j=0; j<anz_lpc[i]; j++) lori[j]='+';
            getNewName( name, "c" );
            n=lcmb_i( name, 0, anz_lpc[i], lori, clines[i] );
            if(n <0 )
            {
              errMsg ("ERROR: lcmb in error, start again\n");
              anz_c=0;
              for (j=0; j<EDGES_PER_SURF; j++) anz_lpc[j]=0;
              return;
            }
            edge[i]=n;
            newc[anz_newc++]=n;
            ctyp[i]='c';
	  }
          else
	  {
            if (printFlag) printf ("use existing lcmb:%s nr:%d from %d\n", lcmb[j].name,j,anzGeo->c );
            edge[i]=j;
            ctyp[i]='c';
          }
        }
      }

      /* try to generate a new surface in any case */
      getNewName( name, "s" );
      entitybuf= surface_i( name, ori[0], shapeNr, anz_c, &cori[0][0], edge, ctyp );

      if (entitybuf>-1)
      {
        /* if a surface was pre-selected to be updated delete the new one and continue */
        if(surfNr>-1)
        {
          delSurf( 1, &entitybuf );
          /* surfNr will be used as a name for the next generated surface and the original surface will be destroyed. The buffer is then reseted */
          strcpy(name,surf[surfNr].name); ori[0]=surf[surfNr].ori; shapeNr=surf[surfNr].sh; surfNr=-1;
          entitybuf= surface_i( name, ori[0], shapeNr, anz_c, &cori[0][0], edge, ctyp );
          if (entitybuf>-1) printf (" name: %s elty:%d\n", name, surf[entitybuf].etyp );
          else  printf(" failed\n");
        }
        else
	{
          if (entitybuf>-1) printf (" name: %s elty:%d\n", name, surf[entitybuf].etyp );
          else  printf(" failed\n");
        }
      }
      else
      {
        /* delete new lcmbs */
        if(anz_newc) delLcmb( anz_newc, newc );
      }
      /* reset counters */
      anz_c=0;
      for (i=0; i<EDGES_PER_SURF; i++) anz_lpc[i]=0;
    }

    /* look if we have a valid integer < 0, then pick a line */
    type[1]='\0';
    icor=atoi(type)-1;
    if (icor>=EDGES_PER_SURF-1)
    {
      errMsg ("ERROR: Only %d edges per surface possible (not:%s!)\n", EDGES_PER_SURF-1, type);
      return;
    }
    if ( icor>=0 )
    {
      if( icor>=anz_c) anz_c=icor+1;
      type[0]='l';
      goPicking(x,y,type);
      nr=getLineNr(type);
      if ( nr>-1)          /* picking was successfull */
      {
        printf ("sum_edges:%d actual_edge:%d lines:%d Name:%s index:%d\n",
               anz_c, icor+1, anz_lpc[icor], type, clines[icor][anz_lpc[icor]]);

        /* check if this line was already selected */
        flag=0;
        for (i=0; i<anz_c; i++)
          for (k=0; k<anz_lpc[i]; k++)
              if (nr==clines[i][k]) { flag=1; break; }
        if(!flag)
        {
          clines[icor][anz_lpc[icor]]=nr;
          anz_lpc[icor]++;
        }
      }
      else
        errMsg ("WARNING: No line picked, please try again! \n");
    }
  }

  else if (compare(pickfunc, "qbod", 4) == 4)
  {
    /* create a gbod */

    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
    }
    else if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf (" mode:%s\n", mode);
    }

    else if ( type[0] == 'b') /* get a body to replace */
    {
      goPicking(x,y,type);
    }
    else if(type[0]=='g')
    {
      // get rid of double names
      sprintf( name, "-gbodbuf");
      delSet(name);
      nr=pre_seta(name,"i", 0);
      for(i=0;i<anz_c;i++) pre_seta(name,"s", face[i]);
      anz_c=set[nr].anz_s;
      for(i=0;i<anz_c;i++) strcpy(face[i],surf[set[nr].surf[i]].name);
      delSet(name);

      if(bodyNr==-1) getNewName( name, "b" );
      else { strcpy(name,body[bodyNr].name); bodyNr=-1; }
      printf (" generate body:%s\n", name);

      /* look if we have the right amount of faces */
      if( anz_c==2)
      {
        /* try to create the missing surfs */
        entitybuf=body_( name, &face[0][0] );
      }
      else
      {
        strcpy( blend, "NORM");
        for (i=0; i<anz_c; i++) strcpy( &cori[i][0], "+");
        entitybuf=gbod( name, blend, anz_c, &cori[0][0], &face[0][0] );
      }

      /* reset counters */
      anz_c=0;
      pickdata[0]=0;
      pick_buffer=0;
    }

    else if((type[0]=='s')||(type[0]=='S'))
    {
      goPicking(x,y,type);
      if(pickdata[0]>SURFS_PER_BODY)
      {
        printf(" ERROR: too many surfaces selected:%d. SURFS_PER_BODY in cgx.h must be increased.\n", pickdata[0]);
        return;
      }
      for (i=0; i<pickdata[0]; i++)
      {
        nr=getSurfNr(surf[pickdata[i+1]].name);
        if (nr>-1)
        {
          strcpy( face[i+anz_c], surf[pickdata[i+1]].name );
          printf (" sum:%d face:%d surf:%s \n", i+1+anz_c, pickdata[i+1], face[i+anz_c]);
        }
        else
          printf("pick failed, surf not known\n");
      }
      anz_c+=i;
      pickdata[0]=0;
    }
  }

  else if (compare(pickfunc, "qlin", 4) == 4)
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }

    /*  type[0] l:end a line b:start a line c:add centerpoint etc. */
    if(type[0]=='l')
    {
      type[1]=type[0];
      goPicking(x,y,type);
    }
    else if((type[0]=='g')||(type[0]=='b')||(type[0]=='c')||(type[0]=='m')||(type[0]=='t'))
    {
      type[1]=type[0];
      type[0]='p';
      goPicking(x,y,type);
    }
    else if( type[0] == 'p')
    {
        xbuf=x; ybuf=y;
        printf(" Displacement: \n");
        glutKeyboardFunc(defineValue);
        return;
    }
    else if(type[0]=='s')
    {
      qspl_i=0;           /* forget the splitted lines */
      if(getSetNr(specialset->tmp)>-1) delSet(specialset->tmp);
      setNrbuf=pre_seta(specialset->tmp, "i", 0 );
      if(setNrbuf<0) { errMsg(" ERROR: could not create set for qdel\n"); return; } 
      type[1]=type[0];
      type[0]='l';
      goPicking(x,y,type);
    }
    else if(type[0]=='x')
    {
      type[1]=type[0];
      type[0]='l';
      goPicking(x,y,type);
    }
    else if(type[0]=='e')
    {
      i=mode[0];
      mode[0]='a';
      type[1]=type[0];
      type[0]='l';
      goPicking(x,y,type);
      mode[0]=i;
    }
  }

  else if (compare(pickfunc, "qali", 4) == 4)
  {
    if((type[0]=='n')||(type[0]=='p'))
      goPicking(x,y,type);
  }

  else if (compare(pickfunc, "qcut", 4) == 4)
  {
    type[1]=type[0];
    if(type[0]=='v')
    {
      type[0]='n';
    }
    if((type[0]=='n')||(type[0]=='p'))
      goPicking(x,y,type);
  }

  /*
  used
    a,b,c,d,e,f,h,i,1,m,n,(q),s,v.X.*./, ,'0-9H
  unused
    9,j,k,0,p.r,t,u,w,y,z
  */
  else if (compare(pickfunc, "qmsh", 4) == 4)
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }

    if ( type[0] == 'f') /* generate shell element, after node selection with n */
    {
      createElem(-1);
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }
    if ( type[0] == 'v') /* generate volume element, after node selection with n */
    {
      createElem(-2);
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }

    if ( type[0] == 'b') /* break (split) the pre-selected surface */
    {
      if((keybuf!='s')&&(keybuf!='e'))
      {
        /* switch from qlin to qbia mode */
        if(set[qmshbuf].anz_l) qbiaFlag=1;
        else printf(" WARNING: The last selected entity must be the surf to split\n");
        return;
      }
      /* delete the mesh */
      if((!qmshbuf)||((!set[qmshbuf].anz_l)&&(!set[qmshbuf].anz_s))) { printf(" WARNING: Nothing selected\n"); return; }
      sprintf(buffer,"me %s", set[qmshbuf].name);
      pre_del(buffer);
      for(i=0; i<set[qmshbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[qmshbuf].line[i]].name);
        line[set[qmshbuf].line[i]].nn=0;
        line[set[qmshbuf].line[i]].ne=0;
      }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[qmshbuf].surf[i]].name);
        surf[set[qmshbuf].surf[i]].nn=0;
        surf[set[qmshbuf].surf[i]].ne=0;
      }
      /* split the surface */
      if(keybuf=='e') entitybuf=set[qmshbuf].surf[0];
      printf("break (split) surface %s with line:%s\n", surf[entitybuf].name, line[set[qmshbuf].line[set[qmshbuf].anz_l-1]].name);
      i=splitSurf( 1, &set[qmshbuf].line[set[qmshbuf].anz_l-1], entitybuf, 1 );
      if(i<0)
      {
	printf(" ERROR:%d: could not split the selected surface\n",i);
	printf(" -1: malloc failed						\n");
	printf(" -2: no line or surf selected					\n");
	printf(" -3: has holes							\n");
	printf(" -4: the lines form only one loop, so no 2 surfs can be generated\n");
	printf(" -5: other error\n");
      }
      else
      {
        /* mesh again */
        seta(qmshbuf,"s",i);
        for(i=0; i<set[qmshbuf].anz_s; i++)
        {
          printf("mesh %s\n",surf[set[qmshbuf].surf[i]].name);
        }
        completeSet(set[qmshbuf].name, "do");
        pre_mesh( set[qmshbuf].name );
        sprintf(buffer,"n %s",set[qmshbuf].name);
        pre_merge(buffer);
        printf(" %d edges in the model\n", anz->g);
      }
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }

    /* skip if in qbia mode */
    if (( type[0] == 'c')&&(qbiaFlag));
    /* combine pre-selected surfaces */
    else if ( type[0] == 'c')
    {
      if((!qmshbuf)||((!set[qmshbuf].anz_l)&&(!set[qmshbuf].anz_s))) { printf(" WARNING: Nothing selected\n"); return; }
      sprintf(buffer,"me %s", set[qmshbuf].name);
      pre_del(buffer);
      for(i=0; i<set[qmshbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[qmshbuf].line[i]].name);
        line[set[qmshbuf].line[i]].nn=0;
        line[set[qmshbuf].line[i]].ne=0;
      }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[qmshbuf].surf[i]].name);
        surf[set[qmshbuf].surf[i]].nn=0;
        surf[set[qmshbuf].surf[i]].ne=0;
      }
      /* combine the surfaces */
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("combine surfs %s\n", surf[set[qmshbuf].surf[i]].name);
      }
      i=combineSurfs( qmshbuf );
      if(i==-2) printf(" ERROR: could not combine the selected surfaces, they are of type'BLEND' \n");
      else if(i<0) printf(" ERROR: could not combine the selected surfaces\n");
      else
      {
        /* mesh again */
        seta(qmshbuf,"s",i);
        for(i=0; i<set[qmshbuf].anz_s; i++)
        {
          printf("mesh %s\n",surf[set[qmshbuf].surf[i]].name);
        }
        completeSet(set[qmshbuf].name, "do");
        pre_mesh( set[qmshbuf].name );
        sprintf(buffer,"n %s",set[qmshbuf].name);
        pre_merge(buffer);
        printf(" %d edges in the model\n", anz->g);
      }
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }

    if ( type[0] == 'm') /* mesh the pre-selected entities (l,s) */
    {
      if((!qmshbuf)||((!set[qmshbuf].anz_l)&&(!set[qmshbuf].anz_s))) { printf(" WARNING: Nothing selected\n"); return; }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("mesh %s\n",surf[set[qmshbuf].surf[i]].name);
      }
      completeSet(set[qmshbuf].name, "do");
      pre_mesh( set[qmshbuf].name );
      sprintf(buffer,"n %s",set[qmshbuf].name);
      pre_merge(buffer);
        printf(" %d edges in the model\n", anz->g);
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }

    if ( type[0] == 'd') /* delete the mesh of the pre-selected entities */
    {
      if((!qmshbuf)||((!set[qmshbuf].anz_l)&&(!set[qmshbuf].anz_s))) { printf(" WARNING: Nothing selected\n"); return; }
      sprintf(buffer,"me %s", set[qmshbuf].name);
      pre_del(buffer);
      for(i=0; i<set[qmshbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[qmshbuf].line[i]].name);
        line[set[qmshbuf].line[i]].nn=0;
        line[set[qmshbuf].line[i]].ne=0;
      }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[qmshbuf].surf[i]].name);
        surf[set[qmshbuf].surf[i]].nn=0;
        surf[set[qmshbuf].surf[i]].ne=0;
      }
      printf(" selection still active!\n");
      //delSet("-qmsh");
      return;
    }

    if( type[0] == 'x') /* change the mesh attribute of the pre-selected surfaces */
    {
      if((!qmshbuf)||(!set[qmshbuf].anz_s)) { printf(" WARNING: No surfaces selected\n"); return; }
      /* change the element attribute */
      for ( i=0; i<set[qmshbuf].anz_s; i++)
      {
        nr=set[qmshbuf].surf[i];
        if(surf[nr].etyp==0) surf[nr].etyp=8;
        if(surf[nr].eattr==0) surf[nr].eattr=-1; else surf[nr].eattr=0;
      }
      /* delete the mesh */
      sprintf(buffer,"me %s", set[qmshbuf].name);
      pre_del(buffer);
      for(i=0; i<set[qmshbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[qmshbuf].line[i]].name);
        line[set[qmshbuf].line[i]].nn=0;
        line[set[qmshbuf].line[i]].ne=0;
      }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[qmshbuf].surf[i]].name);
        surf[set[qmshbuf].surf[i]].nn=0;
        surf[set[qmshbuf].surf[i]].ne=0;
      }
      /* mesh again */
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("mesh %s\n",surf[set[qmshbuf].surf[i]].name);
      }
      completeSet(set[qmshbuf].name, "do");
      pre_mesh( set[qmshbuf].name );
      sprintf(buffer,"n %s",set[qmshbuf].name);
      pre_merge(buffer);
        printf(" %d edges in the model\n", anz->g);
      // free qmshbuf
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }

    if(( type[0] == 'h')||( type[0] == 't')) /* increase or decrease the mesh-density of the pre-selected surfaces */
    {
      if((!qmshbuf)||(!set[qmshbuf].anz_s)) { printf(" WARNING: No surfaces selected\n"); return; }
      /* change the mesh density attribute */
      if ( type[0] == 'h') dist=0.8; else dist=1.2;
      for ( i=0; i<set[qmshbuf].anz_s; i++)
      {
        nr=set[qmshbuf].surf[i];
        if(surf[nr].eparm!=(char *)NULL) sprintf(blend,"%.3f", atof(surf[nr].eparm) * dist);
        else
	{
          sprintf(blend,"%.3f", dist);
	}
        if((surf[nr].eparm= (char *)realloc((char *)surf[nr].eparm, (strlen(blend)+1)*sizeof(char))) == NULL )
        { printf("ERROR: malloc failed\n\n" ); return; }
        strcpy(surf[nr].eparm, blend);
      }
      /* delete the mesh */
      sprintf(buffer,"me %s", set[qmshbuf].name);
      pre_del(buffer);
      for(i=0; i<set[qmshbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[qmshbuf].line[i]].name);
        line[set[qmshbuf].line[i]].nn=0;
        line[set[qmshbuf].line[i]].ne=0;
      }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[qmshbuf].surf[i]].name);
        surf[set[qmshbuf].surf[i]].nn=0;
        surf[set[qmshbuf].surf[i]].ne=0;
      }
      /* mesh again */
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("mesh %s\n",surf[set[qmshbuf].surf[i]].name);
      }
      completeSet(set[qmshbuf].name, "do");
      pre_mesh( set[qmshbuf].name );
      sprintf(buffer,"n %s",set[qmshbuf].name);
      pre_merge(buffer);
        printf(" %d edges in the model\n", anz->g);
      // free qmshbuf
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }

    /* the following keys '*'|'/' must be used before the 'l' key to indicate that 'qdiv' (TBD) should multiply/divide the divisions */
    /* they will be reseted by ' '|'<nr' */
    if ( type[0] == '*') /* div of later selected lines are multiplied by the typed number */
    {
      qmshOperator='*';
      printf ("mode:%c\n", qmshOperator);
      return;
    }
    else if ( type[0] == '/') /* div of later selected lines are divided by the typed number */
    {
      qmshOperator='/';
      printf ("mode:%c\n", qmshOperator);
      return;
    }

    if (type[0]=='n')
    {
      goPicking(x,y,type);
    }

    if ((type[0]=='l')||(type[0]=='s')||(type[0]=='e'))
    {
      if(getSetNr(specialset->tmp)>-1) delSet( specialset->tmp);
      setNrbuf=pre_seta(specialset->tmp, "i", 0 );
      if(setNrbuf<0) { errMsg(" ERROR: could not create set\n"); return; } 
      goPicking(x,y,type);
      //delSet("-setdiv");
      qmshbuf=pre_seta("-qmsh", "i", 0 );
      setNrdiv=pre_seta("-setdiv", "i", 0 );

      /* if an element was selected search the related surf and forget the element */
      if(type[0]=='e')
      {
        for ( i=0; i<set[setNrbuf].anz_e; i++)
        {
          k=-1;
          for (j=0; j<anzGeo->s; j++) if( surf[j].name != (char *)NULL )
          {
            if(ifind(&surf[j].elem, surf[j].ne, set[setNrbuf].elem[i])>-1)
            {
              k=j;
              break;
            }
          }
          if(k>-1)
	  {
            printf(" +surf:%s\n",surf[j].name);
	    seta( qmshbuf, "s", j );
            free(set[setNrbuf].elem);
            set[setNrbuf].elem=(int *)NULL;
            set[setNrbuf].anz_e = 0;
            break;
	  }
	} 
      }

      /* search the surfaces for all selected lines and store the directly selected lines if div will be changed later */
      /* suche abhaengige lcmbs */
      for ( i=0; i<set[setNrbuf].anz_l; i++)
      {
        printf(" +line:%s\n",line[set[setNrbuf].line[i]].name);
        seta(qmshbuf, "l", set[setNrbuf].line[i] );
        seta(setNrdiv, "l", set[setNrbuf].line[i] );
        for (j=0; j<anzGeo->c; j++)
        {
          if( lcmb[j].name != (char *)NULL )
          for (n=0; n<lcmb[j].nl; n++)
          {
            if( set[setNrbuf].line[i] == lcmb[j].l[n] )
            {
              printf(" +lcmb:%s\n",lcmb[j].name);
              seta( qmshbuf, "c", j );
            }
          }
        }
      }
      for ( i=0; i<set[setNrbuf].anz_c; i++)
      {
        j=set[setNrbuf].lcmb[i];
        for (n=0; n<lcmb[j].nl; n++)
        {
          seta( qmshbuf, "l", lcmb[j].l[n] );
          printf(" +line %s\n",line[lcmb[j].l[n]].name);
        }
        for (j=0; j<anzGeo->s; j++)
        {
          if( surf[j].name != (char *)NULL )
          for (n=0; n<surf[j].nl; n++)
          {
            if(( set[setNrbuf].lcmb[i] == surf[j].l[n] )&&( surf[j].typ[n] == 'c' ))
            {
              seta( qmshbuf, "s", j);
              printf(" +surf %s\n",surf[j].name);
            }
          }
        }
      }
      for ( i=0; i<set[setNrbuf].anz_s; i++)
      {
        printf(" +surf %s\n",surf[set[setNrbuf].surf[i]].name);
        seta(qmshbuf, "s", set[setNrbuf].surf[i] );
      }

      /* now go over qmshbuf */
      for ( i=0; i<set[qmshbuf].anz_c; i++)
      {
        j=set[qmshbuf].lcmb[i];
        for (n=0; n<lcmb[j].nl; n++)
        {
          seta( qmshbuf, "l", lcmb[j].l[n] );
          printf(" +line %s\n",line[lcmb[j].l[n]].name);
        }
        for (j=0; j<anzGeo->s; j++)
        {
          if( surf[j].name != (char *)NULL )
          for (n=0; n<surf[j].nl; n++)
          {
            if(( set[qmshbuf].lcmb[i] == surf[j].l[n] )&&( surf[j].typ[n] == 'c' ))
            {
              seta( qmshbuf, "s", j);
              printf(" +surf %s\n",surf[j].name);
            }
          }
        }
      }
      for ( i=0; i<set[qmshbuf].anz_l; i++)
      {
        for (j=0; j<anzGeo->s; j++)
        {
          if( surf[j].name != (char *)NULL )
          for (n=0; n<surf[j].nl; n++)
          {
            if(( set[qmshbuf].line[i] == surf[j].l[n] )&&( surf[j].typ[n] == 'l' ))
            {
              seta( qmshbuf, "s", j);
              printf(" +surf %s\n",surf[j].name);
            }
          }
        }
      }
      if(getSetNr(specialset->tmp)>-1) delSet( specialset->tmp);
      return;
    }

    /* include the functionallity from qdiv */
    if(((type[0]>=48)&&(type[0]<=57))||(type[0]=='c'))
    {
      if(keybuf!='l') { printf(" WARNING: No line selected\n"); return; }
      if(setNrdiv==-1) return;
      /* delete the mesh */
      if((!qmshbuf)||((!set[qmshbuf].anz_l)&&(!set[qmshbuf].anz_s))) { printf(" WARNING: Nothing selected\n"); return; }
      sprintf(buffer,"me %s", set[qmshbuf].name);
      pre_del(buffer);
      for(i=0; i<set[qmshbuf].anz_l; i++)
      {
        printf("del mesh line %s\n",line[set[qmshbuf].line[i]].name);
        line[set[qmshbuf].line[i]].nn=0;
        line[set[qmshbuf].line[i]].ne=0;
      }
      for(i=0; i<set[qmshbuf].anz_s; i++)
      {
        printf("del mesh surf %s\n",surf[set[qmshbuf].surf[i]].name);
        surf[set[qmshbuf].surf[i]].nn=0;
        surf[set[qmshbuf].surf[i]].ne=0;
      }
      type[1] ='\0';
      icor=atoi(type);
      pickbuf=icor;
      if((qbiaFlag)&&(type[0]=='c'))
      {
        /* qbia change direction */
        for(i=0; i<set[setNrdiv].anz_l; i++)
        {
          line[set[setNrdiv].line[i]].bias=1./line[set[setNrdiv].line[i]].bias;
          repLine(set[setNrdiv].line[i]);
        }
      }
      else if(!qbiaFlag)
      {
        /* qdiv change division */
        type[0]='l';
        strcpy(pickfunc, "qdiv");
        goPicking(x,y,type);
        for(i=0; i<set[setNrdiv].anz_l; i++)
        {
          line[set[setNrdiv].line[i]].div=line[entitybuf].div;
          repLine(set[setNrdiv].line[i]);
        }
      }
      else
      {
        /* qbia */
        type[0]='l';
        strcpy(pickfunc, "qbia");
        goPicking(x,y,type);
        for(i=0; i<set[setNrdiv].anz_l; i++)
        {
          line[set[setNrdiv].line[i]].bias=line[entitybuf].bias;
          repLine(set[setNrdiv].line[i]);
        }
      }
      qmshOperator=0;
      qbiaFlag=0;
      strcpy(pickfunc, "qmsh");
      /* mesh again */
      if(set[qmshbuf].anz_s)
      {
        for(i=0; i<set[qmshbuf].anz_s; i++)
        {
          printf("mesh %s\n",surf[set[qmshbuf].surf[i]].name);
        }
        completeSet(set[qmshbuf].name, "do");
        pre_mesh( set[qmshbuf].name );
        sprintf(buffer,"n %s",set[qmshbuf].name);
        pre_merge(buffer);
        printf(" %d edges in the model\n", anz->g);
      }
      delSet("-qmsh");
      delSet("-setdiv");
      setNrdiv=-1;
      return;
    }
    else if(( type[0] == ' ')||(type[0]=='c'))
    {
      if(keybuf!='l') { printf(" WARNING: No line selected\n"); return; }
      /* change division, del me and remesh are done in defineDiv()  */
      xbuf=x; ybuf=y;
      if(!qbiaFlag) printf(" define 2 digits Div: \n");
      else  printf(" define 2 digits Bias: \n");
      if(setNrdiv>-1) { setNrbuf=qmshbuf; glutKeyboardFunc(defineDiv); }
      return;
    }
  }

  else if (compare(pickfunc, "qmov", 4) == 4)
  {
    if ( type[0] == 'm') /*  move item */
    {
      GLubuf[1]=-1;
      moveSet(GLubuf, x, y);
      entitybuf=-1;      /* if not, undo would delete the last created entity */
    }
    else if((type[0]=='n')||(type[0]=='p'))
    {
      goPicking(x,y,type);
    }
  }

  else if (compare(pickfunc, "qnod", 4) == 4)
  {
    if ( type[0] == 'm') /*  move item */
    {
       moveNode(-1, x, y);
       entitybuf=-1;      /* if not, undo would delete the last created entity */
    }
    else if (type[0]=='p')
    {
      type[1]=type[0];
      type[0]='n';
      goPicking(x,y,type);
    }
  }

  else if (compare(pickfunc, "qtxt", 4) == 4)
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }

    if ( type[0] == 'm') /*  move item */
    {
       moveText(-1, x, y);
       entitybuf=-1;      /* if not, undo would delete the last created entity */
    }
    else if (type[0]=='p')
    {
      moveText(-2, x, y);
      type[0]='t';
      type[1]='p';
      goPicking(x,y,type);
    }
    else if (type[0]=='f')
    {
      type[0]='t';
      type[1]='f';
      goPicking(x,y,type);
    }
    else if (type[0]=='n')
    {
      type[0]='t';
      type[1]='n';
      goPicking(x,y,type);
    }
    else if (type[0]=='v')
    {
      type[0]='t';
      type[1]='v';
      goPicking(x,y,type);
    }
    else if (type[0]=='s')
    {
      type[0]='t';
      type[1]='s';
      goPicking(x,y,type);
    }
    else if (type[0]=='t')
    {
      type[0]='t';
      type[1]='t';
      goPicking(x,y,type);
    }
    else if (type[0]=='b')
    {
      type[0]='t';
      type[1]='b';
      goPicking(x,y,type);
    }
    else if (type[0]=='g')
    {
      type[0]='n';
      type[1]=0;
      goPicking(x,y,type);
    }
    else if (type[0]=='d')
    {
      type[0]='t';
      type[1]='d';
      goPicking(x,y,type);
    }
  }

  else if (compare(pickfunc, "qpnt", 4) == 4)
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'm') /*  move item */
    {
      surfNr=-1;
      GLubuf[1]=-1;
      movePoint(GLubuf, x, y);
      entitybuf=-1;      /* if not, undo would delete the last created entity */
    }
    else if ( type[0] == 'g') /* create a new point */
    {
      surfNr=-1;
      entitybuf=createPoint( x, y );
    }
    else if (type[0]=='n')  /* choose a node as a trgt for point */
    {
      goPicking(x,y,type);
    }
    else if (type[0]=='p')  /* choose a point */
    {
      goPicking(x,y,type);
      if(surfNr>-1) /* move points to a surface */
      {
        GLubuf[1]=-1;
        movePoint(GLubuf, x, y);
      }
      if(nurbsNr>-1) /* move points to a nurbs */
      {
        GLubuf[1]=-1;
        movePoint(GLubuf, x, y);
      }
    }
    else if (type[0]=='s')  /* choose a surface for projection */
    {
      if(getSetNr(specialset->tmp)>-1) delSet(specialset->tmp);
      setNrbuf=pre_seta(specialset->tmp, "i", 0 );
      if(setNrbuf<0) { errMsg(" ERROR: could not create set for qpnt\n"); return; } 
      goPicking(x,y,type);
    }
    else if (type[0]=='S')  /* choose a Nurbs for projection */
    {
      if(getSetNr(specialset->tmp)>-1) delSet(specialset->tmp);
      setNrbuf=pre_seta(specialset->tmp, "i", 0 );
      if(setNrbuf<0) { errMsg(" ERROR: could not create set for qpnt\n"); return; } 
      goPicking(x,y,type);
    }
  }

  else if (compare(pickfunc, "qshp", 4) == 4)
  {
    if ((type[0]=='p')||(type[0]=='n'))
    {
      goPicking(x,y,type);
    }
    else if (type[0]=='c') { nurbsNr=shapeNr=surfNr=-1; }
    else if ( type[0] == 'g') /* create a new shape */
    {
      if(shp_pindx==3)
      {
        shapebuf.type=0;
        if(!strlen(shapebuf.name)) getNewName(shapebuf.name,"sh");
        printf(" shape:%s %s %s %s generated\n",shapebuf.name, point[shapebuf.p[0]].name, point[shapebuf.p[1]].name, point[shapebuf.p[2]].name);
        shapeNr=shape_i(shapebuf.name, 0, shapebuf.p[0], shapebuf.p[1], shapebuf.p[2], 0, 0,0,0);
        shapebuf.name[0]=0;
      }
      else printf(" Not the right number of points for a plane selected:%d (must be 3). Other types can only be created with the command line using shpe.\n",shp_pindx);
      shp_pindx=0;
    }
    else if ( type[0] == 's') /* select a surf and attach a pre-selected shape or nurbs */
    {
      goPicking(x,y,type);
    }
    else if (type[0]=='S') /* pre-select a nurbs */
    {
      type[0]='S';
      type[1]='S';
      goPicking(x,y,type);
    }
    else if(type[0]=='h') /* pre-select a shape */
    {
      type[0]='h';
      type[1]='h';
      goPicking(x,y,type);
    }
    else errMsg(" command not known\n");
  }

  else if (compare(pickfunc, "qint", 4) == 4)
  {
    if (type[0]=='l')
    {
      goPicking(x,y,type);
    }
    else errMsg(" command not known\n");
  }

  else if (compare(pickfunc, "qfil", 4) == 4)
  {
    if (type[0]=='l')
    {
      goPicking(x,y,type);
    }
    else errMsg(" command not known\n");
  }

  else if (compare(pickfunc, "qspl", 4) == 4)
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }
    if( type[0] == 's') /*  split lines */
    {
      qspl_i=0;           /* forget the splitted lines */
      if(getSetNr(specialset->tmp)>-1) delSet(specialset->tmp);
      setNrbuf=pre_seta(specialset->tmp, "i", 0 );
      if(setNrbuf<0) { errMsg(" ERROR: could not create set for qdel\n"); return; } 
      type[1]=type[0];
      type[0]='l';
      goPicking(x,y,type);
    }
    else errMsg(" command not known\n");
  }

  else if (compare(pickfunc, "qseq", 4) == 4)
  {
    if ( type[0] == 'a') /*  pick all in range mode */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'l') /* get line */
    {
      type[1]=type[0];
      goPicking(x,y,type);
    }
    else if ( type[0] == 'g') /* generate spline */
    {
      //printf (" mode:%c type:%c ,%d, pickbuf:%f\n",mode[0],type[0],type[1], pickbuf);
      if(set[setNrbuf].anz_l==0) printf (" ERROR: No line selected\n");
      else if((type[1]==' ')||(type[1]==0))
      {
        for (j=0; j<set[setNrbuf].anz_l; j++ ) convertLine( set[setNrbuf].line[j], (int)pickbuf);
      }
      else                               /* we need an lcmb */
      {
        /* do we already have a suitable lcmb? */
        flag=0;
        for (j=0; j<anzGeo->c; j++ )
        {
          if((lcmb[j].name != (char *)NULL)&&(lcmb[j].nl==set[setNrbuf].anz_l))   /* same amount of lines */
	  {
            if (printFlag) printf ("check lcmb:%s \n", lcmb[j].name);
	     
            if((alreadyChecked=(char *)realloc((char *)alreadyChecked,(set[setNrbuf].anz_l)*sizeof(char)))==NULL)
            { printf(" ERROR: realloc failure in pick()\n\n"); return; }
	      for (k=0; k<set[setNrbuf].anz_l; k++) alreadyChecked[k]=0; /* reset */
            for (n=0; n<lcmb[j].nl; n++)
	    {
              flag=0;
              for (k=0; k<set[setNrbuf].anz_l; k++)
	      {
	        /*
                printf ("c:%s == l:%s\n", line[lcmb[j].l[n]].name, line[set[setNrbuf].line[k]].name);
                */
                if ((lcmb[j].l[n]==set[setNrbuf].line[k])&&(alreadyChecked[k]==0))
                {
                  flag=1;               /* share that line */
                  alreadyChecked[k]=1;
                }
	      }
              if (!flag) goto not_equal_lcmb2;
	    }
	     
            if (printFlag) printf ("equal\n");
	     
            break;
	  }
          else flag=0;
          not_equal_lcmb2:;
	}

        if (!flag)  /* no lcmb was found, so create one */
        {
          if((lori=(char *)realloc((char *)lori, (set[setNrbuf].anz_l)*sizeof(char)) ) == NULL )
          { printf(" ERROR: realloc failure in pick()\n\n"); return; }
          for (j=0; j<set[setNrbuf].anz_l; j++) lori[j]='+';
          getNewName( name, "c" );
          n=lcmb_i( name, 0, set[setNrbuf].anz_l, lori, set[setNrbuf].line );
          if(n <0 )
          {
            errMsg ("ERROR: lcmb in error, start again\n");
            anz_c=0;
            for (j=0; j<EDGES_PER_SURF; j++) anz_lpc[j]=0;
            return;
          }
          nr=n;
        }
        else
        {
          if (printFlag) printf ("use existing lcmb:%s nr:%d from %d\n", lcmb[j].name,j,anzGeo->c );
          nr=j;
        }
        convertLCMB( nr );
      }
      type[1] = ' ';
      set[setNrbuf].anz_l=0;
      return;
    }
    else if( type[0] == ' ')
    {
        xbuf=x; ybuf=y;
        printf(" define 2 digits Div: \n");
        glutKeyboardFunc(defineDiv);
        return;
    }
    else
    {
      /* look if we have a valid integer > 0, then pick a line */
      type[1] = ' ';
      icor=atoi(type);
      if (( icor <1 )||( icor >9 ))
      {
        printf(" Key not known\n");
        return;
      }
      type[0]='l';
      pickbuf=icor;
      goPicking(x,y,type);
    }
  }  
  else if (compare(pickfunc, "qflp", 4) == 4)
  {
    if ( type[0] == 'a') /* orient all connected surfs */
    {
      mode[0]='a';
      printf ("mode:%s\n", mode);
      return;
    }
    if ( type[0] == 'i') /* individual pick mode */
    {
      mode[0]='i';
      type[1]='i';
      dx_cur=PICK_LENGTH; dy_cur=PICK_LENGTH;
      printf ("mode:%s\n", mode);
      return;
    }
    if (( type[0] == 'b')||( type[0] == 's')||( type[0] == 'e')) /* flip surf and elems  or only elem */
    {
      if(mode[0]=='a')
      {
	mode[0]='i'; /* change back to indiv because the surf-connection will identify the others */
        type[1]='a';
      }
      goPicking(x,y,type);
    }
  }
}



/* look whats drawn in the picking-rectangle  */
void goPicking( int x, int y, char *type )
{
  int j;
  GLint  viewport[4];

#if TEST
  printf(" in goPicking\n");
#endif 

  glutSetWindow( w1);
  glGetIntegerv (GL_VIEWPORT, viewport);
  /*   glClearColor ( backgrndcol_rgb[0], backgrndcol_rgb[1], backgrndcol_rgb[2], backgrndcol_rgb[3] ); */
  /*   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  */

  /*    glRenderMode (GL_SELECT);   */
  /*    glInitNames();   */
  /*    glPushName (-1);   */
  glRenderMode (GL_SELECT); 
  glInitNames(); 
  glPushName (-1); 

  glLoadIdentity();
  gluPickMatrix( (GLdouble) x, (GLdouble) (viewport[3]-y), dx_cur,dy_cur, viewport);
  moveModel();

  if (drawMode==4)
  {
    drawSets( PICK );
  }
  else
  {
    if(surfFlag==1) 
    {
      if(type[0]=='e')
      {
        for (j=0; j<anzGeo->psets; j++ )
        {
          if(pset[j].type[0]=='f') drawFaces_plot( set[pset[j].nr].anz_f, set[pset[j].nr].face, node, colNr, face, 2, 'e', 1, PICK );
        }
      }
      else if(type[0]=='f')
      {
        for (j=0; j<anzGeo->psets; j++ )
        {
          if(pset[j].type[0]=='f') drawFaces_plot( set[pset[j].nr].anz_f, set[pset[j].nr].face, node, colNr, face, 2, 0, 1, PICK );
        }
      }
      else if(type[0]=='n')
      {
        for (j=0; j<anzGeo->psets; j++ )
        {
          if(pset[j].type[0]=='f') drawFaceNodes_plot( set[pset[j].nr].anz_f, set[pset[j].nr].face, node, face, 2, 0 );
        }
      }
    }
    else 
    {
      if(type[0]=='e')
      {
        for (j=0; j<anzGeo->psets; j++ )
        {
          if(pset[j].type[0]=='e') drawElements_plot( set[pset[j].nr].anz_e, set[pset[j].nr].elem, node, colNr, e_enqire, 2, 0,1,PICK );
        }
      }
      else if(type[0]=='f')
      {
        for (j=0; j<anzGeo->psets; j++ )
        {
          if(pset[j].type[0]=='f') drawFaces_plot( set[pset[j].nr].anz_f, set[pset[j].nr].face, node, colNr, face, 2, 0,1,PICK );
        }
      }
      else if(type[0]=='n')
      {
        for (j=0; j<anzGeo->psets; j++ )
        {
          if(pset[j].type[0]=='e') drawElemNodes_plot( set[pset[j].nr].anz_e, set[pset[j].nr].elem, node, e_enqire, 2, 0 );
        }
      }
    }
  }
  /*     glutSwapBuffers(); */ 
  /*     glFlush(); */
  glFlush();
  hits = glRenderMode (GL_RENDER);
  if (hits<0)
  {
    errMsg("\nWARNING: Overflow occured, select a smaller region and try again!\n");
  }
  else
  {
    processHits( hits, selectBuf, type, mode, x, y );
  }
}



void qali()
{
  qaliCounter=0;
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qali");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qcut()
{
  qcutCounter=0;
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qcut");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qmsh()
{
  if(getSetNr(specialset->tmp)>-1) delSet( specialset->tmp);
  setNrbuf=pre_seta(specialset->tmp, "i", 0 );
  if(setNrbuf<0) { errMsg(" ERROR: could not create set for qmsh\n"); return; } 
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qmsh");
  qmshbuf=0;
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
  /* switch to one thread */
  anz_threads=anz->threads;
  anz->threads=1;
}
void qmov( char *record)
{
  if(strlen(record)) setNr=getSetNr(record);
  if(setNr<0) { printf(" ERROR: no valid setname given\n"); return; }
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qmov");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qnod()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qnod");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qcnt()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qcnt");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qdis()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qdis");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
    generateSetIndexes();
}
void qenq()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qenq");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
    generateSetIndexes();
}
void qflp()
{
  int i,j;

  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qflp");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
  if ((selem = (int *)malloc( (anz->emax+1)*sizeof(int)) ) == NULL )
  { printf("\n\nERROR: malloc failure in qflp\n"); return; }
  /* generate selem[elem]= index of surface */
  for (i=0; i<=anz->emax; i++) selem[i]=0;
  for (i=0; i<anzGeo->s; i++) if( surf[i].name != (char *)NULL )
  {
    for (j=0; j<surf[i].ne; j++)
    {
      if((surf[i].elem[j]<=anz->emax)&&(surf[i].elem[j]>=anz->emin))
        selem[surf[i].elem[j]]=i;
    }
  }
}
void qbia()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qbia");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qdiv()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qdiv");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qseq(char *record)
{
  char setname[MAX_LINE_LENGTH];
  setname[0]=0;
  sscanf(record,"%s",setname);
  if(strlen(setname))
  {
    setNrbuf=pre_seta( setname, "i", 0);
    set[setNrbuf].type=1;
    //set[setNrbuf].anz_p=0;
    strcpy( pickfunc, "qadd");
  }
  else
  { 
    if(getSetNr(specialset->tmp)>-1) delSet(specialset->tmp);
    setNrbuf=pre_seta(specialset->tmp, "i", 0 );
    if(setNrbuf<0) { errMsg(" ERROR: could not create set for qseq\n"); return; } 
    strcpy( pickfunc, "qseq");
  } 
  pickFlag=1;
  mode[0]='i';
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qshp( char *record)
{
  static char name[MAX_LINE_LENGTH];
  name[0]=0;
  shapebuf.name=(char *)&name[0];
  if(strlen(record)) sscanf(record,"%s",shapebuf.name);
  pickFlag=1;
  mode[0]='i';
  shp_pindx=0;
  strcpy( pickfunc, "qshp");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qadd( char *record)
{
  char setname[MAX_LINE_LENGTH], type[MAX_LINE_LENGTH];
  setname[0]=0;
  type[0]=0;
  sscanf(record,"%s %s",setname,type); 
  if(!strlen(setname)) { errMsg(" ERROR: could not create set for qadd\n"); return; } 
  setNrbuf=pre_seta( setname, "i", 0);
  if(setNrbuf<0) { errMsg(" ERROR: could not create set for qadd\n"); return; } 
  if (type[0]=='s') set[setNrbuf].type=1;
  if (type[0]=='t') qaddTol=atof(&type[1]);
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qadd");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qdel( void )
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qdel");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qrem( char *record)
{
  setNrbuf=pre_seta( record, "i", 0);
  if(setNrbuf<0) { errMsg(" ERROR: could not create set for qrem\n"); return; } 
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qrem");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qpnt( char *record)
{
  if(strlen(record)) pntNr=getPntNr(record); else pntNr=-1;
  pickFlag=1;
  mode[0]='i';
  movedp=pntnr=pnt_undo=-1;
  moveFlag=0;
  strcpy( pickfunc, "qpnt");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qnor()
{
  qnorCounter=0;
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qnor");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qlin( char *record)
{
  if(strlen(record)) lineNr=getLineNr(record);
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qlin");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qsur( char *record)
{
  if(strlen(record)) surfNr=getSurfNr(record);
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qsur");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
  if(pickdata==NULL) /* first time call */
  {
    pick_buffer=DPICK_BUFFER;
    if((pickdata=(int *)realloc((int *)pickdata,(pick_buffer)*sizeof(int)))==NULL)
    { errMsg("\n\nERROR: realloc failure in pickstack\n");
    return; }
    pickdata[0]=0;
  } 
}
void qbod( char *record)
{
  if(strlen(record)) bodyNr=getBodyNr(record);
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qbod");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
  if(pickdata==NULL) /* first time call */
  {
    pick_buffer=DPICK_BUFFER;
    if((pickdata=(int *)realloc((int *)pickdata,(pick_buffer)*sizeof(int)))==NULL)
    { errMsg("\n\nERROR: realloc failure in pickstack\n");
    return; }
    pickdata[0]=0;
  } 
}
void qspl()
{
  pickFlag=1;
  mode[0]='a';
  strcpy( pickfunc, "qspl");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qint()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qint");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qfil( char *record)
{
  if(strlen(record)) filletRadius=atof(record)/scale->w;
  if(filletRadius<=0.)
  {
    printf(" ERROR: No radius was defined:%lf. Start again with a value.\n", filletRadius);
    return;
  }
  pickFlag=1;
  pickbuf=-1.;
  mode[0]='i';
  strcpy( pickfunc, "qfil");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}
void qtxt()
{
  pickFlag=1;
  mode[0]='i';
  strcpy( pickfunc, "qtxt");
  glutSetWindow( w1);
  glutKeyboardFunc ( pick );
  glutSetWindow( activWindow );
}



void moveText(int t, int x, int y)
{
  register int i;
  static int *txtnr=NULL, sum=0, moveFlag=0;
  static int xbuf,ybuf;

  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];

  static GLdouble wx, wy, wz;  /*  returned window x, y, z coords  */
  static int flag;

  //printf("t:%d flag:%d moveFlag:%d sum:%d x,y %d %d xbuf,ybuf:%d %d wxyz:%f %f %f\n", t, flag, moveFlag, sum, x,y, xbuf, ybuf, wx, wy, wz);
  glutSetWindow( w1);
  if (t>-1)
  {
    printf(" text at node %d selected\n", ntext[t].node_nr);
    if(moveFlag)
    {
      moveFlag=0;
      free(txtnr);
      txtnr=NULL;
      sum=0;
    }
    /* mark text as one to be moved later */
    if ( (txtnr = (int *)realloc( txtnr, (sum+1) * sizeof(int))) == NULL )
      printf("\n\n ERROR: realloc failed\n\n") ;
    txtnr[sum++]=t;
    xbuf=x;
    ybuf=y;
  }
  else if ((t==-1)&&(sum>0))
  {
    moveFlag=1;
    /* move node to new location */
    /* remember that the point was moved, a new stroke of "p" will choose a new text(node) to move */
    if(sum==1)
    {
      ntext[txtnr[0]].tx= (double)x;
      ntext[txtnr[0]].ty= (double)y;
    }
    else 
    {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glLoadIdentity();
      moveModel();
      glGetIntegerv (GL_VIEWPORT, viewport);
      glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
      glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

      for(i=0; i<sum; i++)
      {
	//printf("i:%d xy:%f %f xbuf:%d %d x: %d %d\n", txtnr[i], ntext[txtnr[i]].tx, ntext[txtnr[i]].ty, xbuf,ybuf,x,y);
        if(ntext[txtnr[i]].tx==-1) //still glued to node-pos
	{
          flag=gluProject( node[ntext[txtnr[i]].node_nr].nx, node[ntext[txtnr[i]].node_nr].ny, node[ntext[txtnr[i]].node_nr].nz, mvmatrix, projmatrix, viewport,  &wx, &wy, &wz);
          if (flag==GL_FALSE) printf("WARNING: Malfunction, please reselect\n");
          ntext[txtnr[i]].tx= wx+(double)(x-xbuf);
          ntext[txtnr[i]].ty= (viewport[3]-wy)+(double)(y-ybuf);
	  //printf("wxyz:%f %f v:%d %f   %f %f\n",wx, wy, viewport[3], wz, ntext[txtnr[i]].tx, ntext[txtnr[i]].ty); 
	}
        else
	{
          ntext[txtnr[i]].tx+= (double)(x-xbuf);
          ntext[txtnr[i]].ty+= (double)(y-ybuf);
	}
      }
    }
    xbuf=x;
    ybuf=y;
    updateDispLists();
  }
  else if ((t==-2)&&(sum>0))
  {
    /* reset, next call with an entity will start fresh */
    moveFlag=1;
  }
}



void moveSet(GLuint *picbuf, int x, int y)
{
  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];

  static GLdouble wx, wy, wz;  /*  returned window x, y, z coords  */
  static GLdouble nx, ny, nz;  /*  new world x, y, z coords  */
  static double    nbuf[3];     /*  coordinate-buffer for undo */
  static int flag, indx, indxbuf=-1,typ;


  glutSetWindow( w1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  moveModel();
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  indx=picbuf[1];
  typ=picbuf[0];

  if ((indx<-1)&&(indxbuf==-1))
  {
    /* undo */
    sprintf(buffer,"%s tra %f %f %f",set[setNr].name,-nx,-ny,-nz);
    pre_move( buffer );
    redraw();
    nx=ny=nz=0.;
  }
  else if((indx!=-1)&&(indxbuf==-1))
  {
    /* mark node as a reference point and determine screen coordinates */
    indxbuf=indx;
    nx=ny=nz=0.;
    if(typ=='n')
    {
      nbuf[0]=node[indx].nx;
      nbuf[1]=node[indx].ny;
      nbuf[2]=node[indx].nz;
    }
    if(typ=='p')
    {
      nbuf[0]=point[indx].px;
      nbuf[1]=point[indx].py;
      nbuf[2]=point[indx].pz;
    }
    flag=gluProject( nbuf[0],nbuf[1],nbuf[2], mvmatrix, projmatrix,
       viewport,  &wx, &wy, &wz);
    if (flag==GL_FALSE)
      printf("WARNING: Malfunction in moveNode(), please reselect\n");
    /* printf (" node:%d x=%lf y=%lf z=%lf \n", n, node[n].nx,node[n].ny,node[n].nz );
    printf (" Win coords are %d (%lf, %lf, %lf)\n",  flag,  wx   ,  wy   , wz     ); */
  }
  else if ((indx==-1)&&(indxbuf>-1))
  {
    /* move set to new location */

    /* remember that the point was moved, a new stroke of "p" will choose a new point to move */
    moveFlag=1;
    wx=(GLdouble)x; wy=(GLdouble)(viewport[3]-y);
    flag=gluUnProject ( wx, wy, wz, mvmatrix, projmatrix, viewport, &nx, &ny, &nz);
    /* printf ("World coords are %d (%lf, %lf, %lf)\n", flag, nx, ny, nz); */
    if (flag==GL_TRUE)
    {
      nx-=nbuf[0];
      ny-=nbuf[1];
      nz-=nbuf[2];
      nx=nx*scale->w;
      ny=ny*scale->w;
      nz=nz*scale->w;
      sprintf(buffer,"%s tra %f %f %f",set[setNr].name,nx,ny,nz);
      pre_move( buffer );
      redraw();
      indxbuf=-1;
    }
    else
      printf("WARNING: Malfunction in moveNode(), please repeat\n");
  }
  else if ((indx>-1)&&(indxbuf>-1))
  {
    /* move the set to the actual node n or point p */
    if(typ=='n')
    {
      nx=node[indx].nx-nbuf[0];
      ny=node[indx].ny-nbuf[1];
      nz=node[indx].nz-nbuf[2];
    }
    if(typ=='p')
    {
      nx=point[indx].px-nbuf[0];
      ny=point[indx].py-nbuf[1];
      nz=point[indx].pz-nbuf[2];
    }
    nx=nx*scale->w;
    ny=ny*scale->w;
    nz=nz*scale->w;
    sprintf(buffer,"%s tra %f %f %f",set[setNr].name,nx,ny,nz);
    printf("%s\n",buffer);
    pre_move( buffer );
    redraw();
    indxbuf=-1;
  }
}



void moveNode(int n, int x, int y)
{
  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];

  static GLdouble wx, wy, wz;  /*  returned window x, y, z coords  */
  static GLdouble nx, ny, nz;  /*  new world x, y, z coords  */
  static double    nbuf[3];     /*  coordinate-buffer for undo */
  static int flag, nodnr=-1, nod_undo=-1, moveFlag=0;


  glutSetWindow( w1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  moveModel();
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  if (((n>-1)&&(nodnr==-1))||((n>-1)&&(moveFlag==1)))
  {
    /* mark node as one to be moved later and determine screen coordinates */
    nodnr=n;
    moveFlag=0;
    nod_undo=n;
    nbuf[0]=node[n].nx;
    nbuf[1]=node[n].ny;
    nbuf[2]=node[n].nz;
    
    flag=gluProject( node[n].nx, node[n].ny, node[n].nz, mvmatrix, projmatrix,
       viewport,  &wx, &wy, &wz);
    if (flag==GL_FALSE)
      printf("WARNING: Malfunction in moveNode(), please reselect\n");
    /* printf (" node:%d x=%lf y=%lf z=%lf \n", n, node[n].nx,node[n].ny,node[n].nz );
    printf (" Win coords are %d (%lf, %lf, %lf)\n",  flag,  wx   ,  wy   , wz     ); */
  }
  else if ((n==-1)&&(nodnr>-1))
  {
    /* move node to new location */

    /* remember that the point was moved, a new stroke of "p" will choose a new point to move */
    moveFlag=1;
    wx=(GLdouble)x; wy=(GLdouble)(viewport[3]-y);
    flag=gluUnProject ( wx, wy, wz, mvmatrix, projmatrix, viewport, &nx, &ny, &nz);
    /* printf ("World coords are %d (%lf, %lf, %lf)\n", flag, nx, ny, nz); */
    if (flag==GL_TRUE)
    {
      node[nodnr].nx=nx;
      node[nodnr].ny=ny;
      node[nodnr].nz=nz;
      getElemNormalen( e_enqire, node, anz->e );
      getFaceNormalen( face, node, anz );
      updateDispLists();
    }
    else
      printf("WARNING: Malfunction in moveNode(), please reselect\n");
  }
  else if ((n>-1)&&(nodnr>-1))
  {
    /* move the node nodnr to the actual node n */
    node[nodnr].nx=node[n].nx;
    node[nodnr].ny=node[n].ny;
    node[nodnr].nz=node[n].nz;
    getElemNormalen( e_enqire, node, anz->e );
    getFaceNormalen( face, node, anz );
    updateDispLists();
    nodnr=-1;
  }
  else if ((n<-1)&&(nod_undo>-1))
  {
    /* undo */
    node[nod_undo].nx=nbuf[0];
    node[nod_undo].ny=nbuf[1];
    node[nod_undo].nz=nbuf[2];
    getElemNormalen( e_enqire, node, anz->e );
    getFaceNormalen( face, node, anz );
    updateDispLists();
    nodnr=nod_undo=-1;
  }
}


void movePoint(GLuint *picbuf, int x, int y)
{
  int i,j,nr,p,typ;
  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];
  char     name[MAX_LINE_LENGTH];
  
  static GLdouble wx, wy, wz;  /*  returned window x, y, z coords  */
  static GLdouble nx, ny, nz;  /*  new world x, y, z coords  */
  static double    pbuf[3];     /*  coordinate-buffer for undo */
  static int flag;
  
  glutSetWindow( w1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  moveModel();
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  p=picbuf[1];
  if(p<0) typ=0; else typ=picbuf[0];

  /* undo */
  if ((p==-2)&&(pnt_undo>-1))
  {
    pntnr=pnt_undo;
    /* pnt_undo=-1; */
    point[pntnr].px=pbuf[0];
    point[pntnr].py=pbuf[1];
    point[pntnr].pz=pbuf[2];
    movedp=pntnr;
    pntnr=pnt_undo=-1;
  }
  else if (surfNr>-1)
  {
    /* project p to surf */
    if(shape[surf[surfNr].sh].type==4)
    {
      projSetToNurbs(shape[surf[surfNr].sh].p[0],set,setNrbuf,point,&node);
      movedp=-1;
    
      /* correct the shape of the lines */
      for (i=0; i<anzGeo->l; i++)
      {
        repLine(i); 
      }
    }
  }
  else if (nurbsNr>-1)
  {
    /* project p to surf */
    projSetToNurbs(nurbsNr,set,setNrbuf,point,&node);
    movedp=-1;
    
    /* correct the shape of the lines */
    for (i=0; i<anzGeo->l; i++)
    {
      repLine(i); 
    }
  }
  else if ((((p>-1)&&(pntnr==-1))||((p>-1)&&(moveFlag==1)))&&(typ!='n'))
  {
    /* mark this point as one to be moved later and determine screen coordinates */  
    pntnr=pnt_undo=p;
    moveFlag=0;
    movedp=-1;

    pbuf[0]=point[p].px;
    pbuf[1]=point[p].py;
    pbuf[2]=point[p].pz;
    
    /* necessary for plotting the moved point */
    flag=gluProject( point[p].px, point[p].py, point[p].pz, mvmatrix, projmatrix,
       viewport,  &wx, &wy, &wz);
    if (flag==GL_FALSE)
      printf("WARNING: Malfunction in movePoint(), please reselect\n");
    
    /* printf (" pnt:%d x=%lf y=%lf z=%lf \n", n, point[p].px, point[p].py, point[p].pz );
    printf (" Win coords are %d (%lf, %lf, %lf)\n",  flag,  wx   ,  wy   , wz     ); */
  }
  else if ((p==-1)&&(pntnr>-1))
  {
    /* move the point to the window-coordinates */

    /* remember that the point was moved, a new stroke of "p" will choose a new point to move */
    moveFlag=1;
    wx=(GLdouble)x; wy=(GLdouble)(viewport[3]-y);
    flag=gluUnProject ( wx, wy, wz, mvmatrix, projmatrix, viewport, &nx, &ny, &nz);
    /* printf ("new World coords are %d (%lf, %lf, %lf)\n", flag, nx, ny, nz); */
    if (flag==GL_TRUE)
    {
      point[pntnr].px=nx;
      point[pntnr].py=ny;
      point[pntnr].pz=nz;
      movedp=pntnr;
    }
    else
      printf("WARNING: Malfunction in movePoint(), please reselect\n");
  }
  else if (((p>-1)&&(pntnr>-1))&&typ=='p')
  {
    /* move the point pntnr to the actual point p */

    point[pntnr].px=point[p].px;
    point[pntnr].py=point[p].py;
    point[pntnr].pz=point[p].pz;
    movedp=pntnr;

    /* in case a line exists between the two points then delete this line and merge the two points */
    delSet("-movePoint");
    i=pre_seta( "-movePoint", "p", point[pntnr].name );
    seta(i,"p",p);
    completeSet("-movePoint", "up");
    pre_del("l0 -movePoint" );
    delSet("-movePoint");

    pntnr=-1;
  }
  else if (((p>-1)&&(pntnr>-1))&&typ=='n')
  {
    /* move the point pntnr to the actual node p */

    point[pntnr].px=node[p].nx;
    point[pntnr].py=node[p].ny;
    point[pntnr].pz=node[p].nz;
    movedp=pntnr;
    pntnr=-1;
  }
  else if (((p>-1)&&(pntnr==-1))&&typ=='n')
  {
    /* create a new point at the actual node p */

    getNewName( name, "p" );
    if(printFlag) printf (" pnt=%s x=%lf y=%lf z=%lf\n",  name, node[p].nx* scale->w+scale->x, node[p].ny* scale->w+scale->y, node[p].nz* scale->w+scale->z);
    if( pnt( name
	       , node[p].nx
	       , node[p].ny
	       , node[p].nz, 0 ) <0) printf("ERROR: point could not be created\n");
  }

  if(movedp>-1)
  {
    /* correct the shape of the connected lines */
    for (i=0; i<anzGeo->l; i++)
    {
      if(line[i].typ=='a')
      {
        if((line[i].p1==movedp)||(line[i].p2==movedp)||(line[i].trk==movedp)) repLine(i);
      }
      else if(line[i].typ=='s')
      {
        for (j=0; j<set[line[i].trk].anz_p; j++)
        {
          if(set[line[i].trk].pnt[j]==movedp)
          {
            repLine(i);
            break;
          }
        }
      }
      else
      {
        if((line[i].p1==movedp)||(line[i].p2==movedp)) repLine(i);
      } 
    }
    /* correct the shape of the connected NURBS */
    for (nr=0; nr<anzGeo->nurs; nr++)
    {
      if( nurbs[nr].name != (char *)NULL )
      {
        for (i=0; i<nurbs[nr].u_npnt; i++)
          for (j=0; j<nurbs[nr].v_npnt; j++)
            if (nurbs[nr].ctlpnt[i][j] == movedp) repNurs(nr);
      }
    }
  }
}



void flipSurfori(int s)
{
  int i,j,k,n,setNr=-1;
  int nbuf[26];
  double *buf;

  /* switch signum of surf */
  if (surf[s].ori=='-') surf[s].ori='+';
  else                  surf[s].ori='-';
  n=0;
  while((surf[s].npgn-n)>0)
  {
    n++; /* jump over the polygon token (ie.GL_POLYGON_TOKEN) */
    j=surf[s].pgn[n++];
    surf[s].pgn[n]*=-1;
    surf[s].pgn[n+1]*=-1;
    surf[s].pgn[n+2]*=-1;
    n+=3;
    if ((buf = (double *)malloc((j*3)*sizeof(double)) ) == NULL )
    { printf("\n\nERROR: realloc failure in flip\n\n"); return; }
    for(k=0; k<j; k++)
    {
      buf[j*3-k*3-3]=surf[s].pgn[n];
      buf[j*3-k*3-2]=surf[s].pgn[n+1];
      buf[j*3-k*3-1]=surf[s].pgn[n+2];
      n+=3;
    }
    n-=3*j;
    for(k=0; k<j; k++)
    {
      surf[s].pgn[n]  =buf[k*3] ; 
      surf[s].pgn[n+1]=buf[k*3+1];
      surf[s].pgn[n+2]=buf[k*3+2];
      n+=3;
    }
    free(buf);
  }

  /* switch orient of all embedded elements */
  if(surf[s].ne)
  {
    delSet(specialset->tmp);
    setNr=pre_seta(specialset->tmp,"i",0);
  }
  for (i=0; i<surf[s].ne; i++)
  {
    if (e_enqire[surf[s].elem[i]].type == 7)
    {
      nbuf[1]=e_enqire[surf[s].elem[i]].nod[1];
      nbuf[2]=e_enqire[surf[s].elem[i]].nod[2];
      e_enqire[surf[s].elem[i]].nod[2]=nbuf[1];
      e_enqire[surf[s].elem[i]].nod[1]=nbuf[2];
      seta(setNr,"e",surf[s].elem[i]);
    }  
    else if (e_enqire[surf[s].elem[i]].type == 8)
    {
      nbuf[0]=e_enqire[surf[s].elem[i]].nod[0];
      nbuf[1]=e_enqire[surf[s].elem[i]].nod[1];
      nbuf[4]=e_enqire[surf[s].elem[i]].nod[4];
      nbuf[5]=e_enqire[surf[s].elem[i]].nod[5];
      e_enqire[surf[s].elem[i]].nod[0]=nbuf[1];
      e_enqire[surf[s].elem[i]].nod[1]=nbuf[0];
      e_enqire[surf[s].elem[i]].nod[4]=nbuf[5];
      e_enqire[surf[s].elem[i]].nod[5]=nbuf[4];
      seta(setNr,"e",surf[s].elem[i]);
    }  
    else if (e_enqire[surf[s].elem[i]].type == 9)
    {
      nbuf[1]=e_enqire[surf[s].elem[i]].nod[1];
      nbuf[3]=e_enqire[surf[s].elem[i]].nod[3];
      e_enqire[surf[s].elem[i]].nod[3]=nbuf[1];
      e_enqire[surf[s].elem[i]].nod[1]=nbuf[3];
      seta(setNr,"e",surf[s].elem[i]);
    }  
    else if (e_enqire[surf[s].elem[i]].type == 10)
    {
      nbuf[1]=e_enqire[surf[s].elem[i]].nod[1];
      nbuf[3]=e_enqire[surf[s].elem[i]].nod[3];
      nbuf[4]=e_enqire[surf[s].elem[i]].nod[4];
      nbuf[5]=e_enqire[surf[s].elem[i]].nod[5];
      nbuf[6]=e_enqire[surf[s].elem[i]].nod[6];
      nbuf[7]=e_enqire[surf[s].elem[i]].nod[7];
      e_enqire[surf[s].elem[i]].nod[3]=nbuf[1];
      e_enqire[surf[s].elem[i]].nod[1]=nbuf[3];
      e_enqire[surf[s].elem[i]].nod[7]=nbuf[4];
      e_enqire[surf[s].elem[i]].nod[6]=nbuf[5];
      e_enqire[surf[s].elem[i]].nod[5]=nbuf[6];
      e_enqire[surf[s].elem[i]].nod[4]=nbuf[7];
      seta(setNr,"e",surf[s].elem[i]);
    }
    else  
    {
      printf(" ERROR: element type %d not known\n", e_enqire[surf[s].elem[i]].type);
    }
  }

  if(setNr!=-1)
  {
    for (i=0; i<anz->f; i++)
    { 
      if(ifind(&set[setNr].elem, set[setNr].anz_e, face[i].elem_nr)>-1)
      {
        if (face[i].type == 7)
        {
          nbuf[1]=face[i].nod[1];
          nbuf[2]=face[i].nod[2];
          face[i].nod[2]=nbuf[1];
          face[i].nod[1]=nbuf[2];
        }  
        else if (face[i].type == 8)
        {
          nbuf[0]=face[i].nod[0];
          nbuf[1]=face[i].nod[1];
          nbuf[4]=face[i].nod[4];
          nbuf[5]=face[i].nod[5];
          face[i].nod[0]=nbuf[1];
          face[i].nod[1]=nbuf[0];
          face[i].nod[4]=nbuf[5];
          face[i].nod[5]=nbuf[4];
        }  
        else if (face[i].type == 9)
        {
          nbuf[1]=face[i].nod[1];
          nbuf[3]=face[i].nod[3];
          face[i].nod[3]=nbuf[1];
          face[i].nod[1]=nbuf[3];
        }  
        else if (face[i].type == 10)
        {
          nbuf[1]=face[i].nod[1];
          nbuf[3]=face[i].nod[3];
          nbuf[4]=face[i].nod[4];
          nbuf[5]=face[i].nod[5];
          nbuf[6]=face[i].nod[6];
          nbuf[7]=face[i].nod[7];
          face[i].nod[3]=nbuf[1];
          face[i].nod[1]=nbuf[3];
          face[i].nod[7]=nbuf[4];
          face[i].nod[6]=nbuf[5];
          face[i].nod[5]=nbuf[6];
          face[i].nod[4]=nbuf[7];
        }
        else if ((face[i].type == 11)||(face[i].type == 12))
        {
          nbuf[1]=face[i].nod[0];
          nbuf[2]=face[i].nod[1];
          face[i].nod[1]=nbuf[1];
          face[i].nod[0]=nbuf[2];
	  j=face[i].elem_nr;
	  if((face[i].nod[0]==e_enqire[j].nod[0])&&(face[i].nod[1]==e_enqire[j].nod[1])) face[i].nr=2;
	  else if((face[i].nod[0]==e_enqire[j].nod[1])&&(face[i].nod[1]==e_enqire[j].nod[2])) face[i].nr=3;
	  else if((face[i].nod[0]==e_enqire[j].nod[2])&&(face[i].nod[1]==e_enqire[j].nod[3])) face[i].nr=4;
	  else if((face[i].nod[0]==e_enqire[j].nod[3])&&(face[i].nod[1]==e_enqire[j].nod[0])) face[i].nr=5;
        }
      }
    }
    delSet(specialset->tmp);
  }
}



int createElem(int n)
{
  static int nodes[26];
  static int nr;

  if(n>-1)
  {
    nodes[nr++]=n;
  }
  else if(n==-1)
  {
    printf(" create shell element:%d\n",anz->emax+1 );
    if     (nr==3) elem_define(anz,&e_enqire, anz->enext++, 7, nodes, 1, 0 );
    else if(nr==4) elem_define(anz,&e_enqire, anz->enext++, 9, nodes, 1, 0 );
    else if(nr==6) elem_define(anz,&e_enqire, anz->enext++, 8, nodes, 1, 0 );
    else if(nr==8) elem_define(anz,&e_enqire, anz->enext++, 10, nodes, 1, 0 );
    else
    {
      printf(" Warning, wrong number of nodes:%d, no element created, start again\n",nr);
      nr=0;
      return(-1);
    }
    nr=0;

    makeSurfaces();
    getElemNormalen( e_enqire, node, anz->e );
    realloc_colNr();
    return(anz->emax);
  }
  else if(n==-2)
  {
    printf(" create volume element:%d\n",anz->emax+1);
    if     (nr==8)  elem_define(anz,&e_enqire, anz->enext++, 1, nodes, 1, 0 );
    else if(nr==20) elem_define(anz,&e_enqire, anz->enext++, 4, nodes, 1, 0 );
    else
    {
      printf(" Warning, wrong number of nodes:%d, no element created, start again\n",nr);
      nr=0;
      return(-1);
    }
    nr=0;

    makeSurfaces();
    getElemNormalen( e_enqire, node, anz->e );
    realloc_colNr();
    return(anz->emax);
  }
  return(-1);
     
}



int createText(int nodenr, int x, int y )
{
  int i, t=-1;

  for(i=0; i<anz->t; i++)
  {
    if(ntext[i].node_nr==nodenr)
    {
      if(ntext[i].text) free(ntext[i].text);
      t=i; break;
    }
  }
  if(t==-1) for(i=0; i<anz->t; i++)
  {
    if(!ntext[i].node_nr) { t=i; break; }
  }
  if(t==-1)
  {
    if ((ntext = (Texts *)realloc( (Texts *)ntext, (anz->t+1)*sizeof(Texts)) ) == NULL )
    { printf("\n\nERROR: realloc failure in createText\n\n"); return(-1); }
    t=anz->t;
    ntext[t].tx= -1;
    ntext[t].ty= 0;
    anz->t++;
  }
  if(y>-1)
  {
    ntext[t].tx= (double)x;
    ntext[t].ty= (double)y;
  }
  ntext[t].nFlag= 1;
  ntext[t].vFlag= 1;
  ntext[t].tFlag= 0;
  ntext[t].fFlag= 0;
  ntext[t].pFlag= 0;
  ntext[t].text=NULL;
  ntext[t].node_nr=nodenr;
  updateDispLists();
  return(t);
}



int createPoint( int x, int y )
{
  /* new point is located in an area parallel to the screen and running through centerPnt */
  char     name[MAX_LINE_LENGTH];
  int      nr=-1;
  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];
  
  static GLdouble wx, wy, wz;  /*  returned world x, y, z coords  */
  static GLdouble nx, ny, nz;  /*  new world x, y, z coords  */
  static int flag;

  printf("createPoint\n");

  glutSetWindow( w1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  moveModel();
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);
  
  flag=gluProject( centerPnt[0], centerPnt[1], centerPnt[2], mvmatrix, projmatrix,
     viewport,  &wx, &wy, &wz);
  if (flag==GL_FALSE)
    printf("WARNING: Malfunction in createPoint(), please reselect\n");
  
  wx=(GLdouble)x; wy=(GLdouble)(viewport[3]-y);
  //printf (" x,y: %d %d Win coords are %d (%lf, %lf, %lf)\n",x,y,  flag,  wx   ,  wy   , wz     );
  flag=gluUnProject ( wx, wy, wz, mvmatrix, projmatrix, viewport, &nx, &ny, &nz);
  //printf ("new World coords are %d (%lf, %lf, %lf)\n", flag, nx, ny, nz);
  if (flag==GL_TRUE)
  {
    if(pntNr==-1) getNewName( name, "p" );
    else { strcpy(name,point[pntNr].name); pntNr=-1; }
    printf(" create point:%s %lf %lf %lf\n", name, nx, ny, nz );
    nr  = pnt( name, nx, ny, nz, 0);
  }
  else
    printf("WARNING: Malfunction in createPoint(), please reselect\n");
  
  updateDispLists();
  return(nr);
}


/* pnts: 3 points defining a plane, 1st pnt is startpnt of line, l:length of line */
/* return point-index of line-end-point p2 or -1 if failed */
int normalLine(char *name, int *pnts, double l)
{
  int i,p1,p2;
  double v0[3], v1[3], v2[3];
  char buf[MAX_LINE_LENGTH];
  double vpnt[3][3];

  if (!l) return(-1);
  for(i=0; i<3; i++)
  {
    //printf("normalLine %d pnt:%s\n", i, point[pnts[i]].name);
    vpnt[i][0]=point[pnts[i]].px* scale->w+scale->x;
    vpnt[i][1]=point[pnts[i]].py* scale->w+scale->y;
    vpnt[i][2]=point[pnts[i]].pz* scale->w+scale->z;
  }
  p1=pnts[0];

  for(i=0; i<3; i++) { v0[i]=vpnt[1][i]-vpnt[2][i]; v1[i]=vpnt[2][i]-vpnt[0][i]; }
  v_prod(v0,v1,v2);
  v_norm(v2,v2);
  v_scal(&l,v2,v2);
  for(i=0; i<3; i++) { v2[i]+=vpnt[0][i]; }
  
  /* create the normal point */
  p2= getNewName( buf, "p" );
  p2=pnt( buf, v2[0], v2[1], v2[2], 1 );
  if ( p2 >-1 )
  {
    /* generate line */
    if(name) line_i( name, p1, p2, -1, ddiv, 1, 0 );
  }
  return(p2);
}



void moveLineEndPoint(int lineNr, int pntNr, double llength)
{
  int p1,p2, flag=0;
  double P1[3], P2[3], u, eva[3], va[3], p0p1_2[3];
  
    p1=line[lineNr].p1;
    p2=line[lineNr].p2;

    /* determine which side of the line has to be moved */
    if(pntNr==p1) flag=-1;
    else if(pntNr==p2) flag=1;
    else
    {
      printf("ERROR: selected point:%s is no line endpoint\n", point[pntNr].name);
      return;
    }
    u=flag*llength;
    u/=scale->w;

    /* calc direction */
    if(line[lineNr].typ=='s')
    {
      if(flag==-1) p2=set[line[lineNr].trk].pnt[1];
      else  p1=set[line[lineNr].trk].pnt[set[line[lineNr].trk].anz_p-2];
    }
    P1[0]=point[p1].px;
    P1[1]=point[p1].py;
    P1[2]=point[p1].pz;
    P2[0]=point[p2].px;
    P2[1]=point[p2].py;
    P2[2]=point[p2].pz;
    v_result( P1, P2, p0p1_2);
    v_norm(p0p1_2,eva);
    v_scal(&u,eva, va); 
    point[pntNr].px+=va[0];
    point[pntNr].py+=va[1];
    point[pntNr].pz+=va[2];
    printf("moved by dxyz= %lf %lf %lf\n",
      (va[0]* scale->w),
      (va[1]* scale->w),
      (va[2]* scale->w));
}


int createLine( char *apnt, int flag )
{
  int i,j=0, nr=-1;
  int imax=0,jmax=0;
  static int px, p1, p2, pc=-1, pm=-1, ps[1000], seq=1;
  double P1P2[3], evamax,evbmax,a,b,c,d;
  double P1[3], P2[3], Pc[3], Pm[3], pbuf[3], u;
  char name[MAX_LINE_LENGTH], setname[MAX_LINE_LENGTH];
  double pmp1[3], pmp2[3], pmp1_2[3], pmp2_2[3], nm12[3];
  double eva[3], evb[3], va[3], vb[3], p0p1_2[3], p0p2_2[3], vr[3];


  if (flag==0) p1  = getPntNr( apnt );
  else if (flag==1)
  {
    p2  = getPntNr( apnt );
    if (p1!=p2)
    {
      if (pc>-1)
      {
        P1[0]=point[p1].px;
        P1[1]=point[p1].py;
        P1[2]=point[p1].pz;
        Pc[0]=point[pc].px;
        Pc[1]=point[pc].py;
        Pc[2]=point[pc].pz;
        /* radius berechnen */
        v_result(Pc, P1, vr);        
        printf(" R:%lf\n", v_betrag(vr)*scale->w);

        if(lineNr==-1) getNewName( name, "l" );
        else { strcpy(name,line[lineNr].name); lineNr=-1; }
        printf(" create line:%s %s %s %s %d\n", name, point[p1].name, point[p2].name, point[pc].name, ddiv );
        nr= line_i( name, p1, p2, pc, ddiv, 1, 'a' );
        pc=-1;
      }
      else if (pm>-1)
      {
        P1[0]=point[p1].px;
        P1[1]=point[p1].py;
        P1[2]=point[p1].pz;
        P2[0]=point[p2].px;
        P2[1]=point[p2].py;
        P2[2]=point[p2].pz;
        Pm[0]=point[pm].px;
        Pm[1]=point[pm].py;
        Pm[2]=point[pm].pz;
    
        /* vprod nm12 = pmp1 x pmp2 */
        v_result( Pm, P1, pmp1 );
        v_result( Pm, P2, pmp2 );
        v_prod( pmp1, pmp2, nm12 );
        if(v_betrag(nm12)==0.) return(-1);
    
        /* Vector halfway between pm and p1 or p2 */
        u=0.5;
        v_scal( &u, pmp1, pmp1_2);
        v_scal( &u, pmp2, pmp2_2);
        v_add( Pm, pmp1_2, p0p1_2); 
        v_add( Pm, pmp2_2, p0p2_2); 


        /* Vector in direction to PC, vprod va = nm12 x pmp1 and vb =  pmp2 x nm12 */
        v_prod( nm12, pmp1, va );
        v_prod( pmp2, nm12, vb );
        v_norm( va, eva );
        v_norm( vb, evb );
    
        // calculation of the intersection between eva and evb
	// two methods, method 1
	evbmax=evamax=0.;
        for (i=0; i<3; i++) if(abs(evb[i])>abs(evbmax)) { evbmax=evb[i]; imax=i; }
        for (j=0; j<3; j++) if((abs(eva[j])>abs(evamax))&&(imax!=j)) { evamax=eva[j]; jmax=j; }
	if(eva[jmax]==0.)
	{
	  evbmax=evamax=0.;
          for (i=0; i<3; i++) if((abs(evb[i])>abs(evbmax))&&(jmax!=i)) { evbmax=evb[i]; imax=i; }
          for (j=0; j<3; j++) if((abs(eva[j])>abs(evamax))&&(imax!=j)) { evamax=eva[j]; jmax=j; }
	}
	if((eva[jmax]==0.)||(evb[imax]==0.))
	{
	  printf(" ERROR: could not create the cener point\n");
  	  return(-1);
	}
	i=imax;
	j=jmax;
	a=eva[i]/evb[i];
	c=evb[j]/eva[j];
	if(a*c!=1.) goto calcPC;
	printf(" ERROR: could not create the cener point\n");
	return(-1);
      calcPC:;
	v_result(pmp1_2,pmp2_2,P1P2);
	b=P1P2[j]/eva[j];
	d=P1P2[i]/evb[i];
	u=(a*b-d)/(1.-a*c);

	/* method 2
        // determine abs-max-component of nm12
        u=0.;

        for (i=0; i<3; i++) if (nm12[i]*nm12[i]>u) { u=nm12[i]*nm12[i]; j=i; }

        // calculation of the intersection between eva and evb
        if (j==0)
	{
          u=(p0p2_2[1]-p0p1_2[1]-eva[1]*(p0p2_2[2]-p0p1_2[2])/eva[2])
           /(eva[1]*evb[2]/eva[2] - evb[1]);
        }
        else if (j==1)
	{
          u=(p0p2_2[0]-p0p1_2[0]-eva[0]*(p0p2_2[2]-p0p1_2[2])/eva[2])
           /(eva[0]*evb[2]/eva[2] - evb[0]);
        }
        else if (j==2)
	{
          u=(p0p2_2[0]-p0p1_2[0]-eva[0]*(p0p2_2[1]-p0p1_2[1])/eva[1])
           /(eva[0]*evb[1]/eva[1] - evb[0]);
        }
        else 
	{
          printf(" ERROR: in createLine, nm12 in error:%d\n", j);
          return(-1);
        }
	*/

        /*  centerpoint Pc = p0p2_2+ evb*u */
        v_scal( &u, evb, pbuf );
        v_add( pbuf, p0p2_2, Pc );

        /* radius berechnen */
        v_result(Pc, P1, vr);        
        printf(" R:%lf\n", v_betrag(vr)*scale->w);

	/*
        printf( "Pc:%lf %lf %lf\n", Pc[0], Pc[1], Pc[2] ); 
        printf(" P1:%s P2:%s PM:%s\n", point[p1].name, point[p2].name,  point[pm].name);
        printf(" u:%lf v:%lf\n", u*scale->w, v*scale->w);
        printf("nm12: %lf %lf %lf\n", nm12[0], nm12[1], nm12[2]);    
        printf("eva: %lf %lf %lf\n", eva[0], eva[1], eva[2]);    
        printf("evb: %lf %lf %lf\n", evb[0], evb[1], evb[2]);    
	*/


        /* centerpoint anlegen pc */
        getNewName( name, "p" );
        printf(" create point:%s %lf %lf %lf\n", name, Pc[0], Pc[1], Pc[2] );
        pc= pnt( name, Pc[0], Pc[1], Pc[2], 0 );

        /* kreisbogen definieren */
        if(lineNr==-1) getNewName( name, "l" );
        else { strcpy(name,line[lineNr].name); lineNr=-1; }
        printf(" create line:%s %s %s %s %d\n" , name, point[p1].name, point[p2].name, point[pc].name, ddiv  );
        nr= line_i( name, p1, p2, pc, ddiv, 1, 'a' );

        pc=-1;
        pm=-1;
      }
      else if (seq>1)
      {
        getNewName( setname, "se" );
        ps[0]=p1;
        ps[seq]=p2;
        for (i=0; i<=seq; i++)
	{
          pre_seta( setname, "ps", point[ps[i]].name);
        }
        nr=getSetNr( setname); 
        /* set[nr].type=1; */
        if(lineNr==-1) getNewName( name, "l" );
        else { strcpy(name,line[lineNr].name); lineNr=-1; }
        printf(" create line:%s %s %s %s %d\n"
        , name, point[p1].name, point[p2].name, setname, ddiv );
        nr= line_i( name, p1, p2, nr, ddiv, 1, 's' );
	seq=1;
      }
      else
      {
        if(lineNr==-1) getNewName( name, "l" );
        else { strcpy(name,line[lineNr].name); lineNr=-1; }
        printf(" create line:%s %s %s %d\n", name, point[p1].name, point[p2].name, ddiv );
        nr= line_i( name, p1, p2, 0, ddiv, 1, 0 );
      }
      p1=p2;
      }
    else errMsg(" p1==p2, try again\n");
    }
  else if (flag==2) /* center point defined */
  {
    pc  = getPntNr( apnt );
  }
  else if (flag==3) /* mit-point defined, create a centerpoint */
  {
    pm  = getPntNr( apnt );
  }
  else if (flag==4) /* seq-point defined, create a set */
  {
    ps[seq]  = getPntNr( apnt );
    if (seq<1000) seq++;
  }
  else if (flag==5)
  {
    if(lineNr==-1)
    {
      printf("ERROR: select line with key l first\n");
      return(-1);
    }
    px  = getPntNr( apnt );
    moveLineEndPoint( lineNr, px, pickbuf);
    for (i=0; i<anzGeo->l; i++) repLine(i);    

    lineNr=-1; 
  }
  return (nr);
}



void movePntAlongNurbs( int nurbsNr, int lmaster)
{
  int i,n;
  double vl, pm[3], pmpn[3], pnp1[3], pnp2[3], pnt[5], vnn[3], vns[3], vnc[3], vn[3];
  Points pntbuf, pntvn;

  // determine the correction direction

  // get uv coordinates of point
  n=set[line[lmaster].trk].pnt[1];
  pnt[0]=point[n].px;
  pnt[1]=point[n].py;
  pnt[2]=point[n].pz;
  proj1PntToNurbs( nurbsNr, pnt);
  // get the normal of the nurbs at point
  evalNurbsWithNormalVector( nurbsNr, 1, &pnt[3], &pnt[4], &pntbuf, &pntvn);
  v_norm(&pntvn.px, vnn);
  // get the normal of the spline
  v_result( &point[n].px,&point[line[lmaster].p1].px,pnp1);
  v_result( &point[n].px,&point[line[lmaster].p2].px,pnp2);
  v_prod(pnp1,pnp2,vns);
  v_norm(vns, vns);
  // get the normal on both (correction direction)
  v_prod(vnn,vns,vn);
  v_norm(vn, vn);
  // get the correction direction
  v_prod(vnn,vn,vnc);
  v_norm(vnc, vn);
  // projection of pmpn on that direction gives the correction
  v_add(&point[line[lmaster].p1].px, &point[line[lmaster].p2].px, pm);
  for(i=0; i<3; i++) pm[i]*=0.5;
  v_result( pm,&point[n].px,pmpn);
  vl=v_sprod( vn, pmpn)*-1.;
  v_scal(&vl,vn,vn);

  v_add(pnt, vn, &point[n].px);
}


// return 1 if the location is inside
int check_lineLoop(double *cg, int np, double *pnt_u, double *pnt_v, int e, int *tri)
{
  int i,j,n;
  Points *ppre=NULL;
  int    inull, iplus, iminus;
  double  vN0N1[3], vN1N2[3], vN2N0[3], vN0N[3], vN1N[3],  vN2N[3],
          vprod0[3], vprod1[3], vprod2[3], vsprod[3], vnorm[3];
  
  if( (ppre = (Points *)realloc(ppre, (np+1)*sizeof(Points) )) == NULL )
  { printf(" ERROR: malloc failure\n\n");  return(-1); }
  for(i=0; i<np; i++)
  {
    ppre[i].px=pnt_u[i];
    ppre[i].py=pnt_v[i];
    ppre[i].pz=0.;
  }

  /* decrement the pointer, tri[0] points to ppre[1] */
  ppre--;

  // go over all tri and determine if cg is inside
  j=0;
  for(n=0; n<e; n++)
  {
    /* berechne die vektoren die von den elem-nodes zum projezierten punkt laufen */
    v_result(&ppre[tri[j]].px,&ppre[tri[j+1]].px, vN0N1);
    v_result(&ppre[tri[j+1]].px,&ppre[tri[j+2]].px, vN1N2);
    v_result(&ppre[tri[j+2]].px,&ppre[tri[j]].px, vN2N0);

    v_result( &ppre[tri[j]].px, cg, vN0N);
    v_result( &ppre[tri[j+1]].px, cg, vN1N);
    v_result( &ppre[tri[j+2]].px, cg, vN2N);

    /* berechne den Normalenvektor auf der flaeche */
    v_prod( vN0N1, vN1N2, vnorm );
    j+=3;
    if(!v_norm( vnorm, vnorm)) continue;

    v_prod( vN0N1, vN0N, vprod0 );
    v_prod( vN1N2, vN1N, vprod1 );
    v_prod( vN2N0, vN2N, vprod2 );

    vsprod[0]= v_sprod( vnorm, vprod0);
    vsprod[1]= v_sprod( vnorm, vprod1);
    vsprod[2]= v_sprod( vnorm, vprod2);
   
    /*   kontrolle ob ein vsprod=0 ist. Dann liegt naemlich                */
    /*   der zu kontrollierende node genau ueber einer der umrandungen     */
    /*   der kontroll-flaeche (vprod ist dann senkrecht auf der normalen). */
    /*   Dann werden vsprod mit den gleichen vorzeichen zusammen-          */
    /*   gezaelt, dabei gelten die nuller als joker                        */

    inull=0;
    for (i=0; i<3; i++) if (vsprod[i] == 0.) inull++;
   
    iplus=inull;
    for (i=0; i<3; i++) if (vsprod[i] > 0.) iplus++;
   
    iminus=inull;
    for (i=0; i<3; i++) if (vsprod[i] < 0.) iminus++;

    //printf (" %d vsprod: %lf %lf %lf inull:%d iplus:%d iminus:%d \n",n,vsprod[0],vsprod[1],vsprod[2], inull, iplus, iminus );
   
    /*   wenn alle vectorprodukte das gleiche vorzeichen haben             */
    /*   dann ist cg innerhalb der flaeche                           */
    if ((iplus == 3)||(iminus == 3)) { ppre++; free(ppre);  return(1); }
  }
  ppre++; free(ppre); return(0);
}



#define TESTX1 0
int reduceOuterLoop(int smaster, int nls, int *ls, char *os, char *typs, int *nsurfs, int **ptr_innerSurf)
{
  int i,j=0, n=0, snew;
  int     flag1, flag2, err=0, lcount=0, lastlcount, closedLoopFlag;
  int     p1,p2, pnt=0, sl, cl, l0;
  int     line0_p1, line0_p2;
  int     nc=0;             /* number of curves */
  int     *cnr=NULL;        /* number of lines in each curve */
  int     *sortLines=NULL;
  char    *sortOri=NULL, *tmpOri=NULL, *sortTyp=NULL;
  char    name[MAX_LINE_LENGTH];
  int     *innerSurf;
  innerSurf=*ptr_innerSurf;
#if TESTX1    
  for (i=0; i<nls; i++)
  {
    printf(" l:%d %s\n",i,line[ls[i]].name);
  }
#endif
  if(nls<1) return(-1);

  if((sortOri = (char *)realloc((char *)sortOri, (nls)*sizeof(char)) ) == NULL )
  { printf(" ERROR: realloc failure in orientSurf()\n\n"); return(-1); }
  if((sortTyp = (char *)realloc((char *)sortTyp, (nls)*sizeof(char)) ) == NULL )
  { printf(" ERROR: realloc failure in orientSurf()\n\n"); return(-1); }
  if((sortLines=(int *)realloc((int *)sortLines, (nls)*sizeof(int) ) )==NULL)
  { printf(" ERROR: realloc failure in orientSurf()\n\n"); return(-1); }
  if((tmpOri = (char *)realloc((char *)tmpOri, (nls)*sizeof(char)) ) == NULL )
  { printf(" ERROR: realloc failure in orientSurf()\n\n"); return(-1); }

  /* loesche alle linien orientierungen */
  for (i=0; i<nls; i++)
  {
    tmpOri[i]=os[i];
    os[i]=0;
  }

  do
  {
    lastlcount=lcount;
#if TESTX1    
    printf("\n curve starts with n:%d lcount:%d \n",n,lcount);
#endif

    /* suche den punkt einer linie die mit der ersten linie verbunden ist */
    /* und orientiere die erste linie */
    l0=ls[n];
    line0_p1=line[l0].p1;
    line0_p2=line[l0].p2;

    /* gehe durch alle Linien und suche die zweite Linie */
    flag1=flag2=0;
    for (i=n+1; i<nls; i++)
    {
      p1=line[ls[i]].p1;
      p2=line[ls[i]].p2;
#if TESTX1    
      printf ("l:%s  p1:%s p2:%s  \n", line[ls[i]].name, point[p1].name, point[p2].name );
#endif
  
      if(tmpOri[n]=='+')
      {
        if      ( line0_p2 == p2 ) flag2=2;
        else if ( line0_p2 == p1 ) flag2=1;
        else if ( line0_p1 == p2 ) flag1=2;
        else if ( line0_p1 == p1 ) flag1=1;
      }
      else
      {
        if      ( line0_p1 == p2 ) flag1=2;
        else if ( line0_p1 == p1 ) flag1=1;
        else if ( line0_p2 == p2 ) flag2=2;
        else if ( line0_p2 == p1 ) flag2=1;
      }

      if (flag1==1)
      {
          os[n]='-';
          os[i]='+';
          pnt=p2;
          break;
      }
      else if (flag1==2)
      {
          os[n]='-';
          os[i]='-';
          pnt=p1;
          break;
      }
      else if (flag2==1)
      {
          os[n]='+';
          os[i]='+';
          pnt=p2;
          break;
      }
      else if (flag2==2)
      {
          os[n]='+';
          os[i]='-';
          pnt=p1;
          break;
      }
      else pnt=0;
    }

    if (i==nls)
    {
      if(typs[n]=='l') errMsg("WARNING: line[%d]:%s not connected\n",n,line[ls[n]].name);
      //typs[n]='d';
      if(typs[n]=='l')  typs[n]='L'; else typs[n]='C';
      n++;
      if(n==nls)
      {
        errMsg("WARNING: surf not closed, no orientation possible\n");
        for (i=0; i<nls; i++)
        {
          os[i]='+';
        }
	err=-1;
        goto orientSurfError;
      }
      else continue;
    }
    else
    {
      sortLines[lcount]=ls[n];      /* merke dir die reihenfolge in der die linien kommen */
      sortOri[lcount]=  os[n];      /* merke dir die dazugehoerigen orientationen */
      sortTyp[lcount]=  typs[n];
      lcount++;
      sortLines[lcount]=ls[i];      /* merke dir die reihenfolge in der die linien kommen */
      sortOri[lcount]=  os[i];      /* merke dir die dazugehoerigen orientationen */
      sortTyp[lcount]=  typs[i];
      lcount++;
    }

#if TESTX1    
    printf ("Curve:%d Start-line:%s ori:%c 2line:%s ori:%c \n",nc,line[ls[n]].name, os[n], line[ls[i]].name, os[i] );
#endif
    
    closedLoopFlag=0;
    for (i=0; i<nls; i++) /* suche die anschliessende durch vergleich (with inner second loop!)*/
    {
      for (j=n+1; j<nls; j++) /* suche die anschliessende durch vergleich */
      {
        if (!os[j])              /* wenn ori noch == 0 */
        {
#if TESTX1    
          printf ("check pnt:%s pnt:%d j:%d\n", point[pnt].name, pnt, j);
#endif
          if (typs[j]=='c')
          {
            sl=ls[j];
            cl=lcmb[sl].nl-1;
            if(lcmb[sl].o[0]=='+')     p1=line[lcmb[sl].l[0]].p1;
            else                       p1=line[lcmb[sl].l[0]].p2;
            if(lcmb[sl].o[cl]=='+')   p2=line[lcmb[sl].l[cl]].p2;
            else                      p2=line[lcmb[sl].l[cl]].p1;
#if TESTX1    
          printf ("compare i:%d j:%d lcmb:%s pnt:%d p1:%d p2:%d\n", i,j, lcmb[ls[j]].name, pnt, p1, p2 );
#endif
          }
          else
          {
            p1=line[ls[j]].p1;
            p2=line[ls[j]].p2;
#if TESTX1    
          printf ("compare i:%d j:%d line:%s ori:%c pnt:%d p1:%d p2:%d\n", i,j, line[ls[j]].name, os[j], pnt, p1, p2 );
#endif
          }

          if ( p1 == pnt )
          {
	    pnt=p2; os[j]='+';
#if TESTX1    
          printf ("found i:%d j:%d p1:%d p2:%d  next pnt:%d\n", i,j, p1, p2, pnt );
#endif
            sortLines[lcount]=ls[j]; sortOri[lcount]=os[j];
            sortTyp[lcount]=typs[j]; lcount++;
          }
          else if ( p2 ==  pnt )
          {
	    pnt=p1; os[j]='-';
#if TESTX1    
          printf ("found i:%d j:%d  p1:%d p2:%d  next pnt:%d\n", i,j, p1, p2, pnt );
#endif
            sortLines[lcount]=ls[j]; sortOri[lcount]=os[j];
            sortTyp[lcount]=typs[j]; lcount++;
          }
        }
      }

#if TESTX1    
      printf ("check connected lines\n");
#endif

      /* check if the curve is closed  */
      if(typs[n]=='l')
      {
        if((pnt==line[ls[n]].p1)||(pnt==line[ls[n]].p2))
        {
          /* number of lines in each curve */
          if((cnr=(int *)realloc((int *)cnr, (nc+1)*sizeof(int)) )==NULL)
          { printf(" ERROR: realloc failure\n"); exit(-1); }
          cnr[nc]=lcount-lastlcount;
          nc++;
	  closedLoopFlag=1;
#if TESTX1    
          printf ("line, closed line loop found, n:%d\n",n);
#endif
        }
      }
      if(closedLoopFlag==1) break;
    }
    if(closedLoopFlag==0) lcount=lastlcount;

    
    /* check for unoriented lines */
    err=0;
    for (n=0; n<nls; n++)
    {
      //if( (!os[n])&&(typs[n]!='d'))
       if( (!os[n]) && ((typs[n]!='L')&&(typs[n]!='C')) )
       {
        err--;
        break;
      }
    }
  }while(err<0);

#if TESTX1    
  printf("lcount:%d nls:%d\n",lcount,nls);
#endif

  for(j=0; j<cnr[0]; j++)
  {
#if TESTX1    
    printf("c:%d l:%d %s %c %c\n",i, j,line[sortLines[j]].name, sortOri[j], sortTyp[j]);
#endif
    ls[j]=sortLines[j];
    os[j]=sortOri[j];
    typs[j]=sortTyp[j];
  }
  lcount=cnr[0];

  if((innerSurf= (int *)realloc((int *)innerSurf, (nc+1)*sizeof(int))) == NULL )
      { printf("ERROR: realloc failure\n"); return(-1); }
  for(i=1; i<nc; i++)
  {
#if TESTX1    
    for(j=lcount; j<lcount+cnr[i]; j++)
    {
      printf("c:%d l:%d %s %c %c\n",i, j,line[sortLines[j]].name, sortOri[j], sortTyp[j]);
    }
#endif
    getNewName( name, "s" );
    snew= surface_i( name, surf[smaster].ori, surf[smaster].sh, cnr[i], &sortOri[lcount], &sortLines[lcount], &sortTyp[lcount] );
    if(snew<0) return(-5);
    surf[snew].etyp=surf[smaster].etyp;
    surf[snew].eattr=surf[smaster].eattr;
    surf[snew].elock=surf[smaster].elock;
    if((surf[smaster].eparm!=NULL)&&(strlen(surf[smaster].eparm)))
    {
      if((surf[snew].eparm= (char *)realloc((char *)surf[snew].eparm, (strlen(surf[smaster].eparm)+1)*sizeof(char))) == NULL )
      { printf("ERROR: realloc failure\n"); return(-1); }
      strcpy(surf[snew].eparm,surf[smaster].eparm);
    }
    innerSurf[*nsurfs]=snew;
    (*nsurfs)++;
    printf("   surfnr:%d %s created\n",*nsurfs, name);
    lcount+=cnr[i];
  }
  *ptr_innerSurf=innerSurf;
  err=cnr[0];
 orientSurfError:;
  free(sortLines);
  free(sortOri);
  free(tmpOri);
  free(sortTyp);
  free(cnr);
  return(err);
}



/*
return(new surf)
if error
return
-1: malloc failed
-2: no line or surf selected
-4: the lines form only one loop, so no 2 surfs can be generated
-5: other error
*/
int splitSurf( int nl, int *ptr_lmaster, int smaster, int projFlag )
{
  int i,j,k,l,lmaster,n,s;
  int nlines, nsurfs=0, nls[2], linoffs;
  int snew=-1, sbuf, Stmp, sh_buf, nurbsnr, patch;
  char repFlag=0;
  static int *ls[2]={NULL,NULL};
  static char *os[2]={NULL,NULL};
  static char *typs[2]={NULL,NULL};
  static int *lines=NULL;
  static char *ori=NULL;
  static char *typ=NULL;
  static int *innerSurf=NULL;
  char name[MAX_LINE_LENGTH];
  
  double p1p2[3], p1pn[3],vn[3],vl,pnt[5],ubuf,vbuf;
  double lineLoop_uvw[1000][3];

  int np, e,p, ii, ll, npc[1], loribuf[2], lbuf[2], ibuf, obuf, secondLoop, splitFlag;
  int *lloop=NULL;
  double pl[2][3], p01[3], es1s0[3], el1l0[3], sp1;
  double *pnt_u=NULL, *pnt_v=NULL;
  int    *pnt_flag=NULL;
  int    *tri=NULL;
  int *usedLoops=NULL;

  printf(" split surf: %s\n", surf[smaster].name);
  for(i=0; i<nl; i++) printf("split-lines:%s\n",line[ptr_lmaster[i]].name);
  
  if(surf[smaster].npgn) repFlag=1;

  /* make sure the master surface is oriented */
  orientSurf( smaster );

  /* proj line to surf */
  if((projFlag==1)&&(surf[smaster].sh>-1))
  {
    for(l=0; l<nl; l++)
    {
      lmaster=ptr_lmaster[l];
      if(( line[lmaster].typ != 's')||(set[line[lmaster].trk].anz_p==3))
      {
        // 1. make spline with just one point
        // loop several times:
        // 2. proj to nurbs
        // 3. set trk pnt back to node pos
        convertLine( lmaster, 1);
        if(shape[surf[smaster].sh].type!=4)
        {
  	  sbuf= s= surfToNurs(smaster);
          if(s<0) return(-5);
        }
        else sbuf=shape[surf[smaster].sh].p[0];
  
        // project the point
        ubuf=vbuf=-1.;
        for( i=0; i<60; i++)
        {
          n=set[line[lmaster].trk].pnt[1];
  	  pnt[0]=point[n].px;
  	  pnt[1]=point[n].py;
  	  pnt[2]=point[n].pz;
          proj1PntToNurbs(sbuf, pnt);
  	  if( (abs(pnt[3]-ubuf)/nurbs[sbuf].uknt[nurbs[sbuf].u_nknt-1]) +
  	    (abs(pnt[4]-vbuf)/nurbs[sbuf].vknt[nurbs[sbuf].v_nknt-1]) < 1e-6)
  	  break;
  	
          repLine(lmaster);
          // move pnt back to line mid position
          j=(line[lmaster].nip-1)/2;
          point[n].px=line[lmaster].ip[j-1]  ;
          point[n].py=line[lmaster].ip[j];
          point[n].pz=line[lmaster].ip[j+1];
  
          // if the point is not at the line between line[lmaster].p1 and line[lmaster].p2
          // correct the projected point 
          v_result( &point[line[lmaster].p1].px,&point[line[lmaster].p2].px,p1p2);
          v_result( &point[line[lmaster].p1].px,&point[n].px,p1pn);
          v_prod( p1p2, p1pn, vn );
          vl=v_betrag(vn);
          if(vl!=0.)
          {
  	    pnt[0]=point[n].px;
  	    pnt[1]=point[n].py;
  	    pnt[2]=point[n].pz;
            proj1PntToNurbs(sbuf, pnt);
            movePntAlongNurbs( sbuf, lmaster);
            repLine(lmaster);
            // move pnt back to line mid position
            j=(line[lmaster].nip-1)/2;
            point[n].px=line[lmaster].ip[j-1]  ;
            point[n].py=line[lmaster].ip[j];
            point[n].pz=line[lmaster].ip[j+1];
          }
  	  ubuf=pnt[3];
  	  vbuf=pnt[4];
        }
        if( ((abs(pnt[3]-ubuf)/nurbs[sbuf].uknt[nurbs[sbuf].u_nknt-1]) >1e-6) &&
  	  ((abs(pnt[4]-vbuf)/nurbs[sbuf].vknt[nurbs[sbuf].v_nknt-1]) >1e-6) ) return(-5);
        
        // 4. make spline with div-1 points
        // 2. proj to nurbs
        convertLine( lmaster, line[lmaster].div-1);
        delPnt( 1, &n );
        // proj normal to smaster
        if(shape[surf[smaster].sh].type==4) projSetToNurbs(shape[surf[smaster].sh].p[0], set, line[lmaster].trk, point,&node);
        else
        {
          if(s>-1)
          {
            projSetToNurbs( s, set, line[lmaster].trk, point,&node);
            delNurs( 1, &s );
          }
        }
      }
      else // line was a spline with more than 3 pnts and trk pnts have to be kept!
      {
        // proj normal to smaster
        if(shape[surf[smaster].sh].type==4) projSetToNurbs(shape[surf[smaster].sh].p[0], set, line[lmaster].trk, point,&node);
        else
        {
          s= surfToNurs(smaster);
          if(s>-1)
          {
            projSetToNurbs( s, set, line[lmaster].trk, point,&node);
            delNurs( 1, &s );
          }
        }
      }
      repLine(lmaster);
    }
  }

  /* copy all lines (also from lcmbs) to an array */
  l=0;
  for(i=0; i<surf[smaster].c[0]; i++)
  {
    if(surf[smaster].typ[i]=='l')
    {
      if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      typ[l]='l';
      ori[l]=surf[smaster].o[i];
      //printf("add %s\n", line[surf[smaster].l[i]].name);
      lines[l++]=surf[smaster].l[i];
    }
    else
    {
      if ((lines = (int *)realloc(lines, (l+lcmb[surf[smaster].l[i]].nl)*sizeof(int)) ) == NULL )
        { printf("ERROR: malloc failure\n"); return(-1); }
      if ((ori = (char *)realloc(ori, (l+lcmb[surf[smaster].l[i]].nl)*sizeof(char)) ) == NULL )
        { printf("ERROR: malloc failure\n"); return(-1); }
      if ((typ = (char *)realloc(typ, (l+lcmb[surf[smaster].l[i]].nl)*sizeof(char)) ) == NULL )
        { printf("ERROR: malloc failure\n"); return(-1); }

      if(surf[smaster].o[i]=='+')
      { 
        for (n=0; n<lcmb[surf[smaster].l[i]].nl; n++)
        {
          typ[l]='l';
          ori[l]=lcmb[surf[smaster].l[i]].o[n];
	  //printf("add %s\n", line[lcmb[surf[smaster].l[i]].l[n]].name);
          lines[l++]=lcmb[surf[smaster].l[i]].l[n];
	}
      }
      else
      { 
        for (n=lcmb[surf[smaster].l[i]].nl-1; n>=0; n--)
        {
          typ[l]='l';
          if(lcmb[surf[smaster].l[i]].o[n]=='+') ori[l]='-'; else ori[l]='+';
	  //printf("add %s\n", line[lcmb[surf[smaster].l[i]].l[n]].name);
          lines[l++]=lcmb[surf[smaster].l[i]].l[n];
	}
      }
    }
  }
  nlines=l;  // number of lines in the outer loop


  /* search the line loops starting from the 1st line */
  l=0;
  s=0;
  nls[0]=nls[1]=0;
  for(i=0; i<nlines; i++)
  {
    if ((ls[s] = (int *)realloc(ls[s], (nls[s]+1)*sizeof(int)) ) == NULL )
    { printf("ERROR: malloc failure\n"); return(-1); }
    if ((os[s] = (char *)realloc(os[s], (nls[s]+1)*sizeof(char)) ) == NULL )
    { printf("ERROR: malloc failure\n"); return(-1); }
    if ((typs[s] = (char *)realloc(typs[s], (nls[s]+1)*sizeof(char)) ) == NULL )
    { printf("ERROR: malloc failure\n"); return(-1); }
    ls[s][nls[s]] = lines[i];
    os[s][nls[s]] = ori[i];
    typs[s][nls[s]] = typ[i];
    //    printf("%d found l:%s\n", i, line[ls[s][nls[s]]].name);
    nls[s]++;

    /* Add the splitting lines to the current surface s when a point of lmaster is identical to the end-point of the actual line */
    /* then change the target surface s until all outer lines are evaluated */
    if(ori[i]=='+')
    {
      for(j=0; j<nl; j++)
      {
        lmaster=ptr_lmaster[j];
        if((line[lines[i]].p2==line[lmaster].p1)||(line[lines[i]].p2==line[lmaster].p2))
	{
	  loribuf[s]=ori[i];
	  lbuf[s]=lines[i];
	  // calc the vector in direction from s1 to s0
	  if(s==0)
	  {
	    p=line[lines[i]].nip-3;
	    pl[0][0]=line[lines[i]].ip[p]  ;
            pl[0][1]=line[lines[i]].ip[p+1];
            pl[0][2]=line[lines[i]].ip[p+2];
	    p-=3;
	    pl[1][0]=line[lines[i]].ip[p]  ;
            pl[1][1]=line[lines[i]].ip[p+1];
            pl[1][2]=line[lines[i]].ip[p+2];
            v_result( &pl[0][0], &pl[1][0], p01 );
            v_norm( p01, es1s0 );
	  }
	  s=!s;
          break;
	}
      }
    }
    else
    {
      for(j=0; j<nl; j++)
      {
        lmaster=ptr_lmaster[j];
        if((line[lines[i]].p1==line[lmaster].p1)||(line[lines[i]].p1==line[lmaster].p2))
	{
	  loribuf[s]=ori[i];
	  lbuf[s]=lines[i];
	  if(s==0)
	  {
	    p=0;
	    pl[0][0]=line[lines[i]].ip[p]  ;
            pl[0][1]=line[lines[i]].ip[p+1];
            pl[0][2]=line[lines[i]].ip[p+2];
	    p=3;
	    pl[1][0]=line[lines[i]].ip[p]  ;
            pl[1][1]=line[lines[i]].ip[p+1];
            pl[1][2]=line[lines[i]].ip[p+2];
            v_result( &pl[0][0], &pl[1][0], p01 );
            v_norm( p01, es1s0 );
	  }
	  s=!s;
	  break;
	}
      }
    }
  }
  if(nls[1]<1) return(-4);

  // add the splitting lines to both surfaces
  for(s=0; s<2; s++)
  {
    for(j=0; j<nl; j++)
    {
      lmaster=ptr_lmaster[j];
      if ((ls[s] = (int *)realloc(ls[s], (nls[s]+1)*sizeof(int)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      if ((os[s] = (char *)realloc(os[s], (nls[s]+1)*sizeof(char)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      if ((typs[s] = (char *)realloc(typs[s], (nls[s]+1)*sizeof(char)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      ls[s][nls[s]] = lmaster;
      os[s][nls[s]] = '+';
      typs[s][nls[s]] = 'l';
      nls[s]++;
    }
  }
  
  /* add the inner line-loops to the surf-lines ls[s][nls[s]] */

  /* get the uv average of the inner line loops and later compare with each inner tri of the 2 prelim-surfs. 
     use the average uv of the splitting line points to decide which surf */ 
  if(surf[smaster].nc>1)
  {
    // search the inner loops between the splitting lines
    // remember and skip them later
    if ( (usedLoops = (int *)calloc( (surf[smaster].nc+1), sizeof(int))) == NULL )
      printf("ERROR: calloc failed \n\n" );
    if ( (lloop = (int *)calloc( (anzGeo->l+1), sizeof(int))) == NULL )
      printf("ERROR: calloc failed \n\n" );

    // -1 Take a certain line loop
    //    go over all loop lines and compare ptr_lmaster[j] with the first point of the loop lines in loop direction
    // -2 if identical switch s to 1, inc j++, and store this line and all following lines in lloop[line]=s(1) until
    //    next ptr_lmaster[j] is  matched, switch s to 0 and keep j until
    //    - either ptr_lmaster[j] is matched again
    //      goto -2
    //    - or if the current loop line l gives lloop[l]=1. Then go to -1 
    
    // over all splitlines j up to nl-1
    j=0;
    do
    {
      lmaster=ptr_lmaster[j];
      linoffs=0;
      for(i=1; i<surf[smaster].nc; i++)
      {
	//printf("loop:%d\n", i);
        linoffs+=surf[smaster].c[i-1];
	if(usedLoops[i]) continue;
        s=0;
	splitFlag=0;
	secondLoop=0;
	// reset lloop[l]
	for(l=0; l<anzGeo->l; l++) { lloop[l]=0; }
	k=linoffs;
	do
        {
          //compare the line with lmaster
          l=surf[smaster].l[k];
     	  if(surf[smaster].typ[k]=='l')
          {
	    //printf("loop %d check line %s with %s flag:%d\n",i, line[l].name, line[lmaster].name, lloop[l]);
  	    if(lloop[l]==1) break;
            ibuf=0;
            if(surf[smaster].o[k]=='+')
	    {
	      obuf=1;
      	      if(line[l].p1==line[lmaster].p1) ibuf++;
      	      if(line[l].p1==line[lmaster].p2) ibuf++;
	    }
	    else
	    {
	      obuf=0;
      	      if(line[l].p2==line[lmaster].p1) ibuf++;
      	      if(line[l].p2==line[lmaster].p2) ibuf++;
	    }
	    // if a point matches switch s and calc the vector of the 1st line in loop[1] if s==1
            if(ibuf==1)
	    {
	      // found the loop, now add all loop-lines until the next split-line matches a loop-point
	      s++;
	      //if(s==2) goto foundLines;
	      if(s==2) { s=0; splitFlag=1; }
	      if(s==1)
	      {
	        j++;
	        lmaster=ptr_lmaster[j];
		if(obuf==1)
		{
	          p=0;
	          pl[0][0]=line[l].ip[p]  ;
                  pl[0][1]=line[l].ip[p+1];
                  pl[0][2]=line[l].ip[p+2];
	          p=3;
	          pl[1][0]=line[l].ip[p]  ;
                  pl[1][1]=line[l].ip[p+1];
                  pl[1][2]=line[l].ip[p+2];
                  v_result( &pl[0][0], &pl[1][0], p01 );
                  v_norm( p01, el1l0 );
		}
		else
		{
	          p=line[l].nip-3;
	          pl[0][0]=line[l].ip[p]  ;
                  pl[0][1]=line[l].ip[p+1];
                  pl[0][2]=line[l].ip[p+2];
	          p-=3;
	          pl[1][0]=line[l].ip[p]  ;
                  pl[1][1]=line[l].ip[p+1];
                  pl[1][2]=line[l].ip[p+2];
                  v_result( &pl[0][0], &pl[1][0], p01 );
                  v_norm( p01, el1l0 );
		}
	      }
	    }
            lloop[l]=s;
	    //printf("s:%d l:%s\n",s,line[l].name);
     	  }
          else
      	  {
	    //printf("lcmb %s\n", lcmb[l].name); 
            if(surf[smaster].o[k]=='+')
	    {
              for (ii=0; ii<lcmb[l].nl; ii++ )
              {
                ll=lcmb[l].l[ii];
                ibuf=0;
	        //printf("loop %d check line %s with %s flag:%d\n",i, line[ll].name, line[lmaster].name, lloop[ll]);
  	        if(lloop[ll]==1) break;
		if (lcmb[l].o[ii]=='+')
		{
	          obuf=1;
      	          if(line[ll].p1==line[lmaster].p1) ibuf++;
      	          if(line[ll].p1==line[lmaster].p2) ibuf++;
	        }
	        else
	        {
	          obuf=0;
      	          if(line[ll].p2==line[lmaster].p1) ibuf++;
      	          if(line[ll].p2==line[lmaster].p2) ibuf++;
	        }
                if(ibuf==1)
	        {
 	          // found the loop, now add all loop-lines until the next split-line matches a loop-point
	          s++;
	          //if(s==2) goto foundLines;
	          if(s==2) { s=0; splitFlag=1; }
	          if(s==1)
	          {
 		    j++;
		    lmaster=ptr_lmaster[j];
		    if(obuf==1)
		    {
	              p=0;
	              pl[0][0]=line[ll].ip[p]  ;
                      pl[0][1]=line[ll].ip[p+1];
                      pl[0][2]=line[ll].ip[p+2];
	              p=3;
	              pl[1][0]=line[ll].ip[p]  ;
                      pl[1][1]=line[ll].ip[p+1];
                      pl[1][2]=line[ll].ip[p+2];
                      v_result( &pl[0][0], &pl[1][0], p01 );
                      v_norm( p01, el1l0 );
		    }
		    else
		    {
	              p=line[ll].nip-3;
	              pl[0][0]=line[ll].ip[p]  ;
                      pl[0][1]=line[ll].ip[p+1];
                      pl[0][2]=line[ll].ip[p+2];
	              p-=3;
	              pl[1][0]=line[ll].ip[p]  ;
                      pl[1][1]=line[ll].ip[p+1];
                      pl[1][2]=line[ll].ip[p+2];
                      v_result( &pl[0][0], &pl[1][0], p01 );
                      v_norm( p01, el1l0 );
		    }
	          }
      	        }
		lloop[ll]=s;
		//printf("s:%d l:%s\n",s,line[ll].name);
	      }
	    }
            if(surf[smaster].o[k]=='-')
	    {
              for (ii=lcmb[l].nl-1; ii>=0; ii-- )
              {
                ll=lcmb[l].l[ii];
                ibuf=0;
	        //printf("loop %d check line %s with %s flag:%d\n",i, line[ll].name, line[lmaster].name, lloop[ll]);
  	        if(lloop[ll]==1) break;
		if (lcmb[l].o[ii]=='-')
		{
	          obuf=1;
      	          if(line[ll].p1==line[lmaster].p1) ibuf++;
      	          if(line[ll].p1==line[lmaster].p2) ibuf++;
	        }
	        else
	        {
	          obuf=0;
      	          if(line[ll].p2==line[lmaster].p1) ibuf++;
      	          if(line[ll].p2==line[lmaster].p2) ibuf++;
	        }
                if(ibuf==1)
	        {
 		  // found the loop, now add all loop-lines until the next split-line matches a loop-point
	          s++;
	          //if(s==2) goto foundLines;
	          if(s==2) { s=0; splitFlag=1; }
	          if(s==1)
	          {
 		    j++;
		    lmaster=ptr_lmaster[j];
		    if(obuf==1)
		    {
	              p=0;
	              pl[0][0]=line[ll].ip[p]  ;
                      pl[0][1]=line[ll].ip[p+1];
                      pl[0][2]=line[ll].ip[p+2];
	              p=3;
	              pl[1][0]=line[ll].ip[p]  ;
                      pl[1][1]=line[ll].ip[p+1];
                      pl[1][2]=line[ll].ip[p+2];
                      v_result( &pl[0][0], &pl[1][0], p01 );
                      v_norm( p01, el1l0 );
		    }
		    else
		    {
	              p=line[ll].nip-3;
	              pl[0][0]=line[ll].ip[p]  ;
                      pl[0][1]=line[ll].ip[p+1];
                      pl[0][2]=line[ll].ip[p+2];
	              p-=3;
	              pl[1][0]=line[ll].ip[p]  ;
                      pl[1][1]=line[ll].ip[p+1];
                      pl[1][2]=line[ll].ip[p+2];
                      v_result( &pl[0][0], &pl[1][0], p01 );
                      v_norm( p01, el1l0 );
		    }
	          }
		}
		lloop[ll]=s;
		//printf("s:%d l:%s\n",s,line[ll].name);
	      }
	    }
	  }
	  if(k<surf[smaster].c[i]+linoffs-1) k++;
	  else
	  {
	    k=linoffs;
	    if(secondLoop) break;
	    secondLoop=1;
	  }
	}while(1);
	//foundLines:;
	//return(-5);
	
        // determine which loop belongs to which surf and add them there (use vector direction vs es1s0)
	//if(s==2)
	if(splitFlag==1)
	{
	  //printf("found the splitting lines\n");
	  usedLoops[i]=1;
	  //v_print("es1s0 ", es1s0);
	  //v_print("el1l0 ", el1l0);
	  sp1=v_sprod(es1s0,el1l0);
	  //printf("sp:%f\n", sp1);
	  for(k=linoffs; k<surf[smaster].c[i]+linoffs; k++)
	  {
	    l=surf[smaster].l[k];
	    if(surf[smaster].typ[k]=='l')
	    {
	      //printf("%d add line %s to %d\n",k, line[l].name, lloop[l]);
	      s=lloop[l];
	      if(sp1>0.) s=!s;
              if ((ls[s] = (int *)realloc(ls[s], (nls[s]+1)*sizeof(int)) ) == NULL )
              { printf("ERROR: malloc failure\n"); return(-1); }
              if ((os[s] = (char *)realloc(os[s], (nls[s]+1)*sizeof(char)) ) == NULL )
              { printf("ERROR: malloc failure\n"); return(-1); }
              if ((typs[s] = (char *)realloc(typs[s], (nls[s]+1)*sizeof(char)) ) == NULL )
              { printf("ERROR: malloc failure\n"); return(-1); }
              ls[s][nls[s]] = l;
              os[s][nls[s]] = '+';
              typs[s][nls[s]] = 'l';
              nls[s]++;
            }
	    else
	    {
	      //printf("lcmb\n");
	      for (ii=0; ii<lcmb[l].nl; ii++)
	      {
		ll=lcmb[l].l[ii];
	        //printf("%d add line %s to %d\n",k, line[ll].name, lloop[ll]);
	        s=lloop[ll];
	        if(sp1>0.) s=!s;
                if ((ls[s] = (int *)realloc(ls[s], (nls[s]+1)*sizeof(int)) ) == NULL )
                { printf("ERROR: malloc failure\n"); return(-1); }
                if ((os[s] = (char *)realloc(os[s], (nls[s]+1)*sizeof(char)) ) == NULL )
                { printf("ERROR: malloc failure\n"); return(-1); }
                if ((typs[s] = (char *)realloc(typs[s], (nls[s]+1)*sizeof(char)) ) == NULL )
                { printf("ERROR: malloc failure\n"); return(-1); }
                ls[s][nls[s]] = ll;
                os[s][nls[s]] = '+';
                typs[s][nls[s]] = 'l';
                nls[s]++;
	      }
            }
	  }
	    
	  //for(s=0; s<2; s++) for(l=0; l<nls[s]; l++) printf("s:%d l:%d %s\n",s,l,line[ls[s][l]].name);
 	}
      }	      	      
      //printf("j:%d nl:%d\n", j,nl);
    }while(j<nl-1);
      
    free(lloop);

    // generate uv, copy from repSurf()
    if(( surf[smaster].name != (char *)NULL )&&(surf[smaster].sh>-1)&&(surf[smaster].nc>0))
    {
      
      // WARNING: The new line loops ls[s][nls[s]] for the splitted surfaces might not be just one closed loop but
      // several ones. The following code can not deal with that in a propper way. Instead the line loops are separated
      // into closed loops were the one which uses lines of the outer loop of smaster is kept and the other ones are
      // directly used to create new surfaces. This is dangerous because inner holes can not be distributed in a correct
      // manner with the below code any more since the new surfs are not regarded in the distribution of the holes!
      nsurfs=0;
      for(s=0; s<2; s++) nls[s]=reduceOuterLoop(smaster, nls[s], &ls[s][0], &os[s][0], &typs[s][0], &nsurfs, &innerSurf);
      if(nsurfs) printf(" nsurfs:%d innerSurf[0]:%s\n",nsurfs,surf[innerSurf[0]].name);
      else printf(" no inner surfs\n");

      /* if shape, generate prelim nurbs */
      sh_buf=-1;
      if(shape[surf[smaster].sh].type>=0) Stmp= surfToNurs(smaster);
      if(Stmp>-1)
      {
        sh_buf=surf[smaster].sh;
        surf[smaster].sh=shape_i( nurbs[Stmp].name, 4, Stmp, 0, 0, 0, 0,0,0);
        if(printFlag) printf (" interior changed to Nurbs: %s\n", nurbs[Stmp].name );
      }
      
      if(shape[surf[smaster].sh].type==4)
      {
        nurbsnr=shape[surf[smaster].sh].p[0];
  
        /* calculate the xyz and uv values and orientation of the trimming loops of the surface (patch) (might be included in orientSurf() later) */
        i=calcTrimLoops(nurbsnr, smaster);
        if(i!=0)
        {
          if(i>0)
  	  {
            printf(" ERROR: surf:%s could not be trimmed. All points are located on ambiguous edges. Try to fix the geometry manually.\n",surf[smaster].name);
            surf[smaster].fail=2;
  	  }
          goto nurbsCouldNotBeTrimmed;
        }
      }
      else goto nurbsCouldNotBeTrimmed;
  
      // evaluate the uv coordinates and determine the location of each inner line loop
      patch=surf[smaster].patch;

      // determine the uv coordinate of the cg of the inner loops in uv 
      for(i=1; i<nurbs[nurbsnr].nc[patch]; i++)
      {
	if(usedLoops[i]==1) continue;
	printf("check line loop %d \n",i);
	
        /* get the uv-values of the single curves which trim a surface */
        n=0;
        ubuf=vbuf=0.;
        //for(j=0; j<nurbs[nurbsnr].np[patch][i]/2; j++)
        for(j=0; j<nurbs[nurbsnr].np[patch][i]-1; j++)
        {
          ubuf+=nurbs[nurbsnr].uv[patch][i][n++];
          vbuf+=nurbs[nurbsnr].uv[patch][i][n++];
        }
        lineLoop_uvw[i][0]=ubuf*2./n;
	lineLoop_uvw[i][1]=vbuf*2./n;
        lineLoop_uvw[i][2]=0.;
      }
	  
      // To which surf belong the lineLoops?	  
      // first create splitted second surf w/o inner loops and fill them internaly with trias
      getNewName( name, "s" );
      snew= surface_i( name, surf[smaster].ori, surf[smaster].sh, nls[1], os[1], ls[1], typs[1] );
      printf("surf:%s created\n", name);
      if(snew<0) return(-5);
      surf[snew].etyp=surf[smaster].etyp;
      surf[snew].eattr=surf[smaster].eattr;
      surf[snew].elock=surf[smaster].elock;
      if((surf[smaster].eparm!=NULL)&&(strlen(surf[smaster].eparm)))
      {
        if((surf[snew].eparm= (char *)realloc((char *)surf[snew].eparm, (strlen(surf[smaster].eparm)+1)*sizeof(char))) == NULL )
        { printf("ERROR: realloc failure\n"); return(-1); }
        strcpy(surf[snew].eparm,surf[smaster].eparm);
      }
      repSurf(snew,1);
  
      // mesh the surf with the uv coordinates 
      patch=surf[snew].patch;
      np=n=0;
      k=1;
      npc[0]=nurbs[nurbsnr].np[patch][0]-1;
      if( (pnt_u = (double *)malloc( (npc[0]+1)*sizeof(double) )) == NULL )
      { printf(" ERROR: malloc failure\n\n");
        return(-1); }
      if( (pnt_v = (double *)malloc( (npc[0]+1)*sizeof(double) )) == NULL )
      { printf(" ERROR: malloc failure\n\n");
        return(-1); }
      for(j=0; j<npc[0]; j++)
      {
        pnt_u[np]=nurbs[nurbsnr].uv[patch][0][n++];
        pnt_v[np]=nurbs[nurbsnr].uv[patch][0][n++];
	np++;
      }
      e=mesh2d(&np, k, &npc[0], &pnt_u, &pnt_v, &pnt_flag, &tri, 0.4, 0.4, 1 );

      // search the inner triangle of the two surfaces in which the lineLoop_uvw[i] lies, then
      // add the line loops to the line buffers of this surface
      patch=surf[smaster].patch;
      linoffs=0;  
      for(i=1; i<nurbs[nurbsnr].nc[patch]; i++)
      {
        linoffs+=surf[smaster].c[i-1];
	if(usedLoops[i]==0) printf("check loop %d \n",i);
	else { printf("skip loop %d \n",i); continue; }
        if(!check_lineLoop(&lineLoop_uvw[i][0], np, pnt_u, pnt_v, e, tri))
        {
          //printf("i belongs to first (orig) surf nls[0], os[0], ls[0], typs[0]\n");
          // i belongs to first (orig) surf nls[0], os[0], ls[0], typs[0]
          s=0;
          for(j=0; j<surf[smaster].c[i]; j++)
          {
            if ((ls[s] = (int *)realloc(ls[s], (nls[s]+1)*sizeof(int)) ) == NULL )
            { printf("ERROR: malloc failure\n"); return(-1); }
            if ((os[s] = (char *)realloc(os[s], (nls[s]+1)*sizeof(char)) ) == NULL )
            { printf("ERROR: malloc failure\n"); return(-1); }
            if ((typs[s] = (char *)realloc(typs[s], (nls[s]+1)*sizeof(char)) ) == NULL )
            { printf("ERROR: malloc failure\n"); return(-1); }
            ls[s][nls[s]] = surf[smaster].l[j+linoffs];
            os[s][nls[s]] = surf[smaster].o[j+linoffs];
            typs[s][nls[s]] = surf[smaster].typ[j+linoffs];
            nls[s]++;
          }
        }
        else
        {
          //printf(" i belongs to second surf nls[1], os[1], ls[1], typs[1]\n");
          // i belongs to second surf nls[1], os[1], ls[1], typs[1]
          s=1;
          for(j=0; j<surf[smaster].c[i]; j++)
          {
            if ((ls[s] = (int *)realloc(ls[s], (nls[s]+1)*sizeof(int)) ) == NULL )
            { printf("ERROR: malloc failure\n"); return(-1); }
            if ((os[s] = (char *)realloc(os[s], (nls[s]+1)*sizeof(char)) ) == NULL )
            { printf("ERROR: malloc failure\n"); return(-1); }
            if ((typs[s] = (char *)realloc(typs[s], (nls[s]+1)*sizeof(char)) ) == NULL )
            { printf("ERROR: malloc failure\n"); return(-1); }
            ls[s][nls[s]] = surf[smaster].l[j+linoffs];
            os[s][nls[s]] = surf[smaster].o[j+linoffs];
            typs[s][nls[s]] = surf[smaster].typ[j+linoffs];
            nls[s]++;
          }
        }
      }
  
      /* redefine selected surf and define new surf */
      surface_i( surf[snew].name, surf[snew].ori, surf[snew].sh, nls[1], os[1], ls[1], typs[1] );
      surface_i( surf[smaster].name, surf[smaster].ori, surf[smaster].sh, nls[0], os[0], ls[0], typs[0] );

      free(pnt_u);
      free(pnt_v);
      free(pnt_flag);
      free(tri);
    nurbsCouldNotBeTrimmed:;

      /* restore the pointer to the shape */
      if(sh_buf>-1)
      {
        surf[smaster].sh=sh_buf;
        surf[snew].sh=sh_buf;
        for (i=0; i<nsurfs; i++) surf[innerSurf[i]].sh=sh_buf;
        for (i=0; i<nurbs[Stmp].u_npnt; i++) delPnt( nurbs[Stmp].v_npnt, nurbs[Stmp].ctlpnt[i] );
        delNurs( 1, &Stmp ); 
      }
    }
    else
    {
      return(-5);
    }
  }
  else
  {
    // create splitted surfs w/o inner loops and fill them internaly with trias
    getNewName( name, "s" );
    snew= surface_i( name, surf[smaster].ori, surf[smaster].sh, nls[1], os[1], ls[1], typs[1] );
    printf("surf:%s created\n", name);
    if(snew<0) return(-5);
    surf[snew].etyp=surf[smaster].etyp;
    surf[snew].eattr=surf[smaster].eattr;
    surf[snew].elock=surf[smaster].elock;
    if((surf[smaster].eparm!=NULL)&&(strlen(surf[smaster].eparm)))
    {
      if((surf[snew].eparm= (char *)realloc((char *)surf[snew].eparm, (strlen(surf[smaster].eparm)+1)*sizeof(char))) == NULL )
      { printf("ERROR: realloc failure\n"); return(-1); }
      strcpy(surf[snew].eparm,surf[smaster].eparm);
    }
    surface_i( surf[smaster].name, surf[smaster].ori, surf[smaster].sh, nls[0], os[0], ls[0], typs[0] );
  }

  /* correct the surface orientation if different to initial orientation of smaster */
  for(k=0; k<2; k++)
  {
    if(k) s=snew;
    else s=smaster;
    //printf("ori surf:%d %s\n",k,surf[s].name);
    for(i=0; i<surf[s].nl; i++)
    {
      if(surf[s].typ[i]=='l')
      {
	if(lbuf[k]==surf[s].l[i])
	{
	  if(loribuf[k]!=surf[s].o[i])
	    flipSurfori(s);
	}
      }
      else
      {
        for(n=0; n<lcmb[surf[s].l[i]].nl; n++)
        {
	  if(lbuf[k]==lcmb[surf[s].l[i]].l[n])
	  {
	    if(loribuf[k]!=lcmb[surf[s].l[i]].o[n])
	      flipSurfori(s);
	  }
	}
      }
    }
  }

  // TBD:
  for(k=0; k<nsurfs; k++)
  {
    s=innerSurf[k];
    //printf("ori surf:%d %s\n",k,surf[s].name);
    for(i=0; i<surf[s].nl; i++)
    {
      if(surf[s].typ[i]=='l')
      {
        for(l=0; l<nl; l++)
        {
  	  if(ptr_lmaster[l]==surf[s].l[i])
	  {
	    if(surf[s].o[i]!='+')
	      flipSurfori(s);
	  }
	}
      }
    }
  }
      
  /* add snew to bodies */
  for(i=0; i<anzGeo->b; i++) if( body[i].name != (char *)NULL )
  {
    for(j=0; j<body[i].ns; j++)
    {
      if(body[i].s[j]==smaster) { n=j; break; }
    }
    if(j<body[i].ns)
    {
      /* body must be redefined */
      if ((body[i].s = (int *)realloc(body[i].s, (body[i].ns+1)*sizeof(int)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      if ((body[i].o = (char *)realloc(body[i].o, (body[i].ns+1)*sizeof(char)) ) == NULL )
      { printf("ERROR: malloc failure\n"); return(-1); }
      body[i].o[body[i].ns]=body[i].o[n];
      body[i].s[body[i].ns]=snew;
      body[i].ns++;
    }
  }

  /* add snew to sets which contain smaster */
  for(i=0; i<anz->sets; i++)
  {
    if(!set[i].type)
    {
      if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&(i!=setall)&&( getIndex(&set[i].surf,set[i].anz_s,smaster) >-1))
      {
	seta(i,"s",snew);
        for(k=0; k<nsurfs; k++) seta(i,"s",innerSurf[k]);
      }
    }
  }

  if(repFlag)
  {
    repSurf(smaster,1);
    repSurf(snew,1);
    for(k=0; k<nsurfs; k++) repSurf(innerSurf[k],1);
  }
  return(snew);
}


/*
return(surf)
if error
return
-1: malloc failed
-2: no 2 surfs selected
-3: the master surf is of type "BLEND"
-4: trim failed
*/
int combineSurfs( int setNr )
{
  int i,ii,j,jj,l,n,s;
  int smaster=-1,sslave, snew, skipFlag;
  static int *lines=NULL;
  static char *ori=NULL;
  static char *typ=NULL;
  double lmax=-MAX_FLOAT, Ll=0.;
  char name[MAX_LINE_LENGTH];
  int  sslaveEvaluated=0;

  /* determine the master surf (the "greater" one) by measure of their line length */
  /* it is supposed that the related nurbs (if any) covers also the smaller surf */
  /* It was successfull if a trimming attempt of the combined surf verifies that */

  /* for the moment just two surfs can be combined */
  if(set[setNr].anz_s!=2) return(-2);

  /* length of the lines of the surfs */
  for (s=0; s<set[setNr].anz_s; s++)
  {
    jj=set[setNr].surf[s];
    if( surf[jj].name != (char *)NULL )
    Ll=0.;
    for (n=0; n<surf[jj].nl; n++)
    {
      if( surf[jj].typ[n] == 'l' )
      {
        Ll+=calcLineLength(surf[jj].l[n]);
      }
      if( surf[jj].typ[n] == 'c' )
      {
        for(l=0; l<lcmb[surf[jj].l[n]].nl; l++)
	{
          Ll+=calcLineLength(lcmb[surf[jj].l[n]].l[l]);
	}
      }
    }
    //printf("surf:%s LENGTH:%lf\n", surf[jj].name, Ll);
    if(Ll>lmax) { lmax=Ll; smaster=s; }
  }
  //printf("msurf:%s LENGTH:%lf sh:%d\n", surf[set[setNr].surf[smaster]].name, lmax, surf[set[setNr].surf[smaster]].sh);

  /* generate a new surface */
  /* for the moment just two surfs can be combined */
  if(smaster==0) sslave=1; else sslave=0;

  smaster=set[setNr].surf[smaster];
  sslave=set[setNr].surf[sslave];

  /* the master surf must not be of type "BLEND" */
  if(surf[smaster].sh<0) return(-3);

  l=0;
  for(i=0; i<surf[smaster].nl; i++)
  {
    skipFlag=0;
    for(j=0; j<surf[sslave].nl; j++)
    {
      // might be necessary to regard lcmb cases later
      if(surf[smaster].l[i]==surf[sslave].l[j]) skipFlag=1;
     
      // if i==1: The lines from the second surf have to follow now since orient() must find now the lines from 2nd surf
      if((i==1)&&(skipFlag))
      {
        // code copy from below if(sslaveEvaluated) 
        for(ii=0; ii<surf[sslave].nl; ii++)
        {
          skipFlag=0;
          for(jj=0; jj<surf[smaster].nl; jj++)
          {
            if(surf[smaster].l[jj]==surf[sslave].l[ii]) { skipFlag=1; break; }
          }
          if(!skipFlag)
          {
            if(surf[sslave].typ[ii]=='l')
            {
              if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
      	      { printf("ERROR: malloc failure\n"); return(-1); }
              if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
      	      { printf("ERROR: malloc failure\n"); return(-1); }
              if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
      	      { printf("ERROR: malloc failure\n"); return(-1); }
              typ[l]='l';
              ori[l]=surf[sslave].o[ii];
              lines[l++]=surf[sslave].l[ii];
            }
            else
            {
              for (n=0; n<lcmb[surf[sslave].l[ii]].nl; n++)
              {
                if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
      	        { printf("ERROR: malloc failure\n"); return(-1); }
                if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
      	        { printf("ERROR: malloc failure\n"); return(-1); }
                if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
      	        { printf("ERROR: malloc failure\n"); return(-1); }
                typ[l]='l';
                ori[l]=lcmb[surf[sslave].l[ii]].o[n];
                if(surf[sslave].o[ii]=='-') { if(ori[l]=='+') ori[l]='-'; else ori[l]='-'; }
                lines[l++]=lcmb[surf[sslave].l[ii]].l[n];
              }
            }
          }
        }
        skipFlag=sslaveEvaluated=1;
        break;
      }
    }
    if(!skipFlag)
    {
      if(surf[smaster].typ[i]=='l')
      {
        if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
	{ printf("ERROR: malloc failure\n"); return(-1); }
        if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
	{ printf("ERROR: malloc failure\n"); return(-1); }
        if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
	{ printf("ERROR: malloc failure\n"); return(-1); }
        typ[l]='l';
        ori[l]=surf[smaster].o[i];
	//printf("l %d add %s\n", l, line[surf[smaster].l[i]].name);
        lines[l++]=surf[smaster].l[i];
      }
      else
      {
        for (n=0; n<lcmb[surf[sslave].l[i]].nl; n++)
        {
          if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
	  { printf("ERROR: malloc failure\n"); return(-1); }
          if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
	  { printf("ERROR: malloc failure\n"); return(-1); }
          if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
	  { printf("ERROR: malloc failure\n"); return(-1); }
          typ[l]='l';
          ori[l]=lcmb[surf[smaster].l[i]].o[n];
          if(surf[smaster].o[i]=='-') { if(ori[l]=='+') ori[l]='-'; else ori[l]='-'; }
	  //printf("c %d add %s\n",l, line[lcmb[surf[smaster].l[i]].l[n]].name);
          lines[l++]=lcmb[surf[smaster].l[i]].l[n];
	}
      }
    }
  }
  if(!sslaveEvaluated)
  for(i=0; i<surf[sslave].nl; i++)
  {
    skipFlag=0;
    for(j=0; j<surf[smaster].nl; j++)
    {
      if(surf[smaster].l[j]==surf[sslave].l[i]) skipFlag=1;
    }
    if(!skipFlag)
    {
      if(surf[sslave].typ[i]=='l')
      {
        if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
	{ printf("ERROR: malloc failure\n"); return(-1); }
        if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
	{ printf("ERROR: malloc failure\n"); return(-1); }
        if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
	{ printf("ERROR: malloc failure\n"); return(-1); }
        typ[l]='l';
        ori[l]=surf[sslave].o[i];
        lines[l++]=surf[sslave].l[i];
      }
      else
      {
        for (n=0; n<lcmb[surf[sslave].l[i]].nl; n++)
        {
          if ((lines = (int *)realloc(lines, (l+1)*sizeof(int)) ) == NULL )
	  { printf("ERROR: malloc failure\n"); return(-1); }
          if ((ori = (char *)realloc(ori, (l+1)*sizeof(char)) ) == NULL )
	  { printf("ERROR: malloc failure\n"); return(-1); }
          if ((typ = (char *)realloc(typ, (l+1)*sizeof(char)) ) == NULL )
	  { printf("ERROR: malloc failure\n"); return(-1); }
          typ[l]='l';
          ori[l]=lcmb[surf[sslave].l[i]].o[n];
          if(surf[sslave].o[i]=='-') { if(ori[l]=='+') ori[l]='-'; else ori[l]='-'; }
          lines[l++]=lcmb[surf[sslave].l[i]].l[n];
	}
      }
    }
  }

  getNewName( name, "s" );
  snew= surface_i( name, surf[smaster].ori, surf[smaster].sh, l, ori, lines, typ );
  surf[snew].etyp=surf[smaster].etyp;
  surf[snew].eattr=surf[smaster].eattr;

  printf(" surfaces combined\n");

  /* update the master surf */
  surface_i( surf[smaster].name, surf[smaster].ori, surf[smaster].sh, surf[snew].nl, surf[snew].o, surf[snew].l, surf[snew].typ );
  surf[smaster].etyp=surf[snew].etyp;
  surf[smaster].eattr=surf[snew].eattr;

  delSurf( 1, &snew );

  /* remove sslave from bodies */
  for(i=0; i<anzGeo->b; i++) if( body[i].name != (char *)NULL )
  {
    for(j=0; j<body[i].ns; j++)
    {
      if(body[i].s[j]==sslave) break;
    }
    if(j<body[i].ns)
    {
      /* body must be redefined */
      n=0;
      for(j=0; j<body[i].ns; j++)
      {
        if(body[i].s[j]!=sslave)
	{
          body[i].o[n]=body[i].o[j];
          body[i].s[n++]=body[i].s[j];
	}
      }
      body[i].ns--;
    }
  }
  delSurf( 1, &sslave );

  return(smaster);
}


/* makes a seqence-line from an lcmb */
int convertLCMB( int c )
{
  int i,j,k,l,n,p, setNr, div=0;
  double  bias;
  char name[MAX_LINE_LENGTH];

  if(lcmb[c].name==(char *)NULL) return(-1); 

  /* define a sequence */
  if( getNewName( name, "se" ) == -1 ) 
  {
    printf("ERROR: set could not be created\n");
    return(-1);
  }
  if( (setNr=pre_seta( name, "is", 0)) <0 ) return(-1);
  bias=1;

  /* create spline points */
  /* and define a sequence */
  for (i=0; i<lcmb[c].nl; i++) 
  {
    l=lcmb[c].l[i];
    repLine(l);
    div+=line[l].div;
    if(lcmb[c].o[i]=='+')
    {
      if(i==0) seta( setNr, "p", line[l].p1 ); 
      for (n=3; n<line[l].nip-3; n+=3)
      {
        if( getNewName( name, "p" ) == -1 ) 
        {
          printf("ERROR: point could not be created\n");
          return(-1);
        }
        p=pnt( name, line[l].ip[n], line[l].ip[n+1], line[l].ip[n+2], 0);
        if( p<0) printf("ERROR: point could not be created\n");
        seta(setNr,"p",p);
      }
      seta( setNr, "p", line[l].p2 ); 
    }
    else
    {
      if(i==0) seta( setNr, "p", line[l].p2 ); 
      for (n=line[l].nip-6; n>=3; n-=3)
      {
        if( getNewName( name, "p" ) == -1 ) 
        {
          printf("ERROR: point could not be created\n");
          return(-1);
        }
        p=pnt( name, line[l].ip[n], line[l].ip[n+1], line[l].ip[n+2], 0);
        if( p<0) printf("ERROR: point could not be created\n");
        seta(setNr,"p",p);
      }
      seta( setNr, "p", line[l].p1 ); 
    }
  }

  if( getNewName( name, "l" ) == -1 ) 
  {
    printf("ERROR: line could not be created\n");
    return(-1);
  }
  l=line_i( name, lcmb[c].p1, lcmb[c].p2, setNr, div, bias, 's' );

  /* replace the lcmb by the new line in all surfaces */
  p=0;
  for (i=0; i<anzGeo->s; i++) if( surf[i].name != (char *)NULL )
  {
    for (j=0; j<surf[i].nl; j++)
    {
      if(( c == surf[i].l[j] )&&( surf[i].typ[j]=='c' ))
      {
        p=1;
        printf ("realloc surf:%s and replace lcmb:%s with line:%s\n",
          surf[i].name, lcmb[c].name, line[l].name );
        surf[i].l[j]=l;
        surf[i].typ[j]='l';
      }
    }
  }
/* if the lcmb was not part of the surface then replace the single lines by the new line */
  if(!p)
  {
    for (i=0; i<anzGeo->s; i++) if( surf[i].name != (char *)NULL )
    {
      p=0;
      for (j=0; j<surf[i].nl; j++)
      {
        for (k=0; k<lcmb[c].nl; k++)
        {
          if(( lcmb[c].l[k] == surf[i].l[j] )&&( surf[i].typ[j]=='l' ))
	  {
            if(!p)
            {
              p=1;
              surf[i].l[j]=l;
              surf[i].typ[j]='l';
            }
            else
	    {
              /* just remove the line from the surf */
              //printf ("remove line:%s \n", line[surf[i].l[j]].name );
              for (n=j; n<surf[i].nl-1; n++)
              {
                //printf (" replace line:%s with line:%s\n", line[surf[i].l[n]].name, line[surf[i].1[n+1]].name );
                surf[i].l[n]=surf[i].l[n+1];
                surf[i].o[n]=surf[i].o[n+1];
                surf[i].typ[n]=surf[i].typ[n+1];
	      }
              surf[i].nl--;
              j--;
            }
          }
	}
      }
    }
  }

  /* replace the lcmb by the new line in all sets */
  for (i=1; i<anz->sets; i++)
  {
    if((set[i].name!=(char *)NULL)&&(set[i].name[0]!='-')&&( getIndex(&set[i].lcmb,set[i].anz_c, c) >-1))
      seta( i, "l", l ); 
  }

  delLine( lcmb[c].nl, lcmb[c].l );
  delLcmb( 1, &c );
  updateDispLists();
  return(l);
}


/* makes a seqence-line from an arc or straight line */
void convertLine( int l, int div )
{
  int k, p, setNr;
  double pn[3], bias;
  char name[MAX_LINE_LENGTH];

  /* define a sequence */
  if( getNewName( name, "se" ) == -1 ) 
  {
    printf("ERROR: point could not be created\n");
    return;
  }
  if( (setNr=pre_seta( name, "is", 0)) <0 ) return;
  seta( setNr, "p", line[l].p1 ); 
  bias=line[l].bias;
  line[l].bias=1.;

  /* create spline points */
  /* and define a sequence */
  for (k=0; k<div; k++)
  {
    if (line[l].typ=='a') arcNodes( l, k, div+1, pn );
    else if (line[l].typ=='n') nurlNodes( l, k, div+1, pn );
    else if (line[l].typ=='s') splineNodes( l, k, div+1, pn );
    else straightNodes( l, k, div+1, pn );
    if( getNewName( name, "p" ) == -1 ) 
    {
      printf("ERROR: point could not be created\n");
      return;
    }
    if(printFlag) printf (" pnt=%s x=%lf y=%lf z=%lf\n",  name, pn[0], pn[1], pn[2]);
    p=pnt( name, pn[0], pn[1], pn[2], 0 );
    if( p<0) printf("ERROR: point could not be created\n");
    seta( setNr, "p", p ); 
  }
  seta( setNr, "p", line[l].p2 );
  line_i( line[l].name, line[l].p1, line[l].p2, setNr, line[l].div, bias, 's' );
  return;
}



double AsplitL_intersect( double *b, double *eu, double *ev, double *eg, double *en )
/******************************************************************/
/*   Grade (eg) schneidet Ebene (eu,ev) return g                  */
/*   determinante dritter Ordnung                                 */
/*   b= Abstand zwischen den Aufpunkten der Linie und Ebene       */
/*   b=eu*u + ev*v + eg *g  (e: Einheitsvektoren )                */
/* return MAX_FLOAT if eg parallel to euXev                       */
/******************************************************************/
{
  double g, D, Dg, a, c;

  a = eu[0]*ev[1]*eg[2]+ ev[0]*eg[1]*eu[2]+ eg[0]*eu[1]*ev[2];
  c = eg[0]*ev[1]*eu[2]+ eu[0]*eg[1]*ev[2]+ ev[0]*eu[1]*eg[2];
  D = a-c;

  a = eu[0]*ev[1]* b[2]+ ev[0]* b[1]*eu[2]+  b[0]*eu[1]*ev[2];
  c =  b[0]*ev[1]*eu[2]+ eu[0]* b[1]*ev[2]+ ev[0]*eu[1]* b[2];
  Dg= a-c;
  g = Dg / D;
  return (g);
}




/* creates the intersection point of two straight lines */
/* not impl:(returns the mismatch in the intersection) */
int intersectionPoint( double *l1p1, double *l1p2, double *l2p1, double *l2p2, double *ps)
{
  double l1_[3], l2_[3], ln_[3], en[3], eg[3], l2l1_p1[3], l2ps[3], l1ps[3];
  double eu[3], ev[3], g;
  double ps1[3], ps2[3], dist, dx,dy,dz;

  /* calc the intersection point 2 times */
  /* first in the direction of line 1 then in the direction of line 2 */
  /* if the difference between the two intersection points is bigger than gtol return(0) else (1) */

  /* ps on vector l2 */
    /* berechne den zweiten Einheitsvektor der Schnittebene, ev==Peilstrahl x Linie */
    v_result( l2p1, l2p2, l1_ );
    v_norm( l1_, eu );
    v_result( l1p1, l1p2, l2_ );
    v_prod( l1_, l2_, ln_ );
    v_norm( ln_, ev );
  
    /* berechne den Normalenvektor der Schnittebene */
    v_prod( eu, ev, en );
    v_norm( en, en );
  
    /* berechne den Einheitsvektor der zu splitenden Linie */
    v_norm( l2_, eg );
  
    /* bestimme den Abstand zwischen den Aufpunkten der Linie und Ebene  */
    v_result( l1p1, l2p1, l2l1_p1 );
  
    /* berechne die Konstante g zur berechnung von ps (Schnittpunkt) ps=p0+eg*g  */
    g = AsplitL_intersect( l2l1_p1, eu, ev, eg, en );
  
    /* erzeuge den Schnittpunkt */
    v_scal( &g, eg, l1ps );
    v_add( l1p1, l1ps, ps1 );

  /* ps on vector l2 */
    /* berechne den zweiten Einheitsvektor der Schnittebene, ev==Peilstrahl x Linie */
    v_result( l1p1, l1p2, l1_ );
    v_norm( l1_, eu );
    v_result( l2p1, l2p2, l2_ );
    v_prod( l1_, l2_, ln_ );
    v_norm( ln_, ev );
  
    /* berechne den Normalenvektor der Schnittebene */
    v_prod( eu, ev, en );
    v_norm( en, en );
  
    /* berechne den Einheitsvektor der zu splitenden Linie */
    v_norm( l2_, eg );
  
    /* bestimme den Abstand zwischen den Aufpunkten der Linie und Ebene  */
    v_result( l2p1, l1p1, l2l1_p1 );
  
    /* berechne die Konstante g zur berechnung von ps (Schnittpunkt) ps=p0+eg*g  */
    g = AsplitL_intersect( l2l1_p1, eu, ev, eg, en );
  
    /* erzeuge den Schnittpunkt */
    v_scal( &g, eg, l2ps );
    v_add( l2p1, l2ps, ps2 );

    /* calc the distance   */
    dx=ps1[0]-ps2[0];
    dy=ps1[1]-ps2[1];
    dz=ps1[2]-ps2[2];
    dist=dx*dx+dy*dy+dz*dz;
    //printf("dist:%lf %lf\n", sqrt(dist), sqrt(dist)*scale->w);
    if(dist<gtol/scale->w)
    {
      ps[0]=ps1[0]-dx*.5;
      ps[1]=ps1[1]-dy*.5;
      ps[2]=ps1[2]-dz*.5;
      return(1);
    }
  return(0);
}


/* return 1 or -1 if failed */
int createFilletCenterPoint( double filletRadius, double *l1p1, double *l1p2, double *l2p1, double *l2p2, double *p1, double *p2, double *ps)
{
  double l1[3], l2[3], p1ps[3], p2ps[3];
  double eu[3], ev[3], en[3], ln[3];
  double l1p1b[3], l1p2b[3], l2p1b[3], l2p2b[3];
  double b[3],eg[3],g, lg[3]; 

  /* at first create the intersection if any */
  /* then create two new artificial lines running away from the intersection */

  if(intersectionPoint( l1p1, l1p2, l2p1, l2p2, ps) == 0) return(-1);

  /* two parallel lines are created. Offset is the radius. Were this lines intersect is the centerpnt ps */
  /* From the ps to the lines are the endpoints of the arc. (p1 p2) */

  v_result(l1p1,l1p2, l1);
  v_result(l2p1,l2p2, l2);

  /* determine the radius-vectors */
  v_prod( l1, l2, ln );
  v_prod( ln, l1, ev );
  v_norm( ev, eu );
  v_scal( &filletRadius, eu, p1ps );

  v_prod( l2, ln, ev );
  v_norm( ev, eu );
  v_scal( &filletRadius, eu, p2ps );
  
  /* new offseted lines which intersect to give ps */
  v_add(l1p1, p1ps, l1p1b);
  v_add(l1p2, p1ps, l1p2b);
  v_add(l2p1, p2ps, l2p1b);
  v_add(l2p2, p2ps, l2p2b);

  if(intersectionPoint( l1p1b, l1p2b, l2p1b, l2p2b, ps)==0) return(-1);

  /* create the endpoints of the arc on line l1 and l2 */

  v_norm( ln, ev );

  v_result(l1p1, ps, b);
  v_norm(p1ps , eu );
  v_norm( l1, eg );
  v_prod(eu,ev,en);
  v_norm( en, en );
  g= AsplitL( b, eu, ev, eg, en );
  v_scal( &g, eg, lg );
  v_add(l1p1, lg, p1); 

  v_result(l2p1, ps, b);
  v_norm(p2ps , eu );
  v_norm( l2, eg );
  v_prod(eu,ev,en);
  v_norm( en, en );
  g= AsplitL( b, eu, ev, eg, en );
  v_scal( &g, eg, lg );
  v_add(l2p1, lg, p2); 

  //for (i=0; i<3; i++) printf("p1:%f p2:%f\n", p1[i], p2[i]);
  //v_result(p1ps, ps, p1); 
  //v_result(p2ps, ps, p2); 

  return(1);
}



int createFillet(int lin, double filletRadius)
{
  int i,j,n,p;
  static int linbuf[2];
  int dirl[2], p1nr, p2nr, psnr, l, l_nr[2], flipflop;
  double p1[3], p2[3], pint[2][3], ps[3], p1p2[3], p1ps[3];
  double lp1p2[2], lp1ps[2], lp1p2_buf[2], lp1ps_buf[2], dv[4][3], dist, min=MAX_FLOAT;
  double lps[2], lbez;
  int lp_ptr[2], lp_ptr_buf[2], lflag[2], icase=0, trkNew, trk, divbuffer[2];

  char name[MAX_LINE_LENGTH];

  typedef struct{
    int np;
    double **pnt;
  }Ltmp;
  Ltmp ltmp[2];

  if((line[lin].typ==' ')||( line[lin].typ=='s'))
  {
    if((intersectFlag)&&( lin==linbuf[0]))
    {
      printf("ERROR: second line:%s not selected because its the same as the first selected.\n", line[lin].name);
      return(-1);
    }
    linbuf[intersectFlag]=lin;
  }
  else
  {
    printf("ERROR: line:%s not selected. Only straight lines and splines are supported\n", line[lin].name);
    return(-1);
  }

  if((intersectFlag)&&( lin==linbuf[0]))
  {
    printf("ERROR: second line:%s not selected because its the same as the first selected.\n", line[lin].name);
    return(-1);
  }
  linbuf[intersectFlag]=lin;

  if(intersectFlag)
  {
    /* determine the closest combi of line-end-points determine the orientation of the lines relative to each other */
    v_result( &point[line[linbuf[0]].p1].px, &point[line[linbuf[1]].p1].px, dv[0] );
    v_result( &point[line[linbuf[0]].p2].px, &point[line[linbuf[1]].p2].px, dv[1] );
    v_result( &point[line[linbuf[0]].p1].px, &point[line[linbuf[1]].p2].px, dv[2] );
    v_result( &point[line[linbuf[0]].p2].px, &point[line[linbuf[1]].p1].px, dv[3] );
    for(i=0; i<4; i++)
    {
      dist=v_betrag(dv[i]);
      if(min>dist) { min=dist; icase=i; }
    }
    switch(icase)
    {
      case 0:
	dirl[0]=1;
	dirl[1]=1;
	break;
      case 1:
	dirl[0]=2;
	dirl[1]=2;
	break;
      case 2:
	dirl[0]=1;
	dirl[1]=2;
	break;
      case 3:
	dirl[0]=2;
	dirl[1]=1;
	break;
    }
 
    /* go over all divisions of the lines and determine the split-points */
    /* create two arrays ltmp[] of points based on the lines regarding their relative orientation */

    for(i=0; i<2; i++)
    {
      if(dirl[i]==1)  // line.p1 is startpnt
      {
        if(line[linbuf[i]].typ==' ')
	{
          ltmp[i].np=2;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          ltmp[i].pnt[0][0]=point[line[linbuf[i]].p1].px; 
          ltmp[i].pnt[0][1]=point[line[linbuf[i]].p1].py; 
          ltmp[i].pnt[0][2]=point[line[linbuf[i]].p1].pz; 
          ltmp[i].pnt[1][0]=point[line[linbuf[i]].p2].px; 
          ltmp[i].pnt[1][1]=point[line[linbuf[i]].p2].py; 
          ltmp[i].pnt[1][2]=point[line[linbuf[i]].p2].pz; 
	}
	else
	{
          /* temporary increase the division */
	  divbuffer[i]=line[linbuf[i]].div;
          line[linbuf[i]].div=98;
          repLine(linbuf[i]);

          ltmp[i].np=line[linbuf[i]].nip/3;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          p=0;
          for(j=0; j<line[linbuf[i]].nip; j+=3)
          {
            ltmp[i].pnt[p][0]=line[linbuf[i]].ip[j];
            ltmp[i].pnt[p][1]=line[linbuf[i]].ip[j+1];
            ltmp[i].pnt[p][2]=line[linbuf[i]].ip[j+2];
            p++;
	  }
	}
      }
      else // line.p2 is startpnt
      {
        if(line[linbuf[i]].typ==' ')
	{
          ltmp[i].np=2;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          ltmp[i].pnt[0][0]=point[line[linbuf[i]].p2].px; 
          ltmp[i].pnt[0][1]=point[line[linbuf[i]].p2].py; 
          ltmp[i].pnt[0][2]=point[line[linbuf[i]].p2].pz; 
          ltmp[i].pnt[1][0]=point[line[linbuf[i]].p1].px; 
          ltmp[i].pnt[1][1]=point[line[linbuf[i]].p1].py; 
          ltmp[i].pnt[1][2]=point[line[linbuf[i]].p1].pz; 
	}
	else
	{
          /* temporary increase the division */
	  divbuffer[i]=line[linbuf[i]].div;
          line[linbuf[i]].div=98;
          repLine(linbuf[i]);

          ltmp[i].np=line[linbuf[i]].nip/3;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          p=0;
          for(j=line[linbuf[i]].nip-3; j>0; j-=3)
          {
            ltmp[i].pnt[p][0]=line[linbuf[i]].ip[j];
            ltmp[i].pnt[p][1]=line[linbuf[i]].ip[j+1];
            ltmp[i].pnt[p][2]=line[linbuf[i]].ip[j+2];
            p++;
	  }
	}
      }
    }

    lflag[0]=1;
    lflag[1]=1;
    /* these pnt-indexes have to be changed due to the results of "searchIntersectionPoints()" */
    lp_ptr_buf[0]=lp_ptr[0]=0;
    lp_ptr_buf[1]=lp_ptr[1]=0;
    flipflop=0;
    do
    {
      if(createFilletCenterPoint( filletRadius, ltmp[0].pnt[lp_ptr[0]], ltmp[0].pnt[lp_ptr[0]+1], ltmp[1].pnt[lp_ptr[1]], ltmp[1].pnt[lp_ptr[1]+1], pint[0],pint[1],ps)==-1)
      {
        /* failed! free ltmp */
        for(i=0; i<2; i++)
        {
          /* restore the division */
          if(line[linbuf[i]].typ!=' ')
	  {
            line[linbuf[i]].div=divbuffer[i];
            repLine(linbuf[i]);
	  }
          for(j=0; j<ltmp[i].np; j++) free(ltmp[i].pnt[j]); free(ltmp[i].pnt);
        }
        printf(" ERROR: lines %s %s do not intersect based on a tolerance of:%lf\n", line[linbuf[0]].name, line[linbuf[1]].name, gtol);
        return(-1);
      } 

      /* check if the intersection of the lines (p1,p2) is outside of the interval (ie. l1ptr, l1ptr+1) */
      for(i=0; i<2; i++)
      {
        v_result(ltmp[i].pnt[lp_ptr[i]], ltmp[i].pnt[lp_ptr[i]+1], p1p2);
        v_result(ltmp[i].pnt[lp_ptr[i]], pint[i], p1ps);
        lp1p2[i]=v_betrag(p1p2);
        lp1ps[i]=v_betrag(p1ps);
        if(v_sprod(p1ps,p1p2)<0.) lp1ps[i]*=-1;

        /* then adjust the interval */

        /* if lp1ps is <0 then decrease lp_ptr[i] */
        /* warning gtol might be set to a very high value to enable intersection, so better not use this */
        //if((lp1ps[i]<(-gtol/scale->w))&&(lp_ptr[i]>0)) lp_ptr[i]--;
        if((lp1ps[i]<(-1e-6))&&(lp_ptr[i]>0)) lp_ptr[i]--;
        /* if lp1ps is >lp1p2 then increase lp_ptr[i] */
        //else if((lp1ps[i]>(lp1p2[i]+gtol/scale->w))&&(lp_ptr[i]<ltmp[i].np-2)) lp_ptr[i]++;
        else if((lp1ps[i]>(lp1p2[i]+1e-6))&&(lp_ptr[i]<ltmp[i].np-2)) lp_ptr[i]++;
        /* if lp1ps is <=lp1p2 then break */
        else lflag[i]=0;
      }
      if(printFlag) printf("lflag %d %d\n", lflag[0],lflag[1]);

      /* escape endless loops */
      if(flipflop)
      {
        if((lp_ptr[0]==lp_ptr_buf[0])&&(lp_ptr[1]==lp_ptr_buf[1]))
        {
          for(i=0; i<2; i++)
          {
	    lp_ptr[i]=lp_ptr_buf[i];
	    lp1ps[i]=lp1ps_buf[i];
	    lp1p2[i]=lp1p2_buf[i];
	  }
          break;
        }
        for(i=0; i<2; i++)
        {
          lp_ptr_buf[i]=lp_ptr[i];
	  lp1ps_buf[i]=lp1ps[i];
	  lp1p2_buf[i]=lp1p2[i];
	}
      }
      flipflop=!flipflop;
    }while( lflag[0] || lflag[1] );

    /* determine the line-length up to the split point, necessary for line-splitting */
    for(i=0; i<2; i++)
    {
      /* check if the split-point is outside */
      if(printFlag) printf("check line[%d]:%s of type:%c\n", i,line[linbuf[i]].name, line[linbuf[i]].typ);
      if(line[linbuf[i]].typ!='s') continue;

      /* restore the division */
      line[linbuf[i]].div=divbuffer[i];
      repLine(linbuf[i]);

      if(printFlag) printf("redefine the trkset of line[%d]:%s\n", i,line[linbuf[i]].name);

      /* determine the line-length up to the splitpoint */
      lps[i]=0.;
      for(j=0; j<lp_ptr[i]; j++)
      {
        v_result(ltmp[i].pnt[j], ltmp[i].pnt[j+1], p1p2);
        lps[i]+=v_betrag(p1p2);
      }
      lps[i]+=lp1ps[i];

      /* redefine the track-set */
      lbez=0.;
      getNewName( name, "se" );
      trkNew=pre_seta(name,"is",0);
      trk=line[linbuf[i]].trk;
      if(dirl[i]==1) /* remove the leading points and add the split-point */
      {
        seta( trkNew, "ps", set[trk].pnt[0] );  // will be redefined later
        for (j=0; j<set[trk].anz_p-1; j++)
        {
          v_result( &point[set[trk].pnt[j]].px, &point[set[trk].pnt[j+1]].px, p1p2 );
          lbez+=v_betrag( p1p2 );
          if(lbez>lps[i]) seta( trkNew, "ps", set[trk].pnt[j+1] );
        }
      }
      if(dirl[i]==2) /* remove the trailing points and add the split-point */
      {
        lps[i]=calcLineLength(linbuf[i])/scale->w-lps[i];
        seta( trkNew, "ps", set[trk].pnt[0] );
        for (j=0; j<set[trk].anz_p-1; j++)
        {
          v_result( &point[set[trk].pnt[j]].px, &point[set[trk].pnt[j+1]].px, p1p2 );
          lbez+=v_betrag( p1p2 );
	  //printf("p %s %s lbez:%f lps[i]:%f\n", point[set[trk].pnt[j]].name,point[set[trk].pnt[j+1]].name, lbez, lps[i]);
          if(lbez<lps[i]) seta( trkNew, "ps", set[trk].pnt[j+1] );
          else break;
        }
        seta( trkNew, "ps", set[trk].pnt[set[trk].anz_p-1] );  // will be redefined later
      }
      line[linbuf[i]].trk=trkNew;
      delSet(set[trk].name);
    }

    /* free ltmp */
    for(i=0; i<2; i++) { for(j=0; j<ltmp[i].np; j++) free(ltmp[i].pnt[j]); free(ltmp[i].pnt); }

    /* create center point */
    p1[0]=pint[0][0];
    p1[1]=pint[0][1];
    p1[2]=pint[0][2];
    p2[0]=pint[1][0];
    p2[1]=pint[1][1];
    p2[2]=pint[1][2];
    getNewName( name, "p" );
    if(printFlag) printf(" create center point:%s %lf %lf %lf\n", name, ps[0], ps[1], ps[2] );
    psnr=pnt( name, ps[0], ps[1], ps[2], 0 );

    /* move one point of line 1 to p1 and one of line 2 to p2 */
    /* but if line 1&2 use a common point then replace the point from line 2 with a new point */

    if(dirl[0]==2) /* line points to the intersection */
    {
      p1nr=line[linbuf[0]].p2;
      point[line[linbuf[0]].p2].px=p1[0];
      point[line[linbuf[0]].p2].py=p1[1];
      point[line[linbuf[0]].p2].pz=p1[2];
    }
    else 
    {
      p1nr=line[linbuf[0]].p1;
      point[line[linbuf[0]].p1].px=p1[0];
      point[line[linbuf[0]].p1].py=p1[1];
      point[line[linbuf[0]].p1].pz=p1[2];
    }
    if(dirl[1]==2) /* line points to the intersection */
    {
      if(p1nr==line[linbuf[1]].p2)
      {
        n= getNewName( name, "p" );
        printf(" create point:%s %lf %lf %lf\n", name, p2[0], p2[1], p2[2] );
        p2nr=pnt( name, p2[0], p2[1], p2[2], 0 );
        line[linbuf[1]].p2=p2nr;
        if(line[linbuf[1]].typ=='s') set[line[linbuf[1]].trk].pnt[set[line[linbuf[1]].trk].anz_p-1]=p2nr;
      }
      else
      {
        p2nr=line[linbuf[1]].p2;
        point[line[linbuf[1]].p2].px=p2[0];
        point[line[linbuf[1]].p2].py=p2[1];
        point[line[linbuf[1]].p2].pz=p2[2];
      }
    }
    else 
    {
      if(p1nr==line[linbuf[1]].p1)
      {
        n= getNewName( name, "p" );
        printf(" create point:%s %lf %lf %lf\n", name, p2[0], p2[1], p2[2] );
        p2nr=pnt( name, p2[0], p2[1], p2[2], 0 );
        line[linbuf[1]].p1=p2nr;
        if(line[linbuf[1]].typ=='s') set[line[linbuf[1]].trk].pnt[0]=p2nr;
      }
      else
      {
        p2nr=line[linbuf[1]].p1;
        point[line[linbuf[1]].p1].px=p2[0];
        point[line[linbuf[1]].p1].py=p2[1];
        point[line[linbuf[1]].p1].pz=p2[2];
      }
    }
    
    /* create the arc */
    l= getNewName( name, "l" );
    if ( l == -1 )
    { printf("copy: could not create new line\n"); return(-1); }
    l_nr[1]=line_i( name, p1nr, p2nr, psnr, 0, 1, 'a' );

    /* add the arc to all higher entities which use line 2 */
    /* code from qsplitLine() */
    l_nr[0]=l=linbuf[1];

    /* untersuche alle lcmbs ob linbuf ein Mitglied ist */
    for (i=0; i<anzGeo->c; i++) if( lcmb[i].name != (char *)NULL )
    {
      for (j=0; j<lcmb[i].nl; j++)
      {
        if( l == lcmb[i].l[j] )
        {
          printf (" realloc lcmb:%s and replace line:%s with %s and %s \n",
             lcmb[i].name, line[l].name,line[l_nr[0]].name, line[l_nr[1]].name );
          if ((lcmb[i].o = (char *)realloc( (char *)lcmb[i].o, (lcmb[i].nl+1)*sizeof(char)) ) == NULL )
          { printf("\n\n ERROR: realloc failure in qspl, lcmb.o:%s not changed\n\n",lcmb[i].name ); return(-1); }
          if ((lcmb[i].l = (int *)realloc( (int *)lcmb[i].l, (lcmb[i].nl+1)*sizeof(int)) ) == NULL )
          { printf("\n\n ERROR: realloc failure in qspl, lcmb.l:%s not changed\n\n", lcmb[i].name); return(-1); }
  	  /* umspeichern der linien beginnend bei der letzten bis einschlieslich j */
          for (n=lcmb[i].nl; n>j; n--)
  	  {
            lcmb[i].o[n]=lcmb[i].o[n-1];
            lcmb[i].l[n]=lcmb[i].l[n-1];
          }
          /* Auffuellen der j, j+1 pos. mit l1, l2 mit der Orientierung der gesplitteten linie */
          lcmb[i].o[j]=lcmb[i].o[j+1];
          lcmb[i].l[j]=l_nr[0];
          lcmb[i].o[j+1]=lcmb[i].o[j];
          lcmb[i].l[j+1]=l_nr[1];
          lcmb[i].nl++;
        }
      }
    }
  
    /* untersuche alle surfs ob linbuf ein Mitglied ist und ersetze sie durch eine lcmb */
    /* kontrolliere ob nicht schon eine geeignete lcmb existiert */
    for (i=0; i<anzGeo->s; i++) if( surf[i].name != (char *)NULL )
    {
      for (j=0; j<surf[i].nl; j++)
      {
        if(( l == surf[i].l[j] )&&( surf[i].typ[j]=='l' ))
        {
  
          /* do we already have a suitable lcmb? */
          for (n=0; n<anzGeo->c; n++ ) if( lcmb[n].name != (char *)NULL )
  	  {
            if (lcmb[n].nl==2)   /* same amount of lines */
  	    {
  	    /*
              printf ("checke lcmb:%s \n", lcmb[n].name);
  	    */
              if (((lcmb[n].l[0]==l_nr[0])||(lcmb[n].l[0]==l_nr[1])) && ((lcmb[n].l[1]==l_nr[0])||(lcmb[n].l[1]==l_nr[1])))
  	      {
  	      /*
                printf ("equal:%s\n",lcmb[n].name);
  	      */
                break;
  	      }
  	    }
          }
          if (n>=anzGeo->c)  /* no lcmb was found, so create one */
          {
            /* create lcmb */
            if ( getNewName( name, "c" ) == -1 )
            { printf("Type c not known, lcmb can not be created\n"); exit(-1); }
            lcmb_i( name, (int)0, (int)2, "++", l_nr );
            n=getLcmbNr( name );
          }
          printf ("realloc surf:%s and replace line:%s with lcmb:%s made of %s and %s \n",
            surf[i].name, line[l].name, name, line[l_nr[0]].name, line[l_nr[1]].name );
          if (n>-1) { surf[i].l[j]=n; surf[i].o[j]='+'; surf[i].typ[j]='c'; }
          else { errMsg("lcmb not known, surface could not be changed \n"); return(-1); }
        }
      }
    }
    updateDispLists();
  }

  intersectFlag=!intersectFlag;
  return(intersectFlag);
}


int intersect(int lin)
{
  int i,j,p;
  static int linbuf[2];
  int dirl[2], p1nr, p2nr, l, flipflop, breakflag=0;
  double ps[3], p1p2[3], p1ps[3];
  double lp1p2[2], lp1ps[2], lp1p2_buf[2], lp1ps_buf[2], dv[4][3], dist, min=MAX_FLOAT;
  double lps[2], lbez=0;
  int lp_ptr[2], lp_ptr_buf[2], lflag[2], icase=0, trkNew, trk, divbuffer[2];

  char name[MAX_LINE_LENGTH];

  typedef struct{
    int np;
    double **pnt;
  }Ltmp;
  Ltmp ltmp[2];

  if((line[lin].typ==' ')||( line[lin].typ=='s'))
  {
    if((intersectFlag)&&( lin==linbuf[0]))
    {
      printf("ERROR: second line:%s not selected because its the same as the first selected.\n", line[lin].name);
      return(-1);
    }
    linbuf[intersectFlag]=lin;
  }
  else
  {
    printf("ERROR: line:%s not selected. Only straight lines and splines are supported\n", line[lin].name);
    return(-1);
  }

  if((intersectFlag)&&( lin==linbuf[0]))
  {
    printf("ERROR: second line:%s not selected because its the same as the first selected.\n", line[lin].name);
    return(-1);
  }
  linbuf[intersectFlag]=lin;

  if(intersectFlag)
  {
    /* determine the closest combi of line-end-points determine the orientation of the lines relative to each other */
    v_result( &point[line[linbuf[0]].p1].px, &point[line[linbuf[1]].p1].px, dv[0] );
    v_result( &point[line[linbuf[0]].p2].px, &point[line[linbuf[1]].p2].px, dv[1] );
    v_result( &point[line[linbuf[0]].p1].px, &point[line[linbuf[1]].p2].px, dv[2] );
    v_result( &point[line[linbuf[0]].p2].px, &point[line[linbuf[1]].p1].px, dv[3] );
    icase=4;
    for(i=0; i<4; i++)
    {
      dist=v_betrag(dv[i]);
      if(min>dist) { min=dist; icase=i; }
    }
    switch(icase)
    {
      case 0:
	dirl[0]=1;
	dirl[1]=1;
	break;
      case 1:
	dirl[0]=2;
	dirl[1]=2;
	break;
      case 2:
	dirl[0]=1;
	dirl[1]=2;
	break;
      case 3:
	dirl[0]=2;
	dirl[1]=1;
	break;
	printf("ERROR: icase:%d not known, talk to the programmer\n", icase);
        return(0);
    }
 
    /* go over all divisions of the lines and determine the split-points */
    /* create two arrays ltmp[] of points based on the lines regarding their relative orientation */

    for(i=0; i<2; i++)
    {
      if(dirl[i]==1)  // line.p1 is startpnt
      {
        if(line[linbuf[i]].typ==' ')
	{
          ltmp[i].np=2;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          ltmp[i].pnt[0][0]=point[line[linbuf[i]].p1].px; 
          ltmp[i].pnt[0][1]=point[line[linbuf[i]].p1].py; 
          ltmp[i].pnt[0][2]=point[line[linbuf[i]].p1].pz; 
          ltmp[i].pnt[1][0]=point[line[linbuf[i]].p2].px; 
          ltmp[i].pnt[1][1]=point[line[linbuf[i]].p2].py; 
          ltmp[i].pnt[1][2]=point[line[linbuf[i]].p2].pz; 
	}
	else
	{
          /* temporary increase the division */
	  divbuffer[i]=line[linbuf[i]].div;
          line[linbuf[i]].div=98;
          repLine(linbuf[i]);

          ltmp[i].np=line[linbuf[i]].nip/3;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          p=0;
          for(j=0; j<line[linbuf[i]].nip; j+=3)
          {
            ltmp[i].pnt[p][0]=line[linbuf[i]].ip[j];
            ltmp[i].pnt[p][1]=line[linbuf[i]].ip[j+1];
            ltmp[i].pnt[p][2]=line[linbuf[i]].ip[j+2];
            p++;
	  }
	}
      }
      else // line.p2 is startpnt
      {
        if(line[linbuf[i]].typ==' ')
	{
          ltmp[i].np=2;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          ltmp[i].pnt[0][0]=point[line[linbuf[i]].p2].px; 
          ltmp[i].pnt[0][1]=point[line[linbuf[i]].p2].py; 
          ltmp[i].pnt[0][2]=point[line[linbuf[i]].p2].pz; 
          ltmp[i].pnt[1][0]=point[line[linbuf[i]].p1].px; 
          ltmp[i].pnt[1][1]=point[line[linbuf[i]].p1].py; 
          ltmp[i].pnt[1][2]=point[line[linbuf[i]].p1].pz; 
	}
	else
	{
          /* temporary increase the division */
	  divbuffer[i]=line[linbuf[i]].div;
          line[linbuf[i]].div=98;
          repLine(linbuf[i]);

          ltmp[i].np=line[linbuf[i]].nip/3;
          if ((ltmp[i].pnt = (double **)malloc( (ltmp[i].np)*sizeof(double *)) ) == NULL )
          { printf("\n\nERROR: malloc failure\n"); return(-1); }
          for(j=0; j<ltmp[i].np; j++)
	  {
            if ((ltmp[i].pnt[j] = (double *)malloc( (3)*sizeof(double)) ) == NULL )
            { printf("\n\nERROR: malloc failure\n"); return(-1); }
	  }
          p=0;
          for(j=line[linbuf[i]].nip-3; j>0; j-=3)
          {
            ltmp[i].pnt[p][0]=line[linbuf[i]].ip[j];
            ltmp[i].pnt[p][1]=line[linbuf[i]].ip[j+1];
            ltmp[i].pnt[p][2]=line[linbuf[i]].ip[j+2];
            p++;
	  }
	}
      }
    }

    lflag[0]=1;
    lflag[1]=1;
    /* these pnt-indexes have to be changed due to the results of "searchIntersectionPoints()" */
    lp_ptr_buf[0]=lp_ptr[0]=0;
    lp_ptr_buf[1]=lp_ptr[1]=0;
    flipflop=0;
    do
    {
      if(intersectionPoint( ltmp[0].pnt[lp_ptr[0]], ltmp[0].pnt[lp_ptr[0]+1], ltmp[1].pnt[lp_ptr[1]], ltmp[1].pnt[lp_ptr[1]+1], ps) == 0)
      {
        /* failed! free ltmp */
        for(i=0; i<2; i++)
        {
          /* restore the division */
          if(line[linbuf[i]].typ!=' ')
	  {
            line[linbuf[i]].div=divbuffer[i];
            repLine(linbuf[i]);
	  }
          for(j=0; j<ltmp[i].np; j++) free(ltmp[i].pnt[j]); free(ltmp[i].pnt);
        }
        printf(" ERROR: lines %s %s do not intersect based on a tolerance of:%lf\n", line[linbuf[0]].name, line[linbuf[1]].name, gtol);
        return(-1);
      } 

      /* check if the intersection of the lines (ps) is outside of the interval (ie. l1ptr, l1ptr+1) */
      for(i=0; i<2; i++)
      {
        v_result(ltmp[i].pnt[lp_ptr[i]], ltmp[i].pnt[lp_ptr[i]+1], p1p2);
        v_result(ltmp[i].pnt[lp_ptr[i]], ps, p1ps);
        lp1p2[i]=v_betrag(p1p2);
        lp1ps[i]=v_betrag(p1ps);
        if(v_sprod(p1ps,p1p2)<0.) lp1ps[i]*=-1;

        /* then adjust the interval */

        /* if lp1ps is <0 then decrease lp_ptr[i] */
        /* warning gtol might be set to a very high value to enable intersection, so better not use this */
        //if((lp1ps[i]<(-gtol/scale->w))&&(lp_ptr[i]>0)) lp_ptr[i]--;
        if((lp1ps[i]<(-1e-6))&&(lp_ptr[i]>0)) lp_ptr[i]--;
        /* if lp1ps is >lp1p2 then increase lp_ptr[i] */
        //else if((lp1ps[i]>(lp1p2[i]+gtol/scale->w))&&(lp_ptr[i]<ltmp[i].np-2)) lp_ptr[i]++;
        else if((lp1ps[i]>(lp1p2[i]+1e-6))&&(lp_ptr[i]<ltmp[i].np-2)) lp_ptr[i]++;
        /* if lp1ps is <=lp1p2 then break */
        else lflag[i]=0;
      }
      if(printFlag) printf("lflag %d %d\n", lflag[0],lflag[1]);

      /* escape endless loops */
      if(flipflop)
      {
        if((lp_ptr[0]==lp_ptr_buf[0])&&(lp_ptr[1]==lp_ptr_buf[1]))
        {
          breakflag=1;
          for(i=0; i<2; i++)
          {
	    lp_ptr[i]=lp_ptr_buf[i];
	    lp1ps[i]=lp1ps_buf[i];
	    lp1p2[i]=lp1p2_buf[i];
	  }
          break;
        }
        for(i=0; i<2; i++)
        {
          lp_ptr_buf[i]=lp_ptr[i];
	  lp1ps_buf[i]=lp1ps[i];
	  lp1p2_buf[i]=lp1p2[i];
	}
      }
      flipflop=!flipflop;
    }while( lflag[0] || lflag[1] );

    /* determine the line-length up to the split point, necessary for line-splitting */
    for(i=0; i<2; i++)
    {
      /* check if the split-point is outside */
      if(printFlag) printf("check line[%d]:%s\n", i,line[linbuf[i]].name);
      if(line[linbuf[i]].typ!='s') continue;

      /* restore the division */
      line[linbuf[i]].div=divbuffer[i];
      repLine(linbuf[i]);

      if(!breakflag) if((lp1ps[i]<=0)||(lp1ps[i]>=lp1p2[i])) continue;

      if(printFlag) printf("redefine the trkset of line[%d]:%s\n", i,line[linbuf[i]].name);

      /* determine the line-length up to the splitpoint */
      lps[i]=0.;
      for(j=0; j<lp_ptr[i]; j++)
      {
        v_result(ltmp[i].pnt[j], ltmp[i].pnt[j+1], p1p2);
        lps[i]+=v_betrag(p1p2);
      }
      lps[i]+=lp1ps[i];


      /* redefine the track-set */
      lbez=0.;
      getNewName( name, "se" );
      trkNew=pre_seta(name,"is",0);
      trk=line[linbuf[i]].trk;
      if(dirl[i]==1) /* remove the leading points and add the split-point */
      {
        seta( trkNew, "ps", set[trk].pnt[0] );  // will be redefined later
        for (j=0; j<set[trk].anz_p-1; j++)
        {
          v_result( &point[set[trk].pnt[j]].px, &point[set[trk].pnt[j+1]].px, p1p2 );
          lbez+=v_betrag( p1p2 );
          if(lbez>lps[i]) seta( trkNew, "ps", set[trk].pnt[j+1] );
        }
      }
      if(dirl[i]==2) /* remove the trailing points and add the split-point */
      {
        lps[i]=calcLineLength(linbuf[i])/scale->w-lps[i];
        seta( trkNew, "ps", set[trk].pnt[0] );
        for (j=0; j<set[trk].anz_p-1; j++)
        {
          v_result( &point[set[trk].pnt[j]].px, &point[set[trk].pnt[j+1]].px, p1p2 );
          lbez+=v_betrag( p1p2 );
          if(lbez<lps[i]) seta( trkNew, "ps", set[trk].pnt[j+1] );
          else break;
        }
        seta( trkNew, "ps", set[trk].pnt[set[trk].anz_p-1] );  // will be redefined later
      }
      line[linbuf[i]].trk=trkNew;
      delSet(set[trk].name);
    }

    /* free ltmp */
    for(i=0; i<2; i++) { for(j=0; j<ltmp[i].np; j++) free(ltmp[i].pnt[j]); free(ltmp[i].pnt); }


    if(dirl[0]==2) /* line points to the intersection */
    {
      p1nr=line[linbuf[0]].p2;
      point[line[linbuf[0]].p2].px=ps[0];
      point[line[linbuf[0]].p2].py=ps[1];
      point[line[linbuf[0]].p2].pz=ps[2];
    }
    else 
    {
      p1nr=line[linbuf[0]].p1;
      point[line[linbuf[0]].p1].px=ps[0];
      point[line[linbuf[0]].p1].py=ps[1];
      point[line[linbuf[0]].p1].pz=ps[2];
    }
    if(dirl[1]==2) /* line points to the intersection */
    {
        p2nr=line[linbuf[1]].p2;
        point[line[linbuf[1]].p2].px=ps[0];
        point[line[linbuf[1]].p2].py=ps[1];
        point[line[linbuf[1]].p2].pz=ps[2];
    }
    else 
    {
        p2nr=line[linbuf[1]].p1;
        point[line[linbuf[1]].p1].px=ps[0];
        point[line[linbuf[1]].p1].py=ps[1];
        point[line[linbuf[1]].p1].pz=ps[2];
    }

    /* check if a line is running from ps_line1 to ps_line2 and delete this (will have 0 length) */
    for (j=0; j<set[setall].anz_l; j++)
    {
      l=set[setall].line[j];
      if( ((line[l].p1==p1nr)||(line[l].p1==p2nr)) && ((line[l].p2==p1nr)||(line[l].p2==p2nr)) )
      {
        printf(" line:%s was detected between the intersecting lines and is deleted\n", line[l].name);
        pre_seta(specialset->zap, "l", line[l].name);
        zap(specialset->zap);
        //delLine( 1, &l );
      }
    }
    updateDispLists();
  }

  intersectFlag=!intersectFlag;
  return(intersectFlag);
}



#define DS 10000
#define MIN_SPROD 1.e30
int qsplitLine( int l, int x, int y )
{
  int k;
  double p0[3], p1[3], p01[3], pm0[3], pm1[3], pm01[3], pm02[3], p0pm0[3], ps[3], p0ps[3];
  double eu[3], ev[3], en[3], eg[3], g;
  int   p1_nr, p2_nr, ps_nr, l_nr[2], line_bias;
  char name[MAX_LINE_LENGTH], trk[2][MAX_LINE_LENGTH];
  double sprod_pk[DS], ek[2][3], pm00[3], pm0k[3];
  double pk[2][3], ps_lbez=0;
  double min_sprod, sprod_euek;

  double pbuf[3];

  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];

  static GLdouble wx, wy, wz;  /*  returned world x, y, z coords  */
  static GLdouble nx, ny, nz;  /*  new world x, y, z coords  */
  static int flag;


  /* erzeugen von zwei Raumpunkten auf dem Mauspeilstrahl (MausPos.-Z-Richtung) */

  /* first create a temporary point in window-z direction at cursor-x,y pos.  */
  glutSetWindow( w1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  moveModel();
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  flag=gluProject( centerPnt[0], centerPnt[1], centerPnt[2], mvmatrix, projmatrix,
     viewport,  &wx, &wy, &wz);
  if (flag==GL_FALSE)
  {
    printf(" ERROR: Malfunction in qsplitLine(), please reselect\n");
    return(-1);
  }
  /* printf (" Win coords are %d (%lf, %lf, %lf)\n",  flag,  wx   ,  wy   , wz     ); */
  
  wx=(GLdouble)x; wy=(GLdouble)(viewport[3]-y);
  flag=gluUnProject ( wx, wy, wz, mvmatrix, projmatrix, viewport, &nx, &ny, &nz);
  /* printf ("new World coords are %d (%lf, %lf, %lf)\n", flag, nx, ny, nz); */
  if (flag==GL_TRUE)
  {
    pm0[0]=nx;
    pm0[1]=ny;
    pm0[2]=nz;
  }
  else
  {
    printf(" ERROR: Malfunction in qsplitLine(), please reselect\n");
    return(-1);
  }
  
  wz++;
  flag=gluUnProject ( wx, wy, wz, mvmatrix, projmatrix, viewport, &nx, &ny, &nz);
  /* printf ("new World coords are %d (%lf, %lf, %lf)\n", flag, nx, ny, nz); */
  if (flag==GL_TRUE)
  {
    pm1[0]=nx;
    pm1[1]=ny;
    pm1[2]=nz;
  }
  else
  {
    printf(" ERROR: Malfunction in qsplitLine(), please reselect\n");
    return(-1);
  }

  /* split the first line (l) */ 

  /* berechne den ersten Einheitsvektor der Schnittebene, eu==Peilstrahlrichtung */
  v_result( pm0, pm1, pm01 );
  v_norm( pm01, eu );

  /* bestimme die Koordinaten der Linienendpunkte */
  p1_nr= line[l].p1;
  p0[0] = point[p1_nr].px;
  p0[1] = point[p1_nr].py;
  p0[2] = point[p1_nr].pz;
  p2_nr= line[l].p2;
  p1[0] = point[p2_nr].px;
  p1[1] = point[p2_nr].py;
  p1[2] = point[p2_nr].pz;
  /*
  printf("endpnt P1:%s %lf %lf %lf\n", point[line[l].p1].name, p0[0],p0[1],p0[2]);
  printf("endpnt P2:%s %lf %lf %lf\n", point[line[l].p2].name, p1[0],p1[1],p1[2]);
  */

  /* berechne den schnittpunkt der zu splittenden linie mit der ebene a(eu,ev) */
  ps_nr=-1;
  if( line[l].typ==' ')
  {
    /* berechne den zweiten Einheitsvektor der Schnittebene, ev==Peilstrahl x Linie */
    v_result( p0, p1, p01 );
    v_prod( pm01, p01, pm02 );
    v_norm( pm02, ev );

    /* berechne den Normalenvektor der Schnittebene */
    v_prod( eu, ev, en );
    v_norm( en, en );

    /* berechne den Einheitsvektor der zu splitenden Linie */
    v_norm( p01, eg );

    /* bestimme den Abstand zwischen den Aufpunkten der Linie und Ebene  */
    v_result( p0, pm0, p0pm0 );

    /* berechne die Konstante g zur berechnung von ps (Schnittpunkt) ps=p0+eg*g  */
    g = AsplitL( p0pm0, eu, ev, eg, en );

    v_scal( &g, eg, p0ps );
    v_add( p0, p0ps, ps );
    strcpy(trk[0], " ");
    strcpy(trk[1], " ");

    /* erzeuge den Punkt */
    getNewName( name, "p" );
    printf(" create point:%s %lf %lf %lf\n", name, ps[0], ps[1], ps[2] );
    ps_nr  = pnt( name, ps[0], ps[1], ps[2], 0 );
  }
  else
  {
    /* Zerlege die linie in DS Teile und bestimme die durchtritte durch die Ebene */
    /* ein Durchtritt erfolgt, wenn das vorzeichen des skalarproduckts sich aendert */
    /* speichere den linienpunkt vor und hinter dem durchtritt fuer die spaetere */
    /* interpolation. Bestimme das dem peilstrahl am naechsten liegende punktepaar */
    /* das naechstliegende punktepaar ist das mit dem kleinsten skalarprodukt */

    /* berechne den Einheitsvektor zum 1. kontrollpunkt */
    pk[0][0]= p0[0];
    pk[0][1]= p0[1];
    pk[0][2]= p0[2];
    v_result( pm1, pk[0], pm00 );
    v_norm( pm00, ek[0] );

    flag=0;
    min_sprod=MIN_SPROD;

    /* set the line-bias to 1 for the search of the split-point */
    /* reset it after this operation to the original one */
    line_bias=line[l].bias;
    line[l].bias=1;
    for (k=0; k<DS; k++)
    {
      if (line[l].typ=='a')
      {
        if(arcNodes( l, k,(int) DS, pbuf )==-1)      { printf("ERROR in spliting\n"); }
      }
      else if (line[l].typ=='s')
      {
        if(splineNodes( l, k,(int) DS, pbuf )==-1)   { printf("ERROR in spliting\n"); }
      }
      pk[1][0]= pbuf[0];
      pk[1][1]= pbuf[1];
      pk[1][2]= pbuf[2];
      

      /* die richtung der schnittebene wird fuer jeden kontrollpunkt neu berechnet */
      /* berechne den zweiten Einheitsvektor der Schnittebene */
      v_result( pk[0], pk[1], p01 );
      v_prod( pm01, p01, pm02 );
      v_norm( pm02, ev );

      /* berechne den Normalenvektor der Schnittebene */
      v_prod( eu, ev, en );

      /* winkel zum ersten kontrollpunkt */
      sprod_pk[0]=v_sprod( en, ek[0]);

      /* berechne den Einheitsvektor zu diesem kontrollpunkt */
      v_result( pm1, pk[1], pm0k );
      v_norm( pm0k, ek[1] );
      sprod_pk[1]=v_sprod( en, ek[1]);

      /* haben wir einen durchstosspunkt? (vorzeichenwechsel von sprod)  */
      if( (sprod_pk[0]*sprod_pk[1]) <=0. )
      {
        flag=1;
        /* durchstoss gefunden, ist er der bisher naechstliegende? */
        sprod_euek=v_sprod( eu, ek[0])+v_sprod( eu, ek[1]);
	/*
        printf("sprod_euek:%lf min_sprod:%lf k:%d DS:%d\n", sprod_euek, min_sprod,k,DS); 
        printf(" pnt0:%lf %lf %lf sprod_pn:%lf\n", pk[0][0], pk[0][1], pk[0][2], sprod_pk[0]);
        printf(" pnt1:%lf %lf %lf sprod_pn:%lf\n", pk[1][0], pk[1][1], pk[1][2], sprod_pk[1]);
	*/
	
        if( sprod_euek<min_sprod )
        {
          min_sprod=sprod_euek;
          p0[0] =pk[0][0] ;
          p0[1] =pk[0][1] ;
          p0[2] =pk[0][2] ;
          p1[0] =pk[1][0] ;
          p1[1] =pk[1][1] ;
          p1[2] =pk[1][2] ;
          /*
          ps_lbez=(double)k/(double)DS;           // position, bezogene linienlaenge
          printf(" min_sprod:%lf ps_lbez:%lf \n", min_sprod, ps_lbez );
          */
	}
      }
      pk[0][0]=pk[1][0] ;
      pk[0][1]=pk[1][1] ;
      pk[0][2]=pk[1][2] ;
      ek[0][0]=ek[1][0] ;
      ek[0][1]=ek[1][1] ;
      ek[0][2]=ek[1][2] ;
    }
    /* reset the line-bias */
    line[l].bias=line_bias;
    
    /* uebernehme einen der durchstoss-punkte als splitpunkt */
    if(flag)
    {  
      ps[0]=p0[0];
      ps[1]=p0[1];
      ps[2]=p0[2];
    }
    else
    {
      errMsg(" ERROR: Malfunction in qsplitLine(), please reselect\n");
      return(-1);
    }

    /* erzeuge den splitpunkt */
    getNewName( name, "p" );
    printf(" create point:%s %lf %lf %lf\n", name, ps[0], ps[1], ps[2] );
    ps_nr  = pnt( name, (double)ps[0], (double)ps[1], (double)ps[2], 0 );
  }
  if(ps_nr==-1) 
  {
    errMsg(" ERROR: Malfunction in qsplitLine(), please reselect\n");
    return(-1);
  }

  splitLine(l, l_nr, ps_nr);

  /* loesche basislinie */
  delLine( 1, &l );

  updateDispLists();
  return(ps_nr);
}



void pre_align( int nr, int flag )
{
  static double p[3][3], p02[3], p01[3], p03[3], p03xz[3];
  double  ay, axz;
  double  p03_b, p03xz_b, scalar, sig;

  GLint    viewport[4];
  GLdouble mvmatrix[16], projmatrix[16];

  static GLdouble wx, wy, wz;  /*  returned world x, y, z coords  */
  static int flag2;

  glutSetWindow( w1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  moveModel();
  glGetIntegerv (GL_VIEWPORT, viewport);
  glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

  if (flag)
  {
    flag2=gluProject( point[nr].px, point[nr].py, point[nr].pz, mvmatrix, projmatrix,
      viewport,  &wx, &wy, &wz);
    if (flag2==GL_FALSE)
      printf("WARNING: Malfunction in movePoint(), please reselect\n");
  }
  else
  {
    flag2=gluProject( node[nr].nx, node[nr].ny, node[nr].nz, mvmatrix, projmatrix,
      viewport,  &wx, &wy, &wz);
    if (flag2==GL_FALSE)
      printf("WARNING: Malfunction in movePoint(), please reselect\n");
  }

  p[qaliCounter][0] = wx;
  p[qaliCounter][1] = wy;
  p[qaliCounter][2] = -wz*height_w1/ds*2.;
  qaliCounter++;

  if (qaliCounter == 1)
  {
    /* 1. Punkt gewaehlt == neuer Drehpunkt und neue Arbeitsebene */
  if (flag==0)
        center( node[nr].nx, node[nr].ny, node[nr].nz);
  if (flag==1)
        center( point[nr].px, point[nr].py, point[nr].pz);
  }
  if (qaliCounter == 3)
  {
    qaliCounter=0;
    /* 3 Punkte gewaehlt, berechne die Normale auf der Ebene */
    v_result( &p[0][0], &p[1][0],  p01);
    v_result( &p[0][0], &p[2][0],  p02);
    v_prod( p01, p02, p03 );

    /* Richtungsumkehr wenn die z-Komponente der Normalen nach hinten zeigt (-) */
    if (p03[2]<0)
    {
       p03[0]*=-1;
       p03[1]*=-1;
       p03[2]*=-1;
    }
    /* return wenn die ebene bereits zu der bildschirmebene paralel ist */
    if ((!p03[0])||(!p03[1])) return;

    /* drehe den Ort des Betrachters  */
    ay=atan(p03[0]/p03[2]);

    p03xz[0]=p03[0];
    p03xz[1]=0.;
    p03xz[2]=p03[2];
    scalar = v_sprod( p03xz, p03 );
    p03xz_b= v_betrag( p03xz );
    p03_b  = v_betrag( p03 );
    axz =acos( scalar/p03xz_b/p03_b ) * 180./PI;

    /* vorzeichen von der y-Komponente der Normalen wird zum vorzeichen von axz */
    scalar= sqrt(p03[1]*p03[1]);
    sig = p03[1]/scalar;
    axz = axz*sig;

    ay=ay*180./PI;

    rot_r( ay );
    rot_u( axz );
  }
}



double *createSplitPoint( double *ns, double *n1, double *n2, double *pn, double *eu, double *ev )
{
  double eg[3], en[3], n1n2[3], n1pn[3], n1ns[3];
  double g;

  v_result( n1, n2, n1n2 );

  /* berechne den Einheitsvektor der zu splitenden Linie */
  v_norm( n1n2, eg );

  /* bestimme den Abstand zwischen den Aufpunkten der Linie und Ebene  */
  v_result( n1, pn, n1pn );

  /* berechne die Konstante g zur berechnung von ps (Schnittpunkt) ps=p0+eg*g  */
  v_prod(eu,ev,en);
  v_norm( en, en );
  g = AsplitL( n1pn, eu, ev, eg, en );

  /* pos von ns bezogen auf |n1n2| */
  ns[3] = g/v_betrag(n1n2);
  if(ns[3]>1.000001) return(NULL);
  if(ns[3]<-0.000001) return(NULL);
  v_scal( &g, eg, n1ns );
  v_add( n1, n1ns, ns );
  return(ns);

  /* winkel zw n1ns und n1n2 bestimmen. wenn >0 dann ist ns zwischen n1 und n2 */
  //if((double)((int)(g*1000)) < (g-1.)) return(NULL); 
  //else if((v_sprod(n1ns,n1n2)>=0)&&((int)(g*1000)<=(int)(v_betrag(n1n2)*1000))) return(ns);
  //else return(NULL);
}


/* calculate values for all Datasets on new nodes */
void updLcase(int lc, int setNr)
{
  int i,j,k,n, compareChars;
  char buffer[2][MAX_LINE_LENGTH];
  double s[6], p[3];
  if(!anz->l) return;

 next_lc:; 

  /* check if the data of the specified lcase (Dataset) are already available */
  if (!lcase[lc].loaded)
  {
    if( pre_readfrdblock(copiedNodeSets , lc, anz, node, lcase )==-1) 
    {
      printf("ERROR in updLcase: Could not read data for Dataset:%d\n", lc+1); 
      return;
    }
    calcDatasets( lc, anz, node, lcase );
    recompileEntitiesInMenu(lc);
  }

  for(i=0; i<lcase[lc].ncomps; i++)
  {
    if ( (lcase[lc].dat[i] = (float *)realloc(lcase[lc].dat[i], (anz->nmax+1) * sizeof(float))) == NULL )
        printf("\n\n ERROR: realloc failure updLcase\n\n" );	               
    for (n=0; n<set[setNr].anz_n; n++)
    {
      lcase[lc].dat[i][qcut_nod[n].nr] = lcase[lc].dat[i][qcut_nod[n].n1]*(1.-qcut_nod[n].val) + lcase[lc].dat[i][qcut_nod[n].n2]*qcut_nod[n].val;
      //printf("val:%lf n0:%lf n1:%lf n2:%lf\n", val, lcase[lc].dat[i][ns], lcase[lc].dat[i][n1], lcase[lc].dat[i][n2]);
    }
  }

  /* normal stress */
  if((lcase[lc].ncomps>=25)&&(lcase[lc].ictype[5]=4)&&(lcase[lc].ictype[6]!=4))
  {
    for (n=0; n<set[setNr].anz_n; n++)
    {
        s[0]=lcase[lc].dat[0][qcut_nod[n].nr];
        s[1]=lcase[lc].dat[1][qcut_nod[n].nr];
        s[2]=lcase[lc].dat[2][qcut_nod[n].nr];
        s[3]=lcase[lc].dat[3][qcut_nod[n].nr];
        s[4]=lcase[lc].dat[4][qcut_nod[n].nr];
        s[5]=lcase[lc].dat[5][qcut_nod[n].nr];
        if(v_betrag(node[qcut_nod[n].nr].nv))
        {
          /* calculate the stress in normal direction */
          p[0]=
          s[0]*node[qcut_nod[n].nr].nv[0]+
          s[3]*node[qcut_nod[n].nr].nv[1]+
          s[5]*node[qcut_nod[n].nr].nv[2];
          p[1]=
          s[3]*node[qcut_nod[n].nr].nv[0]+
          s[1]*node[qcut_nod[n].nr].nv[1]+
          s[4]*node[qcut_nod[n].nr].nv[2];
          p[2]=
          s[5]*node[qcut_nod[n].nr].nv[0]+
          s[4]*node[qcut_nod[n].nr].nv[1]+
          s[2]*node[qcut_nod[n].nr].nv[2];
          lcase[lc].dat[24][qcut_nod[n].nr]= v_sprod(node[qcut_nod[n].nr].nv, p);
        }
      }
    }

  /* update related DISP Dataset */
  /* if the selected lc is not a disp lc, search a related disp lc */
  if(compare(lcase[lc].name, "DISP", 4)!=4)
  {
    /* since real and imaginary part use different names since ccx_2.9 it is necessary to compare the 
       names only for the length excluding the last char if its a 'I' */
    strcpy(buffer[0],"DISP ");
    strcpy(buffer[1],"DISPI");
    compareChars=strlen(lcase[lc].name)-1;
    for(k=compareChars;k>0; k--) if(lcase[lc].name[k]!=' ') break;
    compareChars=k+1;
    if(lcase[lc].name[compareChars-1]=='I') j=1; else j=0;;

    if(lc) { for (i=lc-1; i>=0; i--) { if(lcase[i].step_number!=lcase[lc].step_number) break; } i++; }
    else i=1;
    while((i<anz->l)&&(lcase[i].step_number==lcase[lc].step_number))
    {
      if((compare(lcase[i].name, buffer[j], 5)==5)&&(lcase[i].ictype[0]!= 12))
      {
	//printf("lcase[i].name:%s lcase[lc].name:%s compareChars:%d\n", lcase[i].name,lcase[lc].name,compareChars);
        lc=i;
        goto next_lc;
      }
      i++;
    }
  }
}



void cutElement(int elnr, int nset, int eset, double *en, double *pn, double *eu, double *ev, CutParam *cutParam)
{
  int i,j,k,n1,n2,nm,nn1,nn2,anz_n;
  int   nvalid[CUTS],ipre[CUTS], midnodeFlag=0;
  double v02[3], v01[3], ns[4], cg[3];
  double local_gtol=1.e-9; /* dist to merge new nodes (might be set to gtol) */
  double lengthTolMax, lengthTolMin, refTolMax, refTolMin;
  Rsort *rsort=NULL;

  /* mapping of nodeindexes from frd to pre_cut */
  int nhe20[]  = {8,9,10,11, 12,13,14,15, 16,17,18,19};
  int ntet10[] = {4,5,6, 7,8,9};

  double eg[3], n1n2[3], n1pn[3], n1ns[3];
  double g;

  cutParam->anz_n=anz_n=0;
  for(i=0; i<CUTS; i++) cutParam->nref[i]=0;

  //printf("e:%d t:%d ptr:%x\n", elnr,e_enqire[elnr].type,cutParam);

  // based on try and error!
  lengthTolMax=1.00001;
  lengthTolMin=-1e-10;
  refTolMax=0.99;
  refTolMin=0.01;

  if ((e_enqire[elnr].type == 1)||(e_enqire[elnr].type == 4))  /* HEXA8,20 */
  {
    for(n1=0; n1<4; n1++)
    {
      if(n1==3) n2=0;
      else n2=n1+1;

      v_result( &node[e_enqire[elnr].nod[n1]].nx, &node[e_enqire[elnr].nod[n2]].nx, n1n2 );
      v_result( &node[e_enqire[elnr].nod[n1]].nx, pn, n1pn );
      v_norm( n1n2, eg );
      g = AsplitL( n1pn, eu, ev, eg, en );
      if(g==MAX_FLOAT)
      {
	if(e_enqire[elnr].type == 4) midnodeFlag=1;
        continue;
      }
      ns[3] = g/v_betrag(n1n2);
      if(ns[3]>lengthTolMax) continue;
      if(ns[3]<lengthTolMin) continue;
      if(ns[3]>=refTolMax) cutParam->nref[anz_n]=e_enqire[elnr].nod[n2];
      if(ns[3]<refTolMin) cutParam->nref[anz_n]=e_enqire[elnr].nod[n1];
      if(e_enqire[elnr].type == 4)
      {
        midnodeFlag=1;
	// take the midside node into account
	if(ns[3]<0.5)
	{
	  nn1=n1;
	  nn2=nhe20[n1];
	  ns[3]*=2.;
	}
	else
	{
	  nn1=nhe20[n1];
	  nn2=n2;
	  ns[3]=(ns[3]-0.5)*2.;
	}	  
        ns[0]= node[e_enqire[elnr].nod[nn1]].nx*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nx*ns[3];
        ns[1]= node[e_enqire[elnr].nod[nn1]].ny*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].ny*ns[3];
        ns[2]= node[e_enqire[elnr].nod[nn1]].nz*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nz*ns[3];
        cutParam->npre[anz_n][4]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=e_enqire[elnr].nod[n2];
        cutParam->npre[anz_n][6]=e_enqire[elnr].nod[nn1];
        cutParam->npre[anz_n][7]=e_enqire[elnr].nod[nn2];
      }
      else
      {
        v_scal( &g, eg, n1ns );
        v_add( &node[e_enqire[elnr].nod[n1]].nx, n1ns, ns );
        cutParam->npre[anz_n][4]=cutParam->npre[anz_n][6]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=cutParam->npre[anz_n][7]=e_enqire[elnr].nod[n2];
      }
      cutParam->npre[anz_n][0]=ns[0];
      cutParam->npre[anz_n][1]=ns[1];
      cutParam->npre[anz_n][2]=ns[2];
      cutParam->npre[anz_n][3]=ns[3];
      nvalid[anz_n]=1;
      anz_n++;
      if(anz_n>=CUTS) { printf(" ERROR in pre_cut (talk to the programmer), increase CUTS at least to %d\n", anz_n+1); exit(0); }
    }
    for(n1=0; n1<4; n1++)
    {
      n2=n1+4;
      v_result( &node[e_enqire[elnr].nod[n1]].nx, &node[e_enqire[elnr].nod[n2]].nx, n1n2 );
      v_result( &node[e_enqire[elnr].nod[n1]].nx, pn, n1pn );
      v_norm( n1n2, eg );
      g = AsplitL( n1pn, eu, ev, eg, en );
      if(g==MAX_FLOAT)
      {
	if(e_enqire[elnr].type == 4) midnodeFlag=1;
        continue;
      }
      ns[3] = g/v_betrag(n1n2);
      if(ns[3]>lengthTolMax) continue;
      if(ns[3]<lengthTolMin) continue;
      if(ns[3]>=refTolMax) cutParam->nref[anz_n]=e_enqire[elnr].nod[n2];
      if(ns[3]<refTolMin) cutParam->nref[anz_n]=e_enqire[elnr].nod[n1];
      if(e_enqire[elnr].type == 4)
      {
        midnodeFlag=1;
	// take the midside node into account
	if(ns[3]<0.5)
	{
	  nn1=n1;
	  nn2=nhe20[n1+4];
	  ns[3]*=2.;
	}
	else
	{
	  nn1=nhe20[n1+4];
	  nn2=n2;
	  ns[3]=(ns[3]-0.5)*2.;
	}	  
        ns[0]= node[e_enqire[elnr].nod[nn1]].nx*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nx*ns[3];
        ns[1]= node[e_enqire[elnr].nod[nn1]].ny*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].ny*ns[3];
        ns[2]= node[e_enqire[elnr].nod[nn1]].nz*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nz*ns[3];
        cutParam->npre[anz_n][4]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=e_enqire[elnr].nod[n2];
        cutParam->npre[anz_n][6]=e_enqire[elnr].nod[nn1];
        cutParam->npre[anz_n][7]=e_enqire[elnr].nod[nn2];
      }
      else
      {
        v_scal( &g, eg, n1ns );
        v_add( &node[e_enqire[elnr].nod[n1]].nx, n1ns, ns );
        cutParam->npre[anz_n][4]=cutParam->npre[anz_n][6]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=cutParam->npre[anz_n][7]=e_enqire[elnr].nod[n2];
      }
      cutParam->npre[anz_n][0]=ns[0];
      cutParam->npre[anz_n][1]=ns[1];
      cutParam->npre[anz_n][2]=ns[2];
      cutParam->npre[anz_n][3]=ns[3];
      nvalid[anz_n]=1;
      anz_n++;
      if(anz_n>=CUTS) { printf(" ERROR in pre_cut (talk to the programmer), increase CUTS at least to %d\n", anz_n+1); exit(0); }
    }
    for(n1=4; n1<8; n1++)
    {
      if(n1==7) n2=4;
      else n2=n1+1;
      v_result( &node[e_enqire[elnr].nod[n1]].nx, &node[e_enqire[elnr].nod[n2]].nx, n1n2 );
      v_result( &node[e_enqire[elnr].nod[n1]].nx, pn, n1pn );
      v_norm( n1n2, eg );
      g = AsplitL( n1pn, eu, ev, eg, en );
      if(g==MAX_FLOAT)
      {
        if(e_enqire[elnr].type == 4) midnodeFlag=1;
        continue;
      }
      ns[3] = g/v_betrag(n1n2);
      if(ns[3]>lengthTolMax) continue;
      if(ns[3]<lengthTolMin) continue;
      if(ns[3]>=refTolMax) cutParam->nref[anz_n]=e_enqire[elnr].nod[n2];
      if(ns[3]<refTolMin) cutParam->nref[anz_n]=e_enqire[elnr].nod[n1];
      if(e_enqire[elnr].type == 4)
      {
        midnodeFlag=1;
	// take the midside node into account
	if(ns[3]<0.5)
	{
	  nn1=n1;
	  nn2=nhe20[n1+4];
	  ns[3]*=2.;
	}
	else
	{
	  nn1=nhe20[n1+4];
	  nn2=n2;
	  ns[3]=(ns[3]-0.5)*2.;
	}	  
        ns[0]= node[e_enqire[elnr].nod[nn1]].nx*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nx*ns[3];
        ns[1]= node[e_enqire[elnr].nod[nn1]].ny*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].ny*ns[3];
        ns[2]= node[e_enqire[elnr].nod[nn1]].nz*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nz*ns[3];
        cutParam->npre[anz_n][4]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=e_enqire[elnr].nod[n2];
        cutParam->npre[anz_n][6]=e_enqire[elnr].nod[nn1];
        cutParam->npre[anz_n][7]=e_enqire[elnr].nod[nn2];
      }
      else
      {
        v_scal( &g, eg, n1ns );
        v_add( &node[e_enqire[elnr].nod[n1]].nx, n1ns, ns );
        cutParam->npre[anz_n][4]=cutParam->npre[anz_n][6]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=cutParam->npre[anz_n][7]=e_enqire[elnr].nod[n2];
      }
      cutParam->npre[anz_n][0]=ns[0];
      cutParam->npre[anz_n][1]=ns[1];
      cutParam->npre[anz_n][2]=ns[2];
      cutParam->npre[anz_n][3]=ns[3];
      nvalid[anz_n]=1;
      anz_n++;
      if(anz_n>=CUTS) { printf(" ERROR in pre_cut (talk to the programmer), increase CUTS at least to %d\n", anz_n+1); exit(0); }
    }
    // if an edge was cutted then store the midnodes of all edges
    if(midnodeFlag)
    {
      for(n1=0; n1<4; n1++)
      {
        if(n1==3) n2=0;
        else n2=n1+1;
	// create the midnode reference
        nm=e_enqire[elnr].nod[nhe20[n1]];
        nn1=e_enqire[elnr].nod[n1];
        nn2=e_enqire[elnr].nod[n2];
    sem_wait(&sem_cute);
        if ( (n1nm[nn1].n2 = (int *)realloc( n1nm[nn1].n2, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }    
        if ( (n1nm[nn1].nm = (int *)realloc( n1nm[nn1].nm, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }
        n1nm[nn1].n2[n1nm[nn1].sum]=nn2;
        n1nm[nn1].nm[n1nm[nn1].sum]=nm;
        n1nm[nn1].sum++;
    sem_post(&sem_cute);
      }
      for(n1=0; n1<4; n1++)
      {
        n2=n1+4;
	// create the midnode reference
        nm=e_enqire[elnr].nod[nhe20[n1+4]];
        nn1=e_enqire[elnr].nod[n1];
        nn2=e_enqire[elnr].nod[n2];
    sem_wait(&sem_cute);
        if ( (n1nm[nn1].n2 = (int *)realloc( n1nm[nn1].n2, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }    
        if ( (n1nm[nn1].nm = (int *)realloc( n1nm[nn1].nm, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }
        n1nm[nn1].n2[n1nm[nn1].sum]=nn2;
        n1nm[nn1].nm[n1nm[nn1].sum]=nm;
        n1nm[nn1].sum++;
    sem_post(&sem_cute);
      }
      for(n1=4; n1<8; n1++)
      {
        if(n1==7) n2=4;
        else n2=n1+1;
	// create the midnode reference
        nm=e_enqire[elnr].nod[nhe20[n1+4]];
        nn1=e_enqire[elnr].nod[n1];
        nn2=e_enqire[elnr].nod[n2];
    sem_wait(&sem_cute);
        if ( (n1nm[nn1].n2 = (int *)realloc( n1nm[nn1].n2, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }    
        if ( (n1nm[nn1].nm = (int *)realloc( n1nm[nn1].nm, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }
        n1nm[nn1].n2[n1nm[nn1].sum]=nn2;
        n1nm[nn1].nm[n1nm[nn1].sum]=nm;
        n1nm[nn1].sum++;
    sem_post(&sem_cute);
      }
    }
      
  } /* end he8,he20 */
  if ((e_enqire[elnr].type == 3)||(e_enqire[elnr].type == 6))  /* TET4,10 */
  {
    for(n1=0; n1<3; n1++)
    {
      if(n1==2) n2=0;
      else n2=n1+1;
      v_result( &node[e_enqire[elnr].nod[n1]].nx, &node[e_enqire[elnr].nod[n2]].nx, n1n2 );
      v_result( &node[e_enqire[elnr].nod[n1]].nx, pn, n1pn );
      v_norm( n1n2, eg );
      g = AsplitL( n1pn, eu, ev, eg, en );
      //if(elnr==107052) printf("1n %d %d g:%e sp:%e\n",e_enqire[elnr].nod[n1],e_enqire[elnr].nod[n2],g,v_sprod(n1n2,en));
      if(g==MAX_FLOAT)
      {
	if(e_enqire[elnr].type == 6) midnodeFlag=1;
        continue;
      }
      ns[3] = g/v_betrag(n1n2);
      if(ns[3]>lengthTolMax) continue;
      if(ns[3]<lengthTolMin) continue;
      if(ns[3]>=refTolMax) cutParam->nref[anz_n]=e_enqire[elnr].nod[n2];
      if(ns[3]<refTolMin) cutParam->nref[anz_n]=e_enqire[elnr].nod[n1];
      if(e_enqire[elnr].type == 6)
      {
        midnodeFlag=1;
	// take the midside node into account
	if(ns[3]<0.5)
	{
	  nn1=n1;
	  nn2=ntet10[n1];
	  ns[3]*=2.;
	}
	else
	{
	  nn1=ntet10[n1];
	  nn2=n2;
	  ns[3]=(ns[3]-0.5)*2.;
	}	  
        ns[0]= node[e_enqire[elnr].nod[nn1]].nx*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nx*ns[3];
        ns[1]= node[e_enqire[elnr].nod[nn1]].ny*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].ny*ns[3];
        ns[2]= node[e_enqire[elnr].nod[nn1]].nz*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nz*ns[3];
        cutParam->npre[anz_n][4]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=e_enqire[elnr].nod[n2];
        cutParam->npre[anz_n][6]=e_enqire[elnr].nod[nn1];
        cutParam->npre[anz_n][7]=e_enqire[elnr].nod[nn2];
      }
      else
      {
        v_scal( &g, eg, n1ns );
        v_add( &node[e_enqire[elnr].nod[n1]].nx, n1ns, ns );
        cutParam->npre[anz_n][4]=cutParam->npre[anz_n][6]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=cutParam->npre[anz_n][7]=e_enqire[elnr].nod[n2];
      }
      cutParam->npre[anz_n][0]=ns[0];
      cutParam->npre[anz_n][1]=ns[1];
      cutParam->npre[anz_n][2]=ns[2];
      cutParam->npre[anz_n][3]=ns[3];
      nvalid[anz_n]=1;
      anz_n++;
      if(anz_n>=CUTS) { printf(" ERROR in pre_cut (talk to the programmer), increase CUTS at least to %d\n", anz_n+1); exit(0); }
    }
    for(n1=0; n1<3; n1++)
    {
      n2=3;
      v_result( &node[e_enqire[elnr].nod[n1]].nx, &node[e_enqire[elnr].nod[n2]].nx, n1n2 );
      v_result( &node[e_enqire[elnr].nod[n1]].nx, pn, n1pn );
      v_norm( n1n2, eg );
      g = AsplitL( n1pn, eu, ev, eg, en );
      //if(elnr==107052 printf("n %d %d g:%e sp:%e\n",e_enqire[elnr].nod[n1],e_enqire[elnr].nod[n2],g,v_sprod(n1n2,en));
      if(g==MAX_FLOAT)
      {
        if(e_enqire[elnr].type == 6) midnodeFlag=1;
        continue;
      }
      ns[3] = g/v_betrag(n1n2);
      if(ns[3]>lengthTolMax) continue;
      if(ns[3]<lengthTolMin) continue;
      if(ns[3]>=refTolMax) cutParam->nref[anz_n]=e_enqire[elnr].nod[n2];
      if(ns[3]<refTolMin) cutParam->nref[anz_n]=e_enqire[elnr].nod[n1];
      if(e_enqire[elnr].type == 6)
      {
	midnodeFlag=1;
	// take the midside node into account
	if(ns[3]<0.5)
	{
	  nn1=n1;
	  nn2=ntet10[n1+3];
	  ns[3]*=2.;
	}
	else
	{
	  nn1=ntet10[n1+3];
	  nn2=n2;
	  ns[3]=(ns[3]-0.5)*2.;
	}	  
        ns[0]= node[e_enqire[elnr].nod[nn1]].nx*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nx*ns[3];
        ns[1]= node[e_enqire[elnr].nod[nn1]].ny*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].ny*ns[3];
        ns[2]= node[e_enqire[elnr].nod[nn1]].nz*(1.-ns[3]) + node[e_enqire[elnr].nod[nn2]].nz*ns[3];
        cutParam->npre[anz_n][4]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=e_enqire[elnr].nod[n2];
        cutParam->npre[anz_n][6]=e_enqire[elnr].nod[nn1];
        cutParam->npre[anz_n][7]=e_enqire[elnr].nod[nn2];
      }
      else
      {
        v_scal( &g, eg, n1ns );
        v_add( &node[e_enqire[elnr].nod[n1]].nx, n1ns, ns );
        cutParam->npre[anz_n][4]=cutParam->npre[anz_n][6]=e_enqire[elnr].nod[n1];
        cutParam->npre[anz_n][5]=cutParam->npre[anz_n][7]=e_enqire[elnr].nod[n2];
      }
      cutParam->npre[anz_n][0]=ns[0];
      cutParam->npre[anz_n][1]=ns[1];
      cutParam->npre[anz_n][2]=ns[2];
      cutParam->npre[anz_n][3]=ns[3];
      nvalid[anz_n]=1;
      anz_n++;
      if(anz_n>=CUTS) { printf(" ERROR in pre_cut (talk to the programmer), increase CUTS at least to %d\n", anz_n+1); exit(0); }
    }

    // if an edge was cutted then store the midnodes of all edges
    if(midnodeFlag)
    {
      for(n1=0; n1<3; n1++)
      {
        if(n1==2) n2=0;
        else n2=n1+1;
	// create the midnode reference
        nm=e_enqire[elnr].nod[ntet10[n1]];
        nn1=e_enqire[elnr].nod[n1];
        nn2=e_enqire[elnr].nod[n2];
    sem_wait(&sem_cute);
        if ( (n1nm[nn1].n2 = (int *)realloc( n1nm[nn1].n2, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }    
        if ( (n1nm[nn1].nm = (int *)realloc( n1nm[nn1].nm, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }
        n1nm[nn1].n2[n1nm[nn1].sum]=nn2;
        n1nm[nn1].nm[n1nm[nn1].sum]=nm;
        n1nm[nn1].sum++;
    sem_post(&sem_cute);
      }
      for(n1=0; n1<3; n1++)
      {
        n2=3;
	// create the midnode reference
        nm=e_enqire[elnr].nod[ntet10[n1+3]];
        nn1=e_enqire[elnr].nod[n1];
        nn2=e_enqire[elnr].nod[n2];
    sem_wait(&sem_cute);
        if ( (n1nm[nn1].n2 = (int *)realloc( n1nm[nn1].n2, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }    
        if ( (n1nm[nn1].nm = (int *)realloc( n1nm[nn1].nm, (n1nm[nn1].sum+1) * sizeof(int))) == NULL )
        { printf("\n\n ERROR in mids: realloc\n\n") ; exit(-1); }
        n1nm[nn1].n2[n1nm[nn1].sum]=nn2;
        n1nm[nn1].nm[n1nm[nn1].sum]=nm;
        n1nm[nn1].sum++;
    sem_post(&sem_cute);
      }
    }
  } /* end tet4, tet10 */

  //if(elnr==107052) printf(" create elem on base of elem:%d with nodes:%d\n",elnr, anz_n);


  if(anz_n>0)
  {
    /* aussortieren von nodes zwischen zwei element nodes (nref[j])  */
    for(j=0; j<anz_n; j++)
    {
      if(nvalid[j])
      {
        i=0;
        for(k=0; k<anz_n; k++)
        {
	  if(k==j) continue;
          if(nvalid[k])
          {
            if((cutParam->nref[k]==cutParam->npre[j][4])||(cutParam->nref[k]==cutParam->npre[j][5])) i++;
	  }
	}
	if(i==2) nvalid[j]=0;
      }
    }
    
    /* aussortieren der doppelten nodes */
    for(j=0; j<anz_n; j++)
    {
      if(nvalid[j])
      {
        /* gehe ueber alle nodes und deaktiviere nahe nodes */
        for(k=j+1; k<anz_n; k++)
        {
         if(nvalid[k])
         {
          if(dabs(cutParam->npre[j][0]-cutParam->npre[k][0])>local_gtol) goto nexti; 
          if(dabs(cutParam->npre[j][1]-cutParam->npre[k][1])>local_gtol) goto nexti;
          if(dabs(cutParam->npre[j][2]-cutParam->npre[k][2])>local_gtol) goto nexti;
          nvalid[k]=0;
          nexti:;
         }
        }
	//printf("anz_n:%d store index %d\n", cutParam->anz_n, j); 
        ipre[cutParam->anz_n++]=j;
      }
    }

    /* sortiere die nodes nach winkel zum ersten node, masterrichtung ist en, center ist cg */
    cg[0]=0.;
    cg[1]=0.;
    cg[2]=0.;
    for(j=0; j<cutParam->anz_n; j++)
    {
      cg[0]+=cutParam->npre[ipre[j]][0];
      cg[1]+=cutParam->npre[ipre[j]][1];
      cg[2]+=cutParam->npre[ipre[j]][2];
    }
    cg[0]/=cutParam->anz_n;
    cg[1]/=cutParam->anz_n;
    cg[2]/=cutParam->anz_n;

    if ( (rsort = (Rsort *)malloc( (cutParam->anz_n+1) * sizeof(Rsort))) == NULL )
      printf("ERROR: realloc failed: rsort\n\n" ); 
    v_result( cg, &cutParam->npre[ipre[0]][0],  v01);
    rsort[0].r=0.;
    rsort[0].i=ipre[0];
    k=1;
    for(j=1; j<cutParam->anz_n; j++)
    {
      v_result( cg, &cutParam->npre[ipre[j]][0],  v02);
      rsort[k].r=v_angle_ref(v01,v02,en);
      rsort[k].i=ipre[j];
      k++;
    }
    qsort( rsort, cutParam->anz_n, sizeof(Rsort), (void *)compareRsort );
    k=0;
    for(j=0; j<cutParam->anz_n; j++)
    {
      //printf("anz_n:%d store sorted index %d\n", cutParam->anz_n, rsort[k].i); 
      cutParam->ipre[j]=rsort[k].i;
      //printf(" noderef:%f %f\n",cutParam->npre[cutParam->ipre[j]][4],cutParam->npre[cutParam->ipre[j]][5]);
      k++;
    }
    free(rsort);
  }
}



void *thread_cutElement( void *vargp)
{
  int elem;
  typedef struct {
    int thread;
    int anz_e;
    int nset;
    int eset;
    double *en;
    double *pn;
    double *eu;
    double *ev;
    CutParam *cutParam;
  } Threadargs;
  Threadargs *param;

  param=(Threadargs *)vargp;
  while(1)
  {
    sem_wait(&sem_cute);
    elem=glob_cutElement++;
    sem_post(&sem_cute);
    if(elem>=param->anz_e) break;
    //printf("thread:%d elem:%d\n",param->thread, elem);
    //printf("elem:%d ptr:%x\n",e_enqire[elem].nr,&(param->cutParam[e_enqire[elem].nr]));
    cutElement(e_enqire[elem].nr, param->nset, param->eset, param->en, param->pn, param->eu, param->ev, &(param->cutParam[e_enqire[elem].nr]));
  }
  return((void *)1);
}



void pre_cut( int nr, int flag )
{
  int i,ii,j,jj,k,e,ee,nn,nset,eset,anz_n,anz_nmax,anz_e,anz_nset,nel[3];
  int anz_nini, n1,n2,nm;
  static int areanode[3];
  static double p[3][3];
  double v02[3], v01[3], en[3], vn[3], pn[3], eu[3], ev[3], dl, max_dl=0.;
  char  addDispFlagLocal=0;
  CutParam *cutParam;
  int foundNode,elnod,newnodeFlag=0;
  int *elNodeRef;
  static Rsort *rsort=NULL;
  extern int compareRsort();
  int nearParams, nearParam[CUTS],ncorner,nelem;
  double cg[3];

  int anz_cutElems=0;
  typedef struct {
    int enr;
    double x,y,z;
  }CutElem;
  static CutElem *cutElem=NULL;
  static double *orig_x=NULL, *orig_y=NULL, *orig_z=NULL, *sort_x=NULL, *sort_y=NULL, *sort_z=NULL;
  static int *sort_nx=NULL, *sort_ny=NULL, *sort_nz=NULL;

  
  typedef struct {
    int thread;
    int anz_e;
    int nset;
    int eset;
    double *en;
    double *pn;
    double *eu;
    double *ev;
    CutParam *cutParam;
  } Threadargs;
  Threadargs *targ=NULL;
  pthread_t *tid=NULL;
  int nlocalThreads;
  int threads=NTHREADS_MAX;

  if(flag=='v')
  {
    if(anz->l<1)   { printf(" ERROR: no datasets available. This function needs a vector result.\n"); return; }
    if(lcase[cur_lc].ictype[cur_entity]==2)
    {
      /* search all necessary entities */
      if (lcase[cur_lc].icind1[cur_entity]==1)      { entity_v[0]=cur_entity; entity_v[1]=cur_entity+1; entity_v[2]=cur_entity+2; entity_v[3]=-1; v_dim=3; }         
      else if (lcase[cur_lc].icind1[cur_entity]==2) { entity_v[0]=cur_entity-1; entity_v[1]=cur_entity; entity_v[2]=cur_entity+1; entity_v[3]=-1; v_dim=3; }         
      else if (lcase[cur_lc].icind1[cur_entity]==3) { entity_v[0]=cur_entity-2; entity_v[1]=cur_entity-1; entity_v[2]=cur_entity; entity_v[3]=-1; v_dim=3; }         
      else if (lcase[cur_lc].icind1[cur_entity]==0) { entity_v[0]=cur_entity-3; entity_v[1]=cur_entity-2; entity_v[2]=cur_entity-1; entity_v[3]=cur_entity; v_dim=4; }
      else
      {
        /* vector-components not located before the vector-value */
        printf(" ERROR: no vector result was selected\n");
        return;
      }
    }
    else
    {
      printf(" ERROR: no vector result was selected\n");
      return;
    }
  }

  /* neuer Drehpunkt */
  if ((flag=='n')||(flag=='v'))
  {
    center( node[nr].nx, node[nr].ny, node[nr].nz);
    p[qcutCounter][0] = node[nr].nx;
    p[qcutCounter][1] = node[nr].ny; 
    p[qcutCounter][2] = node[nr].nz;
    areanode[qcutCounter]=nr;
  }
  else if (flag=='p')
  {
    center( point[nr].px, point[nr].py, point[nr].pz);
    p[qcutCounter][0] = point[nr].px;
    p[qcutCounter][1] = point[nr].py; 
    p[qcutCounter][2] = point[nr].pz;
    areanode[qcutCounter]=0;
  }
  qcutCounter++;

  if ((qcutCounter == 3)||(flag=='v'))
  {
    qcutCounter=0;
    pre_view("elem off");

    /* when node coordinates were changed to the deformed ones then switch back before they are copied and then switch again */ 
    if(addDispFlag)
    {
      addDispToCoordinates(node);
      // remember to switch back
      addDispFlagLocal=2;
    }

    if (flag=='v')
    {
      /* 1 Node gewaehlt, bestimme die Normale auf der Ebene auf Basis des knotenwertes */
      for(i=0; i<3; i++)
      {
        vn[i]=lcase[cur_lc].dat[entity_v[i]][nr];
        //printf("n:%d entity:%d val:%f\n", nr, entity_v[i],vn[i] ); 
      }
      v_norm( vn, en );
      v02[0]=en[2];
      v02[1]=en[0];
      v02[2]=en[1];
      v_prod( v02, vn, v01 );
      pn[0] = node[nr].nx;
      pn[1] = node[nr].ny; 
      pn[2] = node[nr].nz;
    }
    else
    {
      /* 3 Punkte gewaehlt, berechne die Normale auf der Ebene */
      v_result( &p[0][0], &p[1][0],  v01);
      v_result( &p[0][0], &p[2][0],  v02);
      v_prod( v01,v02, vn );
      v_norm( vn, en );
  
      /* punkt auf der ebene der nicht mit nodes identisch ist */
      for(i=0; i<3; i++) pn[i]=0.;
      for(i=0; i<3; i++)
      {
        pn[0]+= p[i][0];
        pn[1]+= p[i][1];
        pn[2]+= p[i][2];
      }
      for(i=0; i<3; i++) pn[i]/=3.;
      //center( pn[0], pn[1], pn[2]);
  
      /* eu ev der ebene ausgehend von pn */
      v_result( pn, &p[1][0],  v01);
    }

    v_prod( v01, vn, v02 );
    v_norm( v01, eu );
    v_norm( v02, ev );

    /* open a temp set */
    zap("-qcut");
    if( (nset=pre_seta( "-qcut", "i", 0)) <0 ) return;

    /* schneide alle elemente mit der ebene */
    eset= pre_seta( specialset->etmp, "i", 0 );

    /* daten der geschnittenen elemente */
    if ((cutParam=(CutParam *)calloc((anz->emax+1),sizeof(CutParam)) ) == NULL ) { printf(" ERROR: calloc failure\n\n"); return;}
    if ((elNodeRef=(int *)calloc((anz->nmax+1), sizeof(int)) ) == NULL ) { printf(" ERROR: calloc failure\n\n"); return;}

    anz_nmax=anz->nmax;
    /* create a table for all nodes which points to already created midside nodes */
    if ( (n1nm = (N1nm *)malloc( (anz_nmax+1) * sizeof(N1nm))) == NULL )
    { printf("\n\n ERROR in triSplitter: malloc\n\n") ; return; }    
    for (ii=0; ii<=anz->nmax; ii++) n1nm[ii].sum=0;
    for (ii=0; ii<=anz->nmax; ii++) n1nm[ii].n2=n1nm[ii].nm=NULL;
    
#ifdef SEMINIT
    if(sem_init(&sem_cute, 0, 1) < 0) printf("Error in:sem_init\n");
#else
    cptr = sem_open("/sem_cute", O_CREAT, 0644, 1);
    sem_unlink("/sem_cute");
    if((cptr == SEM_FAILED)) printf("Error in:sem_open\n");
#endif    
    
    if(threads>anz->e) { nlocalThreads=anz->e; }
    else nlocalThreads=threads;
    //nlocalThreads=1;
    if ((tid=(pthread_t *)realloc((pthread_t *)tid, nlocalThreads*sizeof(pthread_t)) ) == NULL ) { printf(" ERROR: malloc failure\n\n"); return;}
    if ((targ=(Threadargs *)realloc((Threadargs *)targ, nlocalThreads*sizeof(Threadargs)) ) == NULL ) { printf(" ERROR: malloc failure\n\n"); return;}
    glob_cutElement=0;
    anz_e=anz->e;
    // Warning: The following printf's are needed to prevent a crash for successive cuts!
    printf("\nplease wait, new elements in this section will be created\n");
    for(i=0; i<nlocalThreads; i++)
    {
      targ[i].thread=i;
      targ[i].anz_e=anz_e;
      targ[i].nset=nset;
      targ[i].eset=eset;
      targ[i].en=en;
      targ[i].pn=pn;
      targ[i].eu=eu;
      targ[i].ev=ev;
      targ[i].cutParam=cutParam;
      pthread_create(&tid[i],NULL,thread_cutElement,(void *)&targ[i]);
    }
    printf("\n");
    for(i=0; i<nlocalThreads; i++)
    {
      pthread_join(tid[i], NULL);
    }
    printf("\n");
    glob_cutElement=0;
    free(tid); tid=NULL;
    free(targ); targ=NULL;
    printf("\n");
#ifdef SEMINIT
    if(sem_destroy(&sem_cute) < 0) printf("Error in:sem_init\n");
#else
    if(sem_close(cptr) < 0) printf("Error in:sem_close\n");
#endif    
    printf("\n");


    anz_nini=anz->nmax;
    if( (checknode = (Checknode *)malloc( (anz_nini+1) * sizeof(Checknode))) == NULL )
      printf("\n\n ERROR: malloc failure\n\n");
    for(i=0; i<=anz_nini; i++)
    {
      checknode[i].n2=NULL;
      checknode[i].ns=NULL;
      checknode[i].nn=0;
    }
    for(i=1; i<=anz->emax; i++) if(cutParam[i].anz_n<3) cutParam[i].skip=1;

    /* generiere neue nodes */
    if ( (rsort = (Rsort *)malloc( (CUTS) * sizeof(Rsort))) == NULL )
      printf("ERROR: realloc failed: Rsort\n\n" ); 
    for(e=0; e<anz_e; e++)
    {
      i=e_enqire[e].nr;
      //printf("e:%d anz_n:%d\n",e_enqire[e].nr, cutParam[i].anz_n); 
      if(!cutParam[i].skip)
      {
        for(j=0; j<cutParam[i].anz_n; j++)
        {
	  // generate a new node only if the referenced corner nodes were not used before
	  foundNode=0;
	  // printf(" npre %f %f\n",cutParam[i].npre[cutParam[i].ipre[j]][4],cutParam[i].npre[cutParam[i].ipre[j]][5]);
	  if(cutParam[i].npre[cutParam[i].ipre[j]][4]>cutParam[i].npre[cutParam[i].ipre[j]][5])
	  {
	    n1=cutParam[i].npre[cutParam[i].ipre[j]][4];
	    n2=cutParam[i].npre[cutParam[i].ipre[j]][5];
	  }
	  else
	  {
	    n2=cutParam[i].npre[cutParam[i].ipre[j]][4];
	    n1=cutParam[i].npre[cutParam[i].ipre[j]][5];
	  }
	  for(k=0; k<checknode[n1].nn; k++)
	  {
	    if(checknode[n1].n2[k]== n2)
	    {
	      //printf("k:%d n1:%d n2:%d\n",k,n1,n2);
	      foundNode=checknode[n1].ns[k];
	      break;
	    }
	  }
	  //printf("foundNode:%d\n",foundNode);
	  if(!foundNode)
	  {
	    // do not create a new node if cutParam[i].nref[cutParam[i].ipre[j]] references a valid element node
	    elnod=cutParam[i].nref[cutParam[i].ipre[j]];
	    if(elnod==0)  // no referenced element corner node
	    {
              nn=nod( anz, &node, 0, anz->nnext++, cutParam[i].npre[cutParam[i].ipre[j]][0], cutParam[i].npre[cutParam[i].ipre[j]][1], cutParam[i].npre[cutParam[i].ipre[j]][2], 0 );
	      newnodeFlag=1;
	    }
	    else if(elNodeRef[elnod]==0)
	    {
              nn=nod( anz, &node, 0, anz->nnext++, cutParam[i].npre[cutParam[i].ipre[j]][0], cutParam[i].npre[cutParam[i].ipre[j]][1], cutParam[i].npre[cutParam[i].ipre[j]][2], 0 );
	      elNodeRef[elnod]=nn;
	      //printf("node %d will reference new node:%d\n",elnod, node[nn].nr);
	      newnodeFlag=1;
	    }
	    else
	    {
	      nn=elNodeRef[elnod]; newnodeFlag=0;
	      //printf("node %d references new node:%d\n",elnod, node[nn].nr);
	    }
	    
	    if(newnodeFlag)
	    {
  	      //printf("j:%d make node:%d from index %d\n",j, node[nn].nr, cutParam[i].ipre[j]); 
              seta(setall, "n", node[nn].nr);
              seta(nset, "n", node[nn].nr);
              anz_nset=set[nset].anz_n-1;
  
              /* store the informations to interpolate node-values for the new nodes */
              if ( (qcut_nod = (Qcut_nodes *)realloc(qcut_nod, (anz_nset+1) * sizeof(Qcut_nodes))) == NULL )
                printf("\n\n ERROR: realloc failure qcut_nod \n\n" );	               
  	      //printf("qcut_nod:%d nod:%d\n",anz_nset, cutParam[i].ipre[j]); 
              qcut_nod[anz_nset].nr=node[nn].nr;
              qcut_nod[anz_nset].n1=cutParam[i].npre[cutParam[i].ipre[j]][6];
              qcut_nod[anz_nset].n2=cutParam[i].npre[cutParam[i].ipre[j]][7];
              qcut_nod[anz_nset].val=cutParam[i].npre[cutParam[i].ipre[j]][3];
	    }
	    
            // use of nnew
            cutParam[i].nnew[j]=node[nn].nr;

            if( (checknode[n1].n2 = (int *)realloc( checknode[n1].n2, (checknode[n1].nn+1) * sizeof(int))) == NULL )
              printf("\n\n ERROR: realloc failure\n\n");
            if( (checknode[n1].ns = (int *)realloc( checknode[n1].ns, (checknode[n1].nn+1) * sizeof(int))) == NULL )
              printf("\n\n ERROR: realloc failure\n\n");
	    checknode[n1].n2[checknode[n1].nn]=n2;
	    checknode[n1].ns[checknode[n1].nn]=node[nn].nr;
	    checknode[n1].nn++;

	    // distance between the corner nodes as a measure for the merging tolerance
	    v_result(&node[n1].nx,&node[n2].nx,v01);
	    dl=v_betrag(v01);
	    if(dl>max_dl) max_dl=dl;
	  }
	  else cutParam[i].nnew[j]=foundNode;
        }

	// avoid identical nodes
	// do not change the order in cutParam[i].nnew
	for( j=0; j<cutParam[i].anz_n; j++)
	{
          rsort[j].r=cutParam[i].nnew[j];
          rsort[j].i=j;
        }
        qsort( rsort, cutParam[i].anz_n, sizeof(Rsort), (void *)compareRsort );
	for( j=0; j<cutParam[i].anz_n-1; j++)
	{
          if(rsort[j].r==rsort[j+1].r) cutParam[i].nnew[rsort[j].i]=0;
	}

	anz_n=0;
	for( j=0; j<cutParam[i].anz_n; j++)
	{
          if(cutParam[i].nnew[j]==0) continue;
	  cutParam[i].nnew[anz_n++]=cutParam[i].nnew[j];
	}
	cutParam[i].anz_n=anz_n;
	if(cutParam[i].anz_n<3)  { cutParam[i].skip=1; continue; }

        if( (cutElem = (CutElem *)realloc((CutElem *)cutElem, (anz_cutElems+1) * sizeof(CutElem))) == NULL )
          printf("ERROR: realloc failed in areampc\n\n" ); 
        cutElem[anz_cutElems].x=cutElem[anz_cutElems].y=cutElem[anz_cutElems].z=0.;
	cutElem[anz_cutElems].enr=i;
	for( j=0; j<cutParam[i].anz_n; j++)
	{
          cutElem[anz_cutElems].x+=node[cutParam[i].nnew[j]].nx;	
          cutElem[anz_cutElems].y+=node[cutParam[i].nnew[j]].ny;	
          cutElem[anz_cutElems].z+=node[cutParam[i].nnew[j]].nz;
	}
        cutElem[anz_cutElems].x/=cutParam[i].anz_n;	
        cutElem[anz_cutElems].y/=cutParam[i].anz_n;	
        cutElem[anz_cutElems].z/=cutParam[i].anz_n;
	cutParam[i].cgx=cutElem[anz_cutElems].x;
	cutParam[i].cgy=cutElem[anz_cutElems].y;
	cutParam[i].cgz=cutElem[anz_cutElems].z;
	anz_cutElems++;
      }
    }

    /* search the closest cutParam */
    if ( (rsort = (Rsort *)malloc( (anz_cutElems+1) * sizeof(Rsort))) == NULL )
      printf("ERROR: realloc failed: Rsort\n\n" ); 
    if ( (orig_x = (double *)malloc( (anz_cutElems+1) * sizeof(double))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (orig_y = (double *)malloc( (anz_cutElems+1) * sizeof(double))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (orig_z = (double *)malloc( (anz_cutElems+1) * sizeof(double))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (sort_x = (double *)malloc( (anz_cutElems+1) * sizeof(double))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (sort_y = (double *)malloc( (anz_cutElems+1) * sizeof(double))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (sort_z = (double *)malloc( (anz_cutElems+1) * sizeof(double))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (sort_nx = (int *)malloc( (anz_cutElems+1) * sizeof(int))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (sort_ny = (int *)malloc( (anz_cutElems+1) * sizeof(int))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    if ( (sort_nz = (int *)malloc( (anz_cutElems+1) * sizeof(int))) == NULL )
      printf("ERROR: realloc failed in areampc\n\n" ); 
    for(i=0; i<anz_cutElems; i++)
    {
      rsort[i].r=orig_x[i]=cutElem[i].x;
      rsort[i].i=i;
    }
    qsort( rsort, anz_cutElems, sizeof(Rsort), (void *)compareRsort );
    for(i=0; i<anz_cutElems; i++)
    {
      sort_x[i]=rsort[i].r;
      sort_nx[i]=rsort[i].i;
    }
    for(i=0; i<anz_cutElems; i++)
    {
      rsort[i].r=orig_y[i]=cutElem[i].y;
      rsort[i].i=i;
    }
    qsort( rsort, anz_cutElems, sizeof(Rsort), (void *)compareRsort );
    for(i=0; i<anz_cutElems; i++)
    {
      sort_y[i]=rsort[i].r;
      sort_ny[i]=rsort[i].i;
    }
    for(i=0; i<anz_cutElems; i++)
    {
      rsort[i].r=orig_z[i]=cutElem[i].z;
      rsort[i].i=i;
    }
    qsort( rsort, anz_cutElems, sizeof(Rsort), (void *)compareRsort );
    for(i=0; i<anz_cutElems; i++)
    {
      sort_z[i]=rsort[i].r;
      sort_nz[i]=rsort[i].i;
    }
    if(anz_cutElems<CUTS) nearParams=anz_cutElems; else nearParams=CUTS;
	
    /* check if double element-faces exist at the same spot */
    for(e=0; e<anz_e; e++)
    {
      i=e_enqire[e].nr;
      //if(i==107052) printf(" create elem on base of elem:%d with nodes:%d\n",i, cutParam[i].anz_n);
      if(!cutParam[i].skip)
      {
	// first check if a second cutParam references at least the same nodes, then set skip=1 and continue
        near3d(orig_x,orig_y,orig_z,sort_x,sort_y,sort_z,sort_nx,sort_ny,sort_nz, cutParam[i].cgx, cutParam[i].cgy,
            cutParam[i].cgz, anz_cutElems, &nearParam[0], nearParams);
        for(ee=0; ee<nearParams; ee++)
        {
          ii=cutElem[nearParam[ee]].enr;
	  if(ii==i) continue;
	  if((!cutParam[ii].skip)&&(e_enqire[i].type == e_enqire[ii].type))
	  {
	    if(cutParam[ii].anz_n>=cutParam[i].anz_n)
	    {
	      nn=0;
	      for(j=0; j<cutParam[i].anz_n; j++)
	      {
		for(k=0; k<cutParam[ii].anz_n; k++) if(cutParam[i].nnew[j]==cutParam[ii].nnew[k]) nn++;
	      }
	      if(nn==cutParam[ii].anz_n)
	      {
		cutParam[ii].skip=1;
		//printf("found cutParam with the same nodes\n");
		break;
	      }
	      else if(nn>2)
	      {
		//printf("found cutParam with %d matching nodes\n",nn);
	        // look if the two cutParam add up to a full element face
	        // this is the case if all corner nodes of a face are included in the combined nodes of both cutParam
		// for simplification it is sufficient if the number of face nodes are matched, the single faces are not checked.
		nelem=ncorner=0;
                if ((e_enqire[ii].type == 1)||(e_enqire[ii].type == 4)) { ncorner=4; nelem=8; } /* HEXA8,20 */
		else if ((e_enqire[ii].type == 3)||(e_enqire[ii].type == 6)) { ncorner=3; nelem=4; } /* TET4,10 */
		nn=0;
                for(j=0; j<nelem; j++)
		{
		  for(k=0; k<cutParam[ii].anz_n; k++) if(e_enqire[ii].nod[j]==cutParam[ii].nref[k]) nn++;
		  for(k=0; k<cutParam[i].anz_n; k++) if(e_enqire[i].nod[j]==cutParam[i].nref[k]) nn++;
		}
		if(nn>=ncorner-1)
		{
		  cutParam[ii].skip=1;
		  j=cutParam[i].anz_n;
		  for(k=0; k<cutParam[ii].anz_n; k++) cutParam[i].nnew[k+j]=cutParam[ii].nnew[k];
		  cutParam[i].anz_n+=cutParam[ii].anz_n;

                  /* sortiere die nodes nach winkel zum ersten node, masterrichtung ist en, center ist cg */
                  cg[0]=0.;
                  cg[1]=0.;
                  cg[2]=0.;
                  for(j=0; j<cutParam[i].anz_n; j++)
                  {
                    cg[0]+=node[cutParam[i].nnew[j]].nx;
                    cg[1]+=node[cutParam[i].nnew[j]].ny;
                    cg[2]+=node[cutParam[i].nnew[j]].nz;
                  }
                  cg[0]/=cutParam[i].anz_n;
                  cg[1]/=cutParam[i].anz_n;
                  cg[2]/=cutParam[i].anz_n;
                
                  if ( (rsort = (Rsort *)malloc( (cutParam[i].anz_n+1) * sizeof(Rsort))) == NULL )
                    printf("ERROR: realloc failed: rsort\n\n" ); 
                  v_result( cg, &node[cutParam[i].nnew[0]].nx,  v01);
                  rsort[0].r=0.;
                  rsort[0].i=cutParam[i].nnew[0];
                  k=1;
                  for(j=1; j<cutParam[i].anz_n; j++)
                  {
                    v_result( cg, &node[cutParam[i].nnew[j]].nx,  v02);
                    rsort[k].r=v_angle_ref(v01,v02,en);
                    rsort[k].i=cutParam[i].nnew[j];
                    k++;
                  }
                  qsort( rsort, cutParam[i].anz_n, sizeof(Rsort), (void *)compareRsort );
                  k=0;
		  cutParam[i].nnew[0]=rsort[0].i;
                  for(j=1; j<cutParam[i].anz_n; j++)
                  {
                    if(cutParam[i].nnew[k]!=rsort[j].i)
		    {
                      k++;
		      cutParam[i].nnew[k]=rsort[j].i;
		    }
                  }
		  cutParam[i].anz_n=k+1;
                  free(rsort);
		  //printf("found two cutParam belonging to the same elem face\n");
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
		  
    /* generiere neue elemente */
    for(e=0; e<anz_e; e++)
    {
      i=e_enqire[e].nr;
      //if(i==107052) printf(" create elem on base of elem:%d with nodes:%d\n",i, cutParam[i].anz_n);
      if(!cutParam[i].skip)
      {
	// add the areanods to the face
        if ((flag=='n')||(flag=='v'))
        {
          // check if the nm matches one areanode, then generate a node there
	  ee=0;
          k=cutParam[i].anz_n;
          for(j=0; j<k; j++)
          {
            n1=cutParam[i].nref[cutParam[i].ipre[j]];   // correlates with cutParam[i].nnew[j]
            if(j==cutParam[i].anz_n-1) n2=cutParam[i].nref[cutParam[i].ipre[0]];
            else n2=cutParam[i].nref[cutParam[i].ipre[j+1]];
            nm=-1;
            for(ii=0; ii<n1nm[n1].sum; ii++) if(n1nm[n1].n2[ii]==n2) nm=n1nm[n1].nm[ii];
            for(ii=0; ii<n1nm[n2].sum; ii++) if(n1nm[n2].n2[ii]==n1) nm=n1nm[n2].nm[ii];
            if(nm!=-1)
            {
              for(ii=0; ii<3; ii++)
              {
	        if(nm==areanode[ii])
	        {
		  //printf(" add aditional node to -qcut based on areanode:%d\n",areanode[ii]);
                  for(jj=cutParam[i].anz_n-1; jj>j+ee; jj--) cutParam[i].nnew[jj+1]=cutParam[i].nnew[jj];
                  nn=nod( anz, &node, 0, anz->nnext++, node[nm].nx, node[nm].ny, node[nm].nz, 0 );
	            elNodeRef[elnod]=nn;
                  seta(setall, "n", node[nn].nr);
                  seta(nset, "n", node[nn].nr);
                  anz_nset=set[nset].anz_n-1;
		    ee++;

                  /* store the informations to interpolate node-values for the new nodes */
                  if ( (qcut_nod = (Qcut_nodes *)realloc(qcut_nod, (anz_nset+1) * sizeof(Qcut_nodes))) == NULL )
                    printf("\n\n ERROR: realloc failure qcut_nod \n\n" );	               
	            //printf("qcut_nod:%d nod:%d\n",anz_nset, cutParam[i].ipre[j]); 
                  qcut_nod[anz_nset].nr=node[nn].nr;
                  qcut_nod[anz_nset].n1=nm;
                  qcut_nod[anz_nset].n2=nm;
                  qcut_nod[anz_nset].val=1.;
                  cutParam[i].nnew[j+ee]=node[nn].nr;
		  cutParam[i].anz_n++;
	          //for(jj=0; jj<cutParam[i].anz_n; jj++) printf("%d ",cutParam[i].nnew[e]); printf("\n");
	        }
	      }
            }
	  }
        }
	  
        for(j=1; j<cutParam[i].anz_n-1; j++)
        {
          nel[0]=cutParam[i].nnew[0];
          nel[1]=cutParam[i].nnew[j];
          nel[2]=cutParam[i].nnew[j+1];
	  // create only elements with an area >0
	  v_result( &node[nel[0]].nx, &node[nel[1]].nx, v01);
	  v_result( &node[nel[0]].nx, &node[nel[2]].nx, v02);
	  v_prod( v01, v02, vn );
	  if(abs(v_betrag(vn)) > abs( (v_betrag(v01)*v_betrag(v02)*1e-10) ) )
	  {
	    // printf("j:%d make elem:%d\n",j, anz->enext); 
            elem_define(anz,&e_enqire, anz->enext++, 7, nel, 0, 0 );
            seta(setall, "e", anz->enext-1);
            seta(nset, "e", anz->enext-1);
            seta(eset, "e", anz->enext-1);
	  }
        }
      }
    }

    for (i=0; i<=anz_nmax; i++) { free(n1nm[i].n2); free(n1nm[i].nm); }
    free(n1nm);
    for (i=0; i<3; i++) areanode[i]=0;
    free(cutParam);
    free(elNodeRef);
    for(i=0; i<=anz_nini; i++)
    {
      free(checknode[i].n2);
      free(checknode[i].ns);
    }
    free(checknode);
    checknode=NULL;

    // remove double nodes and elements
    //sprintf(buffer,"n %s %e",set[nset].name, max_dl*CUTELEMTOL*2.);
    //printf("merg n %s %e",set[nset].name, max_dl*CUTELEMTOL*2.);
    //pre_merge(buffer);
    mergeElem(nset,0);
    illumResultFlag=0;

    /* update Dataset */
    updLcase(cur_lc, nset);

    /* zeige neue elemente */
    makeSurfaces();
    realloc_colNr();        
    updateDispLists(); 

    /* when node coordinates were changed to the deformed ones then switch back  */ 
    if(addDispFlagLocal==2)
    {
      addDispToCoordinates(node);
    }

    if(vectorFlag) pre_view("vector off");
    if(anz->l) plot( "ev -qcut" );
    else plot( "ei -qcut" );
  }
}


void oriAllSurfs(int sur)
{
  int i,j,n,l,ll,lll,cl,s, prod1, prod2,oriflag,counter=0,surl;
  int **ltos, *sori;

  /* go over all surfs.*/
  /* check if one neighbour surf is oriented */
  /* then orient the surf */

  /* first go over all lines and determine all related surfs (should be 2) */
  /* store the surfs in an array which points then to the surf */
  /*  this will be set to -1 if the surf is oriented */
  /*  if all surfs are oriented this array contains only -1 */

  /* relate all surfs to its lines */
  if( (ltos=(int **)malloc((anzGeo->l+1)*sizeof(int *) ) )==NULL)
  printf("ERROR malloc failed in oriAllSurfs()\n");
  for(i=0; i<anzGeo->l; i++)
  {
    if( (ltos[i]=(int *)malloc((3)*sizeof(int) ) )==NULL)
    printf("ERROR malloc failed in oriAllSurfs()\n");
     ltos[i][0]=0; for(j=1;j<3;j++) ltos[i][j]=-1;
  }

  for(s=0; s<anzGeo->s; s++)
  {
    if(surf[s].name!=NULL) for(j=0; j<surf[s].nl; j++)
    {
      if(surf[s].typ[j]=='l')
      {
        n=++ltos[surf[s].l[j]][0];
        if(n>2)
        {
          printf("ERROR: to many related surfs(%d) for line:%s\n", n, line[surf[s].l[j]].name);
          //printf("No inner surfaces are permitted. Command could not be executed\n");
          //return;
          //printf("No inner surfaces are permitted. Some surfs can not be oriented\n");
          ltos[surf[s].l[j]][2]=-1;
        }
        else ltos[surf[s].l[j]][n]=s;
      }
      else
      {
        cl=surf[s].l[j];
        for(l=0; l<lcmb[cl].nl; l++)
	{
          n=++ltos[lcmb[cl].l[l]][0];
          if(n>2)
          {
            printf("ERROR: to many related surfs(%d) for line:%s\n", n, line[lcmb[cl].l[l]].name);
            //printf("No inner surfaces are permitted. Command could not be executed\n");
            //return;
            //printf("No inner surfaces are permitted. Some surfs can not be oriented\n");
            ltos[lcmb[cl].l[l]][2]=-1;
          }
          else ltos[lcmb[cl].l[l]][n]=s;
	}
      }
    }
  }


  for(i=0; i<anzGeo->l; i++) if( line[i].name != (char *)NULL )
  {
    printf("l:%s ", line[i].name);
    if(ltos[i][1]>-1) printf("surf:%s ", surf[ltos[i][1]].name);
    if(ltos[i][2]>-1) printf("surf:%s ", surf[ltos[i][2]].name);
    printf("\n ");
  }

  /* create a link between surfindex and surface-ori */
  /* the sori is "0" as long a surf is not oriented */
  if( (sori=(int *)malloc((anzGeo->s+1)*sizeof(int) ) )==NULL)
    printf("ERROR malloc failed in oriAllSurfs()\n");
  for(i=0; i<anzGeo->s; i++) sori[i]=0;
  
  /* the start-surface is already oriented */
  sori[sur]=1; 
 
  /* go over all surfs and look if one has an oriented neighbour */
 more:;
  oriflag=0;
  for(s=0; s<anzGeo->s; s++)
  {
    /* if the surf is valid and not oriented go over all its lines */
    if((surf[s].name!=NULL)&&(sori[s]==0)) for(j=0; j<surf[s].nl; j++)
    {
      oriflag=1;
      if(surf[s].typ[j]=='l')
      {
        /* check the connected surfs based on the common lines if it is an oriented one */
        for(n=1;n<3;n++) if(ltos[surf[s].l[j]][n]>-1) if((ltos[surf[s].l[j]][n]!=s)&&(sori[ltos[surf[s].l[j]][n]]>0))
	{
	  sur=ltos[surf[s].l[j]][n];
          surl= surf[s].l[j];
          printf("surf:%s line:%s oriented surf:%s\n", surf[s].name, line[surl].name, surf[sur].name);

          /* check if the surf must be inverted */
          /* based on the product of orientations of the oriented surf */
          /* determine the index of the connected line in sur */
          if(surf[sur].ori=='+') prod1=1; else prod1=-1;
          for(ll=0; ll<surf[sur].nl; ll++)
	  {
            if(surf[sur].typ[ll]=='l') 
            {
	      if(surl==surf[sur].l[ll])
              {
                if(surf[sur].o[ll]=='+') prod1*=1; else prod1*=-1;
                goto found1;
              }
	    }
	    else 
            {
              for(lll=0; lll<lcmb[surf[sur].l[ll]].nl; lll++) if(surl==lcmb[surf[sur].l[ll]].l[lll])
              {
                if(surf[sur].o[ll]=='+') prod1*=1; else prod1*=-1;
                if(lcmb[surf[sur].l[ll]].o[lll]=='+') prod1*=1; else prod1*=-1;
                goto found1;
	      }
	    }
	  }
	found1:;

          /* product of orientations of the actual surf */
          if(surf[s].ori=='+') prod2=1; else prod2=-1;
          if(surf[s].o[j]=='+') prod2*=1; else prod2*=-1;

          sori[s]=1;
          if(prod2==prod1)
	  {
            flipSurfori(s);
	  }
          goto new_surf;
	}
      }
      else
      {
        cl=surf[s].l[j];
        for(l=0; l<lcmb[cl].nl; l++)
	{
          /* check the connected surfs based on the common lines if it is an oriented one */
          for(n=1;n<3;n++) if(ltos[lcmb[cl].l[l]][n]>-1) if((ltos[lcmb[cl].l[l]][n]!=s)&&(sori[ltos[lcmb[cl].l[l]][n]]>0))
	  {
	    sur=ltos[lcmb[cl].l[l]][n];
            printf("surf:%s lcmb:%s line:%s oriented surf:%s\n", surf[s].name, lcmb[cl].name, line[lcmb[cl].l[l]].name, surf[sur].name);
            surl= lcmb[cl].l[l];

            /* check if the surf must be inverted */
            /* based on the product of orientations of the oriented surf */
            /* determine the index of the connected line in sur */
            if(surf[sur].ori=='+') prod1=1; else prod1=-1;
            for(ll=0; ll<surf[sur].nl; ll++)
            {
              if(surf[sur].typ[ll]=='l') 
              {
		if(surl==surf[sur].l[ll])
                {
                  if(surf[sur].o[ll]=='+') prod1*=1; else prod1*=-1;
                  goto found2;
                }
   	      }
   	      else 
              {
                for(lll=0; lll<lcmb[surf[sur].l[ll]].nl; lll++) 
                if(surl==lcmb[surf[sur].l[ll]].l[lll])
                {
                  if(surf[sur].o[ll]=='+') prod1*=1; else prod1*=-1;
                  if(lcmb[surf[sur].l[ll]].o[lll]=='+') prod1*=1; else prod1*=-1;
                  goto found2;
   	        }
   	      }
	    }
	   found2:;
 
            /* product of orientations of the actual surf */
            if(surf[s].ori=='+') prod2=1; else prod2=-1;
            if(surf[s].o[j]=='+') prod2*=1; else prod2*=-1;
            if(lcmb[cl].o[l]=='+') prod2*=1; else prod2*=-1;

            sori[s]=1;
            if(prod2==prod1)
	    {
              flipSurfori(s);
	    }
            goto new_surf;
	  }
	}
      }
    }
    new_surf:;    
  }
  if(oriflag)
  {
    counter++;
    if(counter<anzGeo->s) goto more;
    else printf(" WARNING: too much loops. Some surfs might be still unoriented.\n");
  }
}


void flip( char *type, int e)
{
  int i,j,s=0;
  int nbuf[20];

  if (type[0]=='b')
  {
    /* switch body */
    if(body[e].ori=='+') body[e].ori='-'; else body[e].ori='+';
  } 
  else if (type[0]=='s')
  {
    /* switch surface and elements and embedded faces */
    /* identify the surf */
    s=e;

    /* check if all connected surfs should be oriented according to the identified one */
    if(type[1]=='a')
    {
      oriAllSurfs(s);
    }
    else flipSurfori(s);
  }
  /*
  else if ((type[0]=='e')&&(anzGeo->s>0)&&(selem!=NULL))
  {
    // warning: selem is allocated in qflp 
    // switch surface and elements and embedded faces 
    // identify the surf 
    s=selem[e];

    // check if all connected surfs should be oriented according to the identified one 
    if(type[1]=='a')
    {
      oriAllSurfs(s);
    }
    else flipSurfori(s);
  }
  */
  else if (type[0]=='e')
  {
    /* switch orient of the single element */
    if (e_enqire[e].type == 7)
    {
      nbuf[1]=e_enqire[e].nod[1];
      nbuf[2]=e_enqire[e].nod[2];
      e_enqire[e].nod[2]=nbuf[1];
      e_enqire[e].nod[1]=nbuf[2];
    }  
    else if (e_enqire[e].type == 8)
    {
      nbuf[0]=e_enqire[e].nod[0];
      nbuf[1]=e_enqire[e].nod[1];
      nbuf[4]=e_enqire[e].nod[4];
      nbuf[5]=e_enqire[e].nod[5];
      e_enqire[e].nod[0]=nbuf[1];
      e_enqire[e].nod[1]=nbuf[0];
      e_enqire[e].nod[4]=nbuf[5];
      e_enqire[e].nod[5]=nbuf[4];
    }  
    else if (e_enqire[e].type == 9)
    {
      nbuf[1]=e_enqire[e].nod[1];
      nbuf[3]=e_enqire[e].nod[3];
      e_enqire[e].nod[3]=nbuf[1];
      e_enqire[e].nod[1]=nbuf[3];
    }  
    else if (e_enqire[e].type == 10)
    {
      nbuf[1]=e_enqire[e].nod[1];
      nbuf[3]=e_enqire[e].nod[3];
      nbuf[4]=e_enqire[e].nod[4];
      nbuf[5]=e_enqire[e].nod[5];
      nbuf[6]=e_enqire[e].nod[6];
      nbuf[7]=e_enqire[e].nod[7];
      e_enqire[e].nod[3]=nbuf[1];
      e_enqire[e].nod[1]=nbuf[3];
      e_enqire[e].nod[7]=nbuf[4];
      e_enqire[e].nod[6]=nbuf[5];
      e_enqire[e].nod[5]=nbuf[6];
      e_enqire[e].nod[4]=nbuf[7];
    }  
    else  
    {
      printf(" ERROR: element type %d not known\n", e_enqire[e].type);
    }

    for (i=0; i<anz->f; i++)
    { 
      if(face[i].elem_nr==e)
      {
        if (face[i].type == 7)
        {
          nbuf[1]=face[i].nod[1];
          nbuf[2]=face[i].nod[2];
          face[i].nod[2]=nbuf[1];
          face[i].nod[1]=nbuf[2];
        }  
        else if (face[i].type == 8)
        {
          nbuf[0]=face[i].nod[0];
          nbuf[1]=face[i].nod[1];
          nbuf[4]=face[i].nod[4];
          nbuf[5]=face[i].nod[5];
          face[i].nod[0]=nbuf[1];
          face[i].nod[1]=nbuf[0];
          face[i].nod[4]=nbuf[5];
          face[i].nod[5]=nbuf[4];
        }  
        else if (face[i].type == 9)
        {
          nbuf[1]=face[i].nod[1];
          nbuf[3]=face[i].nod[3];
          face[i].nod[3]=nbuf[1];
          face[i].nod[1]=nbuf[3];
        }  
        else if (face[i].type == 10)
        {
          nbuf[1]=face[i].nod[1];
          nbuf[3]=face[i].nod[3];
          nbuf[4]=face[i].nod[4];
          nbuf[5]=face[i].nod[5];
          nbuf[6]=face[i].nod[6];
          nbuf[7]=face[i].nod[7];
          face[i].nod[3]=nbuf[1];
          face[i].nod[1]=nbuf[3];
          face[i].nod[7]=nbuf[4];
          face[i].nod[6]=nbuf[5];
          face[i].nod[5]=nbuf[6];
          face[i].nod[4]=nbuf[7];
        }
        if ((face[i].type == 11)||(face[i].type == 12))
        {
          nbuf[1]=face[i].nod[0];
          nbuf[2]=face[i].nod[1];
          face[i].nod[1]=nbuf[1];
          face[i].nod[0]=nbuf[2];
  	  j=face[i].elem_nr;
  	  if((face[i].nod[0]==e_enqire[j].nod[0])&&(face[i].nod[1]==e_enqire[j].nod[1])) face[i].nr=2;
  	  else if((face[i].nod[0]==e_enqire[j].nod[1])&&(face[i].nod[1]==e_enqire[j].nod[2])) face[i].nr=3;
  	  else if((face[i].nod[0]==e_enqire[j].nod[2])&&(face[i].nod[1]==e_enqire[j].nod[3])) face[i].nr=4;
  	  else if((face[i].nod[0]==e_enqire[j].nod[3])&&(face[i].nod[1]==e_enqire[j].nod[0])) face[i].nr=5;
        }
      }
    }
      
  } 
  getFaceNormalen( face, node, anz );
  getElemNormalen( e_enqire, node, anz->e );
  updateDispLists();
}
