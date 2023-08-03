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

#define MAX_GSUR_PARAMERTER 8
extern Scale     scale[1];
extern Elements  *e_enqire;     /* elem-array by elem-number instead of elem[index]... */
extern Summen    anz[1];
extern Edges     *edge;
extern Nodes     *node;
extern Datasets *lcase;
extern Faces     *face;
extern NodeBlocks *nBlock;
extern BGpicture *bgpicture;

extern Meshp meshp;

extern Alias     *alias;
extern Sets      *set;
extern Shapes    *shape;
extern Materials *material; 
extern Psets     *pset;
extern Values    *value;
extern Points    *point;
extern Lines     *line;
extern Lcmb      *lcmb;
extern Gsur      *surf;
extern Gbod      *body;
extern Nurbl     *nurbl;
extern Nurbs     *nurbs;
extern Materials *material; 
extern Amplitudes *amplitude; 
extern SumGeo    anzGeo[1];
extern SumAsci   sumAsci[1];
extern Texts     *ntext;

extern SpecialSet specialset[1];
extern Eqal eqal;
extern int entitycols;
extern Entitycol *entitycol;
extern Psets     *pset;


extern double dtx, dty, dtz, ds;  /*Verschiebungen */
extern double centerPnt[3];                                /* Rotationszentrum */
extern double lastquat[4];                                 /* letzte Matrix aus Trackball*/
extern GLdouble R[4][4];                                   /* Rotationsmatrix */
extern double vmem[4];                  /* kor. bis auswahl eines neuen drehpkts */
extern int width_w0, height_w0;
extern Scale scale[1];
extern int backgrndcol,legend_font,draw_font,menu_font;
extern char modelEdgeFlag,elemEdgeFlag,surfFlag,illumResultFlag;
extern char  inpformat;       /* defines the start-up mode of cgx */
extern int   cur_entity;                                       /* aktive entity (component), entity in menu= cur_entity+1*/
extern int   entity_v[6];                                         /* components of a vector-entity */
extern double v_factor;                                    /* scaling-factor for the vectors in the vector-plot */
extern double v_scale;                                     /* additional scaling-factor for the vectors */
extern int       pre_lc;                                      /* pre-selected Dataset, active after entity is selected */
extern int       cur_lc;                                          /* aktive Dataset */
extern char  drawMode;               /* protokoliert drawFunktion (Load=1,Light=2,Animate=3,preprocessor=4, vector=5)*/
extern int   steps;                                  /* Schrittweite der Farbscala */


int writefbd(char *setname, char *format )
{
  FILE *handle1;
  int  i, j, k, n, l, se, nr;
  char string[MAX_LINE_LENGTH], datout[MAX_LINE_LENGTH], buffer[MAX_LINE_LENGTH];
  int length, params, setNr, setBuf, nsave;
  int bias_fbd;
  double p1[3], p2[3];
  Points *pnt;
  pnt=point;

  setNr=getSetNr(setname);
  if (setNr<0)
  {
    printf (" ERROR: set:%s does not exist\n", setname);
    return (-1);
  }

  length=sword( setname, string);
  string[length]='.';
  string[length+1]='f';
  string[length+2]='b';
  string[length+3]='d';
  string[length+4]='\0';
  strcpy( datout, string);
  /* Open the files and check to see that it was opened correctly */
  handle1 = fopen (datout, "w+b");
  if (handle1==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n", datout); return(-1);}
  else  printf ("\n%s opened\n",datout);
  fprintf (handle1,"# written by cgx_%s\n", VERSION);
  printf ("\n write fbd  \n");
  delSet(specialset->tmp);

  if( (setBuf=pre_seta("-setbuf", "i", 0)) <0) return(-1);
  for (j=0; j<set[setNr].anz_v; j++) seta( setBuf, "v", set[setNr].valu[j] );
  for (j=0; j<set[setNr].anz_n; j++) seta( setBuf, "n", set[setNr].node[j] );
  for (j=0; j<set[setNr].anz_e; j++) seta( setBuf, "e", set[setNr].elem[j] );
  for (j=0; j<set[setNr].anz_f; j++) seta( setBuf, "f", set[setNr].face[j] );
  for (j=0; j<set[setNr].anz_p; j++) seta( setBuf, "p", set[setNr].pnt[j] );
  for (j=0; j<set[setNr].anz_l; j++) seta( setBuf, "l", set[setNr].line[j] );
  for (j=0; j<set[setNr].anz_c; j++) seta( setBuf, "c", set[setNr].lcmb[j] );
  for (j=0; j<set[setNr].anz_s; j++) seta( setBuf, "s", set[setNr].surf[j] );
  for (j=0; j<set[setNr].anz_b; j++) seta( setBuf, "b", set[setNr].body[j] );
  for (j=0; j<set[setNr].anz_nurl; j++) seta( setBuf, "L", set[setNr].nurl[j] );
  for (j=0; j<set[setNr].anz_nurs; j++) seta( setBuf, "S", set[setNr].nurs[j] );
  for (j=0; j<set[setNr].anz_sh; j++) seta( setBuf, "sh", set[setNr].shp[j] );
  for (j=0; j<set[setNr].anz_se; j++) seta( setBuf, "r", set[setNr].set[j] );
  setNr=setBuf;

  /* do not save entities stored in set :NOSAVE */
  /* remove the contents of nsave from setNr */
  nsave=getSetNr(specialset->nsave);
  if(nsave>0)
  {
    for (j=0; j<set[nsave].anz_v; j++) setr( setNr, "v", set[nsave].valu[j] );
    for (j=0; j<set[nsave].anz_n; j++) setr( setNr, "n", set[nsave].node[j] );
    for (j=0; j<set[nsave].anz_e; j++) setr( setNr, "e", set[nsave].elem[j] );
    for (j=0; j<set[nsave].anz_f; j++) setr( setNr, "f", set[nsave].face[j] );
    for (j=0; j<set[nsave].anz_p; j++) setr( setNr, "p", set[nsave].pnt[j] );
    for (j=0; j<set[nsave].anz_l; j++) setr( setNr, "l", set[nsave].line[j] );
    for (j=0; j<set[nsave].anz_c; j++) setr( setNr, "c", set[nsave].lcmb[j] );
    for (j=0; j<set[nsave].anz_s; j++) setr( setNr, "s", set[nsave].surf[j] );
    for (j=0; j<set[nsave].anz_b; j++) setr( setNr, "b", set[nsave].body[j] );
    for (j=0; j<set[nsave].anz_nurl; j++) setr( setNr, "L", set[nsave].nurl[j] );
    for (j=0; j<set[nsave].anz_nurs; j++) setr( setNr, "S", set[nsave].nurs[j] );
    for (j=0; j<set[nsave].anz_sh; j++) setr( setNr, "sh", set[nsave].shp[j] );
    for (j=0; j<set[nsave].anz_se; j++) setr( setNr, "r", set[nsave].set[j] );
  }
  
  if (anzGeo->p)
  {
      for (i=0; i<set[setNr].anz_p; i++)
      {
        j=set[setNr].pnt[i];
        if(compare(format,"lf",2)==2)
          fprintf (handle1, " PNT %s  %.12lf  %.12lf  %.12lf \n", pnt[j].name, pnt[j].px,pnt[j].py,pnt[j].pz);
        else
          fprintf (handle1, " PNT %s  %.12le  %.12le  %.12le \n", pnt[j].name, pnt[j].px,pnt[j].py,pnt[j].pz);
        //fprintf (handle1, " PNT %s  %.12lf  %.12lf  %.12lf \n", pnt[j].name, pnt[j].px,pnt[j].py,pnt[j].pz);
        //fprintf (handle1, " PNT %-5s %14.5lf %14.5lf %14.5lf \n", pnt[j].name, pnt[j].px,pnt[j].py,pnt[j].pz);
      }
  }

  if (anz->sets)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=1; i<anz->sets; i++)
      {
        if(set[i].type) /* SEQ */
        {
          if ( ((compare( set[i].name, setname, length) == length)
          && (length == strlen(setname) ))
          || (( compare( setname, "all", 3) == 3) 
	      &&(set[i].name != (char *)NULL )) )
          {
            params=MAX_PARAM_PER_RECORD/2;
            if(set[i].anz_p)
            {
              strcpy(string,"   ");
              fprintf (handle1, " SEQA %s %s pnt ", set[i].name, string );
              for (j=0; j<set[i].anz_p; j++)
              {
                if(j==params)
        	   {
        	     params+=MAX_PARAM_PER_RECORD/2;
                  strcpy(string,"END");
                  fprintf (handle1, "\n SEQA %s %s pnt ", set[i].name, string );
                }
                fprintf (handle1, " %s", pnt[set[i].pnt[j]].name);
              }
              fprintf (handle1, "\n");
            }
          }
        }
      }
    }
    else
    {
      if ( set[setNr].name != (char *)NULL )
      {
        if(set[setNr].type==1)
	{
          if(set[setNr].anz_n)
          {
            params=MAX_PARAM_PER_RECORD/2;
            strcpy(string,"   ");
            fprintf (handle1, " SEQA %s %s nod ", set[setNr].name, string );
            for (j=0; j<set[setNr].anz_n; j++)
            {
              if(j==params)
              {
          	params+=MAX_PARAM_PER_RECORD/2;
                strcpy(string,"END");
                fprintf (handle1, "\n SEQA %s %s nod ", set[setNr].name, string );
              }
              fprintf(handle1, "%d ",  set[setNr].node[j]);
            }
            fprintf (handle1, "\n");
          }
	}

        for (i=0; i<set[setNr].anz_l; i++)
        {
          l=set[setNr].line[i];

          if (line[l].typ=='s')     
          {
            se=line[l].trk;
            if ( set[se].name != (char *)NULL )
            {
              params=MAX_PARAM_PER_RECORD/2;
              if(set[se].anz_p)
              {
                strcpy(string,"   ");
                fprintf (handle1, " SEQA %s %s pnt ", set[se].name, string );
                for (j=0; j<set[se].anz_p; j++)
                {
                  if(j==params)
          	  {
          	     params+=MAX_PARAM_PER_RECORD/2;
                    strcpy(string,"END");
                    fprintf (handle1, "\n SEQA %s %s pnt ", set[se].name, string );
                  }
                  fprintf (handle1, " %s", pnt[set[se].pnt[j]].name);
                }
                fprintf (handle1, "\n");
              }
            }
          }
        }
      }
    }
  }

  if (anzGeo->l)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->l; i++) if( line[i].name != (char *)NULL )
      {
        /* load bias and division in buffer */
        if(line[i].div<100)
	{
          bias_fbd=getBias_fbd(i,line);
          if(line[i].div<10) sprintf(buffer, "%d0%d", bias_fbd, line[i].div);
          else sprintf(buffer, "%d%d", bias_fbd, line[i].div);
	}
        else sprintf(buffer, "%d %f", line[i].div, line[i].bias);

        if (line[i].typ=='a') /* its an arc-line */
        {
          if( line[i].name != (char *)NULL )
            fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[i].name, pnt[line[i].p1].name,
            pnt[line[i].p2].name, pnt[line[i].trk].name, buffer );
        }
        else if (line[i].typ=='s') /* its a seq-line */
        {
          /* check if the seq is correct defined */
          /* - the endpoints must be included in the seq */
          if( ((line[i].p1!=set[line[i].trk].pnt[0])&&(line[i].p2!=set[line[i].trk].pnt[0]))
            || ((line[i].p1!=set[line[i].trk].pnt[set[line[i].trk].anz_p-1])&&(line[i].p2!=set[line[i].trk].pnt[set[line[i].trk].anz_p-1])) )
	  {
            /* closer check */
            k=0;
            for(j=0; j<set[line[i].trk].anz_p; j++)
	    {
              if(set[line[i].trk].pnt[j]==line[i].p1) k++;
              if(set[line[i].trk].pnt[j]==line[i].p2) k++;
	    }
            if(k!=2)
	    {
              /* add line to the error-set */
              pre_seta( specialset->tmp, "l", line[i].name);
	    }
	  }

          if( line[i].name != (char *)NULL )
            fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[i].name, pnt[line[i].p1].name,
            pnt[line[i].p2].name, set[line[i].trk].name, buffer );
        }
        else
        {
          if( line[i].name != (char *)NULL )
            fprintf (handle1, " LINE %s %s %s %s\n",
            line[i].name, pnt[line[i].p1].name, pnt[line[i].p2].name, buffer );
        }
      }
      for (i=0; i<anzGeo->l; i++)
      {
        if( line[i].name != (char *)NULL )
        {
          if(line[i].etyp>0)
	  {
            fprintf (handle1, " MSHP %s l %d %d %d\n", line[i].name, line[i].etyp, line[i].eattr, line[i].elock );
	  }
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_l; i++)
      {
        j=set[setNr].line[i];

        /* load bias and division in buffer */
        if(line[j].div<100)
	{
          bias_fbd=getBias_fbd(j,line);
          if(line[j].div<10) sprintf(buffer, "%d0%d", bias_fbd, line[j].div);
          else if(line[j].div<100) sprintf(buffer, "%d%d", bias_fbd, line[j].div);
	}
        else sprintf(buffer, "%d %f", line[j].div, line[j].bias);

        if (line[j].typ=='a') /* its an arc-line */
        {
          fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[j].name, pnt[line[j].p1].name,
            pnt[line[j].p2].name, pnt[line[j].trk].name, buffer );
        }
        else if (line[j].typ=='s') /* its an seq-line */
        {
          /* check if the seq is correct defined */
          /* - the endpoints must be included in the seq */
          if( ((line[j].p1!=set[line[j].trk].pnt[0])&&(line[j].p2!=set[line[j].trk].pnt[0]))
            || ((line[j].p1!=set[line[j].trk].pnt[set[line[j].trk].anz_p-1])&&(line[j].p2!=set[line[j].trk].pnt[set[line[j].trk].anz_p-1])) )
	  {
            /* closer check */
            k=0;
            for(n=0; n<set[line[j].trk].anz_p; n++)
	    {
              if(set[line[j].trk].pnt[n]==line[j].p1) k++;
              if(set[line[j].trk].pnt[n]==line[j].p2) k++;
	    }
            if(k!=2)
	    {
              /* add line to the error-set */
              pre_seta( specialset->tmp, "l", line[j].name);
	    }
	  }

          fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[j].name, pnt[line[j].p1].name, 
            pnt[line[j].p2].name, set[line[j].trk].name, buffer );
        }
        else
        {
          fprintf (handle1, " LINE %s %s %s %s\n",
            line[j].name, pnt[line[j].p1].name, pnt[line[j].p2].name, buffer );
        }
      }
      for (i=0; i<set[setNr].anz_l; i++)
      {
        j=set[setNr].line[i];
        if( line[j].name != (char *)NULL )
        {
          if(line[j].etyp>0)
	  {
            fprintf (handle1, " MSHP %s l %d %d %d\n", line[j].name, line[j].etyp, line[j].eattr, line[j].elock );
	  }
        }
      }

    }
  }

  if (anzGeo->c)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->c; i++)
      {
        if( lcmb[i].name != (char *)NULL )
        {
          fprintf (handle1, " LCMB %s ", lcmb[i].name );
          k=0;
          for (j=0; j<lcmb[i].nl; j++)
	  {
            if (k<8) 
	    {
              fprintf (handle1, " %1c %s", lcmb[i].o[j], line[lcmb[i].l[j]].name );
              k++;
            }
            else
	    {
              k=0; fprintf (handle1, "\n LCMB %s ADD ", lcmb[i].name );
              fprintf (handle1, " %1c %s", lcmb[i].o[j], line[lcmb[i].l[j]].name );
              k++;
            }
          }
          fprintf (handle1, " \n");
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_c; i++)
      {
          n=set[setNr].lcmb[i];
          fprintf (handle1, " LCMB %s ", lcmb[n].name );
          k=0;
          for (j=0; j<lcmb[n].nl; j++)
	  {
            if (k<8) 
	    {
              fprintf (handle1, " %1c %s", lcmb[n].o[j], line[lcmb[n].l[j]].name );
              k++;
            }
            else
	    {
              k=0; fprintf (handle1, "\n LCMB %s ADD ", lcmb[n].name );
              fprintf (handle1, " %1c %s", lcmb[n].o[j], line[lcmb[n].l[j]].name );
              k++;
            }
          }
          fprintf (handle1, " \n");
      }
    }
  }

  if (anzGeo->sh)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->sh; i++)
      {
        if( shape[i].name != (char *)NULL )
        {
          switch (shape[i].type)
          {
            case 0:
              fprintf (handle1, " SHPE %s PLN %s %s %s\n", shape[i].name, pnt[shape[i].p[0]].name, pnt[shape[i].p[1]].name, pnt[shape[i].p[2]].name );
            break;
            case 1:
              v_result( &pnt[shape[i].p[0]].px, &pnt[shape[i].p[2]].px, p1  );
              fprintf (handle1, " SHPE %s CYL %s %s %f\n", shape[i].name, pnt[shape[i].p[0]].name, pnt[shape[i].p[1]].name, v_betrag(p1));
            break;
            case 2:
              v_result( &pnt[shape[i].p[0]].px, &pnt[shape[i].p[2]].px, p1  );
              v_result( &pnt[shape[i].p[1]].px, &pnt[shape[i].p[3]].px, p2  );
              fprintf (handle1, " SHPE %s CON %s %s %f %f\n", shape[i].name, pnt[shape[i].p[0]].name, pnt[shape[i].p[1]].name, v_betrag(p1), v_betrag(p2));
            break;
            case 3:
              v_result( &pnt[shape[i].p[0]].px, &pnt[shape[i].p[1]].px, p1  );
              fprintf (handle1, " SHPE %s SPH %s %f\n", shape[i].name, pnt[shape[i].p[0]].name, v_betrag(p1));
            break;
            case 5:
              v_result( &pnt[shape[i].p[0]].px, &pnt[shape[i].p[2]].px, p1  );
              v_result( &pnt[shape[i].p[2]].px, &pnt[shape[i].p[3]].px, p2  );
              fprintf (handle1, " SHPE %s TOR %s %f %s %f\n", shape[i].name, pnt[shape[i].p[0]].name, v_betrag(p1), pnt[shape[i].p[1]].name, v_betrag(p2));
            break;
          }
	}
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_sh; i++)
      {
        n=set[setNr].shp[i];
        switch (shape[n].type)
        {
          case 0:
            fprintf (handle1, " SHPE %s PLN %s %s %s\n", shape[n].name, pnt[shape[n].p[0]].name, pnt[shape[n].p[1]].name, pnt[shape[n].p[2]].name );
          break;
          case 1:
            v_result( &pnt[shape[n].p[0]].px, &pnt[shape[n].p[2]].px, p1  );
            fprintf (handle1, " SHPE %s CYL %s %s %f\n", shape[n].name, pnt[shape[n].p[0]].name, pnt[shape[n].p[1]].name, v_betrag(p1));
          break;
          case 2:
            v_result( &pnt[shape[n].p[0]].px, &pnt[shape[n].p[2]].px, p1  );
            v_result( &pnt[shape[n].p[1]].px, &pnt[shape[n].p[3]].px, p2  );
            fprintf (handle1, " SHPE %s CON %s %s %f %f\n", shape[n].name, pnt[shape[n].p[0]].name, pnt[shape[n].p[1]].name, v_betrag(p1), v_betrag(p2));
          break;
          case 3:
            v_result( &pnt[shape[n].p[0]].px, &pnt[shape[n].p[1]].px, p1  );
            fprintf (handle1, " SHPE %s SPH %s %f\n", shape[n].name, pnt[shape[n].p[0]].name, v_betrag(p1));
          break;
          case 5:
            v_result( &pnt[shape[n].p[0]].px, &pnt[shape[n].p[2]].px, p1  );
            v_result( &pnt[shape[n].p[2]].px, &pnt[shape[n].p[3]].px, p2  );
            fprintf (handle1, " SHPE %s TOR %s %f %s %f\n", shape[n].name, pnt[shape[n].p[0]].name, v_betrag(p1), pnt[shape[n].p[1]].name, v_betrag(p2));
          break;
        }
      }
    }
  }

  if (anzGeo->nurs)
  {
    if( compare( setname, "all", 3)==3)
      for (nr=0; nr<anzGeo->nurs; nr++)
      {
        if( nurbs[nr].name != (char *)NULL )
        {
          for (j=0; j<nurbs[nr].u_npnt; j++)
            for (k=0; k<nurbs[nr].v_npnt; k++)
	    {
              if((nurbs[nr].ctlpnt[j][k]<0) || (getPntNr(pnt[nurbs[nr].ctlpnt[j][k]].name)<0))
	      {
                printf("ERROR: Control-Point:%d,%d in Nurbs:%s is not defined, Nurbs will not be written. Please talk to the programmer.\n", j+1,k+1, nurbs[nr].name );
                goto nextNurbs;
	      }
	    }

          fprintf(handle1," NURS %s DEFINE %4d %4d %4d %4d %4d %4d\n", nurbs[nr].name, nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, nurbs[nr].u_nknt, nurbs[nr].v_nknt);

          for (j=0; j<nurbs[nr].u_npnt; j++)
            for (k=0; k<nurbs[nr].v_npnt; k++)
              fprintf (handle1, " NURS %s CONTROL %4d %4d %s %f\n", nurbs[nr].name, j+1, k+1, pnt[nurbs[nr].ctlpnt[j][k]].name, nurbs[nr].weight[j][k] );

          for (j=0; j<nurbs[nr].u_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  U %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].uknt[j]);
          for (j=0; j<nurbs[nr].v_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  V %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].vknt[j]);

          fprintf(handle1, " NURS %s END\n", nurbs[nr].name);
	nextNurbs:;
        }
      }
    else
    {
      for (i=0; i<set[setNr].anz_nurs; i++)
      {
        nr=set[setNr].nurs[i];

          fprintf(handle1," NURS %s DEFINE %4d %4d %4d %4d %4d %4d\n", nurbs[nr].name, nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, nurbs[nr].u_nknt, nurbs[nr].v_nknt);

          for (j=0; j<nurbs[nr].u_npnt; j++)
            for (k=0; k<nurbs[nr].v_npnt; k++)
              fprintf (handle1, " NURS %s CONTROL %4d %4d %s %f\n", nurbs[nr].name, j+1, k+1, pnt[nurbs[nr].ctlpnt[j][k]].name, nurbs[nr].weight[j][k]);

          for (j=0; j<nurbs[nr].u_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  U %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].uknt[j]);
          for (j=0; j<nurbs[nr].v_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  V %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].vknt[j]);

          fprintf(handle1, " NURS %s END\n", nurbs[nr].name);

      }
    }
  }

  /* must follow the nurs, shapes */
  if (anzGeo->s)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->s; i++)
      {
        if( surf[i].name != (char *)NULL )
        {
	  if(!surf[i].o[0]) { surf[i].o[0]='+'; k=1; } else k=0;
          if(surf[i].sh>-1) fprintf (handle1, " GSUR %s %1c %s", surf[i].name, surf[i].ori, shape[surf[i].sh].name );
          else fprintf (handle1, " GSUR %s %1c BLEND", surf[i].name, surf[i].ori );
          params=0; 
          for (j=0; j<surf[i].nl; j++)
          {
            if (params==MAX_GSUR_PARAMERTER)
            { params=0; fprintf (handle1, "\n GSUR %s ADD", surf[i].name); }
            else  params++; 
            if (surf[i].typ[j]=='l')
            {
              fprintf (handle1, " %1c %s", surf[i].o[j], line[surf[i].l[j]].name );
            }
            else
            {
              fprintf (handle1, " %1c %s", surf[i].o[j], lcmb[surf[i].l[j]].name );
            }
          }
          fprintf (handle1, " \n");
          if(k) { surf[i].o[0]=0; }
        }
      }
      for (i=0; i<anzGeo->s; i++)
      {
        if( surf[i].name != (char *)NULL )
        {
          if(surf[i].etyp>0)
	  {
            if( surf[i].eparm != (char *)NULL )
              fprintf (handle1, " MSHP %s s %d %d %d %s\n", surf[i].name, surf[i].etyp, surf[i].eattr, surf[i].elock, surf[i].eparm );
            else
              fprintf (handle1, " MSHP %s s %d %d %d\n", surf[i].name, surf[i].etyp, surf[i].eattr, surf[i].elock );
	  }
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_s; i++)
      {
        n=set[setNr].surf[i];
	if(!surf[n].o[0]) { surf[n].o[0]='+'; k=1; } else k=0;
        if(surf[n].sh<=-1) fprintf (handle1, " GSUR %s %1c BLEND", surf[n].name, surf[n].ori );
        else if(surf[n].sh>-1) fprintf (handle1, " GSUR %s %1c %s", surf[n].name, surf[n].ori, shape[surf[n].sh].name );
        params=0; 
        for (j=0; j<surf[n].nl; j++)
        {
          if (params==MAX_GSUR_PARAMERTER)
          { params=0; fprintf (handle1, "\n GSUR %s ADD", surf[n].name); }
          else  params++; 
          if (surf[n].typ[j]=='l')
          {
            fprintf (handle1, " %1c %s", surf[n].o[j], line[surf[n].l[j]].name );
          }
          else
          {
            fprintf (handle1, " %1c %s", surf[n].o[j], lcmb[surf[n].l[j]].name );
          }
        }
        fprintf (handle1, " \n");
        if(k) { surf[n].o[0]=0; }
      }
      for (n=0; n<set[setNr].anz_s; n++)
      {
        i=set[setNr].surf[n];
        if( surf[i].name != (char *)NULL )
        {
          if(surf[i].etyp>0)
	  {
            if( surf[i].eparm != (char *)NULL )
              fprintf (handle1, " MSHP %s s %d %d %d %s\n", surf[i].name, surf[i].etyp, surf[i].eattr, surf[i].elock, surf[i].eparm );
            else
              fprintf (handle1, " MSHP %s s %d %d %d\n", surf[i].name, surf[i].etyp, surf[i].eattr, surf[i].elock );
          }
        }
      }
    }
  }

  if (anzGeo->b)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->b; i++)
      {
        if( body[i].name != (char *)NULL )
        {
          fprintf (handle1, " GBOD %s NORM", body[i].name );
          params=0; 
          for (j=0; j<body[i].ns; j++)
          {
            if (params==MAX_GSUR_PARAMERTER)
            { params=1; fprintf (handle1, "\n GBOD %s ADD", body[i].name); }
            else  params++; 
            fprintf (handle1, " %1c %s", body[i].o[j], surf[body[i].s[j]].name );
          }
          fprintf (handle1, " \n");
        }
      }
      for (i=0; i<anzGeo->b; i++)
      {
        if( body[i].name != (char *)NULL )
        {
          if(body[i].etyp>0)
	  {
            if( body[i].eparm != (char *)NULL )
              fprintf (handle1, " MSHP %s b %d %d %d %s\n", body[i].name, body[i].etyp, body[i].eattr, body[i].elock, body[i].eparm );
            else
              fprintf (handle1, " MSHP %s b %d %d %d\n", body[i].name, body[i].etyp, body[i].eattr, body[i].elock );
	  }
        }
      }
    }
    else
    {
      for (n=0; n<set[setNr].anz_b; n++)
      {
        i=set[setNr].body[n];
        if( body[i].name != (char *)NULL )
        {
          fprintf (handle1, " GBOD %s NORM", body[i].name );
          params=0; 
          for (j=0; j<body[i].ns; j++)
          {
            if (params==MAX_GSUR_PARAMERTER)
            { params=1; fprintf (handle1, "\n GBOD %s ADD", body[i].name); }
            else  params++; 
            fprintf (handle1, " %1c %s", body[i].o[j], surf[body[i].s[j]].name );
          }
          fprintf (handle1, " \n");
        }
      }
      for (n=0; n<set[setNr].anz_b; n++)
      {
        i=set[setNr].body[n];
        if( body[i].name != (char *)NULL )
        {
          if(body[i].etyp>0)
	  {
            if( body[i].eparm != (char *)NULL )
              fprintf (handle1, " MSHP %s b %d %d %d %s\n", body[i].name, body[i].etyp, body[i].eattr, body[i].elock, body[i].eparm );
            else
              fprintf (handle1, " MSHP %s b %d %d %d\n", body[i].name, body[i].etyp, body[i].eattr, body[i].elock );
	  }
        }
      }
    }
  }

  if (anz->v)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anz->v; i++)
      {
        if(( value[i].name != (char *)NULL )&& ( !value[i].flag ))
        {
          fprintf (handle1, " VALU %s %s\n", value[i].name, value[i].string );
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_v; i++)
      {
        n=set[setNr].valu[i];
        fprintf (handle1, " VALU %s %s\n", value[n].name, value[n].string );
      }
    }
  }

  /* write the sets */
  for (i=1; i<anz->sets; i++)
  {
    if((!set[i].type)&&( set[i].name!=(char *)NULL)&&( set[i].name[0]!='_')&&( set[i].name[0]!='-')) /* ordinary SET */
    {
      if( compare( setname, "all", 3)==3);
      else
      {
      k=0; for (j=0; j<set[i].anz_n; j++) if( getIndex(&set[setNr].node,set[setNr].anz_n,set[i].node[j]) >-1)
      {
        if(k==0) { fprintf (handle1, " SETA %s n", set[i].name); }
        fprintf (handle1, " %d", set[i].node[j]);  k++;
        if((k==7)&&(j!=set[i].anz_n-1)) { fprintf (handle1, "\n"); k=0; }
      }
      if(k) { fprintf (handle1, "\n"); }
      k=0; for (j=0; j<set[i].anz_e; j++) if( getIndex(&set[setNr].elem,set[setNr].anz_e,set[i].elem[j]) >-1)
      {
        if(k==0) { fprintf (handle1, " SETA %s e", set[i].name); }
        fprintf (handle1, " %d", set[i].elem[j]);   k++;
        if((k==7)&&(j!=set[i].anz_e-1)) { fprintf (handle1, "\n"); k=0; }
      }
      if(k) { fprintf (handle1, "\n"); }
      }

      /*  FAM vertraegt nur einen wert hinter der spec.  */
      for (j=0; j<set[i].anz_p; j++) if( getIndex(&set[setNr].pnt,set[setNr].anz_p,set[i].pnt[j]) >-1)
      {
        if( pnt[set[i].pnt[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s p %s \n", set[i].name, pnt[set[i].pnt[j]].name);
      }
      for (j=0; j<set[i].anz_l; j++) if( getIndex(&set[setNr].line,set[setNr].anz_l,set[i].line[j]) >-1)
      {
        if( line[set[i].line[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s l %s \n", set[i].name, line[set[i].line[j]].name);
      }
      for (j=0; j<set[i].anz_c; j++) if( getIndex(&set[setNr].lcmb,set[setNr].anz_c,set[i].lcmb[j]) >-1)
      {
        if( lcmb[set[i].lcmb[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s c %s \n", set[i].name, lcmb[set[i].lcmb[j]].name);
      }
      for (j=0; j<set[i].anz_s; j++) if( getIndex(&set[setNr].surf,set[setNr].anz_s,set[i].surf[j]) >-1)
      {
        if( surf[set[i].surf[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s s %s \n", set[i].name, surf[set[i].surf[j]].name);
      }
      for (j=0; j<set[i].anz_b; j++) if( getIndex(&set[setNr].body,set[setNr].anz_b,set[i].body[j]) >-1)
      {
        if( body[set[i].body[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s b %s \n", set[i].name, body[set[i].body[j]].name);
      }
      for (j=0; j<set[i].anz_nurs; j++) if( getIndex(&set[setNr].nurs,set[setNr].anz_nurs,set[i].nurs[j]) >-1)
      {
        if( nurbs[set[i].nurs[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s S %s \n", set[i].name, nurbs[set[i].nurs[j]].name);
      }
      for (j=0; j<set[i].anz_sh; j++) if( getIndex(&set[setNr].shp,set[setNr].anz_sh,set[i].shp[j]) >-1)
      {
        if( shape[set[i].shp[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s sh %s \n", set[i].name, shape[set[i].shp[j]].name);
      }
    }
  }

  /* write the quality thresholds */
  if(eqal.aspr) fprintf (handle1, " EQAL aspr %lf\n", eqal.aspr);
  if(eqal.jbir) fprintf (handle1, " EQAL jbir %lf\n", eqal.jbir);
  if(eqal.mca)  fprintf (handle1, " EQAL mca %lf\n", eqal.mca);

  /* check if specialset->tmp was created */
  i= getSetNr(specialset->tmp);
  if(i>-1)
  { 
    printf (" WARNING: lines in set:%s use seqences which do not include the corner-points. Try to merg the points\n", specialset->tmp);
  }

  if( compare( setname, "all", 3)==3)
  {
    /* write the user colors */
    for(i=SET_COLS; i<entitycols; i++)
    {
      fprintf(handle1, " col %s %f %f %f\n",entitycol[i].name, entitycol[i].r, entitycol[i].g, entitycol[i].b  );
    }

    /* write the nr of threads */
    if(anz->threads>1) fprintf (handle1, " asgn thrds %d\n", anz->threads);

    /* define the tetmesher */
    if(meshp.tetmesher==1) fprintf (handle1, " asgn tetgen\n");

    /* write the initialization commands */
    //if(inpformat) fprintf (handle1, " wpos %d %d\n",glutGet(GLUT_WINDOW_X), glutGet(GLUT_WINDOW_Y));
    writeInit(handle1);
    
    /* write the plot,plus commands */
    for (j=0; j<anzGeo->psets; j++ )
    {
      //if((pset[j].type[0]!='n')&&(pset[j].type[0]!='e')&&(pset[j].type[0]!='f'))
      {
        if(j==0) fprintf (handle1, " plot %s %s %s\n", pset[j].type, set[pset[j].nr].name, entitycol[pset[j].col].name);
        else fprintf (handle1, " plus %s %s %s\n", pset[j].type, set[pset[j].nr].name, entitycol[pset[j].col].name);
      }
    }

    /* write the current dataset */
    if((anz->l)&&((drawMode==1)||(drawMode==4)))
    {
      fprintf (handle1, " ds %d e %d\n", cur_lc+1, cur_entity+1);
      fprintf (handle1, " STEPS %d\n", steps);
      if(scale->smaxr==0) fprintf (handle1, " MAX %lf", scale->smax);
      if(scale->smaxr==1) fprintf (handle1, " MAXR %lf", scale->smax);
      if(scale->smaxr==2) fprintf (handle1, " MAXC %lf", scale->smax);
      if(scale->format) fprintf (handle1, " %c", scale->format);
      if(scale->lock) fprintf (handle1, " %c", scale->lock);
      fprintf (handle1, "\n");
      if(scale->sminr==0) fprintf (handle1, " MIN %lf", scale->smin);
      if(scale->sminr==1) fprintf (handle1, " MINR %lf", scale->smin);
      if(scale->sminr==2) fprintf (handle1, " MINC %lf", scale->smin);
      if(scale->format) fprintf (handle1, " %c", scale->format);
      if(scale->lock) fprintf (handle1, " %c", scale->lock);
      fprintf (handle1, "\n");
    }
  }
  else
  {
    fprintf (handle1, " plus p %s\n", setname);
    fprintf (handle1, " plus lp %s\n", setname);
  }
  fclose(handle1);
  delSet("-setbuf");
  return (1);
}


void writeInit(FILE *handle1)
{
  static GLint ipuf[2];

  if(!inpformat) return;
  glGetIntegerv(GL_CULL_FACE_MODE, ipuf );
  fprintf (handle1, " init %d %d %d =\n", width_w0, height_w0, ipuf[0]);
  fprintf (handle1, " %f %f %f %f %f %f %f %f =\n", ds,dtx,dty,dtz, vmem[0], vmem[1], vmem[2], vmem[3]);
  fprintf (handle1, " %f %f %f %f %f %f %f %f =\n", R[0][0], R[0][1], R[0][2], R[0][3], R[1][0], R[1][1], R[1][2], R[1][3]);
  fprintf (handle1, " %f %f %f %f %f %f %f %f =\n", R[2][0], R[2][1], R[2][2], R[2][3], R[3][0], R[3][1], R[3][2], R[3][3]);
  fprintf (handle1, " %f %f %f %f %f %f %f =\n", lastquat[0], lastquat[1], lastquat[2], lastquat[3], centerPnt[0], centerPnt[1], centerPnt[2]);
  fprintf (handle1, " %d %d %d %d %d %d %d %d\n",backgrndcol,legend_font,draw_font,menu_font,modelEdgeFlag,elemEdgeFlag,surfFlag,illumResultFlag);
}
