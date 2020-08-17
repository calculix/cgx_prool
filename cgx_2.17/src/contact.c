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

#define TEST 0
#define     N_CLOSEST_NODES 10

extern SpecialSet specialset[1];
extern char  printFlag;                     /* printf 1:on 0:off */
extern Sets      *set;                               /* must! else pre_seta will change the address */
extern Entitycol *entitycol;  /* predefined colors of entities */


int mergeCluster(Cluster *cluster, int c1, int c2)
{
  int i,n;
  int sum_n=0, *cnode=NULL, mergFlag=0;

  /* add both node-lists */
  if((cnode =(int *)malloc((cluster[c1].nn + cluster[c2].nn+1) *sizeof(int)))==NULL)
  { errMsg("\n\n ERROR: malloc failed for cluster\n" ); return(-1); }
  for(i=0; i<cluster[c1].nn; i++) cnode[sum_n++]=cluster[c1].node[i];
  for(i=0; i<cluster[c2].nn; i++) cnode[sum_n++]=cluster[c2].node[i];

  /* sort both node-lists */
  qsort( cnode, sum_n, sizeof(int), (void *)compareInt );
#if TEST
  for (i=0; i<sum_n; i++)
    printf("1 %d n:%d\n", i, cnode[i]); 
#endif


  /* go over the node-list and look for double nodenr */
  /* if yes, delete second node from list */
  for (i=0; i<sum_n-1; i++)
  {
    if(cnode[i]==cnode[i+1])
    {
      mergFlag=1;
      cnode[i]=0;
    }
  }

  if(mergFlag)
  {
    n=0;
    for (i=0; i<sum_n; i++)
    {
      if(cnode[i])
      {
        cnode[n]=cnode[i];
        n++;
      }
    }
    sum_n=n;
    /* the first original cluster gets the new address of the node-list */
    free(cluster[c1].node);
    cluster[c1].node=cnode;
    cluster[c1].nn=sum_n;

    /* and all faces of the second cluster in addition */
    if((cluster[c1].face =(int *)realloc((int *)cluster[c1].face, (cluster[c1].nf+cluster[c2].nf) *sizeof(int)))==NULL)
    { errMsg("\n\n ERROR: malloc failed for cluster\n" ); return(-1); }
    for(i=0; i<cluster[c2].nf; i++) cluster[c1].face[cluster[c1].nf++]=cluster[c2].face[i];

    /* delete the second original cluster */
    free(cluster[c2].node);
    free(cluster[c2].face);
    cluster[c2].nn=0;
    cluster[c2].nf=0;
    return(1);
  }

  /* if no, return(0) without action */
  free(cnode);
  return(0);
}


/* return sum of generated sets */
int getMeshSections(int setNr, Summen *anz, Faces *face, Nodes *node)
{
  int i,f,k=2;
  int setNrbuf, setTmp, sum_c=0, c1, sum1,sum2;
  int *fUsed;
  char setname[MAX_LINE_LENGTH];
  char buffer[MAX_LINE_LENGTH];

  printf (" please wait for 'ready'\n");


  if (set[setNr].anz_f<2)
  {
    //printf("ERROR: set:%s contains less than 2 faces!\n", set[setNr].name);
    return(0);
  }
  /*----- complete set with all faces and nodes ---*/

  /* mark all faces as unused */
  if((fUsed = (int *)malloc( (int)(anz->f)*sizeof(int)))==NULL)
  { errMsg("\n\n ERROR: malloc failed for elUsed\n" ); return(0); }
  for(i=0; i<anz->f; i++) fUsed[i]=1;
  for(i=0; i<set[setNr].anz_f; i++) fUsed[set[setNr].face[i]]=0;

  /* use "comp set up/do" until no more entities are added */
  f=set[setNr].face[0];
  do
  {
    sprintf( setname, "%s%d", specialset->cf,sum_c+1); 
    delSet(setname);
    c1=pre_seta(setname,"i", 0);
    seta(c1,"f", f);                                            // 1st face in specialset->cf
    fUsed[f]=1; 

    delSet("+buf");
    sprintf(buffer,"%d", f);
    setNrbuf=pre_seta("+buf","f", buffer);
    completeSet_Faces(  setNrbuf, 0, fUsed, 1);                 // down, nodes from new face in +buf
    for(i=0; i<set[setNrbuf].anz_n; i++)
    {
      seta(c1,"n", set[setNrbuf].node[i]);                      // nodes from new face in specialset->cf
    }

    do
    {
      sum1=set[c1].anz_n+set[c1].anz_f;

      delSet("+buf2");
      setTmp=pre_seta("+buf2", "i", 0 );
      completeSet_Faces( setNrbuf, setTmp, fUsed, 0);         // up, only new faces from +buf (setNrbuf) in +buf2 (setTmp)
      completeSet_Faces(  setTmp, 0, fUsed, 1);               // down, nodes from new faces in +buf2 (setTmp)
      for(i=0; i<set[setNrbuf].anz_n; i++)
        setr(setTmp,"n", set[setNrbuf].node[i]);              // remove used nodes in +buf from +buf2
      for(i=0; i<set[setTmp].anz_f; i++)
      {
        seta(c1,"f", set[setTmp].face[i]);                    // new faces in specialset->cf
        fUsed[set[setTmp].face[i]]=1; 
      }

      delSet("+buf");
      setNrbuf=pre_seta("+buf","i", 0);
      for(i=0; i<set[setTmp].anz_n; i++)
      {
        seta(c1,"n", set[setTmp].node[i]);
        seta(setNrbuf,"n", set[setTmp].node[i]);
      }
      sum2=set[c1].anz_n+set[c1].anz_f;
    }while(sum1<sum2);

    //completeSet( setname, "f");                  // faces (draw-nodes are missed)
    completeSet( setname, "do");                  // faces, nodes (to include also the additional midface-nodes for drawing)

    if(k>=SET_COLS) k=BAS_COLS-1;
    printf("set%d: %s col:%s\n",sum_c+1, setname, entitycol[k].name);
    sprintf(buffer,"f %s %s", setname, entitycol[k++].name);
    if(!sum_c) plot(buffer);
    else plus(buffer);

    sum_c++;
    f=0;
    for(i=0; i<set[setNr].anz_f; i++)
    {
      if(fUsed[set[setNr].face[i]]==0)
      {
        f=set[setNr].face[i];
        break;
      }
    }
  }while(f);

  free(fUsed);
  delSet("+buf");
  delSet("+buf2");

  printf (" ready\n");
  return(sum_c);
}


/* search dependent nodes of unconnected sets of faces */
/* return 1 if depnodes were found */

int getDepNodes( Summen *anz, int set1, int set2, Elements *e_enqire, Faces *face, Nodes *node, double tol, int *sum_n, int **nodesets)
{
  int  i, j, n, e, f;
  int  n_closest_nodes, eset;

  int    elbuf=0;
  double dx,dy,dz;
  double xmin,ymin,zmin,xmax,ymax,zmax;
  double offset=0., min_offset, offsetbuf=0.;

  double *orig_x=NULL, *orig_y=NULL, *orig_z=NULL, *sort_x=NULL, *sort_y=NULL, *sort_z=NULL;
  int *sort_nx=NULL, *sort_ny=NULL, *sort_nz=NULL, near_node[N_CLOSEST_NODES];
  Rsort *rsort=NULL;

  int *dep=NULL, *ind=NULL;
  int sum_dep=0, sum_ind=0;

  CTri3     *ctri3=NULL;             /* triangulierte indep-nodes flaeche */
  static int **tri3_index;           /* am jeweiligen node haengenden dreiecke */
  int       *ntri_nodes=NULL;        /* anzahl der anhaengenden dreiecke       */
  int       sum_tri3=0;                /* anzahl der ctri3 */

  /* get elem-faces and a reference from nodes to faces */
  delSet("-eset");
  eset=pre_seta("-eset", "i", 0);
  for (f=0; f<set[set2].anz_f; f++)
  {
    i=set[set2].face[f];
    seta(eset, "e", face[i].elem_nr);
  }

  sum_tri3 = makeTriFromElems(set2, eset, anz->nmax, set, e_enqire, &ctri3, &tri3_index, &ntri_nodes);
  delSet("-eset");

  if(printFlag) printf (" %d ctri3 created from elements\n\n", sum_tri3 );
  if(!sum_tri3)
  {
    printf("ERROR: found no valid element\n\n" );
    return(-1);
  }

  /* get the dimensions of a cubus in which the set2-nodes are located */ 
  xmin=ymin=zmin=MAX_FLOAT;
  xmax=ymax=zmax=-MAX_FLOAT;
  for (j=0; j<set[set2].anz_n; j++)
  {
    if(xmin > node[set[set2].node[j]].nx) xmin = node[set[set2].node[j]].nx;
    if(xmax < node[set[set2].node[j]].nx) xmax = node[set[set2].node[j]].nx;
    if(ymin > node[set[set2].node[j]].ny) ymin = node[set[set2].node[j]].ny;
    if(ymax < node[set[set2].node[j]].ny) ymax = node[set[set2].node[j]].ny;
    if(zmin > node[set[set2].node[j]].nz) zmin = node[set[set2].node[j]].nz;
    if(zmax < node[set[set2].node[j]].nz) zmax = node[set[set2].node[j]].nz;
  }

  xmax+=tol;
  xmin-=tol;
  ymax+=tol;
  ymin-=tol;
  zmax+=tol;
  zmin-=tol;

  if ( (rsort = (Rsort *)malloc( (set[set2].anz_n+1) * sizeof(Rsort))) == NULL )
  {  printf("ERROR: realloc failed: Rsort\n\n" );     return(-1); }
  if ( (orig_x = (double *)malloc( (set[set2].anz_n+1) * sizeof(double))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (orig_y = (double *)malloc( (set[set2].anz_n+1) * sizeof(double))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (orig_z = (double *)malloc( (set[set2].anz_n+1) * sizeof(double))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (sort_x = (double *)malloc( (set[set2].anz_n+1) * sizeof(double))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (sort_y = (double *)malloc( (set[set2].anz_n+1) * sizeof(double))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (sort_z = (double *)malloc( (set[set2].anz_n+1) * sizeof(double))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (sort_nx = (int *)malloc( (set[set2].anz_n+1) * sizeof(int))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (sort_ny = (int *)malloc( (set[set2].anz_n+1) * sizeof(int))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 
  if ( (sort_nz = (int *)malloc( (set[set2].anz_n+1) * sizeof(int))) == NULL )
    printf("ERROR: realloc failed in areampc\n\n" ); 

  for(i=0; i<set[set2].anz_n; i++)
  {
    rsort[i].r=orig_x[i]=node[set[set2].node[i]].nx;
    rsort[i].i=i;
  }
  qsort( rsort, set[set2].anz_n, sizeof(Rsort), (void *)compareRsort );
  for(i=0; i<set[set2].anz_n; i++)
  {
    sort_x[i]=rsort[i].r;
    sort_nx[i]=rsort[i].i;
  }
  for(i=0; i<set[set2].anz_n; i++)
  {
    rsort[i].r=orig_y[i]=node[set[set2].node[i]].ny;
    rsort[i].i=i;
  }
  qsort( rsort, set[set2].anz_n, sizeof(Rsort), (void *)compareRsort );
  for(i=0; i<set[set2].anz_n; i++)
  {
    sort_y[i]=rsort[i].r;
    sort_ny[i]=rsort[i].i;
  }
  for(i=0; i<set[set2].anz_n; i++)
  {
    rsort[i].r=orig_z[i]=node[set[set2].node[i]].nz;
    rsort[i].i=i;
  }
  qsort( rsort, set[set2].anz_n, sizeof(Rsort), (void *)compareRsort );
  for(i=0; i<set[set2].anz_n; i++)
  {
    sort_z[i]=rsort[i].r;
    sort_nz[i]=rsort[i].i;
  }

  if((int)N_CLOSEST_NODES<set[set2].anz_n) n_closest_nodes=(int)N_CLOSEST_NODES; else n_closest_nodes=set[set2].anz_n;

  for (i=0; i<set[set1].anz_n; i++ )
  {
    if(node[set[set1].node[i]].nx > xmax) goto next_node;
    if(node[set[set1].node[i]].nx < xmin) goto next_node;
    if(node[set[set1].node[i]].ny > ymax) goto next_node;
    if(node[set[set1].node[i]].ny < ymin) goto next_node;
    if(node[set[set1].node[i]].nz > zmax) goto next_node;
    if(node[set[set1].node[i]].nz < zmin) goto next_node;

    /* suche die naechst-liegenden indep-nodes  */

    near3d(orig_x,orig_y,orig_z,sort_x,sort_y,sort_z,sort_nx,sort_ny,sort_nz, node[set[set1].node[i]].nx,node[set[set1].node[i]].ny,
           node[set[set1].node[i]].nz, set[set2].anz_n, &near_node[0], n_closest_nodes);
    for (j=0; j<n_closest_nodes; j++) rsort[j].i=set[set2].node[near_node[j]];
    dx= node[rsort[0].i].nx - node[set[set1].node[i]].nx;
    dy= node[rsort[0].i].ny - node[set[set1].node[i]].ny;
    dz= node[rsort[0].i].nz - node[set[set1].node[i]].nz;      
    rsort[0].r=dx*dx + dy*dy + dz*dz;      

#if TEST
    for (j=0; j<n_closest_nodes; j++)
      printf("%d n:%d\n", j, rsort[j].i)); 
#endif
    if(printFlag)
      printf(" node:%d, found close node:%d at dist:%lf\n", set[set1].node[i], rsort[0].i, sqrt(rsort[0].r));

    /* knoten liegt zwischen den indep-nodes, suche das passende element.       */
    /* fange mit den elementen des naechstliegenden nodes an, und kontrolliere   */
    /* mit dem vectorprodukt welche elementflaeche den node umschliest          */

    min_offset=MAX_INTEGER;
    for (n=0; n<n_closest_nodes; n++)
    {
      for (j=0; j<ntri_nodes[ rsort[n].i ]; j++)
      {
        e=tri3_index[ rsort[n].i ][j];
        if(printFlag)
          printf("check n:%d from e:%d  nod:%d %d %d\n",rsort[n].i, e+1,ctri3[e].nod[0],ctri3[e].nod[1],ctri3[e].nod[2] );
        offset=check_tri3( &node[set[set1].node[i]].nx, &node[ctri3[e].nod[0]].nx, &node[ctri3[e].nod[1]].nx, &node[ctri3[e].nod[2]].nx, 1);
        if(printFlag)  printf("offset:%f\n", offset);
        if( offset != MAX_INTEGER )
        {
          if(printFlag) printf(" found enclosing tri3:%d for node:%d, offset:%lf\n", e+1, set[set1].node[i], offset);
          
          /* merke dir das element und den offset */
          if (min_offset > offset*offset ) { min_offset=offset*offset;  elbuf=e; offsetbuf=offset;}
        }
      }
    }
    if(printFlag) 
      printf("min_offset!=(double)MAX_INTEGER:%d tol*tol>=offsetbuf*offsetbuf:%d\n", min_offset!=(double)MAX_INTEGER,tol*tol>=offsetbuf*offsetbuf);
    if((min_offset!=(double)MAX_INTEGER)&&(tol*tol>=offsetbuf*offsetbuf))
    {
      if(printFlag) 
        printf("closest tri3:%d (%d %d %d)for node:%d, offset:%lf\n", elbuf+1, ctri3[elbuf].nod[0],ctri3[elbuf].nod[1],ctri3[elbuf].nod[2], set[set1].node[i], offsetbuf);
      if((dep =(int *)realloc((int *)dep, (sum_dep+2) *sizeof(int)))==NULL)
      { errMsg("\n\n ERROR: realloc failed\n" ); return(-1); }
      dep[sum_dep]=set[set1].node[i];
      sum_dep++;
      if((ind =(int *)realloc((int *)ind, (sum_ind+4) *sizeof(int)))==NULL)
      { errMsg("\n\n ERROR: realloc failed\n" ); return(-1); }
      ind[sum_ind]=ctri3[elbuf].nod[0];
      sum_ind++;
      ind[sum_ind]=ctri3[elbuf].nod[1];
      sum_ind++;
      ind[sum_ind]=ctri3[elbuf].nod[2];
      sum_ind++;
    }
    next_node:;
  }
  
  if (sum_dep)
  {
    sum_n[0]=sum_dep;
    sum_n[1]=sum_ind;
    nodesets[0]=dep;
    nodesets[1]=ind;
  }
  else
  {
    if(dep) free(dep);
    if(ind) free(ind);
    sum_n[0]=0;
    sum_n[1]=0;
    nodesets[0]=NULL;
    nodesets[1]=NULL;
  }

  
  if(ntri_nodes) free(ntri_nodes);
  if(ctri3) free(ctri3);
  if(rsort) free(rsort);
  if(orig_x) free(orig_x);
  if(orig_y) free(orig_y);
  if(orig_z) free(orig_z);
  if(sort_x) free(sort_x);
  if(sort_y) free(sort_y);
  if(sort_z) free(sort_z);
  if(sort_nx) free(sort_nx);
  if(sort_ny) free(sort_ny);
  if(sort_nz) free(sort_nz);

  if (sum_dep) return(1);
  else return(-1);
}



/* search face-pairs of unconnected sets of faces */
/* deside which side is the dep side */

int getFacePair( Summen *anz, int set1, int set2, Elements *e_enqire, Faces *face, Nodes *node, double tol, int *mpcset, char *par1)
{
  int i,j,e,f,n,s,sset,mset;
  int snod_indx[10], depflag[4]={0,0,0,0};
  int *snod=NULL;
  int *nodesets[4]={NULL,NULL,NULL,NULL},*facesets[4]={NULL,NULL,NULL,NULL};
  int sum_n[4],sum_f[4];
  char setname[MAX_LINE_LENGTH];

  printf ("getFacePair of %s %s\n", set[set1].name, set[set2].name);

  getDepNodes( anz, set1, set2, e_enqire, face, node, tol, &sum_n[0], &nodesets[0]);
  getDepNodes( anz, set2, set1, e_enqire, face, node, tol, &sum_n[2], &nodesets[2]);
  if(!sum_n[0]&&!sum_n[2])
  {
    for(i=0; i<4; i++) { if(nodesets[i]!=NULL) free(nodesets[i]); }
    return(0);
  }

  // which node-set is suited as dep? (and add all complete defined faces.
  /* try to determine the dependent side */
  for(s=0; s<3; s+=2)
  {
    if(!sum_n[s]) continue;
    sum_f[s]=0;
    if((snod =(int *)realloc((int *)snod, (sum_n[s]+1) *sizeof(int)))==NULL)
    { errMsg("\n\n ERROR: realloc failed\n" ); return(-1); }
    for(i=0; i<sum_n[s]; i++) snod[i]=0;

    if(s==0) mset=set1;
    else mset=set2;
    /* add the face only if the complete face is defined */
    for(f=0; f<set[mset].anz_f; f++)
    {
      if (face[set[mset].face[f]].type == 7)       j = 3;  /* TRI3  */
      else if (face[set[mset].face[f]].type == 8)  j = 6;  /* TRI6  */
      else if (face[set[mset].face[f]].type == 9)  j = 4;  /* QUAD4 */
      else if (face[set[mset].face[f]].type == 10) j = 8;  /* QUAD8 */
      /* add the face only if the complete face is defined */
      e=0;
      for(n=0; n<sum_n[s]; n++)
      {
        for(i=0; i<j; i++)
	{
          if(face[set[mset].face[f]].nod[i]==nodesets[s][n]) { e++; snod_indx[i]=n; }
        }
        if(e==j)
        {
          for(i=0; i<j; i++) snod[snod_indx[i]]=1;
          if((facesets[s] =(int *)realloc((int *)facesets[s], (sum_f[s]+1) *sizeof(int)))==NULL)
          { errMsg("\n\n ERROR: realloc failed\n" ); return(-1); }
          facesets[s][sum_f[s]++]=set[mset].face[f];
	}
      }
    }

    /* if all nodes are marked as used then this set should be used as the dependent set */
    for(n=0; n<sum_n[s]; n++) if(!snod[n]) break;
    if(n==sum_n[s])
    {
      //printf(" set %d is qualified as dependent set\n",s );
      depflag[s]=1;
    }
  }
  free(snod);

  // judge which set is dep and indep
  // if the two potentially depsets (0,2) are both valid or not, then the one with more nodes is dep and the following set is indep (1,3)
  printf("nodes:%d %d %d %d    %d %d\n",sum_n[0],sum_n[1],sum_n[2],sum_n[3], depflag[0], depflag[2]);
  if(depflag[0]==depflag[2])
  {
    /* eventually switch dep and ind */
    if(sum_n[2]>sum_n[0])
    {
      printf("1dep set nodes:%d\n",sum_n[2]);
      sset=set2; mset=set1;
      sprintf( setname, "D%s_%s", &set[sset].name[1], &set[mset].name[1]); 
      delSet(setname);
      mpcset[0]=pre_seta(setname, "i", 0);
      for (n=0; n<sum_n[2]; n++) seta(mpcset[0], "n", nodesets[2][n]);
      for (n=0; n<sum_f[2]; n++) seta(mpcset[0], "f", facesets[2][n]);

      printf("ind set nodes:%d\n",sum_n[3]);
      sprintf( setname, "I%s_%s", &set[mset].name[1], &set[sset].name[1]); 
      delSet(setname);
      mpcset[1]=pre_seta(setname, "i", 0);
      for (n=0; n<sum_n[3]; n++) seta(mpcset[1], "n", nodesets[3][n]);
    }
    else
    {
      printf("2dep set nodes:%d\n",sum_n[0]);
      sset=set1; mset=set2;
      sprintf( setname, "D%s_%s", &set[sset].name[1], &set[mset].name[1]); 
      delSet(setname);
      mpcset[0]=pre_seta(setname, "i", 0);
      for (n=0; n<sum_n[0]; n++) seta(mpcset[0], "n", nodesets[0][n]);
      for (n=0; n<sum_f[0]; n++) seta(mpcset[0], "f", facesets[0][n]);

      printf("ind set nodes:%d\n",sum_n[1]);
      sprintf( setname, "I%s_%s", &set[mset].name[1], &set[sset].name[1]); 
      delSet(setname);
      mpcset[1]=pre_seta(setname, "i", 0);
      for (n=0; n<sum_n[1]; n++) seta(mpcset[1], "n", nodesets[1][n]);
    }
  }
  // which depsets (0,2) is valid?
  else
  {
    if(depflag[0]==1) { s=0;  sset=set1; mset=set2; }
    else { s=2;  sset=set2; mset=set1; }
      printf("3dep set nodes:%d\n",sum_n[s]);
      sprintf( setname, "D%s_%s", &set[sset].name[1], &set[mset].name[1]); 
      delSet(setname);
      mpcset[0]=pre_seta(setname, "i", 0);
      for (n=0; n<sum_n[s]; n++) seta(mpcset[0], "n", nodesets[s][n]);
      for (n=0; n<sum_f[s]; n++) seta(mpcset[0], "f", facesets[s][n]);

    s++;
      printf("ind set nodes:%d\n",sum_n[s]);
      sprintf( setname, "I%s_%s", &set[mset].name[1], &set[sset].name[1]); 
      delSet(setname);
      mpcset[1]=pre_seta(setname, "i", 0);
      for (n=0; n<sum_n[s]; n++) seta(mpcset[1], "n", nodesets[s][n]);
  }
  printf(" found in set:%s %d dependent nodes and %d independent nodes in set:%s\n", set[mpcset[0]].name, set[mpcset[0]].anz_n,set[mpcset[1]].anz_n, set[mpcset[1]].name);

  // add the faces to the indep set
  for(f=0; f<set[mset].anz_f; f++)
  {
    j=0;
    switch(face[set[mset].face[f]].type)
    {
      case 7: j = 3; break;
      case 8: j = 6; break;
      case 9: j = 4; break;
      case 10: j = 8; break;
    }    
    s=0;
    for(i=0; i<j; i++)
    {
      for (n=0; n<set[mpcset[1]].anz_n; n++)
      {
        if(face[set[mset].face[f]].nod[i]==set[mpcset[1]].node[n]) { s++; break; }
      }
      if(s>=3)
      { seta(mpcset[1], "f", set[mset].face[f]); break; }
    }
  }

  for(i=0; i<4; i++) { if(nodesets[i]!=NULL) free(nodesets[i]); if(facesets[i]!=NULL) free(facesets[i]); }
  if((set[mpcset[0]].anz_f)&&(set[mpcset[1]].anz_f)) return(2);
  else if((set[mpcset[0]].anz_n)&&(set[mpcset[1]].anz_f)) return(1);
  return(0);
}
