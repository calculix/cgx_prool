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

extern CopiedNodeSets copiedNodeSets[1];
extern Nodes    *node;                                                  
extern Faces    *face;                                            
extern Elements *e_enqire;                                               
extern Datasets *lcase;
extern Sets     *set;
extern Points   *point;
extern Lines    *line ;
extern Lcmb     *lcmb ;
extern Gsur     *surf ;
extern Gbod     *body ;
extern Nurbl    *nurbl;
extern Nurbs    *nurbs;
extern Shapes   *shape;


      
void dom(int sum_n)
{
  int i,j,f,n1,n2;
  int **edge, setNr, setNr2, setNr3,sumSets, set_maxe=0;
  char buffer[MAX_LINE_LENGTH];
  setNr=0;
  setNr2=pre_seta("-failedElem","i",0);
  setNr3=pre_seta("dom","i",0);

  
  if( (edge = (int **)malloc((sum_n)*sizeof(int *) ))==NULL )
    { printf(" ERROR: realloc failure\n\n"); return; }
  for (i=0; i<sum_n; i++)
  {
    if( (edge[i] = (int *)calloc((sum_n),sizeof(int) ))==NULL )
    { printf(" ERROR: realloc failure\n\n"); return; }
  }
  // search volume elements which are connected to each other by just one or two corner nodes
  for(i=0; i<set[setNr].anz_f; i++)
  {
    f=set[setNr].face[i];
    if((face[f].type == 7)||(face[f].type == 8))
    {
      n1=node[face[f].nod[0]].indx;
      n2=node[face[f].nod[1]].indx;
      edge[n1][n2]++;
      //printf("e:%d n %d %d ed:%d\n", face[f].elem_nr,  node[n1].nr,  node[n2].nr, edge[n1][n2]);
      if(edge[n1][n2]>1) seta(setNr2,"e",face[f].elem_nr);
      n1=node[face[f].nod[1]].indx;
      n2=node[face[f].nod[2]].indx;
      edge[n1][n2]++;
      //printf("e:%d n %d %d ed:%d\n", face[f].elem_nr,  node[n1].nr,  node[n2].nr, edge[n1][n2]);
      if(edge[n1][n2]>1) seta(setNr2,"e",face[f].elem_nr);
      n1=node[face[f].nod[2]].indx;
      n2=node[face[f].nod[0]].indx;
      edge[n1][n2]++;
      //printf("e:%d n %d %d ed:%d\n", face[f].elem_nr,  node[n1].nr,  node[n2].nr, edge[n1][n2]);
      if(edge[n1][n2]>1) seta(setNr2,"e",face[f].elem_nr);
    }
  }
  for (i=0; i<sum_n; i++) free(edge[i]); free(edge);

  printf(" found %d edges\n", set[setNr2].anz_e/2);

  completeSet(set[setNr2].name, "do") ;
  completeSet(set[setNr2].name, "up") ;
  for(i=0; i<set[setNr].anz_e; i++) seta(setNr3,"e",set[setNr].elem[i]);
  for(i=0; i<set[setNr2].anz_e; i++) setr(setNr3,"e",set[setNr2].elem[i]);
  sumSets=separateMeshes(set[setNr3].name, "+dom");

  for(i=0; i<sumSets; i++)
  {
    sprintf(buffer,"+dom%d", i+1);
    j=getSetNr(buffer);
    if(j>-1) if(set[j].anz_e>set_maxe) set_maxe=j;
  }
  completeSet(set[set_maxe].name, "do") ;
  printf("new dom in set:%s stored in 'valu dom'\n",set[set_maxe].name);
  sprintf(buffer,"dom %s", set[set_maxe].name);
  pre_value(buffer);

  // remove elems from the domain
  //if(set[setNr2].anz_f) { for(i=0; i<set[setNr2].anz_f; i++) setr(setNr,"e",face[set[setNr2].face[i]].elem_nr); delSet("-failedElem"); }
}

/* --------------------------------------------------------------------  */
/* Userfunctions                                                         */
/* interface to the mesh, geometry and datasets                          */
/*                                                                       */
/* string:  parameter list from command line                             */
/* sum:     mesh related number of entities                              */
/* sumGeo:  geometrie related number of entities                         */
/*                                                                       */
/* REMARK:                                                               */
/* If you intend to create additional elements then be aware that you    */
/* must use the following commands before you do that:                   */
/* //free the additional midside-nodes for higher order elements         */
/*   for(i=anz->orign; i<anz->n; i++) node[node[i].nr].pflag=-1;         */
/*   anz->n= anz->orign;                                                 */
/*   anz->nmax=anz->orignmax;                                            */
/* And after the elements are created:                                   */
/*   adjustDrawNodes(1);                                                 */
/*                                                                       */
/* --------------------------------------------------------------------  */
                                                                
void userFunction(char *string, Summen   *sum, SumGeo   *sumGeo )
{                                                               
  int i,l,n;
  static int setNr=0;
                                    
  FILE *handle;
  char buf1[MAX_LINE_LENGTH], buf2[MAX_LINE_LENGTH], buf3[MAX_LINE_LENGTH], rec_str[MAX_LINE_LENGTH];
  int e, length, sum_lc, nodnr;
  int lmin=0, lmax=MAX_INTEGER;
  double Fres, tref=0.;
  double vx[100], vy[100], val=0.;


  /* list the implemented user-functions if no parameter was provided for the "call" command */
  if(strlen(string)==0)                         
  {
    printf("  abs <dataset-nr> <dataset-nr> // changes the entities to abs(value) between both datasets\n");
    printf("  add <set> <dataset-nr> <entity-nr>  // just sums up all node-vals (ie. RF)\n");
    printf("  dom  // eliminates linear element connections and stores the biggest element block in valu 'dom'\n");
    printf("  hydro // hydrostatic pressure\n");
    printf("  move <set> <file> //moves nodes in set by interpolated values from <file> \n");
    printf("         1st line defines operation (column1-descriptor(xyz) column2-descriptor(xyz) operator(*,+)\n");
    printf("         all other lines: data data\n");
    printf("  bend\n"); 
  }
  else if(compare(string, "abs", 3)==3)                         
  {                                                             
    sscanf(string, "%s %d %d", buf1, &lmin, &lmax);

    for(l=0; l<sum->l; l++) if((l+1>=lmin)&&(l+1<=lmax))                              
    {                                                           
      printf("change value to abs(value) for dataset %d\n", l+1);
      /* check if the data of the specified lcase (Dataset) are already available */
      if (!lcase[l].loaded)
      {
       if( pre_readfrdblock(copiedNodeSets , l, sum, node, lcase )==-1) 
       {
         printf("ERROR in userfunction: Could not read data for Dataset:%d\n", l+1); 
         return;
       }
       calcDatasets( l, sum, node, lcase );
       recompileEntitiesInMenu(l);
      }
      for(e=0; e<lcase[l].ncomps; e++)
      {
	printf("change value to abs(value) for dataset %d e:%d\n", l+1, e+1);
        lcase[l].max[e]=-MAX_INTEGER;
        lcase[l].min[e]=MAX_INTEGER;
        for(n=0; n<sum->n; n++)
        {
          lcase[l].dat[e][node[n].nr]=abs(lcase[l].dat[e][node[n].nr]);
          if (lcase[l].dat[e][node[n].nr] >  lcase[l].max[e])
          {  lcase[l].max[e]=lcase[l].dat[e][node[n].nr]; lcase[l].nmax[e]=node[n].nr;}
          if (lcase[l].dat[e][node[n].nr] <  lcase[l].min[e])
          {  lcase[l].min[e]=lcase[l].dat[e][node[n].nr]; lcase[l].nmin[e]=node[n].nr;}
        }
      }
    }
  }
  else if(compare(string, "add", 3)==3)                         
  {                                                             
    sscanf(string, "%*s %s %d %d", buf1, &l, &e);
    printf("set:%s ds:%d e:%d\n", buf1, l,e);
    l--; e--;


    setNr=getSetNr(buf1);
    if (setNr<0)
    {
      printf (" calcFres: set:%s does not exist\n", buf1);
      return;
    }
    if(l>=sum->l)
    {
      printf (" calcFres: dataset:%d does not exist\n", l+1);
      return;
    }
    if(e>=lcase[l].ncomps)
    {
      printf (" calcFres: entity:%d does not exist\n", e+1);
      return;
    }

    /* check if the data of the specified lcase (Dataset) are already available */
    if (!lcase[l].loaded)
    {
     if( pre_readfrdblock(copiedNodeSets , l, sum, node, lcase )==-1) 
     {
       printf("ERROR in userfunction: Could not read data for Dataset:%d\n", l+1); 
       return;
     }
     calcDatasets( l, sum, node, lcase );
     recompileEntitiesInMenu(l);
    }

    Fres=0.;
    for(i=0; i<set[setNr].anz_n; i++)
    {
        Fres+=lcase[l].dat[e][set[setNr].node[i]];
    }
    printf("\n Fres:%lf\n", Fres);
  }
  else if(compare(string, "dom", 3)==3)                         
  {
    dom(sum->n);                                                             
  }
  else if(compareStrings(string, "hydro")>0)                         
  {                                                             
    /* calculate the hydrostatic stress and the deviator */     
    for(l=0; l<sum->l; l++)                                     
    {                                                           
      /* use only stresses */
      if (compare(lcase[l].name,"STRESS",6)==6)
      {
        printf(" calc hydrostatic stress and deviator for LC[%d]: %s \n", l, lcase[l].name); 

        /* check if the data of the specified lcase (Dataset) are already available */
        if (!lcase[l].loaded)
        {
         if( pre_readfrdblock(copiedNodeSets , l, sum, node, lcase )==-1) 
         {
           printf("ERROR in userfunction: Could not read data for Dataset:%d\n", l+1); 
           return;
         }
         calcDatasets( l, sum, node, lcase );
         recompileEntitiesInMenu(l);
        }
  
        /* create a new dataset */
        generateDataset(sum,&lcase,"dummy",4,lcase[l].value,"",1,1,"");
        sprintf(buf1,"%f",tref);
        createDSparam(sum->l-1,"TREF",buf1);
        sum->l--;

        sprintf( lcase[sum->l].name, "HYD&DEV");
        strcpy( lcase[sum->l].compName[0], "SM");
        strcpy( lcase[sum->l].compName[1], "SXX-SM");
        strcpy( lcase[sum->l].compName[2], "SYY-SM");
        strcpy( lcase[sum->l].compName[3], "SZZ-SM");
  
        for(n=0; n<sum->n; n++)
        {
          /* calc average Princ Stress sigm */
          lcase[sum->l].dat[0][node[n].nr]=
            (lcase[l].dat[8][node[n].nr]+lcase[l].dat[9][node[n].nr]+lcase[l].dat[10][node[n].nr])/3.;
         
          /* calc the deviator */
          for(i=0; i<3; i++)
            lcase[sum->l].dat[i+1][node[n].nr]=lcase[l].dat[i][node[n].nr]-lcase[sum->l].dat[0][node[n].nr];

          for(i=0; i<lcase[sum->l].ncomps; i++)
          {
            if (lcase[sum->l].dat[i][node[n].nr] >  lcase[sum->l].max[i])
            {  lcase[sum->l].max[i]=lcase[sum->l].dat[i][node[n].nr]; lcase[sum->l].nmax[i]=node[n].nr; }
            if (lcase[sum->l].dat[i][node[n].nr] <  lcase[sum->l].min[i])
            {  lcase[sum->l].min[i]=lcase[sum->l].dat[i][node[n].nr]; lcase[sum->l].nmin[i]=node[n].nr; }
          }
        }
        sum->l++;
      }
    }
  }
  else if(compare(string, "move", 4)==4)                         
  {                                                             
    sscanf(string, "%*s %s %s", buf1, buf2);
    printf("set:|%s| file:|%s|\n", buf1, buf2);

    setNr=getSetNr(buf1);
    if (setNr<0)
    {
      printf (" ERROR: set:%s does not exist\n", buf1);
      return;
    }

    handle = fopen (buf2, "r");
    if ( handle== NULL )  { printf ("\nThe input file \"%s\" could not be opened.\n\n", buf2); return; }
    else  printf ("\n%s opened\n\n",buf2);

    length = frecord( handle, rec_str);
    sscanf(rec_str, "%s %s %s", buf1, buf2, buf3);
    printf("   data %s, move in %s by operator:%s \nData:\n", buf1, buf2, buf3);

    i=0;
    do
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
      else rec_str[length] =(char)0;
      if (!length) break;
      //printf ("record:%s\n", rec_str);
      sscanf( rec_str, "%lf %lf", &vx[i], &vy[i]);
      printf( "%lf %lf\n", vx[i], vy[i]);
      i++;
    }while(length);
    sum_lc=i;

    /* go over all nodes in set and move them by operator buf3 */
    for(i=0; i<set[setNr].anz_n; i++)
    {
      nodnr=set[setNr].node[i];

      if(buf1[0]=='x') val= intpol( vx, vy, sum_lc, node[nodnr].nx );
      if(buf1[0]=='y') val= intpol( vx, vy, sum_lc, node[nodnr].ny );
      if(buf1[0]=='z') val= intpol( vx, vy, sum_lc, node[nodnr].nz );
      if(buf2[0]=='x')
      {
        if(buf3[0]=='*') node[nodnr].nx*=val;
        else if(buf3[0]=='+') node[nodnr].nx+=val;
        else { printf("ERROR: Operator:%s not known\n", buf3); return; }
      }
      if(buf2[0]=='y')
      {
        if(buf3[0]=='*') node[nodnr].ny*=val;
        else if(buf3[0]=='+') node[nodnr].ny+=val;
        else { printf("ERROR: Operator:%s not known\n", buf3); return; }
      }
      if(buf2[0]=='z')
      {
        if(buf3[0]=='*') node[nodnr].nz*=val;
        else if(buf3[0]=='+') node[nodnr].nz+=val;
        else { printf("ERROR: Operator:%s not known\n", buf3); return; }
      }
      //printf("n:%d dy:%f y:%f\n",nodnr, val,node[nodnr].ny);
    }
  }
  else printf(" ERROR userfunction:%s not known\n", string);

}
