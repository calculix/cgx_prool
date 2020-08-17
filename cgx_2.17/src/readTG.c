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
Jan 2016 Peter Heppel added an optional parameter to read <file> ng 
 new parameter is ndsb (NoDeleteShellsorBeams) if  you want to keep them.
 (It would be better if this parameter contained vsb  for read volume ,
 shell or beams)

Jan 2016 Peter Heppel added a reader for quad8 elements - 
although I'm guessing that .mat should be 10
*/

#include <cgx.h>

#define TEST     0

#define INI_FIELD_SIZE 100000

extern Sets      *setx;
extern Summen    *anzx;
extern Alias     *alias;
extern SumGeo    anzGeo[1];
extern SumAsci   sumAsci[1];

int readTG( char *datin, Summen *apre, Sets **sptr, Nodes **nptr, Elements **eptr, Datasets **lptr )
{
  FILE *handle;
  int i=0,sp;

  char rec_str[MAX_LINE_LENGTH], buffer[MAX_LINE_LENGTH], name[MAX_LINE_LENGTH];
  int  node_field_size, elem_field_size;
  int  e_nmax=1, e_nmin=1;
  int  length, sum,n;

  Nodes     *node=NULL;
  Elements  *elem=NULL;

  anzx=apre;
  setx=*sptr;

  node_field_size=INI_FIELD_SIZE;
  do
  {
    if ( (node = (Nodes *)realloc( (Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
    {
      printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", node_field_size );
      node_field_size/=2;
    }
    if(node_field_size<100)
    {
      printf("\n\n ERROR: not enough memory in readfrd()\n\n");
      exit(-1);
    }
  }while(!node);

  elem_field_size=INI_FIELD_SIZE;
  do
  {
    if((elem = (Elements *)realloc( (Elements *)elem, (elem_field_size+1) * sizeof(Elements))) == NULL )
    {
      printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", elem_field_size );
      elem_field_size/=2;
    }
    if(elem_field_size<100)
    {
      printf("\n\n ERROR: not enough memory in readfrd()\n\n");
      exit(-1);
    }
  }while(!elem);


  /* Open the files and check to see that it was opened correctly */
  
  printf (" reading Tetgen format\n");

  sp=strlen(datin);
  while((datin[--sp]!='.')&&(sp>0));
  datin[sp]=0;

  strcpy(anzx->model, datin);
  printf (" MODEL NAME:  %s", anzx->model);
  
  /* nodes */
  sprintf(&datin[sp],".node");
  handle = fopen (datin, "r");
  if ( handle== NULL )  { printf ("ERROR: The input file \"%s\" could not be opened.\n\n", datin); return(-1); }
  else  printf (" file:%s opened\n", datin);

  length = frecord( handle, rec_str);
  sscanf(rec_str, "%d", &sum);
  for(i=0; i<sum; i++)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;
    sscanf(rec_str,"%d", &node[anzx->n].nr);
    //printf("node:%d\n",node[anzx->n].nr);
    //node[anzx->n].nr = anzx->n+1;
    if (node[anzx->n].nr>=node_field_size)
    {
      node_field_size=node[anzx->n].nr+100;
      if ( (node = (Nodes *)realloc((Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
      {
        printf("\n\n ERROR: realloc failed, nodenr:%d\n\n", node[anzx->n].nr) ;
        return(-1);
      }
    }
    node[node[anzx->n].nr].indx=anzx->n;

    sscanf(rec_str,"%*d %lf %lf %lf", &node[node[anzx->n].nr].nx,&node[node[anzx->n].nr].ny,&node[node[anzx->n].nr].nz);
    if (node[anzx->n].nr >  anzx->nmax)  anzx->nmax=node[anzx->n].nr;
    if (node[anzx->n].nr <  anzx->nmin)  anzx->nmin=node[anzx->n].nr;
    anzx->n++;
  }

  /* elems */
  sprintf(&datin[sp],".ele");
  handle = fopen (datin, "r");
  if ( handle== NULL )  { printf ("ERROR: The input file \"%s\" could not be opened.\n\n", datin); return(-1); }
  else  printf (" file:%s opened\n", datin);
  
  do{ length = frecord( handle, rec_str); }while(rec_str[0]=='#');
  sscanf(rec_str, "%d %d", &sum, &n );
  for(i=0; i<sum; i++)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;
    sscanf(rec_str,"%d", &elem[anzx->e].nr);
    //printf("elem:%d\n",elem[anzx->e].nr);

    if (anzx->e>=elem_field_size)
    {
      elem_field_size=anzx->e+100;
      if((elem=(Elements *)realloc((Elements *)elem,(elem_field_size+1)*sizeof(Elements)))==NULL)
      {
        printf("\n\n ERROR: realloc failed, elem-index:%d\n\n", anzx->e);
        return(-1);
      }
    }
    elem[anzx->e].group = 1;
    elem[anzx->e].type  = 0;
    if (elem[anzx->e].nr >  anzx->emax)  anzx->emax=elem[anzx->e].nr;
    if (elem[anzx->e].nr <  anzx->emin)  anzx->emin=elem[anzx->e].nr;
    if(n==4)
    {
      elem[anzx->e].mat   = 3;
      elem[anzx->e].type  = 3;
      sscanf(rec_str, "%*d %d %d %d %d", &elem[anzx->e].nod[0], &elem[anzx->e].nod[1], &elem[anzx->e].nod[2], &elem[anzx->e].nod[3] );
    }
    else printf("elem-type not known, nr of nodes:%d\n", n);
    anzx->etype[elem[anzx->e].type]++;
    anzx->e++;
  }


  fclose(handle);

  node_field_size=anzx->nmax+1;
  if((node =(Nodes *)realloc( (Nodes *)node, node_field_size*sizeof(Nodes)))==NULL)
    printf("\n\n ERROR: realloc failed\n\n") ;
  else
    printf ("\n %d nodes reallocated \n",anzx->nmax);

  elem_field_size=anzx->e+1;
  if ( (elem = (Elements *)realloc((Elements *)elem, elem_field_size * sizeof(Elements))) == NULL )
    printf("\n\n ERROR: in readfrd realloc failed\n\n") ;
  else
    printf ("\n %d elements reallocated \n", anzx->e);

  if ( e_nmax > (anzx->nmax) )
  {
    printf ("\nWARNING: element requestes a nodename higher than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }
  if ( e_nmin < 1 )
  {
    printf ("\nWARNING: element requestes a nodename lower than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }


  for (i=0; i<anzx->e; i++)
  {
    sprintf(name,"%d",  elem[i].nr);
    sprintf(buffer, "+set%d", elem[i].mat );
    pre_setax( buffer, "e", name );
  }
  for (i=0; i<anzx->e; i++)
  {
    sprintf(name,"%d", elem[i].nr);
    sprintf(buffer, "+typ%d", elem[i].type );
    pre_setax( buffer, "e", name );
  }

  *sptr = setx;
  *nptr = node;
  *eptr = elem;
  return(1);
}

