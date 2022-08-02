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

#include <extUtil.h>

#define TEST     0

#define INI_FIELD_SIZE 100000

#define UINT8  char
#define UINT16  short int
#define UINT32  int
#define REAL32  float

int readStl( char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr )
{
  FILE *handle;
  int i=0,j;

  char rec_str[MAX_LINE_LENGTH];
  char buffer[MAX_LINE_LENGTH];
  int  node_field_size, elem_field_size;
  int  e_nmax=1, e_nmin=1, sum_tri=0, attr=0;
  int  length;

  float nv[3],np[3][3];

  Nodes     *node=NULL;
  Elements  *elem=NULL;

  int  HEAD_CHARS=80;



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

  anz->n=anz->e=anz->l=-1;
  anz->emax=0;  anz->emin=MAX_INTEGER;
  anz->nmax=0;  anz->nmin=MAX_INTEGER;

  /* Open the files and check to see that it was opened correctly */
  handle = fopen (datin, "r");
  if ( handle== NULL )  { printf ("ERROR: The input file \"%s\" could not be opened.\n\n", datin); return(-1); }
  else  printf (" file:%s opened\n", datin);

  printf (" reading stl format\n");

  length = frecord( handle, rec_str);
  sscanf(rec_str, "%s", buffer);
  for(i=0;i<strlen(buffer); i++) buffer[i]=toupper(buffer[i]);
  if(compare(buffer, "SOLID", 5)==5)
  {
    strcpy(anz->model, rec_str);
    printf (" MODEL NAME:  %s", anz->model);

    while(length)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
  
      sscanf(rec_str, "%s", buffer);
      for(i=0;i<strlen(buffer); i++) buffer[i]=toupper(buffer[i]);
      if(compare(buffer, "FACET", 5)==5)
      {
        /* overread "outer loop" */
        length = frecord( handle, rec_str);
        for(i=0; i<3; i++)
        {
          anz->n++;
          node[anz->n].nr = anz->n+1;
          if (node[anz->n].nr>=node_field_size)
          {
            node_field_size=node[anz->n].nr+100;
            if ( (node = (Nodes *)realloc((Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
            {
              printf("\n\n ERROR: realloc failed, nodenr:%d\n\n", node[anz->n].nr) ;
              return(-1);
            }
          }
          node[node[anz->n].nr].indx=anz->n;
  
          length = frecord( handle, rec_str);
          sscanf(rec_str,"%*s %lf %lf %lf", &node[node[anz->n].nr].nx,&node[node[anz->n].nr].ny,&node[node[anz->n].nr].nz);
          if (node[anz->n].nr >  anz->nmax)  anz->nmax=node[anz->n].nr;
          if (node[anz->n].nr <  anz->nmin)  anz->nmin=node[anz->n].nr;
#if TEST
          printf (" n=%d x=%lf y=%lf z=%lf \n",  node[anz->n].nr,
          node[node[anz->n].nr].nx, node[node[anz->n].nr].ny,
          node[node[anz->n].nr].nz); 
#endif 
        }
        anz->e++;
        if (anz->e>=elem_field_size)
        {
          elem_field_size=anz->e+100;
          if((elem=(Elements *)realloc((Elements *)elem,(elem_field_size+1)*sizeof(Elements)))==NULL)
          {
            printf("\n\n ERROR: realloc failed, elem-index:%d\n\n", anz->e);
            return(-1);
          }
        }
        elem[anz->e].nr    = anz->e+1;
        elem[anz->e].type  = 7;
        elem[anz->e].group = 1;
        elem[anz->e].mat   = 1;
        anz->etype[elem[anz->e].type]++;
        if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
        if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
#if TEST
        printf (" e=%d typ=%d mat=%d \n", elem[anz->e].nr,
        elem[anz->e].type, elem[anz->e].group, elem[anz->e].mat );
#endif
        for (i=0; i<3; i++) elem[anz->e].nod[i]=node[anz->n-2+i].nr;
      }
    }
  }
  else
  {
    printf ("\n\n: No ascii file, trying to read binary file-format.\n\n");
    strcpy(anz->model, "binstl");
    fclose(handle);
    handle = fopen (datin, "rb");

    // HEAD_CHARS chars at the beginning of the binary file
    length=fread(buffer,sizeof(UINT8),HEAD_CHARS,handle);
    // nr of triangles
    length=fread(&sum_tri,sizeof(UINT32),1,handle);
    printf("header:%s\n",buffer);
    printf("triangles:%d\n",sum_tri);

    if ( (node = (Nodes *)realloc((Nodes *)node, (3*sum_tri+1) * sizeof(Nodes))) == NULL )
    {
      printf("\n\n ERROR: realloc failed\n\n") ;
      return(-1);
    }
    if((elem=(Elements *)realloc((Elements *)elem,(sum_tri+1)*sizeof(Elements)))==NULL)
    {
      printf("\n\n ERROR: realloc failed\n\n");
      return(-1);
    }

    for(i=0; i<sum_tri; i++)
    {
      length=fread(nv, sizeof(REAL32), 3, handle);
      if(length<1) break;
      //printf("nv:%f %f %f\n",nv[0],nv[1],nv[2]);
      length=fread(&np[0][0], sizeof(REAL32), 3, handle);
      if(length<1) break;
      //printf("n1:%f %f %f\n",np[0][0],np[0][1],np[0][2]);
      length=fread(&np[1][0], sizeof(REAL32), 3, handle);
      if(length<1) break;
      //printf("n1:%f %f %f\n",np[1][0],np[1][1],np[1][2]);
      length=fread(&np[2][0], sizeof(REAL32), 3, handle);
      if(length<1) break;
      //printf("n1:%f %f %f\n",np[2][0],np[2][1],np[2][2]);
      length=fread(&attr, sizeof(UINT16), 1, handle);
      //printf("attr:%d\n",attr);

      for(j=0; j<3; j++)
      {
        anz->n++;
        node[anz->n].nr = anz->n+1;
        node[node[anz->n].nr].indx=anz->n;
        node[node[anz->n].nr].nx=np[j][0];
        node[node[anz->n].nr].ny=np[j][1];
        node[node[anz->n].nr].nz=np[j][2];
      }
      anz->e++;
      elem[anz->e].nr    = anz->e+1;
      elem[anz->e].type  = 7;
      elem[anz->e].group = 1;
      elem[anz->e].mat   = 1;
      anz->etype[elem[anz->e].type]++;
      for (i=0; i<3; i++) elem[anz->e].nod[i]=node[anz->n-2+i].nr;
    }
    anz->nmax=anz->n+1;
    anz->nmin=1;
    anz->emax=anz->e+1;
    anz->emin=1;
    fclose(handle);
  }

  anz->n++;
  anz->e++;
  anz->l++;
  fclose(handle);

  node_field_size=anz->nmax+1;
  if((node =(Nodes *)realloc( (Nodes *)node, node_field_size*sizeof(Nodes)))==NULL)
    printf("\n\n ERROR: realloc failed\n\n") ;
  else
    printf ("\n %d nodes reallocated \n",anz->nmax);

  elem_field_size=anz->e+1;
  if ( (elem = (Elements *)realloc((Elements *)elem, elem_field_size * sizeof(Elements))) == NULL )
    printf("\n\n ERROR: in readfrd realloc failed\n\n") ;
  else
    printf ("\n %d elements reallocated \n", anz->e);

  if ( e_nmax > (anz->nmax) )
  {
    printf ("\nWARNING: element requestes a nodename higher than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }
  if ( e_nmin < 1 )
  {
    printf ("\nWARNING: element requestes a nodename lower than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }

  *nptr = node; *eptr = elem;
  return(1);
}

