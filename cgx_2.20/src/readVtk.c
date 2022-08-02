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

#define TEST     1

#define INI_FIELD_SIZE 100000

#define UINT8  char
#define UINT16  short int
#define UINT32  int
#define REAL32  float

int readVtk( char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr )
{
  FILE *handle;
  fpos_t *filepntr=NULL;
  int i=0,j,k;

  char n,rec_str[MAX_LINE_LENGTH];
  char buffer[MAX_LINE_LENGTH];
  char format[MAX_LINE_LENGTH];
  char key[MAX_LINE_LENGTH];
  int  node_field_size, elem_field_size;
  int  e_nmax=1, e_nmin=1, sum_n=0, sum_typ=0, sum_e=0;
  int  length;

  double version=0;
  int *eltype=NULL;
  double nv[3],np[3][3];

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

  length = frecord( handle, rec_str);
  for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
  printf("%s\n", buffer);
  sscanf(buffer,"# %*s %*s %s %lf",key,&version);
  if(compare(key, "VERSION", 5)!=5)
  {
    printf (" it seems not to be vtk format, key:%s\n",key);
    return(-1);
  }
  else
    printf (" reading vtk format\n");

  length = frecord( handle, rec_str);
  printf("header %s\n", rec_str);
  strcpy(anz->model, rec_str);

  length = frecord( handle, rec_str);
  for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
  if(compare(rec_str, "ASCII", 5)!=5)
  {
    printf (" it seems not to be ascii format, key:%s\n",key);
    return(-1);
  }

  length = frecord( handle, rec_str);
  for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
  sscanf(buffer,"%*s %s",key);
  if(compare(key, "UNSTRUCTURED_GRID", 17)!=17)
  {
    printf (" it seems not to be an unstructured grid, key:%s\n",key);
    return(-1);
  }
  else
  {
    // first read the file and look for CELL_TYPES
    /* store the beginning of the data-block for later reading */
    if(filepntr==NULL)
    {  if( (filepntr=(fpos_t *)malloc(1*sizeof(fpos_t))) == NULL ) printf(" ERROR: malloc failed\n"); }
    if(fgetpos( handle, (fpos_t *)filepntr)!=0) { printf("error in fgetpos"); return(-1); }
    while(length)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
       printf (" pre read line:%s\n",rec_str);
  
      for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
      if(compare(buffer, "CELL_TYPES", 9)==9)
      {
	buffer[0]=0;
        sscanf(rec_str, "%*s %d", &sum_typ);
	printf("reading %d CELL_TYPES\n",sum_typ);
        if ( (eltype = (int *)malloc( (sum_typ) * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
	
        for(i=0; i<sum_typ; i++)
        {
          fscanf(handle, "%d", &k);
	  if(k==9) eltype[i]=9;
	  if(k==10) eltype[i]=3;
	}
        break;
      }
    }
    if( fsetpos( handle, (fpos_t *)filepntr)!=0) { printf("error in fsetpos"); return(-1); }

    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
       printf (" new line:%s\n",rec_str);
  
      sscanf(rec_str, "%s", buffer);
      for(i=0;i<=strlen(buffer); i++) buffer[i]=toupper(buffer[i]);
      if(compare(buffer, "POINTS", 6)==6)
      {
	buffer[0]=0;
        sscanf(rec_str, "%*s %d %s", &sum_n, format);
	printf("reading %d POINTS in format:%s\n",sum_n, format);
	
        for(i=0; i<sum_n; i++)
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
  
          do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
          k=0; do{ buffer[k++]=n=getc(handle); printf("%c",n); }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	  node[node[anz->n].nr].nx=atof(buffer);
          do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
          k=0; do{ buffer[k++]=n=getc(handle); printf("%c",n); }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	  node[node[anz->n].nr].ny=atof(buffer);
          do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
          k=0; do{ buffer[k++]=n=getc(handle); printf("%c",n); }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	  node[node[anz->n].nr].nz=atof(buffer);
	  
          if (node[anz->n].nr >  anz->nmax)  anz->nmax=node[anz->n].nr;
          if (node[anz->n].nr <  anz->nmin)  anz->nmin=node[anz->n].nr;
#if TEST
          printf ("%d n=%d x=%lf y=%lf z=%lf \n", anz->n,  node[anz->n].nr,
          node[node[anz->n].nr].nx, node[node[anz->n].nr].ny,
          node[node[anz->n].nr].nz); 
#endif 
        }
	printf(" nodes allocated:%d\n",anz->n+1);
      }
      if(compare(buffer, "CELLS", 5)==5)
      {
        sscanf(rec_str, "%*s %d %d", &sum_e, &sum_n);
	printf("reading %d CELLS with data points:%d\n",sum_e, sum_n);

	if(sum_typ!=sum_e)
	{
          if(sum_typ==sum_e-1)
	  {
            length = frecord( handle, rec_str);
            for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
            if(compare(buffer, "OFFSETS", 7)!=7)
            {
  	      printf(" ERROR %d cells but only %d types specified\n",k,sum_e);
	      goto jumpElems;
	    }
            for(i=0; i<sum_e; i++)
            {
              do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
              k=0; do{ buffer[k++]=n=getc(handle); printf("%c",n); }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	    }
            do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
	    printf("\n");
            length = frecord( handle, rec_str);
            for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
            if(compare(buffer, "CONNECTIVITY", 7)!=7)
            {
  	      printf(" ERROR CONNECTIVITY not found\n");
	      goto jumpElems;
	    }
	    sum_e-=1;
	  }
	  printf(" ERROR %d cells but only %d types specified\n",k,sum_e);
          goto jumpElems;
	}
	  
        for(i=0; i<sum_e; i++)
        {
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
          elem[anz->e].type  = eltype[anz->e];
          elem[anz->e].group = 1;
          elem[anz->e].mat   = 1;
          anz->etype[elem[anz->e].type]++;
          if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
          if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
#if TEST
          printf (" e=%d typ=%d mat=%d \n", elem[anz->e].nr,
          elem[anz->e].type, elem[anz->e].group, elem[anz->e].mat );
#endif
          if(version>3.)
	  {
	    if((elem[anz->e].type==3)||(elem[anz->e].type==9))
	    {
	      for(j=0; j<4; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0; do{ buffer[k++]=n=getc(handle); printf("%c",n); }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		elem[anz->e].nod[j]=atoi(buffer)+1;
	      }
	    }	
	  }
	  else
	  {
            length = frecord( handle, rec_str);
            if (rec_str[length] == (char)EOF) break;
	  
	    if((elem[anz->e].type==3)||(elem[anz->e].type==9))
	    {
              sscanf(rec_str, "%*d %d %d %d %d", &elem[anz->e].nod[0], &elem[anz->e].nod[1], &elem[anz->e].nod[2], &elem[anz->e].nod[3]);
             for (i=0; i<4; i++) elem[anz->e].nod[i]+=1;
	    }
	  }
        }
      }
    }
  }
 jumpElems:;
  
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

  free( filepntr);
  *nptr = node; *eptr = elem;
  return(1);
}

