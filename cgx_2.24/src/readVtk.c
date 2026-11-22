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
void generateDataset( Summen *anz, Datasets **ptr_lcase, char *name, int ncomps, double value, char *dataset_text, int analysis_type, int step_number, char *analysis_name);

/* --- ADDED BLOCK FOR WINDOWS COMPATIBILITY --- */
#ifdef WIN32
#include <ctype.h>
#include <string.h>

char *strcasestr(const char *haystack, const char *needle) {
    if (!*needle) return (char *)haystack;
    for (; *haystack; haystack++) {
        if (toupper((unsigned char)*haystack) == toupper((unsigned char)*needle)) {
            const char *h, *n;
            for (h = haystack, n = needle; *h && *n; h++, n++) {
                if (toupper((unsigned char)*h) != toupper((unsigned char)*n)) break;
            }
            if (!*n) return (char *)haystack;
        }
    }
    return NULL;
}
#endif
/* ------------------------------------------------- */

#define TEST     0

#define INI_FIELD_SIZE 100000

#define UINT8  char
#define UINT16  short int
#define UINT32  int
#define REAL32  float

int readVtk( char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr )
{
  FILE *handle;
  fpos_t *filepntr=NULL;
  
  /* REMOVED UNUSED VARIABLES 'l' AND 'nn' TO SILENCE WARNINGS */
  int i=0,j,k,e; 
  char *ptr=NULL;

  char n,rec_str[MAX_LINE_LENGTH];
  char buffer[MAX_LINE_LENGTH];
  char format[MAX_LINE_LENGTH];
  char key[MAX_LINE_LENGTH];
  int  node_field_size, elem_field_size;
  int  e_nmax=1, e_nmin=1, sum_n=0, sum_typ=0, sum_e=0;
  int  length;
  int  nodeMap[100];

  double version=0;
  int *eltype=NULL;

  Nodes     *node=NULL;
  Elements  *elem=NULL;
  Datasets  *lcase=NULL;



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

  anz->n=anz->e=-1;anz->l=0;
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
    printf (" it seems not to be ascii format, key:%s\n",rec_str);
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
    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
#if TEST
       printf (" pre read line:%s\n",rec_str);
#endif  
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
          switch (k) {
            case 12: eltype[i]=1;   break; // Hex
            case 13: eltype[i]=2;   break; // Wedge
            case 10: eltype[i]=3;   break; // Tetra
            case 25: eltype[i]=4;   break; // Quad Hex
            case 26: eltype[i]=5;   break; // VTK_QUADRATIC_WEDGE (15 nodes)
            case 24: eltype[i]=6;   break; // Quad Tetra
	    case 5:  eltype[i]=7;   break; // Triangle
            case 22: eltype[i]=8;   break; // Quad Triangle
	    case 9:  eltype[i]=9;   break; // Quad
            case 23: eltype[i]=10;  break; // Quad Quad
	    case 3:  eltype[i]=11;  break; // Line
            case 21: eltype[i]=12;  break; // Quad Line
	  default: eltype[i] = 0;  break; }
	}
        break;
      }
    }
    if( fsetpos( handle, (fpos_t *)filepntr)!=0) { printf("error in fsetpos"); return(-1); }

    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
#if TEST
       printf (" new line:%s\n",rec_str);
#endif  
  
      if((strcasestr(rec_str, "POINTS")) !=NULL)
      {
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
          k=0;
	  do{
	    buffer[k++]=n=getc(handle);
	    //printf("%c",n);
	  }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	  node[node[anz->n].nr].nx=atof(buffer);
          do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
          k=0;
	  do{
	    buffer[k++]=n=getc(handle);
	    //printf("%c",n);
	  }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	  node[node[anz->n].nr].ny=atof(buffer);
          do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
          k=0;
	  do{
	    buffer[k++]=n=getc(handle);
	    //printf("%c",n);
	  }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
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
      
      if((strcasestr(rec_str, "CELLS")) !=NULL)
      {
        sscanf(rec_str, "%*s %d %d", &sum_e, &sum_n);
	printf("reading %d CELLS with data points:%d, elements of known type:%d\n",sum_e, sum_n, sum_typ);

	if(sum_typ!=sum_e)
	{
          if(sum_typ==sum_e-1)
	  {
            length = frecord( handle, rec_str);
            for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
            if(compare(buffer, "OFFSETS", 7)!=7)
            {
  	      printf(" ERROR %d cells found but only %d types specified\n",k,sum_e);
	      goto jumpElems;
	    }
            for(i=0; i<sum_e; i++)
            {
              do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
              k=0;
	      do{
		buffer[k++]=n=getc(handle);
		//printf("%c",n);
	      }while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
	    }
            do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
	    //printf("\n");
            length = frecord( handle, rec_str);
            for(i=0;i<=strlen(rec_str); i++) buffer[i]=toupper(rec_str[i]);
            if(compare(buffer, "CONNECTIVITY", 7)!=7)
            {
  	      printf(" ERROR CONNECTIVITY not found\n");
	      goto jumpElems;
	    }
	    sum_e-=1;
	  }
	}
	  
        for(e=0; e<sum_e; e++)
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
          printf (" e=%d typ=%d group=%d mat=%d \n", elem[anz->e].nr,
          elem[anz->e].type, elem[anz->e].group, elem[anz->e].mat );
#endif
          //if(version>3.)
	  {

	    if(elem[anz->e].type==1)
	    {
	      for(j=0; j<9; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[j-1]=atoi(buffer)+1;
	      }
	    }	
	    else if(elem[anz->e].type==4)
	    {
	      j=0;
	      for(k=0; k<12; k++) nodeMap[j++]=k;
	      for(k=16; k<20; k++) nodeMap[j++]=k;
	      for(k=12; k<16; k++) nodeMap[j++]=k;
	      for(j=0; j<21; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[nodeMap[j-1]]=atoi(buffer)+1;
		//printf("%d node:%d %d\n",j,nodeMap[j-1],elem[anz->e].nod[nodeMap[j-1]]);
	      }
	    }	
	    if(elem[anz->e].type==6)
	    {
	      for(j=0; j<11; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[j-1]=atoi(buffer)+1;
	      }
	    }	
	    else if(elem[anz->e].type==7)
	    {
	      for(j=0; j<4; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[j-1]=atoi(buffer)+1;
	      }
	    }	
	    else if(elem[anz->e].type==8)
	    {
	      for(j=0; j<7; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[j-1]=atoi(buffer)+1;
	      }
	    }	
	    else if(elem[anz->e].type==10)
	    {
	      for(j=0; j<9; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[j-1]=atoi(buffer)+1;
	      }
	    }	
	    else if((elem[anz->e].type==3)||(elem[anz->e].type==9))
	    {
	      for(j=0; j<5; j++)
	      {
                do{ n=getc(handle); }while((n==' ')||(n=='\n')||(n=='\r')||(n=='\0')||(n==(char)EOF)); ungetc(n, handle);
                k=0;
		do{
		  buffer[k++]=n=getc(handle);
		}while((n!=' ')&&(n!='\n')&&(n!='\r')&&(n!='\0')&&(n!=(char)EOF));
		buffer[k++]=0;
		
		// jump the first number
		if(!j) continue;
		elem[anz->e].nod[j-1]=atoi(buffer)+1;
	      }
	    }	
	  }
#if TEST
	  printf("node: %d %d\n",j,elem[anz->e].nod[j]);
#endif
        }
      }
    }
  }
 jumpElems:;

  
  anz->n++;
  anz->e++;
  fclose(handle);
  free( filepntr);

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
    printf ("\nWARNING: element requests a nodename higher than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }
  if ( e_nmin < 1 )
  {
    printf ("\nWARNING: element requests a nodename lower than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }


  /* read results */
  
  /* Open the files and check to see that it was opened correctly */
  handle = fopen (datin, "r");
  if ( handle== NULL )  { printf ("ERROR: The input file \"%s\" could not be opened.\n\n", datin); return(-1); }
  else  printf (" file:%s opened\n", datin);

  /* default: */
  sum_n=anz->n;
  while(1)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;
#if TEST
    printf (" new line:%s\n",rec_str);
#endif
    
    if((ptr=strcasestr(rec_str, "POINT_DATA")) !=NULL)
    {
      sum_n=atoi((ptr+11));
      printf("data points:%d nodes:%d\n",sum_n, anz->n);
      if(sum_n>anz->n) { printf("ERROR: More point data than nodes\n"); return(-1); }
    }

    /* scalars */
    if((ptr=strcasestr(rec_str, "SCALARS")) !=NULL)
    {
      sscanf(rec_str,"%s %s %s", buffer, key, format);
      printf("found key:%s format:%s\n",key, format);

      /* create a new dataset */
      // Summen Datasets name ncomps value dataset_text analysis_type step_number analysis_name
      generateDataset(anz,&lcase,"dummy",1,0.,"",1,1,"");
      anz->l--;

      sprintf( lcase[anz->l].name, key);
      strcpy( lcase[anz->l].compName[0], key);

      /* further commands might follow */
      do {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        //printf (" new line:%s\n",rec_str);
        sscanf(rec_str,"%s", buffer);
      }while(checkIfNumber(buffer)!=1);

      for(k=0; k<sum_n; k++)
      {
        sscanf(rec_str,"%f", &lcase[anz->l].dat[0][node[k].nr]);
	//printf("%f\n", lcase[anz->l].dat[0][node[k].nr]);
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          lcase[anz->l].ictype[i] = 1;
          lcase[anz->l].icind1[i] = 1+i;
          //printf("comp%d n:%d val:%f min:%f\n",i,node[k].nr,lcase[anz->l].dat[i][node[k].nr],lcase[anz->l].min[i]);
          if (lcase[anz->l].dat[i][node[k].nr] >  lcase[anz->l].max[i])
          {  lcase[anz->l].max[i]=lcase[anz->l].dat[i][node[k].nr]; lcase[anz->l].nmax[i]=node[k].nr; }
          if (lcase[anz->l].dat[i][node[k].nr] <  lcase[anz->l].min[i])
          {  lcase[anz->l].min[i]=lcase[anz->l].dat[i][node[k].nr]; lcase[anz->l].nmin[i]=node[k].nr; }
        }
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        //printf (" new line:%s\n",rec_str);
      }
      anz->l++;
    }

    /* vectors */
    if((ptr=strcasestr(rec_str, "VECTORS")) !=NULL)
    {
      sscanf(rec_str,"%s %s %s", buffer, key, format);
      printf("found key:%s format:%s\n",key, format);

      /* create a new dataset */
      // Summen Datasets name ncomps value dataset_text analysis_type step_number analysis_name
      generateDataset(anz,&lcase,"dummy",3,0.,"",1,1,"");
      anz->l--;

      sprintf( lcase[anz->l].name, key);
      strcpy( lcase[anz->l].compName[0], "DX");
      strcpy( lcase[anz->l].compName[1], "DY");
      strcpy( lcase[anz->l].compName[2], "DZ");

      /* further commands might follow */
      do {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        //printf (" new line:%s\n",rec_str);
        sscanf(rec_str,"%s", buffer);
      }while(checkIfNumber(buffer)!=1);
	
      for(k=0; k<sum_n; k++)
      {
        sscanf(rec_str,"%f %f %f", &lcase[anz->l].dat[0][node[k].nr], &lcase[anz->l].dat[1][node[k].nr], &lcase[anz->l].dat[2][node[k].nr]);
	  //printf("%f %f %f\n", lcase[anz->l].dat[0][node[k].nr], lcase[anz->l].dat[1][node[k].nr], lcase[anz->l].dat[2][node[k].nr]);
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          lcase[anz->l].ictype[i] = 2;
          lcase[anz->l].icind1[i] = 1+i;
          //printf("comp%d n:%d val:%f min:%f\n",i,node[k].nr,lcase[anz->l].dat[i][node[k].nr],lcase[anz->l].min[i]);
          if (lcase[anz->l].dat[i][node[k].nr] >  lcase[anz->l].max[i])
          {  lcase[anz->l].max[i]=lcase[anz->l].dat[i][node[k].nr]; lcase[anz->l].nmax[i]=node[k].nr; }
          if (lcase[anz->l].dat[i][node[k].nr] <  lcase[anz->l].min[i])
          {  lcase[anz->l].min[i]=lcase[anz->l].dat[i][node[k].nr]; lcase[anz->l].nmin[i]=node[k].nr; }
        }
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        //printf (" new line:%s\n",rec_str);
      }
      anz->l++;
    }

    /* tensors */
    if((ptr=strcasestr(rec_str, "TENSORS")) !=NULL)
    {
      sscanf(rec_str,"%s %s %s", buffer, key, format);
      printf("found key:%s format:%s\n",key, format);

      /* create a new dataset */
	// Summen Datasets name ncomps value dataset_text analysis_type step_number analysis_name
      generateDataset(anz,&lcase,"dummy",6,0.,"",1,1,"");
      anz->l--;

      sprintf( lcase[anz->l].name, key);
      strcpy( lcase[anz->l].compName[0], "SXX");
      strcpy( lcase[anz->l].compName[1], "SYY");
      strcpy( lcase[anz->l].compName[2], "SZZ");
      strcpy( lcase[anz->l].compName[3], "SXY");
      strcpy( lcase[anz->l].compName[4], "SYZ");
      strcpy( lcase[anz->l].compName[5], "SZX");

      /* further commands might follow */
      do {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        //printf (" new line:%s\n",rec_str);
        sscanf(rec_str,"%s", buffer);
      }while(checkIfNumber(buffer)!=1);
	
      for(k=0; k<sum_n; k++)
      {
        sscanf(rec_str,"%f %f %f %f %f %f", &lcase[anz->l].dat[0][node[k].nr], &lcase[anz->l].dat[1][node[k].nr], &lcase[anz->l].dat[2][node[k].nr], &lcase[anz->l].dat[3][node[k].nr], &lcase[anz->l].dat[4][node[k].nr], &lcase[anz->l].dat[5][node[k].nr]);
	  //printf("%f %f %f\n", lcase[anz->l].dat[0][node[k].nr], lcase[anz->l].dat[1][node[k].nr], lcase[anz->l].dat[2][node[k].nr]);
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          lcase[anz->l].ictype[i] = 4;
          //printf("comp%d n:%d val:%f min:%f\n",i,node[k].nr,lcase[anz->l].dat[i][node[k].nr],lcase[anz->l].min[i]);
          if (lcase[anz->l].dat[i][node[k].nr] >  lcase[anz->l].max[i])
          {  lcase[anz->l].max[i]=lcase[anz->l].dat[i][node[k].nr]; lcase[anz->l].nmax[i]=node[k].nr; }
          if (lcase[anz->l].dat[i][node[k].nr] <  lcase[anz->l].min[i])
          {  lcase[anz->l].min[i]=lcase[anz->l].dat[i][node[k].nr]; lcase[anz->l].nmin[i]=node[k].nr; }
        }
        for(i=0; i<3; i++) { lcase[anz->l].icind1[i] = 1+i; lcase[anz->l].icind1[i+3] = 1+i; }
        for(i=0; i<3; i++) { lcase[anz->l].icind2[i] = 1+i; }
	  lcase[anz->l].icind2[3] = 2; lcase[anz->l].icind2[4] = 3; lcase[anz->l].icind2[5] = 1;
	  
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        //printf (" new line:%s\n",rec_str);
      }
      anz->l++;
    }
  }
  printf("end\n");

 

  *nptr = node; *eptr = elem; *lptr = lcase;
  return(1);
}

