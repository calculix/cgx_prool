/* ------------------------------------------------------------------------  */
/* write2dyna schreibt das frd-file aus einer struktur 22.04.1999 Wittig     */
/* ------------------------------------------------------------------------  */

#include <extUtil.h>



int write2dyna(char *datout, Summen *anz, Nodes *node, Elements *elem, Datasets *lcase )
{
  FILE *handle1;
  int  i;

  sprintf(datout,"%s.dyn",datout);
  /* Open the files and check to see that it was opened correctly */
  handle1 = fopen (datout, "w+b");
  if (handle1==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf ("\n%s opened\n",datout);


  printf ("\n write ls-dyna data  \n");

  if (anz->n>0)
  {
    fprintf (handle1, "*NODE\n" );
    for (i=0; i<anz->n; i++)
    {
      fprintf (handle1, "%8d%16.3f%16.3f%16.3f\n", node[i].nr, node[node[i].nr].nx, node[node[i].nr].ny, node[node[i].nr].nz );
    }
  }

  if (anz->e>0)
  {
    fprintf (handle1, "*ELEMENT_SOLID\n");
    for (i=0; i<anz->e; i++)
    {
      if (elem[i].type == 1)
      {
        fprintf (handle1, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
        elem[i].nr, elem[i].group,elem[i].nod[0], elem[i].nod[1], elem[i].nod[2], elem[i].nod[3],
        elem[i].nod[4], elem[i].nod[5], elem[i].nod[6], elem[i].nod[7] );
      }
      else if (elem[i].type == 3)
      {
        fprintf (handle1, "%8d%8d%8d%8d%8d%8d\n",
        elem[i].nr, elem[i].group,elem[i].nod[0], elem[i].nod[1], elem[i].nod[2], elem[i].nod[3] );
      }
      else if (elem[i].type == 4)
      {
        fprintf (handle1, "%8d%8d\n",
        elem[i].nr, elem[i].group );
        fprintf (handle1, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  elem[i].nod[0], elem[i].nod[1], elem[i].nod[2], elem[i].nod[3],
          elem[i].nod[4], elem[i].nod[5], elem[i].nod[6], elem[i].nod[7],
      	  elem[i].nod[8], elem[i].nod[9] );
        fprintf (handle1, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
          elem[i].nod[10], elem[i].nod[11], elem[i].nod[16], elem[i].nod[17],
          elem[i].nod[18], elem[i].nod[19], elem[i].nod[12], elem[i].nod[13],
          elem[i].nod[14], elem[i].nod[15] );
      }
      else if (elem[i].type == 6)
      {
        fprintf (handle1, "%8d%8d\n",
        elem[i].nr, elem[i].group );
        fprintf (handle1, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  elem[i].nod[0], elem[i].nod[1], elem[i].nod[2], elem[i].nod[3],
          elem[i].nod[4], elem[i].nod[5], elem[i].nod[6], elem[i].nod[8],
      	  elem[i].nod[9], elem[i].nod[7] );
      }
      else
      {
        printf (" elem not a known type (%d)\n",  elem[i].type);
      }
    }
  }

  fclose(handle1);
  return (1);
}

