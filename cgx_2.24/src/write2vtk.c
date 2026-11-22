/* --------------------------------------------------------------------  */
/* CALCULIX                                     */
/* - GRAPHICAL INTERFACE -                             */
/* */
/* MODIFIED VERSION:                                                 */
/* - RAW OUTPUT (No scaling mm->m, No Temp conversion)               */
/* - FIXED MESH TOPOLOGY (Node ID Mapping included)                  */
/* - FIXED PE15 ELEMENT TYPE (15-node Wedge)                         */
/* --------------------------------------------------------------------  */

#include <extUtil.h>
#include <stdlib.h> 

void recompileEntitiesInMenu(int lc);
void calcDatasets( int lc, Summen *anz, Nodes *node, Datasets *lcase );

int write2vtk( char *datout, Summen *anz, Nodes *node, Elements *elem, Datasets *lcase, char **dat )
{
  FILE *handle1=NULL, *handle2=NULL;
  int i, j, k, lc, lendatout, step_number_buf=-1, memFlag=0;
  int seff=10;       //PS1
  int formatFlag=0;  //v7.0: 0, v7.1: 71
  
  // Tablica mapuj?ca ID w?z?a (z CalculiX) na indeks (dla VTK)
  int *nodeMap = NULL;

  /* version 7.1 */
  if((dat!=0)&&( (compare( dat[0], "v7.1", 4)>= 4)||(compare( dat[1], "v7.1", 4)>= 4) )) formatFlag=71;

  // Informacja dla u?ytkownika - teraz wyniki s? SUROWE
  printf (" Exporting RAW results to VTK (No unit scaling applied).\n");

  lendatout=strlen(datout);

  /* open the mesh file */
  if(1)
  {
    strcpy (&datout[lendatout], ".vtk");
    handle1 = fopen (datout, "w+b");
    if (handle1==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
       datout); return(-1);}
    else  printf (" file %s opened\n",datout);
  
    // Alokacja mapy w?z?ów
    nodeMap = (int*)malloc((anz->nmax + 1) * sizeof(int));
    if (nodeMap == NULL) {
        printf("ERROR: Not enough memory for node map allocation.\n");
        fclose(handle1);
        return -1;
    }
    for(i=0; i <= anz->nmax; i++) nodeMap[i] = -1;

    fprintf (handle1, "# vtk DataFile Version 3.1\n");
    //fprintf (handle1, "TimeStep %lf\n", lcase[lc].value ); //lc undefined here, replacement:
    fprintf (handle1, "vtk output\n");
    fprintf (handle1, "ASCII\n");
    fprintf (handle1, "DATASET UNSTRUCTURED_GRID\n");
    
    // U?ywamy anz->n (ilo?? w?z?ów)
    fprintf (handle1, "POINTS %d double\n", anz->n);

    for (i=0; i<anz->n; i++)
    {
      // USUNI?TO SKALOWANIE * MM_TO_M
      fprintf (handle1, "%.12e %.12e %.12e\n", 
               node[node[i].nr].nx,
               node[node[i].nr].ny, 
               node[node[i].nr].nz);
               
      // Wype?nienie mapy
      if (node[i].nr <= anz->nmax) {
          nodeMap[node[i].nr] = i;
      }
    }

    int total_elements = anz->e;
    int total_size = 0;
    
    // --- KROK 1: Obliczenie sumarycznego rozmiaru sekcji CELLS ---
    for (i = 0; i < total_elements; i++) {
        switch (elem[i].type) {
            case 1:  total_size += (8 + 1); break; // he8
            case 2:  total_size += (6 + 1); break; // pe6
            case 3:  total_size += (4 + 1); break; // te4
            case 4:  total_size += (20 + 1); break; // he20
            case 5:  total_size += (15 + 1); break; // pe15 (15 w?z?ów)
            case 6:  total_size += (10 + 1); break; // te10
            case 7:  total_size += (3 + 1); break; // tr3
            case 8:  total_size += (6 + 1); break; // tr6
            case 9:  total_size += (4 + 1); break; // qu4
            case 10: total_size += (8 + 1); break; // qu8
            case 11: total_size += (2 + 1); break; // be2
            case 12: total_size += (3 + 1); break; // be3
            default: break;
        }
    }
    
    // --- KROK 2: Wypisanie sekcji CELLS ---
    fprintf(handle1, "\nCELLS %d %d\n", total_elements, total_size);
    
    for (i = 0; i < total_elements; i++) {
        if (elem[i].type == 1) { // Hexahedron (8)
            fprintf(handle1, "8 %d %d %d %d %d %d %d %d\n",
            nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]], nodeMap[elem[i].nod[2]], nodeMap[elem[i].nod[3]],
            nodeMap[elem[i].nod[4]], nodeMap[elem[i].nod[5]], nodeMap[elem[i].nod[6]], nodeMap[elem[i].nod[7]]);
        }
        else if (elem[i].type == 2) { // Wedge (6)
            fprintf(handle1, "6 %d %d %d %d %d %d\n", 
            nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]], nodeMap[elem[i].nod[2]], 
            nodeMap[elem[i].nod[3]], nodeMap[elem[i].nod[4]], nodeMap[elem[i].nod[5]]);
        }
        else if (elem[i].type == 5) { // Penta15 / Quadratic Wedge (15)
            fprintf(handle1, "15");
            for(j=0; j<15; j++) fprintf(handle1, " %d", nodeMap[elem[i].nod[j]]);
            fprintf(handle1, "\n");
        }
        else if (elem[i].type == 3 || elem[i].type == 9) { // Tetra / Quad (4 nodes)
            fprintf(handle1, "4 %d %d %d %d\n",
            nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]], nodeMap[elem[i].nod[2]], nodeMap[elem[i].nod[3]]);
        }
        else if (elem[i].type == 4) { // Quadratic Hex (20 nodes)
            fprintf(handle1, "20");
            for(j=0; j<12; j++) fprintf(handle1, " %d", nodeMap[elem[i].nod[j]]);
            fprintf(handle1, " %d %d %d %d %d %d %d %d\n",
            nodeMap[elem[i].nod[16]], nodeMap[elem[i].nod[17]], nodeMap[elem[i].nod[18]], nodeMap[elem[i].nod[19]],
            nodeMap[elem[i].nod[12]], nodeMap[elem[i].nod[13]], nodeMap[elem[i].nod[14]], nodeMap[elem[i].nod[15]]);
        }
        else if (elem[i].type == 6) { // Quadratic Tetra (10 nodes)
            fprintf(handle1, "10");
            for(j=0; j<10; j++) fprintf(handle1, " %d", nodeMap[elem[i].nod[j]]);
            fprintf(handle1, "\n");
        }
        else if (elem[i].type == 7 || elem[i].type == 12) { // Triangle / Quad Line (3 nodes)
            fprintf(handle1, "3 %d %d %d\n",
            nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]], nodeMap[elem[i].nod[2]]);
        }
        else if (elem[i].type == 8) { // Quadratic Triangle (6 nodes)
            fprintf(handle1, "6 %d %d %d %d %d %d\n",
            nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]], nodeMap[elem[i].nod[2]], 
            nodeMap[elem[i].nod[3]], nodeMap[elem[i].nod[4]], nodeMap[elem[i].nod[5]]);
        }
        else if (elem[i].type == 10) { // Quadratic Quad (8 nodes)
            fprintf(handle1, "8 %d %d %d %d %d %d %d %d\n",
            nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]], nodeMap[elem[i].nod[2]], nodeMap[elem[i].nod[3]],
            nodeMap[elem[i].nod[4]], nodeMap[elem[i].nod[5]], nodeMap[elem[i].nod[6]], nodeMap[elem[i].nod[7]]);
        }
        else if (elem[i].type == 11) { // Line (2 nodes)
            fprintf(handle1, "2 %d %d\n", nodeMap[elem[i].nod[0]], nodeMap[elem[i].nod[1]]);
        }
        else {
            printf(" WARNING: elem(%d) unknown type\n", elem[i].nr);
        }
    }
    
    // --- KROK 3: Wypisanie sekcji CELL_TYPES ---
    fprintf(handle1, "\nCELL_TYPES %d\n", total_elements);
    
    for (i = 0; i < total_elements; i++) {
        int vtk_id = 0;
        switch (elem[i].type) {
            case 1:  vtk_id = 12; break; // Hex
            case 2:  vtk_id = 13; break; // Wedge
            case 3:  vtk_id = 10; break; // Tetra
            case 4:  vtk_id = 25; break; // Quad Hex
            case 5:  vtk_id = 26; break; // VTK_QUADRATIC_WEDGE (15 nodes)
            case 6:  vtk_id = 24; break; // Quad Tetra
            case 7:  vtk_id = 5;  break; // Triangle
            case 8:  vtk_id = 22; break; // Quad Triangle
            case 9:  vtk_id = 9;  break; // Quad
            case 10: vtk_id = 23; break; // Quad Quad
            case 11: vtk_id = 3;  break; // Line
            case 12: vtk_id = 21; break; // Quad Line
            default: vtk_id = 0;  break;
        }
        fprintf(handle1, "%d\n", vtk_id);
    }
    
    free(nodeMap);
    fclose(handle1);
  }

  /* open node-value file */
  if((anz->l>0)&&(dat!=0)&&(compare( dat[0], "ds", 2)>= 2))
  {
    strcpy (&datout[lendatout], ".vtk");
    handle2 = fopen (datout, "a");
    if (handle2==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n", datout); return(-1);}
    
    fprintf (handle2, "\nPOINT_DATA %d", anz->n);

    j=0;
    for (lc=0; lc<anz->l; lc++)
    {
      if((strlen(dat[0])>2) && (atoi(&dat[0][2])!=lc+1)) goto next;

      if(compare( lcase[lc].name, "DISP", 4) == 4)
      {
        if (!lcase[lc].loaded)
        {
          if( readfrdblock( lc, anz, node, lcase )==-1) 
          {
            printf("ERROR in nodalDataset: Could not read data for Dataset:%d\n", lc+1); 
            return(-1);
          }
          calcDatasets( lc, anz, node, lcase );
          recompileEntitiesInMenu(lc);
          memFlag=1;
        }
        else memFlag=0;
	
        if(step_number_buf!=lcase[lc].step_number)
	{
          j++;
          step_number_buf=lcase[lc].step_number;
	}

        fprintf (handle2, "\nVECTORS DISP double\n");
        for (i=0; i<anz->n; i++)
        {
          fprintf (handle2, "%12.5f %12.5f %12.5f\n",
	  lcase[lc].dat[0][node[i].nr],lcase[lc].dat[1][node[i].nr],lcase[lc].dat[2][node[i].nr]);
        }

	if(memFlag)
        {
          for(k=0; k<lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
          free(lcase[lc].dat);
          lcase[lc].dat=NULL;
          lcase[lc].loaded=0;
        }
      }
      
      else if((compare( lcase[lc].name, "NDTEMP", 6) == 6)||( compare( lcase[lc].name, "TEMP", 4) == 4)||( compare( lcase[lc].name, "TT3D", 4) == 4))
      {
        if (!lcase[lc].loaded)
        {
          if( readfrdblock( lc, anz, node, lcase )==-1) 
          {
            printf("ERROR in nodalDataset: Could not read data for Dataset:%d\n", lc+1); 
            return(-1);
          }
          calcDatasets( lc, anz, node, lcase );
          recompileEntitiesInMenu(lc);
          memFlag=1; 
        }
        else memFlag=0;
	
        if(step_number_buf!=lcase[lc].step_number)
	{
          j++;
          step_number_buf=lcase[lc].step_number;
	}

        fprintf (handle2, "\nSCALARS TEMP double\n");
        fprintf (handle2, "LOOKUP_TABLE default\n");
        for (i=0; i<anz->n; i++)
        {
          // SUROWE DANE (Bez odejmowania 273.15)
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[0][node[i].nr]);
        }

	if(memFlag)
        {
          for(k=0; k<lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
          free(lcase[lc].dat);
          lcase[lc].dat=NULL;
          lcase[lc].loaded=0;
        }
      }

      else if(compare( lcase[lc].name, "STRESS", 6) == 6)
      {
        if (!lcase[lc].loaded)
        {
          if( readfrdblock( lc, anz, node, lcase )==-1) 
          {
            printf("ERROR in nodalDataset: Could not read data for Dataset:%d\n", lc+1); 
            return(-1);
          }
          calcDatasets( lc, anz, node, lcase );
          recompileEntitiesInMenu(lc);
          memFlag=1;
        }
        else memFlag=0;
	
        if(step_number_buf!=lcase[lc].step_number)
	{
          j++;
          step_number_buf=lcase[lc].step_number;
	}

        fprintf (handle2, "\nTENSORS STRESS double\n");
        for (i=0; i<anz->n; i++)
        {
          fprintf (handle2, "%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",
	  lcase[lc].dat[0][node[i].nr],lcase[lc].dat[3][node[i].nr],lcase[lc].dat[4][node[i].nr],
	  lcase[lc].dat[3][node[i].nr],lcase[lc].dat[1][node[i].nr],lcase[lc].dat[5][node[i].nr],
	  lcase[lc].dat[4][node[i].nr],lcase[lc].dat[5][node[i].nr],lcase[lc].dat[2][node[i].nr]);
        }
	
	fprintf (handle2, "\nSCALARS vMises double\n");
	fprintf (handle2, "LOOKUP_TABLE default\n");
	for (i=0; i<anz->n; i++)
	{
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[6][node[i].nr]);
	}
	
	fprintf (handle2, "\nSCALARS PS1 double\n");
	fprintf (handle2, "LOOKUP_TABLE default\n");
	for (i=0; i<anz->n; i++)
	{
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[10][node[i].nr]);
	}
	
	fprintf (handle2, "\nSCALARS PS2 double\n");
	fprintf (handle2, "LOOKUP_TABLE default\n");
	for (i=0; i<anz->n; i++)
	{
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[14][node[i].nr]);
	}
	
	fprintf (handle2, "\nSCALARS PS3 double\n");
	fprintf (handle2, "LOOKUP_TABLE default\n");
	for (i=0; i<anz->n; i++)
	{
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[18][node[i].nr]);
	}
	
	fprintf (handle2, "\nSCALARS worstPS double\n");
	fprintf (handle2, "LOOKUP_TABLE default\n");
	for (i=0; i<anz->n; i++)
	{
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[22][node[i].nr]);
	}
	
	fprintf (handle2, "\nSCALARS sgMises double\n");
	fprintf (handle2, "LOOKUP_TABLE default\n");
	for (i=0; i<anz->n; i++)
	{
          fprintf (handle2, "%12.5f\n", lcase[lc].dat[25][node[i].nr]);
	}

	if(memFlag)
        {
          for(k=0; k<lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
          free(lcase[lc].dat);
          lcase[lc].dat=NULL;
          lcase[lc].loaded=0;
        }
      }
      
      else if(compare( lcase[lc].name, "TOSTRAIN", 8) == 8)
      {
        if (!lcase[lc].loaded)
        {
          if( readfrdblock( lc, anz, node, lcase )==-1) 
          {
            printf("ERROR in nodalDataset: Could not read data for Dataset:%d\n", lc+1); 
            return(-1);
          }
          calcDatasets( lc, anz, node, lcase );
          recompileEntitiesInMenu(lc);
          memFlag=1;
        }
        else memFlag=0;
	
        if(step_number_buf!=lcase[lc].step_number)
	{
          j++;
          step_number_buf=lcase[lc].step_number;
	}

        fprintf (handle2, "\nTENSORS TOSTRAIN double\n");
        for (i=0; i<anz->n; i++)
        {
          fprintf (handle2, "%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",
	  lcase[lc].dat[0][node[i].nr],lcase[lc].dat[3][node[i].nr],lcase[lc].dat[4][node[i].nr],
	  lcase[lc].dat[3][node[i].nr],lcase[lc].dat[1][node[i].nr],lcase[lc].dat[5][node[i].nr],
	  lcase[lc].dat[4][node[i].nr],lcase[lc].dat[5][node[i].nr],lcase[lc].dat[2][node[i].nr]);
        }

	if(memFlag)
        {
          for(k=0; k<lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
          free(lcase[lc].dat);
          lcase[lc].dat=NULL;
          lcase[lc].loaded=0;
        }
      }
      
      else if(compare( lcase[lc].name, "PE", 2) == 2)
      {
        if (!lcase[lc].loaded)
        {
          if( readfrdblock( lc, anz, node, lcase )==-1) 
          {
            printf("ERROR in nodalDataset: Could not read data for Dataset:%d\n", lc+1); 
            return(-1);
          }
          calcDatasets( lc, anz, node, lcase );
          recompileEntitiesInMenu(lc);
          memFlag=1; 
        }
        else memFlag=0;
	
        if(step_number_buf!=lcase[lc].step_number)
	{
          j++;
          step_number_buf=lcase[lc].step_number;
	}

        fprintf (handle2, "\nSCALARS PE double\n");
        fprintf (handle2, "LOOKUP_TABLE default\n");
        for (i=0; i<anz->n; i++)
        {
          fprintf (handle2, "%12.5f\n",
	  lcase[lc].dat[0][node[i].nr]);
        }

	if(memFlag)
        {
          for(k=0; k<lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
          free(lcase[lc].dat);
          lcase[lc].dat=NULL;
          lcase[lc].loaded=0;
        }
      }
      
      else if(compare( lcase[lc].name, "CONTACT", 7) == 7)
      {
        if (!lcase[lc].loaded)
        {
          if( readfrdblock( lc, anz, node, lcase )==-1) 
          {
            printf("ERROR in nodalDataset: Could not read data for Dataset:%d\n", lc+1); 
            return(-1);
          }
          calcDatasets( lc, anz, node, lcase );
          recompileEntitiesInMenu(lc);
          memFlag=1;
        }
        else memFlag=0;
	
        if(step_number_buf!=lcase[lc].step_number)
	{
          j++;
          step_number_buf=lcase[lc].step_number;
	}

        fprintf (handle2, "\nSCALARS CPRESS double\n");
        fprintf (handle2, "LOOKUP_TABLE default\n");
        for (i=0; i<anz->n; i++)
        {
          fprintf (handle2, "%12.5f\n",
	  lcase[lc].dat[3][node[i].nr]);
        }

	if(memFlag)
        {
          for(k=0; k<lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
          free(lcase[lc].dat);
          lcase[lc].dat=NULL;
          lcase[lc].loaded=0;
        }
      }
      next:;
    }
    fclose(handle2);
  }
  return(1);
}
