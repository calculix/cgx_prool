/*
This function uses the following empirical functions (si units)
    cpGas=1005.+0.0785*(temp-273.15); // sp. heat cap. J/Kg/K
    wl=0.0245*(1.+0.00225*(temp-273.15));  // thermal conductivity
    // sutherland:
    dynVisc=pow((temp/273.),(3./2.)) * (273.+111.)/(temp+111.) * dynVisc0;
*/

#include <extUtil.h>
extern Values   *value;
int getValuNr(char *name);

int write2isaac( char *name, Summen *anz, Nodes *node, Elements *elem, Datasets *lcase, NodeBlocks *nBlock, int bouNr, NodeBlockbou *blockbou )
{
  FILE *handle1, *handle2, *handle3;
  int  i,j,k,l,n,b=0, edges=0, m[3], jetFlag;
  char datout[MAX_LINE_LENGTH];
  int ic,jc,kc, **buf;
  int        pn=0;

  double p1=-1,p2=-1,pt=-1,tt=-1, v,vu,vv;
  char buffer[MAX_LINE_LENGTH];
  double  gasConstant=287.1, gamma, dynVisc, dynVisc0=1.786e-5, temp=0., twall=0., cpGas, wl, pr, re, rho, va, alpha, lref=1., velFactor=1.;

  /* change the name of files from <setname> to isaac */
  strcpy(name, "isaac");

  /* Open the files and check to see that it was opened correctly */
  sprintf(datout, "%s.grd", name);
  handle1 = fopen (datout, "w+b");
  if (handle1==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf (" file %s opened\n",datout);

  sprintf(datout, "%s.dat", name);
  handle2 = fopen (datout, "w+b");
  if (handle2==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf (" file %s opened\n",datout);

  sprintf(datout, "%s.fbd", name);
  handle3 = fopen (datout, "w+b");
  if (handle3==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf (" file %s opened\n",datout);

 
  /*** grid file ***/
  
  printf (" write isaac mesh \n");

  if(anz->b>1) fprintf (handle1, "%d\n", anz->b);
  for (b=0; b<anz->b; b++)
  {
    /* if 2D switch i,j (face normals must point against z direction, in duns in z direction) */
    if(nBlock[b].dim==2) fprintf (handle1, "%d %d\n", nBlock[b].j, nBlock[b].i);
    else                 fprintf (handle1, "%d %d %d\n", nBlock[b].i, nBlock[b].j, nBlock[b].k);
  }

  for (b=0; b<anz->b; b++)
  {
    /* if 2D switch i,j */
    if(nBlock[b].dim==2)
    {  
      if ( (buf = (int **)malloc( ( nBlock[b].j) * sizeof(int *))) == NULL ) printf(" ERROR: malloc failed\n\n");
      for (j=0; j<nBlock[b].j; j++)
        if ( (buf[j] = (int *)malloc(( nBlock[b].i) * sizeof(int))) == NULL ) printf(" ERROR: malloc failed\n\n");

      n=0;
      for (j=0; j<nBlock[b].j; j++)
      {
        for (i=0; i<nBlock[b].i; i++)
        {
          buf[j][i]=nBlock[b].nod[n];
          n++;
        }
      }

      for (l=0; l<nBlock[b].dim; l++)
      {
        n=0;
        for (i=0; i<nBlock[b].i; i++)
        {
          for (j=0; j<nBlock[b].j; j++)
          {
            if(l==0) fprintf (handle1, "%le\n", node[buf[j][i]].nx);
            if(l==1) fprintf (handle1, "%le\n", node[buf[j][i]].ny);
            if(l==2) fprintf (handle1, "%le\n", node[buf[j][i]].nz);
            n++;
          }
        }
      }
      for (j=0; j<nBlock[b].j; j++) free(buf[j]); free(buf);
    }
    else
    {
      for (l=0; l<nBlock[b].dim; l++)
      {
        n=0;
        for (k=0; k<nBlock[b].k; k++)
        {
          for (j=0; j<nBlock[b].j; j++)
          {
            for (i=0; i<nBlock[b].i; i++)
            {
              if(l==0) fprintf (handle1, "%le\n", node[nBlock[b].nod[n]].nx);
              if(l==1) fprintf (handle1, "%le\n", node[nBlock[b].nod[n]].ny);
              if(l==2) fprintf (handle1, "%le\n", node[nBlock[b].nod[n]].nz);
              n++;
            }
  	  }
        }
      }

      /* fbd file */
      fprintf (handle3, "seto b%d\n",b+1); 
      n=0;
      for (k=0; k<nBlock[b].k; k++)
      {
        for (j=0; j<nBlock[b].j; j++)
        {
          for (i=0; i<nBlock[b].i; i++)
          {
            fprintf (handle3, "pnt %d %le %le %le\n",pn, node[nBlock[b].nod[n]].nx, node[nBlock[b].nod[n]].ny, node[nBlock[b].nod[n]].nz);
            pn++; n++;
          }
        }
      }
      fprintf (handle3, "setc b%d\n",b+1); 
    }
  }
  fclose(handle1);
  fclose(handle3);

  // is a JET condition defined?
  jetFlag=0;
  for (b=0; b<anz->b; b++)
  {
    if(nBlock[b].dim==2) edges=4;
    if(nBlock[b].dim==3) edges=6;
    for(i=0; i<edges; i++) if(nBlock[b].map[i][0]<0 )
    {
      if(compare(nBlock[b].bctype[i],"JET", 3)==3) jetFlag=1;
    }
  }

  /*** input file ***/

  // write the gas-properties and initialization block only if the temperature is defined
  sprintf(buffer,"temperature");
  if(getValuNr(buffer)>-1)
  {
    temp=atof(value[getValuNr(buffer)].string);
    sprintf(buffer,"ref-length");
    if(getValuNr(buffer)>-1)
    {
      lref=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"WallTemperature");
    if(getValuNr(buffer)>-1)
    {
      twall=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"stagnation-temperature");
    if(getValuNr(buffer)>-1)
    {
      tt=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"pressure");
    if(getValuNr(buffer)>-1)
    {
      p2=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"stagnation-pressure");
    if(getValuNr(buffer)>-1)
    {
      pt=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"u-velocity");
    if(getValuNr(buffer)>-1)
    {
      vu=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"v-velocity");
    if(getValuNr(buffer)>-1)
    {
      vv=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"ini-velocity-factor");
    if(getValuNr(buffer)>-1)
    {
      velFactor=atof(value[getValuNr(buffer)].string);
    }

    cpGas=1005.+0.0785*(temp-273.15); // sp. heat cap. J/Kg/K
    gamma=1./(1-(gasConstant/cpGas));
    wl=0.0245*(1.+0.00225*(temp-273.15));
    // sutherland:
    dynVisc=pow((temp/273.),(3./2.)) * (273.+111.)/(temp+111.) * dynVisc0;
    pr=dynVisc*cpGas/wl;

    v=sqrt(vu*vu+vv*vv);
    if(vu!=0.) alpha=atan(vv/vu)*180/PI; else alpha=90.;
    rho=p2/(gasConstant*temp);
    re=lref*v*rho/dynVisc;

    //tt=v*v*.5/cpg + ts;
    if((pt!=-1)&&(tt!=-1)) p1=pt*pow( (temp/tt), (gamma/(gamma-1)) );
    else p1=p2;
    
    va=sqrt(gamma*gasConstant*temp);

    fprintf (handle2, "'SIUNIT'    1\n");
    fprintf (handle2, "'TINF'    %f\n", temp );
    fprintf (handle2, "'BACK PRESSURE'    %f\n", p2/p2 );
    if(twall) fprintf (handle2, "'TWALL'    %f\n", twall );
    fprintf (handle2, "'RE'    %f\n", re );
    fprintf (handle2, "'PR'    %f\n", pr );
    if(jetFlag)
    {
      fprintf (handle2, "'MACH'    %f\n", v/va*velFactor );
      fprintf (handle2, "'ALPHA'    %f\n", alpha*velFactor );
      fprintf (handle2, "'JET CONDITIONS'    1\n");
      fprintf (handle2, "'U'    %f\n", vu/va );
      fprintf (handle2, "'V'    %f\n", vv/va );
      fprintf (handle2, "'P'    %f\n", p1/p2 );
      fprintf (handle2, "'T'    %f\n", 1. );
      fprintf (handle2, "'END JET'    1\n");
    }
    else
    {
      fprintf (handle2, "'MACH'    %f\n", v/va );
      fprintf (handle2, "'ALPHA'    %f\n", alpha );
    }
    fprintf (handle2, "'GAMMA'    %f\n", gamma );
    fprintf (handle2, "'REFERENCE AREA'    %f\n", lref );
  }
  
  if(nBlock[0].dim==2)
  {
    /* if 2D switch i,j */
    ic=1;
    jc=0;
    kc=2;
    fprintf (handle2, "'TWOD'       1\n" );
  }
  else
  {
    ic=0;
    jc=1;
    kc=2;
  }
  fprintf (handle2, "'GRID FORMATTED'       1\n" );
  fprintf (handle2, "'GRID'       1\n" );
  fprintf (handle2, "  %s.grd\n", name );

  /* connectivity */
  fprintf (handle2, "'CUT'    %d\n", anz->c/2);
  for (b=0; b<anz->b; b++)
  {
    if(nBlock[b].dim==2) edges=4;
    if(nBlock[b].dim==3) edges=6;
    for(i=0; i<edges; i++) if(nBlock[b].map[i][0]>-1)
    {
      /* block-borders appear twice, write only the couple were the neighbour-nr is higher, may be equal if 'periodic' */
      if((nBlock[b].neighbor[i]>= nBlock[b].geo)&&(nBlock[b].strt2[i][jc]!=0))
      {
        for(j=0; j<3; j++) if( nBlock[b].map[i][j]>3) m[j]=nBlock[b].map[i][j]-3; else m[j]=nBlock[b].map[i][j];

        fprintf(handle2,"  'b%db%d'   %4d   %4d %4d %4d   %4d %4d %4d\n", nBlock[b].geo+1, nBlock[b].neighbor[i]+1, nBlock[b].geo+1, nBlock[b].strt1[i][ic], nBlock[b].strt1[i][jc], nBlock[b].strt1[i][kc], nBlock[b].end_1[i][ic], nBlock[b].end_1[i][jc], nBlock[b].end_1[i][kc]);
        fprintf(handle2,"  'b%db%d'   %4d   %4d %4d %4d   %4d %4d %4d", nBlock[b].geo+1, nBlock[b].neighbor[i]+1, nBlock[b].neighbor[i]+1, nBlock[b].strt2[i][ic], nBlock[b].strt2[i][jc], nBlock[b].strt2[i][kc], nBlock[b].end_2[i][ic], nBlock[b].end_2[i][jc], nBlock[b].end_2[i][kc]);
        fprintf(handle2,"   %d %d %d\n", m[0],m[1],m[2]);
      }
    }
  }
  
  /* boundaries */
  for (b=0; b<anz->b; b++)
  {
    n=0; for(i=0; i<edges; i++) if(nBlock[b].map[i][0]<0 ) n++;

    if(1)
    {
      fprintf (handle2, "'BLOCK'    %d\n", nBlock[b].geo+1);
      fprintf (handle2, "  'DIMENSIONS'    1\n" );
      if(nBlock[b].dim==2) fprintf (handle2, "    %d %d 2\n", nBlock[b].j, nBlock[b].i);
      else                 fprintf (handle2, "    %d %d %d\n", nBlock[b].i, nBlock[b].j, nBlock[b].k);
      fprintf (handle2, "  'FLUX'    1\n" );
      fprintf (handle2, "    'ROE'    1\n" );
      if(nBlock[b].dim==2) fprintf (handle2, "  'VISCOUS'    2\n    1 2  \n" );
      else                 fprintf (handle2, "  'VISCOUS'    3\n    1 2 3\n" );
      if(nBlock[b].dim==2) fprintf (handle2, "  'BC'    %d\n", n+2 );
      else                 fprintf (handle2, "  'BC'    %d\n", n );  
      if(nBlock[b].dim==2) edges=4;
      else                 edges=6;
      for(i=0; i<edges; i++) if(nBlock[b].map[i][0]<0 )
      {
        fprintf(handle2,"    'BC'    %4d %4d %4d   %4d %4d %4d  '%s'\n", nBlock[b].strt1[i][ic], nBlock[b].strt1[i][jc], nBlock[b].strt1[i][kc], nBlock[b].end_1[i][ic], nBlock[b].end_1[i][jc], nBlock[b].end_1[i][kc], nBlock[b].bctype[i]);
      }

      /* if dim==2 set the top and bottom face */
      if(nBlock[b].dim==2)
      {
        fprintf(handle2,"    'BC'    %4d %4d %4d   %4d %4d %4d  'EXTRAPOLATE'\n",1,1,1,nBlock[b].j, nBlock[b].i, 1);
        fprintf(handle2,"    'BC'    %4d %4d %4d   %4d %4d %4d  'EXTRAPOLATE'\n",1,1,2,nBlock[b].j, nBlock[b].i, 2);
      }
      fprintf (handle2, "'END BLOCK'    %d\n", nBlock[b].geo+1 );
    }
  }

  fprintf (handle2, "'END'       0\n" );
  if(temp!=0.)
  {
    fprintf (handle2, "'PINF'    %f\n", p2 );
    fprintf (handle2, "'RG'    %f\n", gasConstant );
    fprintf (handle2, "'CPG'    %f\n", cpGas );
    fprintf (handle2, "'TC'    %f\n", wl );
  }
  fclose(handle2);
  return (1);
}

