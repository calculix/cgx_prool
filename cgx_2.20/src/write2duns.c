
#include <cgx.h>
extern Values   *value;



int write2duns( char *datout, Summen *anz, Nodes *node, Elements *elem, Datasets *lcase, NodeBlocks *nBlock, int bouNr, NodeBlockbou *dunsbou )
{
  FILE *handle1, *handle2, *handle3;
  int  i,j,k,n,b=0,edges=0;
  double dl, p1=-1, pt1=-1,p2=-1, vu,vv,vel,flowangle;
  char buffer[MAX_LINE_LENGTH];
  double  gasConstant=287.1, gamma, dynVisc, dynVisc0=1.786e-5, temp,tt1, cpGas, wl, pr, massflow, velFactor=1.;

  /* Open the files and check to see that it was opened correctly */
  i=strlen(datout);
  strcpy (datout, "duns.grid");
  handle1 = fopen (datout, "w+b");
  if (handle1==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf (" file %s opened\n",datout);

  strcpy (datout, "duns.conn");
  handle2 = fopen (datout, "w+b");
  if (handle2==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf (" file %s opened\n",datout);

  strcpy (datout, "duns.bou");
  handle3 = fopen (datout, "w+b");
  if (handle3==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n",
     datout); return(-1);}
  else  printf (" file %s opened\n",datout);

  printf (" write duns mesh \n");

  /* grid file */
  fprintf (handle1, "%d\n", anz->b);
  for (b=0; b<anz->b; b++)
  {
    if(nBlock[b].dim==2) fprintf (handle1, "%d %d\n", nBlock[b].i, nBlock[b].j);
    else                 fprintf (handle1, "%d %d %d\n", nBlock[b].i, nBlock[b].j, nBlock[b].k);
    n=0;
    for (k=0; k<nBlock[b].k; k++)
    {
      for (j=0; j<nBlock[b].j; j++)
      {
        for (i=0; i<nBlock[b].i; i++)
        {
          if(nBlock[b].dim==2) fprintf (handle1, "%le %le\n", node[nBlock[b].nod[n]].nx, node[nBlock[b].nod[n]].ny);
          else                 fprintf (handle1, "%le %le %le\n", node[nBlock[b].nod[n]].nx, node[nBlock[b].nod[n]].ny, node[nBlock[b].nod[n]].nz);
          n++;
	}
      }
    }
  }
  fclose(handle1);

  /* connectivity file */
  fprintf (handle2, "%3d       blocks\n", anz->b);
  for (b=0; b<anz->b; b++)
  {
    fprintf (handle2,"B %3d",nBlock[b].geo+1-nBlock[0].geo );
    if(nBlock[b].dim==2) edges=4;
    if(nBlock[b].dim==3) edges=6;
    for(i=0; i<edges; i++)
    {
      if(nBlock[b].map[i][0]>-1) fprintf(handle2," b  0%3d %1d%1d%1d",  nBlock[b].neighbor[i]+1-nBlock[0].geo, nBlock[b].map[i][0], nBlock[b].map[i][1], nBlock[b].map[i][2]);
      else  fprintf(handle2," s %2d  0 000",  nBlock[b].neighbor[i]);
    }
    fprintf (handle2,"\n" );
  }
  fclose(handle2);

  // write the gas-properties and initialization block only if the temperature is defined
  sprintf(buffer,"temperature");
  if(getValuNr(buffer)>-1)
  {
    /* boundary file */
    
    temp=atof(value[getValuNr(buffer)].string);
    cpGas=1005.+0.0785*(temp-273.15); // sp. heat cap. J/Kg/K
    gamma=1./(1-(gasConstant/cpGas));
    wl=0.0245*(1.+0.00225*(temp-273.15));
    // sutherland:
    dynVisc=pow((temp/273.),(3./2.)) * (273.+111.)/(temp+111.) * dynVisc0;
    pr=dynVisc*cpGas/wl;
    fprintf(handle3,"  gas-properties\n");
    fprintf(handle3,"      compressible\n");
    fprintf(handle3,"      gas-constant %f\n",gasConstant);
    fprintf(handle3,"      gamma %f\n",gamma);
    fprintf(handle3,"      viscosity %e\n",dynVisc);
    fprintf(handle3,"      prandtl-number %f\n",pr);
    fprintf(handle3,"  turbulence\n");
    fprintf(handle3,"      prandtl-number %f\n",pr*1.33);

    /* initialize */

    fprintf(handle3,"   initialize\n");
    sprintf(buffer,"ini-velocity-factor");
    if(getValuNr(buffer)>-1)
    {
      velFactor=atof(value[getValuNr(buffer)].string);
    }
    sprintf(buffer,"pressure");
    if(getValuNr(buffer)>-1)
    {
      p2=atof(value[getValuNr(buffer)].string);
      sprintf(buffer,"stagnation-pressure");
      if(getValuNr(buffer)>-1)
      {
        pt1=atof(value[getValuNr(buffer)].string);
        fprintf(handle3,"      pressure %f\n", (p2+pt1)*.5);
      }
      else fprintf(handle3,"      pressure %f\n", p2);
    }
    fprintf(handle3,"      temperature %f\n",temp);
    sprintf(buffer,"u-velocity");
    if(getValuNr(buffer)>-1)
    {
      vu=atof(value[getValuNr(buffer)].string);
      fprintf(handle3,"      u-velocity %f\n",vu*velFactor);
    }
    sprintf(buffer,"v-velocity");
    if(getValuNr(buffer)>-1)
    {
      vv=atof(value[getValuNr(buffer)].string);
      fprintf(handle3,"#      v-velocity %f\n",vv*velFactor);
    }
    fprintf(handle3,"      v-velocity 0.\n");
    fprintf(handle3,"      turbulent-intensity 0.06\n");
    fprintf(handle3,"      turbulent-viscosity %e\n", dynVisc*10.);
    fprintf(handle3,"\n");
  }
  vel=sqrt(vu*vu+vv*vv);
  flowangle=asin(-vv/vel);
  
  fprintf(handle3,"  boundary-conditions \n");
  for (b=0; b<bouNr; b++)
  {
    fprintf(handle3,"#\n# set:%s\n", dunsbou[b].name);
    fprintf(handle3,"  surface ");
    for(i=0; i<dunsbou[b].surfs; i++) fprintf(handle3,"%d ", dunsbou[b].surf[i]);
    fprintf(handle3,"\n    %s\n", dunsbou[b].bctype);
    if(compare(dunsbou[b].bctype,"inviscid-wall",11)==11)
    {
      sprintf(buffer,"WallTemperature");
      if(getValuNr(buffer)>-1)
	{  fprintf(handle3,"        temperature %s\n",value[getValuNr(buffer)].string); }
      else fprintf(handle3,"        adiabatic\n");
    }
    if(compare(dunsbou[b].bctype,"viscous-wall",11)==11)
    {
      sprintf(buffer,"WallTemperature");
      if(getValuNr(buffer)>-1)
	{  fprintf(handle3,"        temperature %s\n",value[getValuNr(buffer)].string); }
      else fprintf(handle3,"        adiabatic\n");
    }
    if(compare(dunsbou[b].bctype,"subsonic-outflow",14)==14)
    {
      sprintf(buffer,"pressure");
      if(getValuNr(buffer)>-1) fprintf(handle3,"        pressure %s\n",value[getValuNr(buffer)].string);
    }
    if(compare(dunsbou[b].bctype,"subsonic-inflow",14)==14)
    {
      sprintf(buffer,"stagnation-temperature");
      if(getValuNr(buffer)>-1)
      {
        tt1=atof(value[getValuNr(buffer)].string);
	fprintf(handle3,"        stagnation-temperature %e\n",tt1);
        sprintf(buffer,"stagnation-pressure");
        if(getValuNr(buffer)>-1)
	{
	  pt1=atof(value[getValuNr(buffer)].string);
          fprintf(handle3,"        stagnation-pressure %e\n",pt1);
	}
        fprintf(handle3,"#        velocity %f\n",vel);

        p1=pt1*pow( (temp/tt1), (gamma/(gamma-1)) );
        // specific massflow per area
        massflow= vu * p1/gasConstant/temp;
        fprintf(handle3,"# mass-flow calculation based on p:%e\n",p1);
        fprintf(handle3,"#        mass-flow %e\n",massflow);
        fprintf(handle3,"#        temperature %e\n",temp);
      }
      else
      {
        fprintf(handle3,"        velocity %f\n",vel);
        sprintf(buffer,"temperature");
        if(getValuNr(buffer)>-1) fprintf(handle3,"        temperature %s\n",value[getValuNr(buffer)].string);
      }
      fprintf(handle3,"        flow-angle %f\n",flowangle);
      sprintf(buffer,"temperature");
      if(getValuNr(buffer)>-1)
      {
        fprintf(handle3,"         turbulent-intensity 0.06\n");
        fprintf(handle3,"         turbulent-viscosity %e\n", dynVisc*10.);
      }
    }
  }

  // define precondition only if the pressure is known
  sprintf(buffer,"pressure");
  if(getValuNr(buffer)>-1)
  {
    fprintf(handle3,"  precondition\n");
    fprintf(handle3,"      inviscid        on\n");
    fprintf(handle3,"      viscous         on\n");
    sprintf(buffer,"pressure");
    if(getValuNr(buffer)>-1) fprintf(handle3,"        pressure %s\n",value[getValuNr(buffer)].string);
    if(getValuNr(buffer)>-1) fprintf(handle3,"        velocity %f\n",vel);
  }
  
  fclose(handle3);

  return (1);
}

