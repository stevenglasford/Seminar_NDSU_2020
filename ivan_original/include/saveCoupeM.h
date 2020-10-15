//- coupe[d]=0 ==> on coupe dans la direction d en la valeur x[d]=val[d]
//- [DEC 2016: vtab parameter added]
void HJB::saveCoupeM (const char* file, const int* coupe, const double* val, double* vtab)
{
  if(my_rank==0)
    printf("Save coupe  : [saveCoupeM]   (%s)...\n",file);

  int i, rk, d, *testval;
  bool test;

  if(PRINTBUG && my_rank==0)
    printf("saveCoupeM    : Saving (%19s) on axes [",file);
  testval = new int[dim];
  for(d=0;d<dim;++d){
    if(coupe[d]){
      testval[d] = -1;
      if (PRINTBUG && my_rank==0)
        cout << d << " ";
    }
  }
  int itest=0;
  for(d=0;d<dim;++d){
    itest=max(coupe[d],itest);
  }
  if (itest==0) {
    printf(" all coupe[d]==0 ==> no savings !\n");
    return;
  }

  if(PRINTBUG && my_rank==0)
    printf("] with ");
  for(d=0;d<dim;++d){
    //printf("---mesh->GLOBAL_lowBounds[d]= %8.4f; \n",mesh->GLOBAL_lowBounds[d]);
    //printf("---mesh->GLOBAL_highBounds[d]=%8.4f; \n",mesh->GLOBAL_highBounds[d]);
    //printf("---mesh->lowBounds[d]= %8.4f; \n",mesh->lowBounds[d]);
    //printf("---mesh->highBounds[d]=%8.4f; \n",mesh->highBounds[d]);
    if(!coupe[d]){
      //printf("---val[%i]=%8.5f\n",d,val[d]);
      //printf("---testval[%i]=%5i\n",d,testval[d]);
      //printf("---mesh->inn_gridDim[d]=%5i; \n",mesh->inn_gridDim[d]);
      //printf("---mesh->GLOBAL_lowBounds[d]= %8.4f; \n",mesh->GLOBAL_lowBounds[d]);
      //printf("---mesh->GLOBAL_highBounds[d]=%8.4f; \n",mesh->GLOBAL_highBounds[d]);
      //printf("---mesh->lowBounds[d]= %8.4f; \n",mesh->lowBounds[d]);
      //printf("---mesh->highBounds[d]=%8.4f; \n",mesh->highBounds[d]);
      //(NO) if (val[d]<mesh->GLOBAL_lowBounds[d] || val[d]>mesh->GLOBAL_highBounds[d]){
      if (val[d]<mesh->lowBounds[d] || val[d]>mesh->highBounds[d]){
        printf(" 'val[d] out of boundary' => no savings !\n");
        return;
      }
      else{
        testval[d] = (int)floor((val[d]-mesh->lowBounds[d])*divdx[d]);
        if(PRINTBUG && my_rank==0)
          printf("val(%i)=%8.5f => %8.5f",d+1,val[d],mesh->discretization[d][testval[d]]);
      }
    }
  }
  //exit(1);

  //- xx contains an exact x, where we interpolate vtab[x].
  double xx[dim];
  for(d=0; d<dim; d++) xx[d]=val[d];

  double value;

  FILE * pFile;
  pFile = fopen (file,"w");
  if(pFile==NULL)
    cout << "09.Impossible to write in the given file!" << endl;
  else{
    for(i=0; i<mesh->inn_nbPoints; i++){
      rk=mesh->getOuterRank(i);
      test=true;
      mesh->setrank(rk,iCoord);

      //- TEST MAY 2013
      for(d=0; d<dim; d++){
        if(coupe[d]==1)  //- put in xx[d] the value of mesh coord[d]
          xx[d]=mesh->discretization[d][iCoord[d]];
      }

      for(d=0; d<dim; d++){
        //printf("d=%2i, coupe[d]=%i, testval[d]=%i, iCoord[d]=%i\n",d,coupe[d],testval[d],iCoord[d]);//- debug
        if(!coupe[d] && testval[d] != iCoord[d]){
          test=false;
          break;
        }
      }

      if(test){

        //- normal ascii spacing
        for(d=0; d<dim; d++){
          if(coupe[d]){
            fprintf(pFile,"%d\t", iCoord[d]);
            //fprintf(pFile,"%8.5f  ", xx[d]);
          }
        }
	//value=vtab[rk];                        //- (no interpolation)
	value= (*this.*interpolation)(xx,vtab);  //- (with interpolation)
        fprintf(pFile,"%8.5f",  value);
	// ---- BEGIN SUPPLEMENT: true coordinates shown at right
	int SHOW_COORDINATES=0;
	if (SHOW_COORDINATES){
          fprintf(pFile,"\t");
          for(d=0; d<dim; d++){
            if(coupe[d]){
              fprintf(pFile,"%8.4f  ", xx[d]);
            }
          }
	}
	// ---- END SUPPLEMENT
        fprintf(pFile,"\n");
      }
    }
  }
  fclose (pFile);
  delete[] testval;
  if (PRINTBUG && my_rank==0)
    printf("... end.saveCoupeM\n");
  return;
}

