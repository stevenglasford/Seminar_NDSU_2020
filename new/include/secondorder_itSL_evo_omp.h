//- NEW SL ORDER 2 function: BASIC ORDER 1 SCHEME (OK)
void HJB_SL::secondorder_itSL_evo_omp(double t, double dt)
{
  //- P=PARAMP of the data file.
  int i,j,c,e,d,p,test;

  #pragma omp parallel num_threads(OMP_NUM_THREADS) private(d,i,j,c,p,test,e) shared(t,dt) default(none)
  {


  #pragma omp for
  //for(i=0; i<mesh->inn_nbPoints; i++)
  for(j=0;j<ranksize;j++){
    i   =rank[j];
    vold[i]=v[i];
  }

  #pragma omp for
  //for(i=0; i<mesh->inn_nbPoints; i++){
  for(j=0;j<ranksize;j++){
    i   =rank[j];

    double a=0.0;
    double minA=0.0;
    double coordsomp[dim];
    double rCoordomp[dim];
    double dvectdoubleomp[dim];

    (mesh->*(mesh->setcoords))(i,rCoordomp);

    for(c=0;c<ncall;c++){
      a = 0.0;

      //- this: (p<=P) for use with data/data_SL2_optionput_approach2.h
      for(p=0;p<P;p++){
        for(e=-1;e<=1;e+=2){ //- two points formula
          funcY(rCoordomp, p, e, (*u)[c], t, dt, dvectdoubleomp);
          if(periodic_mesh)
            ComputePeriodic(dvectdoubleomp, coordsomp);
          else
            for(d=0;d<dim;d++)
              coordsomp[d] = dvectdoubleomp[d];
          test=0;
          for(d=0;d<dim;d++){
            if(coordsomp[d] <= lb[d] || coordsomp[d] >= hb[d]){
              test=1;
              break;
            }
          }
          if(test)
            a += (*this.*VbordCompute)(t,coordsomp,vold);
          else
            a += (*this.*interpolation)(coordsomp,vold);
        }
      }
      a  = a / (2*P);

      a += dt * (*distributed_cost)(rCoordomp,(*u)[c],t);

      a *= exp(-funcR(rCoordomp, (*u)[c], t)*dt);

      a = a/(1.0+dt*discount_factor(rCoordomp)); // for + lambda * u parameter (or steady equations)

      if(c==0)
        minA = a;
      else
        minA = (*OPTIM_SL)(minA,a);
    }
    v[i]=minA;
  }

  }// end of pragma parallel

  return;
}

