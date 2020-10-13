//- NEW SL ORDER 2 function: BASIC ORDER 1 SCHEME (OK)
void HJB_SL::secondorder_itSL_evo(double t, double dt)
{
  //- P=PARAMP of the data file.
  double a=0.0;
  double minA=0.0;
  double coords[dim];
  int i,j,c,e,d,p,test;

  //for(i=0; i<mesh->inn_nbPoints; i++){
  for(j=0;j<ranksize;j++){
    i   =rank[j];
    vold[i]=v[i];
  }

  //for(i=0; i<mesh->inn_nbPoints; i++){
  for(j=0;j<ranksize;j++){
    i   =rank[j];

    (mesh->*(mesh->setcoords))(i,rCoord);

    for(c=0;c<ncall;c++){
      a = 0.0;

      //- this: (p<=P) for use with data/data_SL2_optionput_approach2.h
      for(p=0;p<P;p++){
        for(e=-1;e<=1;e+=2){ //- two points formula
          funcY(rCoord, p, e, (*u)[c], t, dt, dvectdouble);
          if(periodic_mesh)
            ComputePeriodic(dvectdouble, coords);
          else
            for(d=0;d<dim;d++)
              coords[d] = dvectdouble[d];
          test=0;
          for(d=0;d<dim;d++){
            //double epsilon_mesh2=1e-12;
            //double epsilon_mesh2=0.0;
            //if(coords[d] <= mesh->lowBounds[d]-epsilon_mesh2 || coords[d] >= mesh->highBounds[d]+epsilon_mesh2){
            if(coords[d] < lb[d] || coords[d] > hb[d]){
              test=1;
              break;
            }
          }
          if(test)
            a += (*this.*VbordCompute)(t,coords,vold);
          else
            a += (*this.*interpolation)(coords,vold);
        }
      }
      a  = a / (2*P);

      a += dt * (*distributed_cost)(rCoord,(*u)[c],t);

      a *= exp(-funcR(rCoord, (*u)[c], t)*dt);

      a = a/(1.0+dt*discount_factor(rCoord)); // for + lambda * u parameter (or steady equations)

      if(c==0)
        minA = a;
      else
        minA = (*OPTIM_SL)(minA,a);
    }
    v[i]=minA;
  }



  return;
}

