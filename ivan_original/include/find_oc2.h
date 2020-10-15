//---------------------------------
//- find the optimal command for each iteration of the compute trajectory algorithm
//-     interpolation on vtab + h distributed_cost
//- [NOV 2016] 2 control case
//---------------------------------
int HJB::find_optimal_control_val2
  //(const double* y, double h, double t, double& val, double *vtab)
  (const double* y, double h, double t, double& val, double *vtab, int& c2opt)
{
  //- returns the index of an optimal control ("c1opt")
  //- returns in "c2opt" the corresponding worst-case index for control u2[c2].
  double val_star,vect[ncall]; 
  int c,d,test;
  int ncall2=u2->ncall;    // - [NOV 2016] 2 controls
  double vect2[u2->ncall]; // - [NOV 2016] 2 controls

  int c_star;
  double ff[dim], ff2[dim], vectdouble[dim];

  //- note: epsOPTIM=1 if OPTIM=MAXIMUM, and epsOPTIM=-1 if OPTIM=MINIMUM
  double eps=epsOPTIM;

  //%OPTIM_PARAM = general_params.OPTIM;
  //switch(OPTIM_PARAM)
  //{
  //  case MAXIMUM:  //- Max HJB equation; scheme= control will try to minimize in the DPP
  //    eps= 1.0;
  //    break;
  //  case MINIMUM:  //- Min HJB equation; scheme= control will try to maximize in the DPP
  //    eps=-1.0;
  //    break;
  //  default:
  //    printf("find_optimal_control_val(): this OPTIM_PARAM not programmed for definition of eps!. Aborting.\n");
  //    exit(1);
  //}

  //- OB:2015
  double INFLOC=INF;
  if (TOPT_TYPE==1) INFLOC=-1e-5;

  val_star=eps*INFLOC; 
  c_star=-1;

/*
  double z, maxA, minA=INF;
  for(int c=0;c<ncall;c++){
    maxA = -INF;
    for(int c2=0;c2<u2->ncall;c2++){
      (*dynamics2)(x,(*u)[c],(*u2)[c2],t,dvectdouble);
      z = 0.;
      for(int d=0;d<dim;d++)
        z += MAX(0.,dvectdouble[d])*dv[2*d] + MIN(0.,dvectdouble[d])*dv[2*d+1];
      z=z+(*distributed_cost2)(x,(*u)[c],(*u2)[c2],t);
      maxA = MAX(maxA,z);
    }
    minA = MIN(maxA, minA);
  }
  return minA;
*/

  // ------------------
  // - OB [NOV 2016]: if eps=1 then look for optimal control u[c] realizing  { min_{u[c]} max_{u2[c2]} vtab[traj] }
  // ------------------
  for(c=0;c<ncall;c++) {

    // - [NOV 2016]  We aim to compute first max_{u2[c2]] topt[ y^{u[c],u2[c2]}_x(h) ]  (if eps=1) : worst adverse control
    //               We stop if there is an adverse control that leads the traj. out of the domain (or if topt=+INF)

    // - initialisation for 2nd control: [NOV 2016]
    double INFLOC2;
    if (TOPT_TYPE==0)
      INFLOC2=-INF; // case INFLOC= INF;
    else
      INFLOC2= INF; // case INFLOC=-1e-5;
    //<> printf("TOPT_TYPE=%i\n", TOPT_TYPE); exit(1);

    double val_star2=INFLOC2; 
    int c_star2=-1;

    
    for(int c2=0;c2<ncall2;c2++){

      // RK2 computation

      (*dynamics2)(y,(*u)[c],(*u2)[c2],t,ff2);
      for(d=0;d<dim;d++)
        vectdouble[d]=y[d]+h*ff2[d];

      (*dynamics2)(vectdouble,(*u)[c],(*u2)[c2],t+h,ff);
      for(d=0;d<dim;d++)
        vectdouble[d]=y[d]+h/2.*(ff[d]+ff2[d]);

      //- Possible periodization:
      //- for periodic variables some components could be corrected.
      if(periodic_mesh)
        periodizePoint(vectdouble);

      //- "test=1"  if next point out of the domain
      test=0;
      for(d=0; d<dim; d++){
        if(vectdouble[d]<=mesh->GLOBAL_lowBounds[d] || vectdouble[d]>=mesh->GLOBAL_highBounds[d]){
          test=1;
          break;
        }
      }



      //- test if next position is in the domain.
      //- if yes : compute the value
      //- if no  : value=INFLOC
      //double Fnorm, fact1=1.0,fact2=0.00, rc,
      //Fnorm=0;
      //for(d=0; d<dim; d++)
      //  Fnorm=max(Fnorm, abs(ff2[d]));
      //fact1=1.000, fact2=0.000; //- DEFAULT VALUES
      //fact1=1.000, fact2=0.1*h;   //- SMALL WEIGHT ON THE NORM OF THE DYNAMICS
      //vect[c] = fact1*rc + fact2*Fnorm;
      //vect[c] = rc;
      //- REM: for a minimal time problem, distributed cost should be defined as 0.
      //-      however, adding a small dist. cost proportional
      //-      to the norm of the dynamics can stabilize the reconstruction
      double v_interpol;

      if(!test){

        if(mesh->MESH_TYPE == 0){
          //- We test if current point out of the interior domain (useful only if MESH=0)
          int test_out=0;
          for(d=0; d<dim; d++){
            if(vectdouble[d]<=  mesh->lowBounds[d] || vectdouble[d]>=mesh->highBounds[d]){
              test_out=1;
              break;
            }
          }
          if(test_out==0)
            v_interpol=(*this.*interpolation)(vectdouble,vtab);  // interpolation inside comoput. domain.
          else
            v_interpol=(*this.*VbordCompute)(t,vectdouble,vtab); // extrapolation outside of comput. domain.
        }
        else
          v_interpol=(*this.*interpolation)(vectdouble,vtab);  //- NORMAL THING

        //vect[c]= v_interpol +  h*(*distributed_cost)(y,(*u)[c],t);
        vect2[c2]= v_interpol +  h*(*distributed_cost2)(y,(*u)[c],(*u2)[c2],t);  // - [NOV 2016] 2 controls

      }
      else
      {
        vect2[c2]= eps*INFLOC2;
      }
      //printf("vtab(interpol)=%8.5f, vect[c]=%8.5f\n",v_interpol,vect[c]);
       //printf("hb=%10.5f, GLOBAL_hb=%10.5f\n",mesh->highBounds[1], mesh->GLOBAL_highBounds[1]);

      // -----------------------------------------------------
      // -- [NOV 2016] 2 controls &  adapted to TOPT_TYPE=1 : exit time / tmax problem
      // -- find MAXimizing (resp. MINimising) optimal control u2[c2]
      // -----------------------------------------------------
      if(TOPT_TYPE==0){
        if(-eps*(vect2[c2]-val_star2)<0.0){
          val_star2=vect2[c2];
          c_star2=c2;
          //printf("yes ! val_star=%5.2f, c=%i,\n",val_star,c);
        }
      }
      else{
        if(-eps*(vect2[c2]-val_star2)>0.0){
          val_star2=vect2[c2];
          c_star2=c2;
          //printf("yes ! val_star=%5.2f, c=%i,\n",val_star,c);
        }
      }

    } // end of c2 loop


    
    //val_star=val_star2;
    vect[c]=val_star2;

    //- 2015 adapted to TOPT_TYPE=1 : exit time / tmax problem
    //- find minimising (resp. maximising) optimal control
    if(TOPT_TYPE==0){
      if(eps*(vect[c]-val_star)<0.0){
        val_star=vect[c];
        c_star=c;
        c2opt=c_star2;
        //printf("yes ! val_star=%5.2f, c=%i,\n",val_star,c);
      }
    }
    else{
      if(eps*(vect[c]-val_star)>0.0){
        val_star=vect[c];
        c_star=c;
        c2opt=c_star2;
        //printf("yes ! val_star=%5.2f, c=%i,\n",val_star,c);
      }
    }

  }

  //- nov 2015
  if(TOPT_TYPE==0){
    if(eps*val_star>=INFLOC){
      //printf("eps=%6.3f, val_star=%8.5f, INFLOC=%8.5f\n",eps,val_star,INFLOC);
      printf("c=%i; val_star>=INFLOC ==> no admissible controls !",c);
      c_star=-1;
    }
  }
  else{ // TOPT_TYPE==1
    if(eps*val_star<=INFLOC){
      printf("c=%i; val_star<=INFLOC ==> no admissible controls !",c);
      c_star=-1;
    }
  }


  val=val_star;
  return(c_star);
}

