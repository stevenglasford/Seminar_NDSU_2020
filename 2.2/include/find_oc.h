//---------------------------------
//- find the optimal command for each iteration of the compute trajectory algorithm
//-     interpolation on vtab + h distributed_cost
//---------------------------------
int HJB::find_optimal_control_val
  (const double* y, double h, double t, double& val, double *vtab)
{
  //- returns the index of an optimal control
  double val_star,vect[ncall];
  int c,d,test;

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


  for(c=0;c<ncall;c++) {

    // RK2 computation
    { 
      (*dynamics)(y,(*u)[c],t,ff2);
      for(d=0;d<dim;d++)
        vectdouble[d]=y[d]+h*ff2[d];

      (*dynamics)(vectdouble,(*u)[c],t+h,ff);
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

    }

    //- test if next position is in the domain.
    //- if yes : compute the value
    //- if no  : value=INFLOC
    //double Fnorm, fact1=1.0,fact2=0.00, rc,
    //- Fnorm=norm(ff2)
    double Fnorm=0.0; 
    for(d=0; d<dim; d++) Fnorm=max(Fnorm, abs(ff2[d]));				// norm-Linf
    //for(d=0; d<dim; d++) Fnorm+=abs(ff2[d]);  				// norm-L1
    //for(d=0; d<dim; d++) Fnorm+=pow(abs(ff2[d]),2.0); Fnorm=sqrt(Fnorm);	// norm-L2
    //vect[c] = fact1*rc + fact2*Fnorm;
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

      vect[c] = v_interpol +  h*(*distributed_cost)(y,(*u)[c],t);

      //- to add a small penalisation with respect to the control, do for instance:  
      //C uu=(*u)[c]; vect[c]+= 0.001*abs(uu[1]);
      //- to add a small penalisation with respect to norm of the dynamics, do for instance:  
      //double fact1=0.001; vect[c]+= fact1*Fnorm;

    }
    else
    {
      vect[c]=eps*INFLOC;
    }
    //printf("vtab(interpol)=%8.5f, vect[c]=%8.5f\n",v_interpol,vect[c]);
    //printf("hb=%10.5f, GLOBAL_hb=%10.5f\n",mesh->highBounds[1], mesh->GLOBAL_highBounds[1]);


    //- 2015 adapted to TOPT_TYPE=1 : exit time / tmax problem
    //- find minimising (resp. maximising) optimal control
    if(TOPT_TYPE==0){
      if(eps*(vect[c]-val_star)<0.0){
        val_star=vect[c];
        c_star=c;
        //printf("yes ! val_star=%5.2f, c=%i,\n",val_star,c);
      }
    }
    else{
      if(eps*(vect[c]-val_star)>0.0){
        val_star=vect[c];
        c_star=c;
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

