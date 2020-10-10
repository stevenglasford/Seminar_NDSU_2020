#include "HJB.h"

//---------------------------------------------------------------------
//- [can be distributed] : contains following elements for file HJB.cpp
//---------------------------------------------------------------------
//- int HJB::compute_optimal_trajectory    (const double* initialpoint, const char* filename, int number, double* terminalpoint, double &terminaltime)
//-    : reconstruction based on on tmin
//- int HJB::compute_optimal_trajectory2   (const double* initialpoint, const char* filename, int number, double* terminalpoint, double &terminaltime, ...)
//-    : reconstruction based on on tmin - for 2 controls
//- int HJB::compute_optimal_trajectory_val(const double* initialpoint, const char* filename, int number, double time_start) 
//-    : reconstr. based on value
//- int HJB::find_optimal_control_val(const double* y, double h, double t, double& val, double *vtab) : find optimal control based on vtab[.]
//-
//- Warning: here only for the case OPTIM=MAXIMUM problems (means "max" in the HJB equation)
//-

#include "compute_ot.h"        //- based on topt, for min or max problems
#include "compute_ot2.h"       //- based on topt, for two player games (min/max or max/min pbs)
#include "compute_otval.h"     //- based on value

//---------------------------------
//- find the optimal control for a given point 
//---------------------------------
#include "find_oc.h"           //- for min or max problems:    function HJB::find_optimal_control_val
#include "find_oc2.h" 	       //- including min/max problems: function HJB::find_optimal_control_val (-> will become val2)


