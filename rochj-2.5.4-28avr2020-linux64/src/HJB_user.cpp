#include "HJB.h"

// Based on topt, for min or max problems
#include "compute_ot.h"
// Based on topt, for two player games (min/max or max/min pbs)
#include "compute_ot2.h"
// Base on value
#include "compute_otval.h"

///////////////////////////////
//Find the optimal control for a given point
///////////////////////////////
//for min or max problems: function HJB::find_optimal_control_val
#include "find_oc.h"
//including min/max problems: function HJB::find_optimal_control_val (-> will become val2)
