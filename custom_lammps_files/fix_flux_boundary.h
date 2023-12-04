/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(flux/boundary,FixFluxBoundary)

#else

#ifndef LMP_FIX_FLUX_BOUNDARY_H
#define LMP_FIX_FLUX_BOUNDARY_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixFluxBoundary : public Fix { 
 public:
  FixFluxBoundary(class LAMMPS *, int, char **);
  ~FixFluxBoundary();
  int setmask();
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  void post_integrate(); 
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);


 protected:
  
  int me; 
  int dim; // x y or z direction to be used for entrance/exit comparisons
  int out; // -1/+1 to indicate which direction brings you out of simulation box, 0 if both directions are valid
  double box_half;
  double **x;

  long filepos;
  
  FILE *fp; // file to output to

};

}

#endif
#endif
