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

FixStyle(passage,FixPassage)

#else

#ifndef LMP_FIX_PASSAGE_H
#define LMP_FIX_PASSAGE_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixPassage : public Fix { 
 public:
  FixPassage(class LAMMPS *, int, char **);
  ~FixPassage();
  int setmask();
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  void end_of_step(); 
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);


 protected:
  
  int me; 
  int dim; // x y or z direction to be used for entrance/exit comparisons
  double entr, exit1, exit2, max_step; // exit1 is aborted, exit2 is successful
  double **x;

  long filepos;
  
  FILE *fp; // file to output to

};

}

#endif
#endif
