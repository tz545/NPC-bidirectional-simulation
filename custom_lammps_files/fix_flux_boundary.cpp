/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

written by Tiantian Zheng

NOTE: this code checks for crossings of periodic boundaries only

------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "fix_flux_boundary.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "domain.h"
#include "region.h"
#include "force.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;



// pass in as ID groupID flux/boundary nevery z out box_half filename 


FixFluxBoundary::FixFluxBoundary(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp, narg, arg),
 me(0), fp(NULL), dim(2), box_half(0), out(0){

  restart_peratom = 1;

  if (narg < 8) error->all(FLERR, "Fix flux/boundary requires 8 arguments");

  nevery = atoi(arg[3]);
  if (nevery != 1)
    error->warning(FLERR, "Not performing fix flux/boundary check every timestep");
  if (nevery <= 0)
    error->all(FLERR, "Illegal fix flux/boundary command");

  if (strcmp(arg[4], "x") == 0) dim = 0;
  else if (strcmp(arg[4], "y") == 0) dim = 1;
  else if (strcmp(arg[4], "z") == 0) dim = 2;
  else error->all(FLERR, "Unrecognized axis specified for fix flux/boundary");

  out = atof(arg[5]); // +/-1 corresponds to record upper/lower boundary crossing
  box_half = atoi(arg[6]); // half of length of simulation box in required dimension


  if (me == 0){
    fp = fopen(arg[7], "a");
    printf("%s\n", arg[7]);
    if (fp == NULL){
      char str[128];
      sprintf(str, "Cannot open fix flux/boundary file %s", arg[7]);
      error->one(FLERR, str);
    }
  }


  if (fp && me == 0) {
    clearerr(fp);
    if (ferror(fp)) error->one(FLERR,"Error writing file header");
 
    filepos = ftell(fp);
  }


  int nmax = atom->nmax;
  int nlocal = atom->nlocal;
  double **x_c = atom->x;

  memory->create(this->x, nmax, 4, "FixFluxBoundary:x");

  // atom->nmax is the total number of local and ghost atoms that can be held on a processor
  atom->add_callback(0); // this registers the fix to have grow_arrays called whenever nmax changes
  atom->add_callback(1); // or whenever a restart is read in. 

  for (int particleInd = 0; particleInd < nlocal; ++particleInd){
    if (atom->mask[particleInd] & groupbit){
      x[particleInd][0] = x_c[particleInd][0];
      x[particleInd][1] = x_c[particleInd][1];
      x[particleInd][2] = x_c[particleInd][2];
      x[particleInd][3] = 0.0;
    }
  }

}
// 
/* ---------------------------------------------------------------------- */

FixFluxBoundary::~FixFluxBoundary()
{

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  memory->destroy(x);

  if (fp && me == 0) fclose(fp);

}

/* ---------------------------------------------------------------------- */

int FixFluxBoundary::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;

  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixFluxBoundary::memory_usage() // check whether we want sum of memory of all arrays
{
  int nmax = atom->nmax;
  double bytes = nmax * 4 * sizeof(double);

  return bytes;

}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixFluxBoundary::grow_arrays(int nmax)
{
  memory->grow(this->x,nmax,4,"FixFluxPlane:x");
  // fprintf(fp, "grow arrays called \n");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixFluxBoundary::copy_arrays(int i, int j, int delflag)
{
  memcpy(this->x[j], this->x[i], sizeof(double) * 4);
  // fprintf(fp, "copy arrays called \n");

}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixFluxBoundary::set_arrays(int i)
{
  memset(this->x[i], 0, sizeof(double) * 4);
  // fprintf(fp, "set arrays called \n");

}

/*------------------------------------------------------------------------
  store atom's data in a buffer
-------------------------------------------------------------------------- */

int FixFluxBoundary::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = x[i][3];

  // fprintf(fp, "pack exchange called \n");

  return m;
}

/*------------------------------------------------------------------------
  retrieve atom's data from a buffer
-------------------------------------------------------------------------- */

int FixFluxBoundary::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  x[nlocal][3] = buf[m++];

  // fprintf(fp, "unpack exchange called \n");

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixFluxBoundary::pack_restart(int i, double *buf)
{
  int m = 1;
  // buf[m++] = 6;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = x[i][3];
  buf[0] = m;

  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixFluxBoundary::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  x[nlocal][0] = extra[nlocal][m++];
  x[nlocal][1] = extra[nlocal][m++];
  x[nlocal][2] = extra[nlocal][m++];
  x[nlocal][3] = extra[nlocal][m++];

}


/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixFluxBoundary::maxsize_restart()
{

  return 5;
}

//  ----------------------------------------------------------------------
//    size of atom nlocal's restart data
// ------------------------------------------------------------------------- 

int FixFluxBoundary::size_restart(int nlocal)
{

  return 5;
}


/*------------------------------------------------------------------------
  calculations, comparisons here
-------------------------------------------------------------------------- */

void FixFluxBoundary::post_integrate()
{
  bigint ntimestep;
  ntimestep = update->ntimestep;

  double **x_c = atom->x; // current position

  // tagint *tag = atom->tag;

  int nlocal = atom->nlocal;


  for (int particleInd = 0; particleInd < nlocal; ++particleInd){

    if (atom->mask[particleInd] & groupbit){

      // check if particle has covered more than half the simulation box
      // (implying it has crossed a periodic boundary)

      if (abs(x_c[particleInd][dim] - x[particleInd][dim]) > box_half){

        // crossed lower boundary
        if (x_c[particleInd][dim] - x[particleInd][dim] > box_half){

          if (out == -1){
            // record lower boundary crossing
            if (fp && me == 0) {   
              clearerr(fp);
              fprintf(fp, "%d %ld\n", atom->type[particleInd], ntimestep); 
              fflush(fp);
            }
          }

        }
        

        // crossed upper boundary
        else {

           if (out == 1){
            if (fp && me == 0) {   
              clearerr(fp);
              fprintf(fp, "%d %ld\n", atom->type[particleInd], ntimestep); 
              fflush(fp);
            }
          }

        }

        
      }

      x[particleInd][0] = x_c[particleInd][0];
      x[particleInd][1] = x_c[particleInd][1];
      x[particleInd][2] = x_c[particleInd][2];
    }
  }
}