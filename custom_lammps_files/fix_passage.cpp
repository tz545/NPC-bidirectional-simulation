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

------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "fix_passage.h"
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



// pass in as ID groupID fix/passage nevery z entr exit filename max_step


FixPassage::FixPassage(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp, narg, arg),
 me(0), fp(NULL), dim(2), entr(0), exit1(0), exit2(0){

  restart_peratom = 1;

  if (narg < 10) error->all(FLERR, "Fix/passage requires 9 arguments");

  nevery = atoi(arg[3]);
  if (nevery != 1)
    error->warning(FLERR, "Not performing fix/passage check every timestep");
  if (nevery <= 0)
    error->all(FLERR, "Illegal fix passage command");

  if (strcmp(arg[4], "x") == 0) dim = 0;
  else if (strcmp(arg[4], "y") == 0) dim = 1;
  else if (strcmp(arg[4], "z") == 0) dim = 2;
  else error->all(FLERR, "Unrecognized axis specified for fix passage");

  entr = atof(arg[5]);
  exit1 = atof(arg[6]);
  exit2 = atof(arg[7]);


  if (me == 0){
    fp = fopen(arg[8], "a");
    printf("%s\n", arg[8]);
    if (fp == NULL){
      char str[128];
      sprintf(str, "Cannot open fix flux/boundary file %s", arg[8]);
      error->one(FLERR, str);
    }
  }


  if (fp && me == 0) {
    clearerr(fp);
    if (ferror(fp)) error->one(FLERR,"Error writing file header");
 
    filepos = ftell(fp);
  }

  // how large to expect max step size, used to determine if particle is entering pore or passing through PBC
  max_step = atof(arg[9]); 

  int nmax = atom->nmax;
  int nlocal = atom->nlocal;
  double **x_c = atom->x;

  memory->create(this->x, nmax, 5, "FixPassage:x");

  // atom->nmax is the total number of local and ghost atoms that can be held on a processor
  atom->add_callback(0); // this registers the fix to have grow_arrays called whenever nmax changes
  atom->add_callback(1); // or whenever a restart is read in. 

  for (int particleInd = 0; particleInd < nlocal; ++particleInd){
    if (atom->mask[particleInd] & groupbit){
      x[particleInd][0] = x_c[particleInd][0];
      x[particleInd][1] = x_c[particleInd][1];
      x[particleInd][2] = x_c[particleInd][2];
      x[particleInd][3] = 0.0;
      x[particleInd][4] = 0.0;
    }
  }

}
// 
/* ---------------------------------------------------------------------- */

FixPassage::~FixPassage()
{

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  memory->destroy(x);

  if (fp && me == 0) fclose(fp);

}

/* ---------------------------------------------------------------------- */

int FixPassage::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;

  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPassage::memory_usage() // check whether we want sum of memory of all arrays
{
  int nmax = atom->nmax;
  double bytes = nmax * 5 * sizeof(double);

  return bytes;

}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPassage::grow_arrays(int nmax)
{
  memory->grow(this->x,nmax,5,"FixPassage:x");
  // fprintf(fp, "grow arrays called \n");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixPassage::copy_arrays(int i, int j, int delflag)
{
  memcpy(this->x[j], this->x[i], sizeof(double) * 5);
  // fprintf(fp, "copy arrays called \n");

}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPassage::set_arrays(int i)
{
  memset(this->x[i], 0, sizeof(double) * 5);
  // fprintf(fp, "set arrays called \n");

}

/*------------------------------------------------------------------------
  store atom's data in a buffer
-------------------------------------------------------------------------- */

int FixPassage::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = x[i][3];
  buf[m++] = x[i][4];

  // fprintf(fp, "pack exchange called \n");

  return m;
}

/*------------------------------------------------------------------------
  retrieve atom's data from a buffer
-------------------------------------------------------------------------- */

int FixPassage::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  x[nlocal][3] = buf[m++];
  x[nlocal][4] = buf[m++];

  // fprintf(fp, "unpack exchange called \n");

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPassage::pack_restart(int i, double *buf)
{
  int m = 1;
  // buf[m++] = 6;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = x[i][3];
  buf[m++] = x[i][4];
  buf[0] = m;


  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPassage::unpack_restart(int nlocal, int nth)
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
  x[nlocal][4] = extra[nlocal][m++];

}


/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPassage::maxsize_restart()
{

  return 6;
}

//  ----------------------------------------------------------------------
//    size of atom nlocal's restart data
// ------------------------------------------------------------------------- 

int FixPassage::size_restart(int nlocal)
{

  return 6;
}


/*------------------------------------------------------------------------
  calculations, comparisons here
-------------------------------------------------------------------------- */

void FixPassage::end_of_step()
{
  bigint ntimestep;
  ntimestep = update->ntimestep;

  double **x_c = atom->x;

  // tagint *tag = atom->tag; //remove

  int nlocal = atom->nlocal;

  // exit1 is aborted, exit2 is successful
  if (exit1 > exit2){
    for (int particleInd = 0; particleInd < nlocal; ++particleInd){

      if (atom->mask[particleInd] & groupbit){

        if (x[particleInd][3] == 1.0){
          // particle in pore, check if it has exited
          
          if (x_c[particleInd][dim] > exit1 && x[particleInd][dim] < exit1){
            // nanoparticle exited through entrance (abortive translocation)

            // write info to file in format | atom type | successful/not | Tfinal | Tinitial
            if (fp && me == 0) {   
              clearerr(fp);
              fprintf(fp, "%d 0 %ld %ld\n", atom->type[particleInd], ntimestep, static_cast<long int>(x[particleInd][4])); 
              fflush(fp);
            }


            // set Si, Ti back to zero
            x[particleInd][3] = 0.0;
            x[particleInd][4] = 0.0;
          }

          else if (x_c[particleInd][dim] < exit2 && x[particleInd][dim] > exit2){
            
            // successful translocation

            // write info to file
            if (fp && me == 0) {   
              clearerr(fp);
              fprintf(fp, "%d 1 %ld %ld\n", atom->type[particleInd], ntimestep, static_cast<long int>(x[particleInd][4])); 
              fflush(fp);
            }
            

            // set Si, Ti back to zero
            x[particleInd][3] = 0.0;
            x[particleInd][4] = 0.0;

          }

        }

        else if (x[particleInd][3] == 0.0){
          // particle not in pore, check if it has entered

          if (x_c[particleInd][dim] < entr && x[particleInd][dim] > entr && abs(x[particleInd][dim]-x_c[particleInd][dim]) < max_step){
            // nanoparticle has entered into pore

            x[particleInd][3] = 1.0;
            x[particleInd][4] = static_cast<double>(ntimestep);

          }
        }

        x[particleInd][0] = x_c[particleInd][0];
        x[particleInd][1] = x_c[particleInd][1];
        x[particleInd][2] = x_c[particleInd][2];
      }
    }
  }

  else if (exit1 < exit2){
    for (int particleInd = 0; particleInd < nlocal; ++particleInd){

      if (atom->mask[particleInd] & groupbit){

        // // TEST! remove
        // if (tag[particleInd] == 100){
        //   // if (fp && me == 0) {   
        //   // clearerr(fp);
        //   fprintf(fp, "%f %f %ld %ld\n", x[particleInd][dim], x_c[particleInd][dim], static_cast<long int>(x[particleInd][3]), static_cast<long int>(ntimestep));
        //   // fflush(fp);
        //   // }
        // }


        if (x[particleInd][3] == 1.0){
          // particle in pore, check if it has exited
          
          if (x_c[particleInd][dim] < exit1 && x[particleInd][dim] > exit1){
            // nanoparticle exited through entrance (abortive translocation)

            // write info to file in format | atom type | successful/not | Tfinal | Tinitial
            if (fp && me == 0) {   
              clearerr(fp);
              fprintf(fp, "%d 0 %ld %ld\n", atom->type[particleInd], ntimestep, static_cast<long int>(x[particleInd][4]));         
              fflush(fp);
            }

            
            // set Si, Ti back to zero
            x[particleInd][3] = 0.0;
            x[particleInd][4] = 0.0;
          }

          else if (x_c[particleInd][dim] > exit2 && x[particleInd][dim] < exit2){
            
            // successful translocation

            // write info to file
            if (fp && me == 0) {   
              clearerr(fp);
              fprintf(fp, "%d 1 %ld %ld\n", atom->type[particleInd], ntimestep, static_cast<long int>(x[particleInd][4])); 
              fflush(fp);
            }
            
            // set Si, Ti back to zero
            x[particleInd][3] = 0.0;
            x[particleInd][4] = 0.0;

          }

        }

        else if (x[particleInd][3] == 0.0){
          // particle not in pore, check if it has entered

          if (x_c[particleInd][dim] > entr && x[particleInd][dim] < entr && abs(x[particleInd][dim]-x_c[particleInd][dim]) < max_step){
            // nanoparticle has entered into pore

            x[particleInd][3] = 1.0;
            x[particleInd][4] = static_cast<double>(ntimestep);
          }
      
        }

        x[particleInd][0] = x_c[particleInd][0];
        x[particleInd][1] = x_c[particleInd][1];
        x[particleInd][2] = x_c[particleInd][2];




      }
    }
  }
}