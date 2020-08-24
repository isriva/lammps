/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_strain_rate.h"
#include "fix_nh_sphere.h"
#include <cmath>
#include <cstring>
#include <mpi.h>
#include "update.h"
#include "domain.h"
#include "math_extra.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "error.h"
#include "modify.h"


#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeStrainRate::ComputeStrainRate(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute strain_rate command");
  if (igroup) error->all(FLERR,"Compute strain_rate must use group all");


  vector_flag = 1;
  extvector = 0;
  size_vector = 7;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

void ComputeStrainRate::init()
{
  ifix = modify->find_fix_by_style("^nph/sphere");
  if (ifix == -1)
    error->all(FLERR,"Fix nph/sphere required for compute strain_rate");
}

/* ---------------------------------------------------------------------- */

ComputeStrainRate::~ComputeStrainRate()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeStrainRate::compute_vector()
{
  // hrate, hinv and velgrad are of the form:
  // [ 0 5 4 ]
  // [ - 1 3 ]
  // [ - - 2 ]

  double velgrad[6];
  double *hrate = ((FixNHSphere *) modify->fix[ifix])->hrate;
  MathExtra::multiply_shape_shape(hrate,domain->h_inv,velgrad);

  for (int i=0;i<6;i++) {
    vector[i] = velgrad[i];
  }

  double d0 = (1.0/3.0)*(velgrad[0] + velgrad[1] + velgrad[2]);
  double sum = 0.0;
  sum += (velgrad[0]-d0)*(velgrad[0]-d0);
  sum += (velgrad[1]-d0)*(velgrad[1]-d0);
  sum += (velgrad[2]-d0)*(velgrad[2]-d0);
  sum += 0.5*velgrad[3]*velgrad[3];
  sum += 0.5*velgrad[4]*velgrad[4];
  sum += 0.5*velgrad[5]*velgrad[5];
  vector[6] = sqrt(0.5*sum);
}
