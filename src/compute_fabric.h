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

#ifdef COMPUTE_CLASS

ComputeStyle(fabric,ComputeFabric)

#else

#ifndef LMP_COMPUTE_FABRIC_H
#define LMP_COMPUTE_FABRIC_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFabric : public Compute {
 public:
  ComputeFabric(class LAMMPS *, int, char **);
  ~ComputeFabric();
  void init();
  void init_list(int, class NeighList *);  
  void compute_vector();
  
 private:
  int ntensors, pstyle, cutstyle;
  int *tensor_style;              // style of each requested tensor
  class NeighList *list;
  
  int cn_flag, br_flag, fn_flag, ft_flag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
