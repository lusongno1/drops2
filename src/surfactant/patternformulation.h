/// \file
/// \brief Discretization for PDEs on an interface.
/// \author LSEC: Song Lu

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "surfactant/heatsolver.h"
#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "levelset/levelsetmapper.h"
#include "out/output.h"
#include "out/vtkOut.h"
#include "misc/dynamicload.h"
#include "misc/funcmap.h"
#include "misc/omp_variable.h"
#include "geom/subtriangulation.h"
#include "num/gradient_recovery.h"
/*
#include "surfphasesep/separation.
*/
#include <cmath>
#include <fstream>
#include <string>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "out/output.h"
#include "geom/geomselect.h"
#include "misc/funcmap.h"

#include "geom/deformation.h"

// include numeric computing!
#include "num/fe.h"
#include "num/krylovsolver.h"
#include "num/MGsolver.h"
#include "poisson/integrTime.h"
#include "num/prolongation.h"

// include problem class
#include "misc/params.h"
#include "poisson/poissonParam.h"      // poissonCoeffCL
#include "poisson/poisson.h"           // setting up the Poisson problem
#include "num/bndData.h"

// include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

//include output
#include "out/ensightOut.h"
#include "out/vtkOut.h"

// include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#endif

// include function container
#include "misc/funcmap.h"
#include "num/poissonsolverfactory.h"
#include "poisson/ale.h"

#include "misc/progressaccu.h"
#include "misc/dynamicload.h"


#ifndef DROPS_PATTERNFORMULATION_H
#define DROPS_PATTERNFORMULATION_H

namespace DROPS
{
//typedef double (*dist_funT) (const DROPS::Point3DCL&, double);
//void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t= 0.);
//void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, instat_vector_fun_ptr LsgVel, double t= 0.);
//SurfactantP1BaseCL* make_surfactant_timedisc( MultiGridCL& mg, LevelsetP2CL& lset,
//        VecDescCL& v, VecDescCL& nd, const BndDataCL<Point3DCL>& Bnd_v,
//        const ParamCL& P, const double & dist=0);
//std::unique_ptr<VTKOutCL> vtkwriter;
//template<class DiscP1FunType>
//double L2_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
//                 const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol);
//double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
//                DROPS::instat_scalar_fun_ptr extsol);
//                template<class DiscP1FunType>
//double H1_error (const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
//                 const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol);

void ContantInit (double uw0, VecDescCL& ic, const MultiGridCL& mg, double t);
class PatternFormulationCL
{
public:
    VecDescCL ic;  ///< concentration on the interface at current time
    VecDescCL icw;  ///< concentration on the interface at current time w.r.t
    DROPS::MultiGridCL& mg;
    DROPS::LevelsetP2CL& lset;
    DROPS::AdapTriangCL& adap;
    VecDescCL lsgradrec;//the recovered gradient of the level set function to calculate curvature
    instat_scalar_fun_ptr the_lset_fun;
    instat_vector_fun_ptr the_normal_fun;
    instat_scalar_fun_ptr the_rhs_fun;
    instat_scalar_fun_ptr the_sol_fun;
    double d1,d2,gamma,a,b,tEnd,delta,epsilon,dT;
    double cur_time=0;
    double dist;
    IdxDescCL idx;
private:
//    SurfactantP1BaseCL& timedisc;
protected:

public:
    PatternFormulationCL (DROPS::MultiGridCL& mg,DROPS::AdapTriangCL& adap,
                          DROPS::LevelsetP2CL& lset,instat_scalar_fun_ptr the_lset_fun,instat_vector_fun_ptr the_normal_fun,
                          instat_scalar_fun_ptr the_rhs_fun,instat_scalar_fun_ptr the_sol_fun);
    //virtual void DoStepRD (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset)
    virtual void DoStepRD ();
    virtual void DoStepHeat();
    virtual void DoStepHeat2();
    void P2ConstantInit (double uw0, VecDescCL& ic, const MultiGridCL& mg, double t);
    void P1ConstantInit (double uw0, VecDescCL& ic, const MultiGridCL& mg, double t);
    void GetGradientOfLevelSet();

};

}

//#include "surfactant/ifacetransp.tpp"

#endif

