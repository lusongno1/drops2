/// \file poissonParam.h
/// \brief PoissonCoeffCL
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Liang Zhang, Thorolf Schulte; SC RWTH Aachen:
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

#ifndef POISSONPARAM_H_
#define POISSONPARAM_H_

#include "misc/container.h"
#include "misc/params.h"
#include "misc/funcmap.h"
#include <sstream>
#include "geom/multigrid.h"

namespace DROPS
{

/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/
static double delta_;
static VecDescCL ic_;
class PoissonCoeffCL
{
private:
    static ParamCL C_;
    static double dx_, dy_;
    static int    nx_, ny_;
    static double dt_;

public:
    //reaction
    static scalar_tetra_function q;
    //diffusion
    static double alpha;
    //source
    static scalar_tetra_function f;
    //initial condition
    static instat_scalar_fun_ptr InitialCondition;
    //solution
    static scalar_tetra_function Solution;
    //static instat_scalar_fun_ptr Solution;
    //velocity
    static vector_tetra_function Vel;
    //Free interface function
    static instat_scalar_fun_ptr interface;

    static double ff(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t)
    {
        //const DROPS::Point3DCL p= DROPS::GetWorldCoord( tet, b);
        NoBndDataCL<> nobnddata;
        LocalP1CL<> icp1(tet, ic_, nobnddata);
        return -delta_*icp1(b);
        //return 0;
        //val(const TetraCL& s, const BaryCoordCL&) const;
        //return (-2.0*std::exp(t)*std::exp(p[0]+p[1]+p[2]));
    }

    PoissonCoeffCL( ParamCL& P)
    {
        C_=P;
        nx_= P.get<int>("Mesh.N1");
        ny_= P.get<int>("Mesh.N2");
        dx_= norm(P.get<DROPS::Point3DCL>("Mesh.E1"));
        dy_= norm(P.get<DROPS::Point3DCL>("Mesh.E2"));
        dt_= P.get<int>("Time.NumSteps")!=0 ? P.get<double>("Time.FinalTime")/P.get<int>("Time.NumSteps") : 0;  //step size used in ALEVelocity
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        DROPS::ScaTetMap & scatet = DROPS::ScaTetMap::getInstance();
        q = scatet[P.get<std::string>("Poisson.Coeff.Reaction")];
        alpha = P.get<double>("Poisson.Coeff.Diffusion");
        f = scatet[P.get<std::string>("Poisson.Coeff.Source")];
        Solution = scatet[P.get<std::string>("Poisson.Solution")];
        //Solution = scamp[P.get<std::string>("Poisson.Solution")];
        InitialCondition = scamap[P.get<std::string>("Poisson.InitialValue")];
        DROPS::VecTetMap & vectet = DROPS::VecTetMap::getInstance();
        Vel = vectet[P.get<std::string>("Poisson.Coeff.Flowfield")];

        interface = scamap[P.get<std::string>("ALE.Interface")];

    }

    PoissonCoeffCL( ParamCL& P, ParamCL& P2)
    {
        C_=P2;
        nx_= P.get<int>("Mesh.N1");
        ny_= P.get<int>("Mesh.N2");
        dx_= norm(P.get<DROPS::Point3DCL>("Mesh.E1"));
        dy_= norm(P.get<DROPS::Point3DCL>("Mesh.E2"));
        dt_= P2.get<int>("Time.NumSteps")!=0 ? P2.get<double>("Time.FinalTime")/P2.get<int>("Time.NumSteps") : 0;  //step size used in ALEVelocity
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        DROPS::ScaTetMap & scatet = DROPS::ScaTetMap::getInstance();
        q = scatet[P2.get<std::string>("Poisson.Coeff.Reaction")];
        alpha = P2.get<double>("Poisson.Coeff.Diffusion");
        f = scatet[P2.get<std::string>("Poisson.Coeff.Source")];
        Solution = scatet[P2.get<std::string>("Poisson.Solution")];
        //Solution = scamp[P.get<std::string>("Poisson.Solution")];
        InitialCondition = scamap[P2.get<std::string>("Poisson.InitialValue")];
        DROPS::VecTetMap & vectet = DROPS::VecTetMap::getInstance();
        Vel = vectet[P2.get<std::string>("Poisson.Coeff.Flowfield")];

        interface = scamap[P2.get<std::string>("ALE.Interface")];

    }


    PoissonCoeffCL( ParamCL& P1,ParamCL& P2,double epsilon,double delta,VecDescCL ic,double dT)
    {
        delta_ = delta;
        ic_ = ic;
        C_=P2;
        nx_= P1.get<int>("Mesh.N1");
        ny_= P1.get<int>("Mesh.N2");
        dx_= norm(P1.get<DROPS::Point3DCL>("Mesh.E1"));
        dy_= norm(P1.get<DROPS::Point3DCL>("Mesh.E2"));
        //dt_= P2.get<int>("Time.NumSteps")!=0 ? dT/P2.get<int>("Time.NumSteps") : 0;  //step size used in ALEVelocity
        dt_= P2.get<int>("Time.NumSteps")!=0 ? dT/P2.get<int>("Time.NumSteps") : 0;  //step size used in ALEVelocity
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        DROPS::ScaTetMap & scatet = DROPS::ScaTetMap::getInstance();
        q = scatet[P2.get<std::string>("Poisson.Coeff.Reaction")];
        //alpha = P2.get<double>("Poisson.Coeff.Diffusion");
        alpha = epsilon;
        //f = scatet[P.get<std::string>("Poisson.Coeff.Source")];
        f = &ff;
        Solution = scatet[P2.get<std::string>("Poisson.Solution")];
        //Solution = scamp[P.get<std::string>("Poisson.Solution")];
        InitialCondition = scamap[P2.get<std::string>("Poisson.InitialValue")];
        DROPS::VecTetMap & vectet = DROPS::VecTetMap::getInstance();
        Vel = vectet[P2.get<std::string>("Poisson.Coeff.Flowfield")];

        interface = scamap[P2.get<std::string>("ALE.Interface")];

    }

    static Point3DCL ALEDeform(const Point3DCL& p, double t)
    {

        DROPS::Point3DCL ret= Point3DCL(0.);
        if(t == -1)
        {
            ret[0] = p[0];
            ret[1] = p[1] *  interface(p, 0.)/dy_;
            ret[2] = p[2];
        }
        else
        {
            ret[0] = p[0];
            ret[1] = p[1] * interface(p, t + dt_)/dy_;
            ret[2] = p[2];
        }
        return ret;
    }
    static Point3DCL ALEVelocity(const DROPS::TetraCL& tet, const DROPS::BaryCoordCL& b, double t)
    {
        double eps =1.0e-7 * dt_;
        DROPS::Point3DCL ret= Vel(tet, b, t);
        const DROPS::Point3DCL& p= DROPS::GetWorldCoord(tet,b);
        ret[1] -= p[1]/interface(p, t)*(interface(p, t+eps)-interface(p, t))/eps;  //y/h(p,t)*h_p'(t)
        return ret;
    }


};

}//end of namespace
#endif
