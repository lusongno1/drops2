// include geometric computing
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
#include "levelset/levelset.h"
#include "misc/problem.h"

DROPS::ParamCL P;   //Parameter class, read in json file in main function
using namespace DROPS;

double the_lset_fun( const DROPS::Point3DCL& p, double t)
{
    return p.norm()-1;
    //return (std::exp(t)*std::exp(p[0]+p[1]+p[2]));
}
void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, double (*d)(const Point3DCL&, double), double t= 0.)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
    ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
    ls.t= t;
}
//DROPS::ParamCL P;
int main (int argc, char** argv)
{
    DROPS::read_parameter_file_from_cmdline( P, argc, argv, "../../param/poisson/cdrdrops/instatpoissonEx3.json");
    P.put_if_unset<std::string>("VTK.TimeFileName",P.get<std::string>("VTK.VTKName"));
    std::cout << P << std::endl;

    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"),
                       P.get<std::vector<std::string> >("General.DynamicLibs") );
    if (P.get<int>("General.ProgressBar"))
        DROPS::ProgressBarTetraAccumulatorCL::Activate();



    DROPS::MultiGridCL* mg= 0;
    DROPS::PoissonBndDataCL* bdata = new DROPS::PoissonBndDataCL(0);
    std::unique_ptr<DROPS::MGBuilderCL> builder( DROPS::make_MGBuilder( P));
    mg = new DROPS::MultiGridCL( *builder);

    for ( int ref=1; ref <= P.get<int>("Mesh.AdaptRef.FinestLevel"); ++ref)
    {
        std::cout << " refine (" << ref << ")\n";
        DROPS::MarkAll( *mg);
        mg->Refine();
    }

                    DROPS::dynamicLoad(P.get<std::string>("General.DynamicLibsPrefix"),
                           P.get<std::vector<std::string> >("General.DynamicLibs") );
    read_BndData( *bdata, *mg, P.get_child( "Poisson.BoundaryData"));
    for (int i=0; i<mg->GetBnd().GetNumBndSeg(); ++i)
        std::cout << i << ": BC = " << bdata->GetBndSeg(i).GetBC() << std::endl;

    //MeshDeformationCL& md = MeshDeformationCL::getInstance();
    //md.Initialize(&mg);
    instat_scalar_fun_ptr sigma( 0);
    SurfaceTensionCL sf( sigma, 0);
    DROPS::LsetBndDataCL lsbnd( 6);
    //DROPS::MultiGridCL mg2( *builder);
    //DROPS::LevelsetP2CL& lset( *LevelsetP2CL::Create( *mg, lsbnd, sf, P2_FE) );
    DROPS::LevelsetP2CL& lset( *LevelsetP2CL::Create( *mg, lsbnd, sf, P.get_child("Levelset")) );
    //DROPS::LevelsetP2CL & lset( * DROPS::LevelsetP2CL::Create( mg, lsbnd, sf,P1_FE) );
    lset.CreateNumbering( mg->GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    LSInit( *mg, lset.Phi, the_lset_fun, 0.);

    DROPS::PoissonP1CL<DROPS::PoissonCoeffCL> *Poisson_point = 0;
    DROPS::SUPGCL supg;
    DROPS::PoissonCoeffCL tmp = DROPS::PoissonCoeffCL( P);
    Poisson_point = new DROPS::PoissonP1CL<DROPS::PoissonCoeffCL>( *mg, tmp, *bdata, supg, P.get<int>("ALE.wavy"));
    //DROPS::Strategy<DROPS::PoissonP1CL<DROPS::PoissonCoeffCL> >(*probP1);
    auto Poisson = *Poisson_point;
    Poisson.idx.SetFE( P2_FE);
    Poisson.CreateNumbering( mg->GetLastLevel(), &Poisson.idx);
    //Poisson.SetupInstatSystem( Poisson.A, Poisson.M, Poisson.x.t);
    Poisson.x.SetIdx( &Poisson.idx);                            // tell x about numbering
    Poisson.Init( Poisson.x, Poisson.Coeff_.InitialCondition, 0.0);//give the initial value of solution




    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.Freq",0))
    {
        vtkwriter = new VTKOutCL(*mg, "DROPS data",
                                 P.get<int>("Time.NumSteps")+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<std::string>("VTK.TimeFileName"),
                                 P.get<int>("VTK.Binary"),
                                 P.get<int>("VTK.UseOnlyP1"),
                                 -1,  /* <- level */
                                 P.get<int>("VTK.ReUseTimeFile"),
                                 P.get<int>("VTK.UseDeformation"));
        vtkwriter->Register( make_VTKScalar(      lset.GetSolution(),              "Levelset") );
        vtkwriter->Register( make_VTKScalar( Poisson.GetSolution(), "ConcenT"));
        vtkwriter->Write( Poisson.x.t);
        //std::cout<<"to here"<<std::endl;
    }

    //vtkwriter->Write( Poisson.x.t);

}
