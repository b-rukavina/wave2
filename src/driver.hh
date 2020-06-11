#ifndef DRIVER_HH
#define DRIVER_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/common/fvector.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/newton/newton.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "bctype.hh"
#include "operator.hh"
#include "norm.hh"


template<typename GV>
void driver(const GV& gv, double dt, double c, double T_end,  std::string name)
{
  using namespace Dune::PDELab;

  const int dim = GV::dimension;
  const int k = 1;
  typedef PkLocalFiniteElementMap<GV,double,double,k>               FEM;
  typedef ConformingDirichletConstraints                            CONSTRAINTS;
  typedef ISTL::VectorBackend<>                                     VBE;
  typedef GridFunctionSpace<GV,FEM,CONSTRAINTS,VBE>                 GFS;
  typedef typename GFS::template ConstraintsContainer<double>::Type CC;
  typedef Backend::Vector<GFS,double>                               U;
  typedef DiscreteGridFunction<GFS,U>                               DGF;
  typedef DiscreteGridFunctionGradient<GFS,U>                       DGG;
  typedef WaveLOP<DGF, DGG, BCType<GV>, FEM>                        LOP;
  typedef ISTL::BCRSMatrixBackend<>                                 MBE;
  typedef GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC>  GO;


  double time = 0.0;

  FEM fem(gv);
  GFS gfs(gv,fem);
  CC cc;
  BCType<GV> bctype(gv);
  constraints(bctype, gfs, cc);
  std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;

  U u(gfs,0.0), u0(gfs,0.0), u1(gfs,0.0), ut(gfs,0.0), u2(gfs,0.0), x(gfs,0.0), e(gfs,0.0);
  DGF udgf(gfs,u), udgf0(gfs,u0), udgf1(gfs,u1), utdgf(gfs, ut), ex(gfs, x), error(gfs,e);
  DGG gradu(gfs, u), gradu0(gfs,u0), gradu1(gfs,u1), gradu2(gfs,u2);
  LOP lop(dt, c, udgf0, gradu0, udgf1, gradu1 , fem);
  MBE mbe(5);
  GO go(gfs,cc,gfs,cc,lop,mbe);

  typedef BCExtension0<GV,double>                BCE0;
  typedef BCExtension1<GV,double>                BCE1;
//  typedef ISTLBackend_SEQ_SuperLU                LS;
  typedef ISTLBackend_SEQ_BCGS_SSOR              LS;
  typedef StationaryLinearProblemSolver<GO,LS,U> PDESOLVER;
//  typedef Newton<GO, LS, U>                     PDESOLVER;


  BCE0 U0(gv);
  U0.setTime(time);
  BCE1 U1(gv);

  //inicijalni uvjet za u
  interpolate(U0, gfs, u0);
  //inicijalni uvjet za u_t
  interpolate(U1, gfs, u1);

  Dune::FieldVector<double, 10000> Energy(0.0);
  Energy[1] = pow(L2norm(gv, udgf1), 2);
  // u u1 se zapravo nalazi u0_t, diferencijom unaprijed dobivamo aproksimaciju u1
  u1 *= dt;
  u1 += u0;

  //u^(1/2)
  u2 = u0;
  u2 += u1;
  u2 *= 0.5;
  Energy[1] += c*c*pow(H1halfnorm(gv, gradu2), 2);



  LS ls(5000,true);
  PDESOLVER pdesolver(go,ls,u,1e-10);

/*  PDESOLVER pdesolver(go, u, ls) ;
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setMaxIterations(25);
  pdesolver.setLineSearchMaxIterations(10);*/

  using VTKW = Dune::SubsamplingVTKWriter<GV>;
  VTKW vtkwriter(gv,Dune::RefinementIntervals{2});
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf,"u"));
  vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(ex,"exact"));
  Dune::VTKSequenceWriter<GV> writer(std::make_shared<VTKW>(vtkwriter), name);

  u=u0;
  Exact<GV, double> exact(gv);
  exact.setTime(time);
  interpolate(exact, gfs, x);
  writer.write(time);
  Dune::FieldVector<double, 10000> Error(0.0);
  Error[1] = 0; //inicijalan uvjet je egzaktan
  time+=dt;

  u=u1;
  exact.setTime(time);
  interpolate(exact, gfs, x);
  e = x;
  e -= u;
  Error[2]= L2norm(gv,error); //L2 norma greške
  writer.write(time);
  time+=dt;

  int i = 2;

  while (time < T_end-1e-5){

    pdesolver.apply();

    //\partial_t u^n
    ut = u;
    ut -= u1;
    ut *= (1/dt);

    //u^(n+1/2)
    u2 = u;
    u2 += u1;
    u2 *= 0.5;

    Energy[i]= pow(L2norm(gv, utdgf), 2) + c*c*pow(H1halfnorm(gv, gradu2), 2);

    u0=u1;
    u1=u;

    exact.setTime(time);
    interpolate(exact, gfs, x);
    writer.write(time);
    e = x;
    e -= u;
    Error[i+1]= L2norm(gv,error);
    time = time + dt;
    i++;
    }


    std::cout << "ENERGIJA" << std::endl;
    for (i = 1; i<T_end/dt; i++) std::cout << Energy[i] << std::endl;
    double max_error = 0;
    std::cout << "GRESKA" << std::endl;{
      for (i = 1; i<T_end/dt; i++) std::cout << Error[i] << std::endl;
      if (Error[i]>max_error) max_error = Error[i];
    }
    std::cout << "gfs global size="  << gfs.globalSize() << std::endl;
    std::cout << "maksimalna L2 norma greške="  << max_error << std::endl;

}

#endif
