#ifndef OPERATOR_SIMPLE_ITER_HH
#define OPERATOR_SIMPLE_ITER_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>


template<typename DGF, typename DGG, typename BCType, typename FEM>
class WaveLOP : //
  public Dune::PDELab::NumericalJacobianApplyVolume  <WaveLOP<DGF, DGG, BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianVolume       <WaveLOP<DGF, DGG, BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<WaveLOP<DGF, DGG, BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianBoundary     <WaveLOP<DGF, DGG, BCType, FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:

  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  WaveLOP(const double& dt_, const double& c_,
        const DGF& udgf0_, const DGG& gradu0_,
        const DGF& udgf1_, const DGG& gradu1_,
        const FEM & fem_, unsigned int intorder_=3) :
    dt(dt_), c(c_), udgf0(udgf0_), gradu_dgf0(gradu0_), udgf1(udgf1_),
        gradu_dgf1(gradu1_), fem(fem_), intorder( intorder_ )
  {}

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    const int dim = EG::Geometry::mydimension;
    using Gradient = Dune::FieldVector<double,dim>;

    auto gt = eg.geometry().type();
    auto & rule = Dune::QuadratureRules<double,dim>::rule(gt, intorder);


    for (auto qpoint : rule ){
      auto const & xi = qpoint.position();
      auto & phihat = cache.evaluateFunction(xi, lfsu.finiteElement().localBasis());
      auto & gradphihat = cache.evaluateJacobian(xi,lfsu.finiteElement().localBasis());

      auto const & jac = eg.geometry().jacobianInverseTransposed(xi);
      std::vector<Gradient> gradphi(lfsu.size());
      for (int i=0; i<lfsu.size();++i)
        jac.mv(gradphihat[i][0],gradphi[i]);

      double u=0.0;
      for (size_t i=0; i<lfsu.size(); i++)
      u += x(lfsu,i)*phihat[i];

      Gradient gradu(0.0);
      for (int i=0; i<lfsu.size();++i)
        gradu.axpy(x(lfsu, i), gradphi[i]);

      auto xi_glob = eg.geometry().global(xi);
      Dune::FieldVector<double,1> u0, u1;
      Dune::FieldVector<double,dim>  gradu0, gradu1;

      udgf0.evaluate(eg.entity(),qpoint.position(),u0);
      udgf1.evaluate(eg.entity(),qpoint.position(),u1);
      gradu_dgf0.evaluate(eg.entity(),qpoint.position(),gradu0);
      gradu_dgf1.evaluate(eg.entity(),qpoint.position(),gradu1);

      auto dx = qpoint.weight()*eg.geometry().integrationElement(xi);

      for (int i=0; i<lfsu.size(); ++i){
      r.accumulate(lfsu, i, (u*phihat[i]+c*c*dt*dt/4*(gradu*gradphi[i]))*dx );
      r.accumulate(lfsu, i, (-2*u1*phihat[i]+c*c*dt*dt/2*(gradu1*gradphi[i]))*dx );
      r.accumulate(lfsu, i, (u0*phihat[i]+c*c*dt*dt/4*(gradu0*gradphi[i]))*dx );
      }
    }
  }


private:
  double const& dt;
  double const& c;
  DGF const & udgf0;
  DGF const & udgf1;
  DGG const & gradu_dgf0;
  DGG const & gradu_dgf1;
  FEM const & fem;
  unsigned int intorder;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
#endif
