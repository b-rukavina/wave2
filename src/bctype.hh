#ifndef BCTYPE_HH
#define BCTYPE_HH

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

int c = 1;

//rubni i inicijalni uvjeti za u
template<typename GV, typename RF>
class BCExtension0
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, BCExtension0<GV,RF> > {
  const GV& gv;
  double time;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
  BCExtension0 (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    auto x = e.geometry().global(xlocal);

    // primjer 1., rješenje y = sin(x)*cos(ct)
    if (time==0) y = sin(x);
    else y = 0;

    //primjer 2.
/*    if (time == 0) y = exp(-25*x*x);
    else y = 0;*/

    //primjer 3.
/*    if (time == 0) y = sin(M_PI*x[0])*sin(M_PI*x[1]);
    else y = 0; */

    //primjer 4.
/*  if (time == 0) y = exp(-25*x[0]*x[0]-25*x[1]*x[1]);
   else y = 0;*/

  return;
  }

  inline const GV& getGridView () {return gv;}
  void setTime (double t) {time = t;}
};


//inicijalni uvjeti za u_t
template<typename GV, typename RF>
class BCExtension1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, BCExtension1<GV,RF> > {
  const GV& gv;
  double time;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
  BCExtension1 (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    //auto x = e.geometry().global(xlocal);
    y = 0;
    return;
  }
  inline const GV& getGridView () {return gv;}
};

//Klasa koja određuje Dirichletove rubne uvjete
template <typename GV>
class BCType : public Dune::PDELab::DirichletConstraintsParameters
{
  double time;
  const GV& gv;
public:
BCType(const GV& gv_) : gv(gv_) {}
public:
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
  return true;
  }
  inline const GV& getGridView () {return gv;}
};



//egzaktno rješenje
template<typename GV, typename RF>
class Exact
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Exact<GV,RF> > {
  const GV& gv;
  double time;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
  Exact (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    auto x = e.geometry().global(xlocal);

    //primjer 1.
    y = sin(x)*cos(c*time);

    //primjer 3.
    //y = sin(M_PI*x[0])*sin(M_PI*x[1])*cos(M_PI*c*time);

    return;
  }

  inline const GV& getGridView () {return gv;}
  void setTime (double t) {time = t;}
};
#endif
