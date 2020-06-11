#pragma once
#include <dune/geometry/quadraturerules.hh>


//L2 norma
template<class GV, class DGF>
double L2norm(GV const & gv, DGF const & dgf, std::size_t intorder = 4)
{
    double integral = 0.0;
    for(auto const & elem : elements(gv))
    {
        const auto & rule = Dune::QuadratureRules<double,GV::dimension>::rule( elem.geometry().type() , intorder );
        for(auto const & qpoint : rule)
        {
            auto const xi = qpoint.position();
            auto const dx = qpoint.weight() * elem.geometry().integrationElement(xi);

            Dune::FieldVector<double,1> f;
            dgf.evaluate(elem, xi, f);

            integral += f * f * dx;
        }
    }
    return std::sqrt(integral);
}

//H1 polunorma
template<class GV, class DGG>
double H1halfnorm(GV const & gv, DGG const & dgg, std::size_t intorder = 4)
{
    double integral = 0.0;
    for(auto const & elem : elements(gv))
    {
        const int dim = GV::dimension;
        const auto & rule = Dune::QuadratureRules<double,GV::dimension>::rule( elem.geometry().type() , intorder );
        for(auto const & qpoint : rule)
        {
            auto const xi = qpoint.position();
            auto const dx = qpoint.weight() * elem.geometry().integrationElement(xi);

            Dune::FieldVector<double,dim> f;
            dgg.evaluate(elem, xi, f);

            integral += f*f*dx;
        }
    }
    return std::sqrt(integral);
}
