// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include<dune/grid/onedgrid.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include"driver.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

 // Pročitaj ulaznu datoteku
    Dune::ParameterTree input_data;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("wave.ini",input_data);
    ptreeparser.readOptions(argc,argv,input_data);

    int refinement    = input_data.get<int>("refinement");             //broj profinjenja
    double dt         = input_data.get<double>("dt");                  //vremenski korak
    double c          = input_data.get<double>("c");                   // brzina vala
    double T_end      = input_data.get<double>("T");                   // konačno vrijeme
    std::string name  = input_data.get<std::string>("filename");

    //2D mreža
  /*  constexpr int dim = 2;  // dimenzija mreže
    using GridType = Dune::UGGrid<dim>;
    std::unique_ptr<GridType> pgrid{Dune::GmshReader<GridType>::read("domena2.msh", true, false)};
    pgrid->globalRefine(refinement);
    using GV =  GridType::LeafGridView;
    const GV &gv = pgrid->leafGridView();*/

    //1D mreža
    constexpr int dim = 1;  // dimenzija mreže
    using GridType = Dune::OneDGrid;
    typedef Dune::OneDGrid::ctype DF;
    int N = 50;
    DF a = 0, b = 4*M_PI;
    GridType grid (N, a, b);
    if(refinement > 0)
         grid.globalRefine(refinement);
    auto gv = grid.leafGridView();


    driver(gv, dt, c, T_end, name);

    return 0;
}
