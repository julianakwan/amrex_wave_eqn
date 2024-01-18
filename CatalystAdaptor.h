// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h


#include <catalyst.hpp>
#include<AMReX_MultiFab.H>


/**
 * The namespace hold wrappers for the three main functions of the catalyst API
 * - catalyst_initialize
 * - catalyst_execute
 * - catalyst_finalize
 * Although not required it often helps with regards to complexity to collect
 * catalyst calls under a class /namespace.
 */
namespace CatalystAdaptor
{

/**
 * In this example, we show how we can use Catalysts's C++
 * wrapper around conduit's C API to create Conduit nodes.
 * This is not required. A C++ adaptor can just as
 * conveniently use the Conduit C API to setup the
 * `conduit_node`. However, this example shows that one can
 * indeed use Catalyst's C++ API, if the developer so chooses.
 */

/*   // Populate the catalyst_initialize argument based on the "initialize" protocol [1]. */
/*   // [1] https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-initialize */
/*   conduit_cpp::Node node; */

  void Initialize(int argc, char* argv[]); 

/*   // Populate the catalyst_execute argument based on the "execute" protocol [3]. */
/*   // [3] https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-execute */

//void Execute(int cycle, double time, Grid& grid, Attributes& attribs)
  void Execute(int cycle, double time, amrex::MultiFab& S); 

// Although no arguments are passed for catalyst_finalize  it is required in
// order to release any resources the ParaViewCatalyst implementation has
// allocated.
  void Finalize(); 
}

#endif
