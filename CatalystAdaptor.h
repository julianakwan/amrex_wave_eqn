// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <AMReX_Conduit_Blueprint.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <catalyst.hpp>

/**
 * The namespace hold wrappers for the three main functions of the catalyst API
 * - catalyst_initialize
 * - catalyst_execute
 * - catalyst_finalize
 * Although not required it often helps with regards to complexity to collect
 * catalyst calls under a class /namespace.
 */
namespace CatalystAdaptor {

/**
 * In this example, we show how we can use Catalysts's C++
 * wrapper around conduit's C API to create Conduit nodes.
 * This is not required. A C++ adaptor can just as
 * conveniently use the Conduit C API to setup the
 * `conduit_node`. However, this example shows that one can
 * indeed use Catalyst's C++ API, if the developer so chooses.
 */

/*   // Populate the catalyst_initialize argument based on the "initialize"
 * protocol [1]. */
/*   // [1]
 * https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-initialize
 */
/*   conduit_cpp::Node node; */

//  void Initialize(int argc, char* argv[]);
void Initialize(std::string filename, std::string catalyst_options,
                std::string paraview_impl_dir);

/*   // Populate the catalyst_execute argument based on the "execute" protocol
 * [3]. */
/*   // [3]
 * https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-execute
 */

void Execute(int verbosity, int cycle, double time, int iteration,
             int output_levs, const amrex::Vector<amrex::Geometry> &geoms,
             const amrex::Vector<amrex::IntVect> &ref_ratios,
             const amrex::Vector<const amrex::MultiFab *> &mfs);


void MultiLevelToParaviewConduitBlueprint(
    int verbosity, int numLevels,
    const amrex::Vector<const amrex::MultiFab *> &mfs,
    const amrex::Vector<std::string> &varnames,
    const amrex::Vector<amrex::Geometry> &geoms, amrex::Real time_value,
    const amrex::Vector<int> &level_steps,
    const amrex::Vector<amrex::IntVect> &ref_ratios, conduit_cpp::Node &res);

void TestMultiLevelToParaviewConduitBlueprint(
    int verbosity, int numLevels,
    const amrex::Vector<const amrex::MultiFab *> &mfs,
    const amrex::Vector<std::string> &varnames,
    const amrex::Vector<amrex::Geometry> &geoms, amrex::Real time_value,
    const amrex::Vector<int> &level_steps,
    const amrex::Vector<amrex::IntVect> &ref_ratios, conduit_cpp::Node &res);

void FabToBlueprintTopology(int verbosity, const amrex::Geometry &geom,
                            const amrex::FArrayBox &fab, int ngrow,
                            conduit_cpp::Node &res);

bool Nestsets(const int level, const int numLevels, const amrex::FArrayBox &fab,
              const amrex::Vector<const amrex::BoxArray *> box_arrays,
              const amrex::Vector<amrex::IntVect> &ref_ratio,
              const amrex::Vector<int> &domain_offsets,
              conduit_cpp::Node &nestset);

// Although no arguments are passed for catalyst_finalize  it is required in
// order to release any resources the ParaViewCatalyst implementation has
// allocated.
void Finalize();
} // namespace CatalystAdaptor

#endif
