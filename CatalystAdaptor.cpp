// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause

#include <catalyst.hpp>


#include <string>
#include<AMReX_Print.H>
#include<AMReX_Geometry.H>
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
void Initialize(int argc, char* argv[])
{
  // Populate the catalyst_initialize argument based on the "initialize" protocol [1].
  // [1] https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-initialize
  conduit_cpp::Node node;

  // Using the arguments given to the driver set the filename for the catalyst
  // script and pass the rest of the arguments as arguments of the script
  // itself. To retrieve these  arguments from the script  use the `get_args()`
  // method of the paraview catalyst module [2]
  // [2] https://kitware.github.io/paraview-docs/latest/python/paraview.catalyst.html
  node["catalyst/scripts/script/filename"].set_string(argv[1]);
  for (int cc = 2; cc < argc; ++cc)
  {
    conduit_cpp::Node list_entry = node["catalyst/scripts/script/args"].append();
    list_entry.set(argv[cc]);
  }

  // For this example we hardcode the implementation name to "paraview" and
  // define the "PARAVIEW_IMPL_DIR" during compilation time (see the
  // accompanying CMakeLists.txt). We could however defined them via
  // environmental variables  see [1].
  node["catalyst_load/implementation"] = "paraview";
  //  node["catalyst_load/search_paths/paraview"] = PARAVIEW_IMPL_DIR;
  node["catalyst_load/search_paths/paraview"] = "/rds/project/rds-YVo7YUJF2mk/shared/paraview/build-v5.11.0/install/lib/catalyst";
  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    amrex::Print() << "Failed to initialize Catalyst: " << err << "\n";
  }
  else
    {
      amrex::Print() << "Initialized Catalyst: " << err << "\n";
    }    
}

//void Execute(int cycle, double time, Grid& grid, Attributes& attribs)
  void Execute(int cycle, double time, amrex::Geometry& geom, amrex::MultiFab& S)
{
  // Populate the catalyst_execute argument based on the "execute" protocol [3].
  // [3] https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-execute

  conduit_cpp::Node exec_params;

  // State: Information about the current iteration. All parameters are
  // optional for catalyst but downstream filters may need them to execute
  // correctly.

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(cycle);
  state["time"].set(time);
  state["multiblock"].set(1); //number of channels

  // Channels: Named data-sources that link the data of the simulation to the
  // analysis pipeline in other words we map the simulation datastructures to
  // the ones expected by ParaView.  In this example we use the Mesh Blueprint
  // to describe data see also bellow.

/*   // Add channels. */
/*   // We only have 1 channel here. Let's name it 'grid'. */
  auto channel = exec_params["catalyst/channels/grid"];

/*   // Since this example is using Conduit Mesh Blueprint to define the mesh, */
/*   // we set the channel's type to "mesh". */

  channel["type"].set("mesh");

  // now create the mesh.
  conduit_cpp::Node mesh = channel["data"]; //TODO: change to "amrdata"

  // set up the fields on the mesh
  conduit_cpp::Node fields = mesh["fields"];

  // cell data corresponding to MPI process id
  conduit_cpp::Node proc_id_field = fields["procid"];
  proc_id_field["association"] = "element";
  proc_id_field["topology"] = "topo";


  int NumberOfAMRLevels = 1;

/*   // populate the data node following the Mesh Blueprint [4] */
/*   // [4] https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html */

/*   // start with coordsets (of course, the sequence is not important, just make */
/*   // it easier to think in this order). */


  amrex::ParallelFor(S,
  [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {

      for (unsigned int level = 0; level < NumberOfAMRLevels; level++)
	{
	  std::string patch_name = "domain_" + std::to_string(bi); //use the box index the label
	  conduit_cpp::Node patch = mesh[patch_name];
	  // add basic state info
	  patch["state/domain_id"] = bi; 
	  patch["state/cycle"] = cycle;
	  patch["state/time"] = time;
	  patch["state/level"] = level;

	  patch["coordsets/coords/type"] = "uniform";

	  const auto problo = geom.ProbLoArray();
	  const auto probhi = geom.ProbHiArray();
	  const auto dx = geom.CellSizeArray();

	  patch["coordsets/coords/dims/i"] = probhi[0] - problo[0] + 1;
	  patch["coordsets/coords/dims/j"] = probhi[1] - problo[1] + 1;
	  patch["coordsets/coords/dims/k"] = probhi[2] - problo[2] + 1;

	  patch["coordsets/coords/spacing/dx"] = dx[0];
	  patch["coordsets/coords/spacing/dy"] = dx[1];
	  patch["coordsets/coords/spacing/dz"] = dx[2];


	  patch["coordsets/coords/origin/x"] = problo[0] + 0.5*dx[0];
	  patch["coordsets/coords/origin/y"] = problo[1] + 0.5*dx[1];
	  patch["coordsets/coords/origin/z"] = problo[2] + 0.5*dx[2];

	  // create a rectilinear topology that refs our coordset
	  patch["topologies/topo/type"] = "uniform";
	  patch["topologies/topo/coordset"] = "coords";

	  // add logical elements origin
	  patch["topologies/topo/elements/origin/i0"] = problo[0];
	  patch["topologies/topo/elements/origin/j0"] = problo[1];
	  patch["topologies/topo/elements/origin/k0"] = problo[2];

	  conduit_cpp::Node nest_set;
	  nest_set["association"] = "element";
	  nest_set["topology"] = "topo";
	  if (level > 0)
	    {
	      // int parent_id = amr.BlockId[level - 1];
	      // std::string parent_name = "windows/window_" + std::to_string(parent_id);
	      // conduit_cpp::Node parent = nest_set[parent_name];
	      // parent["domain_id"] = parent_id;
	      // parent["domain_type"] = "parent";
	      // std::array<int, 6> parentLevelIndices = amr.GetLevelIndices(level - 1);
	      // parent["origin/i"] = levelIndices[0] / 2;
	      // parent["origin/j"] = parentLevelIndices[2];
	      // parent["origin/k"] = parentLevelIndices[4];
	      // parent["dims/i"] = parentLevelIndices[1] - levelIndices[0] / 2 + 1;
	      // parent["dims/j"] = parentLevelIndices[3] - parentLevelIndices[2] + 1;
	      // ;
	      // parent["dims/k"] = parentLevelIndices[5] - parentLevelIndices[4] + 1;
	      // ;
	      // parent["ratio/i"] = 2;
	      // parent["ratio/j"] = 2;
	      // parent["ratio/k"] = 2;
	    }
	  if (level < NumberOfAMRLevels - 1)
	    {
	      // int child_id = amr.BlockId[level];
	      // std::string child_name = "windows/window_" + std::to_string(child_id);
	      // conduit_cpp::Node child = nest_set[child_name];
	      // child["domain_id"] = child_id;
	      // child["domain_type"] = "child";

	      // child["origin/i"] = levelIndices[0];
	      // child["origin/j"] = levelIndices[2];
	      // child["origin/k"] = levelIndices[4];
	      
	      // child["dims/i"] = levelIndices[1] - levelIndices[0] + 1;
	      // child["dims/j"] = levelIndices[3] - levelIndices[2] + 1;
	      // child["dims/k"] = levelIndices[5] - levelIndices[4] + 1;

	      // child["ratio/i"] = 2;
	      // child["ratio/j"] = 2;
	      // child["ratio/k"] = 2;
	    }
	  // patch["nestsets/nest"].set(nest_set);
	}
    });
/*   // Finally, add fields. */

/*   // First component of the path is the name of the field . The rest are described */
/*   // in https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#fields */
/*   // under the Material-Independent Fields section. */
/*   auto fields = mesh["fields"]; */
/*   fields["velocity/association"].set("vertex"); */
/*   fields["velocity/topology"].set("mesh"); */
/*   fields["velocity/volume_dependent"].set("false"); */

/*   // velocity is stored in non-interlaced form (unlike points). */
/*   fields["velocity/values/x"].set_external( */
/*     attribs.GetVelocityArray(), grid.GetNumberOfPoints(), /\*offset=*\/0); */
/*   fields["velocity/values/y"].set_external(attribs.GetVelocityArray(), grid.GetNumberOfPoints(), */
/*     /\*offset=*\/grid.GetNumberOfPoints() * sizeof(double)); */
/*   fields["velocity/values/z"].set_external(attribs.GetVelocityArray(), grid.GetNumberOfPoints(), */
/*     /\*offset=*\/grid.GetNumberOfPoints() * sizeof(double) * 2); */

/*   // pressure is cell-data. */
/*   fields["pressure/association"].set("element"); */
/*   fields["pressure/topology"].set("mesh"); */
/*   fields["pressure/volume_dependent"].set("false"); */
/*   fields["pressure/values"].set_external(attribs.GetPressureArray(), grid.GetNumberOfCells()); */

  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if (err != catalyst_status_ok)
  {
    amrex::Print() << "Failed to execute Catalyst: " << err << "\n";
  }
  else
    {
      amrex::Print() << "Running Catalyst at timestep: " << cycle << "\n";
    }

}

// Although no arguments are passed for catalyst_finalize  it is required in
// order to release any resources the ParaViewCatalyst implementation has
// allocated.
void Finalize()
{
  conduit_cpp::Node node;
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    amrex::Print() << "Failed to finalize Catalyst: " << err << "\n";
  }
  else
    {
      amrex::Print() << "Finalized Catalyst: " << err << "\n";
    }
}
}


