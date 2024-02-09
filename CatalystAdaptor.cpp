// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause



#include"CatalystAdaptor.h"



#include<AMReX.H>
#include<AMReX_Print.H>
#include<AMReX_Geometry.H>
#include<AMReX_MultiFab.H>
#include<AMReX_Conduit_Blueprint.H>


#include <catalyst.hpp>


#include <string>

/**
 * The namespace hold wrappers for the three main functions of the catalyst API
 * - catalyst_initialize
 * - catalyst_execute
 * - catalyst_finalize
 * Although not required it often helps with regards to complexity to collect
 * catalyst calls under a class /namespace.
 */

using namespace amrex;

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
  void Initialize(std::string filename, std::string catalyst_options)
{
  // Populate the catalyst_initialize argument based on the "initialize" protocol [1].
  // [1] https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-initialize
  conduit_cpp::Node node;

  // Using the arguments given to the driver set the filename for the catalyst
  // script and pass the rest of the arguments as arguments of the script
  // itself. To retrieve these  arguments from the script  use the `get_args()`
  // method of the paraview catalyst module [2]
  // [2] https://kitware.github.io/paraview-docs/latest/python/paraview.catalyst.html
  //  node["catalyst/scripts/script/filename"].set_string(argv[1]);

  node["catalyst/scripts/script/filename"] = filename; 


  conduit_cpp::Node list_entry = node["catalyst/scripts/script/args"].append();
  list_entry.set(catalyst_options);
  amrex::Print() << "Read catalyst option: " << catalyst_options << "\n"; 


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

				       

  void Execute(int cycle, double time, int iteration, int output_levs, const amrex::Vector<amrex::Geometry>& geoms, const amrex::Vector<amrex::IntVect>& ref_ratios, const amrex::Vector<const amrex::MultiFab*>& mfs)
{
  // Populate the catalyst_execute argument based on the "execute" protocol [3].
  // [3] https://docs.paraview.org/en/latest/Catalyst/blueprints.html#protocol-execute

  conduit_cpp::Node exec_params;
  //  conduit::Node mesh; //NB: this is not the same as conduit_cpp which is provided by Catalyst!!

  // State: Information about the current iteration. All parameters are
  // optional for catalyst but downstream filters may need them to execute
  // correctly.

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(cycle);
  state["time"].set(time);
  //  state["multiblock"].set(1); //number of channels

  // Channels: Named data-sources that link the data of the simulation to the
  // analysis pipeline in other words we map the simulation datastructures to
  // the ones expected by ParaView.  In this example we use the Mesh Blueprint
  // to describe data see also bellow.

/*   // Add channels. */
/*   // We only have 1 channel here. Let's name it 'input'. */
  auto channel = exec_params["catalyst/channels/input"];

/*   // Since this example is using Conduit Mesh Blueprint to define the mesh, */
/*   // we set the channel's type to "mesh". */

  channel["type"].set("mesh"); //TODO: change to "amrmesh"

  // now create the mesh.
  //when using multimesh, the additional meshes must be named e.g grid/domain for this to be
  //a valid conduit node.

  auto mesh = channel["data"];  


  mesh["coordsets/coords/type"].set_string("uniform");

  mesh["topologies/topo/type"] = "uniform";
  // reference the coordinate set by name
  mesh["topologies/topo/coordset"] = "coords";

  
  amrex::Vector<int> level_steps;  
  
  for(int i = 0; i < output_levs; i++)
    level_steps.push_back(iteration);

  amrex::Vector<std::string> varnames; 
  varnames.push_back("phi0");



  conduit::Node bp_mesh;
  MultiLevelToParaviewConduitBlueprint( output_levs,      //how many levels? 
  					mfs,         //MultiFab object
  					varnames,         //name of fields passed to conduit
  					geoms,       //Simulation geometry 
  					time,             //
  					level_steps,      //these are the iterations
  					ref_ratios,         //ref ratio
  					mesh);     //conduit node object

  //  amrex::WriteBlueprintFiles(bp_mesh, "/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/AMReX/wave/insitu-viz/conduit_example_", cycle);


  exec_params.print();
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



void 
TestMultiLevelToParaviewConduitBlueprint (int n_levels,
				       const amrex::Vector<const amrex::MultiFab*>& mfs,
				       const amrex::Vector<std::string>& varnames,
				       const amrex::Vector<amrex::Geometry>& geoms,
				       amrex::Real time_value,
				       const amrex::Vector<int>& level_steps,
				       const amrex::Vector<amrex::IntVect>& ref_ratios,
				       conduit_cpp::Node &mesh)
{
  const amrex::Geometry &geom = geoms[0];
  const amrex::MultiFab &mf = *mfs[0];
  const amrex::FArrayBox& fab = mf[0];

  int ngrow = mf.nGrow(); //number of ghost cells


  FabToBlueprintTopology(geom, fab, ngrow, mesh);

  // int numPerDim = 128;
  // // create the coordinate set
  // mesh["coordsets/coords/type"] = "uniform";
  // mesh["coordsets/coords/dims/i"] = numPerDim;
  // mesh["coordsets/coords/dims/j"] = numPerDim;
  // mesh["coordsets/coords/dims/k"] = numPerDim;

  // // // add origin and spacing to the coordset (optional)
  // mesh["coordsets/coords/origin/x"] = -10.0;
  // mesh["coordsets/coords/origin/y"] = -10.0;
  // mesh["coordsets/coords/origin/z"] = -10.0;
  // double distancePerStep = 20.0/(numPerDim-1);
  // mesh["coordsets/coords/spacing/dx"] = 0.0078125; //distancePerStep;
  // mesh["coordsets/coords/spacing/dy"] = 0.0078125; //distancePerStep;
  // mesh["coordsets/coords/spacing/dz"] = 0.0078125; //distancePerStep;

  // //  add the topology
  // //  this case is simple b/c it's implicitly derived from the coordinate set
  // mesh["topologies/topo/type"] = "uniform";
  // //  reference the coordinate set by name
  // mesh["topologies/topo/coordset"] = "coords";

  // now extract the FAB details
  const amrex::Box &fab_box = fab.box();


  int i_min = fab_box.smallEnd(0);
  int i_max = fab_box.bigEnd(0);

  int j_min = fab_box.smallEnd(1);
  int j_max = fab_box.bigEnd(1);

  int k_min = fab_box.smallEnd(2);
  int k_max = fab_box.bigEnd(2);

  amrex::Print() << "x : " << i_min << " " << i_max << "\n";
  amrex::Print() << "y : " << j_min << " " << j_max << "\n";
  amrex::Print() << "z : " << k_min << " " << k_max << "\n";

  //  int numPerDim = 128-1; 
  //  int numVertices = (numPerDim-1)*(numPerDim-1)*(numPerDim-1); //(i_max-i_min)*(j_max-j_min)*(k_max-k_min);



  //  amrex::Print()<< "n_vert: " << numVertices << "\n";

  // float *vals = new float[numVertices]; //if association = "element" then each dim is reduced by 1
  // for (int i = 0 ; i < numVertices ; i++)
  //   vals[i] = ( (i%2)==0 ? 0.0 : 1.0);

  //  amrex::Print() << fab_box.numPts() << "\n";

  // create a vertex associated field 
  mesh["fields/phi0/association"] = "element";
  mesh["fields/phi0/topology"] = "topo";
  //  mesh["fields/phi0/values"].set_external(vals, numVertices);


  amrex::Real *data_ptr = const_cast<amrex::Real*>(fab.dataPtr(0)); //0=only use phi0
  mesh["fields/phi0/values"].set_external(data_ptr,fab.box().numPts());




};


void 
MultiLevelToParaviewConduitBlueprint (int n_levels,
				       const amrex::Vector<const amrex::MultiFab*>& mfs,
				       const amrex::Vector<std::string>& varnames,
				       const amrex::Vector<amrex::Geometry>& geoms,
				       amrex::Real time_value,
				       const amrex::Vector<int>& level_steps,
				       const amrex::Vector<amrex::IntVect>& ref_ratios,
				       conduit_cpp::Node &mesh)
{
    BL_PROFILE("MultiLevelToBlueprint()");

    BL_ASSERT(n_levels <= mfs.size());
    BL_ASSERT(n_levels <= geoms.size());
    BL_ASSERT(n_levels <= ref_ratio.size()+1);
    BL_ASSERT(n_levels <= level_steps.size());
    BL_ASSERT(mfs[0]->nComp() == varnames.size());

    // get mpi rank and # of tasks
    int rank   = amrex::ParallelDescriptor::MyProc();
    int ntasks = amrex::ParallelDescriptor::NProcs();

    // get global domains already present in node
    long domain_offset = (long)mesh.number_of_children();
    amrex::ParallelDescriptor::ReduceLongSum(domain_offset);



    amrex::Vector<const amrex::BoxArray*> box_arrays;
    amrex::Vector<int> box_offsets;

    box_offsets.resize(n_levels);

    for(int i = 0; i < n_levels; i++)
    {

      const amrex::BoxArray &boxs = mfs[i]->boxArray();
      box_arrays.push_back(&boxs);
      if(i == 0)
      {
        box_offsets[i] = 0;
      }
      else
      {
        box_offsets[i] = box_offsets[i-1] + mfs[i]->size();
      }

    }


    conduit_cpp::Node nest_set;

    nest_set["association"] = "element";
    nest_set["topology"] = "mesh";


    int num_domains = 0;
    for(int level = 0; level < n_levels; level++)
    {
        //
        // Geometry represents the physical and logical space of an entire level.
        //
        // Multifab contains the patches or blocks for this level.
        // In Blueprint speak, Each Multifab contains several domains and
        // each fab has "components" which map to a Blueprint field.



      const amrex::Geometry &geom = geoms[level];
      const amrex::MultiFab &mf = *mfs[level];
      auto ref_ratio = ref_ratios[level];

        // ngrow tells us how many layers of ghosts
        int ngrow = mf.nGrow();

        // mfiter allows us to iterate over local patches
        for(amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            // domain_id is mfi.index + all patches on lower levels
	  // int domain_id = mfi.index() + num_domains + domain_offset;
            // const std::string& patch_name = amrex::Concatenate("domain_",
            //                                                    domain_id,
            //                                                    6);

	    std::string patch_name = "domain_" + std::to_string(level + n_levels * rank);


            const amrex::FArrayBox& fab = mf[mfi];
	    const amrex::Box &box = fab.box();

            // create coordset and topo
	    FabToBlueprintTopology(geom, fab, ngrow, mesh);


            // add the nesting relationship
            if(n_levels > 1)
            {

	      for (int i = 0; i < n_levels; i++)
		amrex::Print() << "is box_array ok: " <<  box_arrays[i]->ok() << " " << "is box ok: " << box.ok() << "\n"; 
	      
	      Nestsets(level, n_levels, fab, box_arrays, ref_ratios, box_offsets, nest_set);


	    	// if (level > 0)
	    	//   {

	    	//     std::vector<std::pair<int,Box> > isects
	    	//       = box_arrays[level-1]->intersections(amrex::coarsen(box, ref_ratio[level-1]));


		
	    	//     for (int b = 0; b < isects.size();++b){

	    	//       amrex::Print() << "currently on parent box " << b << "\n"; 
		      
	    	//       int parent_id = isects[b].first+box_offsets[level-1];
	    	//       std::string parent_name = "windows/window_" + std::to_string(parent_id);
	    	//       conduit_cpp::Node parent = nest_set[parent_name];

	    	//       amrex::Box parent_box  = amrex::refine(isects[b].second, ref_ratio[level-1]);
	    	//       amrex::Box overlap = box & parent_box; 

	    	//       parent["domain_id"] = parent_id;
	    	//       parent["domain_type"] = "parent";
	    	//       //		    std::array<int, 6> parentLevelIndices = amr.GetLevelIndices(level - 1);
	    	//       parent["origin/i"] = overlap.smallEnd()[0]-parent_box.smallEnd()[0];
	    	//       parent["origin/j"] = overlap.smallEnd()[1]-parent_box.smallEnd()[1];
	    	//       parent["origin/k"] = overlap.smallEnd()[2]-parent_box.smallEnd()[2];

	    	//       parent["dims/i"] = overlap.size()[0];
	    	//       parent["dims/j"] = overlap.size()[1];
	    	//       parent["dims/k"] = overlap.size()[2];
  	      
	    	//       parent["ratio/i"] = ref_ratios[level-1][0];
	    	//       parent["ratio/j"] = ref_ratios[level-1][1];
	    	//       parent["ratio/k"] = ref_ratios[level-1][2];
	    	//     }

	    	//   }
	    	// if (level < n_levels - 1)
	    	//   {
	    	//     std::vector<std::pair<int,Box> > isects
	    	//       = box_arrays[level+1]->intersections(amrex::refine(box, ref_ratio[level]));

	    	//     for (int b = 0; b < isects.size();++b){

	    	//       amrex::Print() << "currently on child box " << b << "\n"; 

	    	//       int child_id = isects[b].first + box_offsets[level+1];
	    	//       std::string child_name = "windows/window_" + std::to_string(child_id);
	    	//       amrex::Box child_box = amrex::coarsen(isects[b].second, ref_ratio[level]);
	    	//       amrex::Box overlap = box & child_box;

	    	//       conduit_cpp::Node child = nest_set[child_name];
	    	//       child["domain_id"] = child_id;
	    	//       child["domain_type"] = "child";

	    	//       child["origin/i"] = overlap.smallEnd()[0]-box.smallEnd()[0];
	    	//       child["origin/j"] = overlap.smallEnd()[1]-box.smallEnd()[1];
	    	//       child["origin/k"] = overlap.smallEnd()[2]-box.smallEnd()[2];
	      
	    	//       child["dims/i"] = overlap.size()[0];
	    	//       child["dims/j"] = overlap.size()[1];
	    	//       child["dims/k"] = overlap.size()[2];

	    	//       child["ratio/i"] = ref_ratios[level][0];
	    	//       child["ratio/j"] = ref_ratios[level][1];
	    	//       child["ratio/k"] = ref_ratios[level][2];
	    	//     }
	    	//   }

	      mesh["nestsets/nest"].set(nest_set);

	    }


            // add fields
	    // set up the fields on the mesh
	    conduit_cpp::Node fields = mesh["fields"];


	    // cell data corresponding to MPI process id
	    conduit_cpp::Node proc_id_field = fields[varnames[0]]; //TODO: check that varnames is correct, not ID
	    // proc_id_field["association"] = "element";
	    // proc_id_field["topology"] = "mesh";
	    // proc_id_field["volume_dependent"] = "false";

	    amrex::Real *data_ptr = const_cast<amrex::Real*>(fab.dataPtr(0)); //0=only use phi0
	    proc_id_field["values"].set_external(data_ptr,fab.box().numPts());

	    proc_id_field["association"] = "element";
	    proc_id_field["topology"] = "mesh";	    
	    

	    // make sure we are not asking for more components than exist.
	    BL_ASSERT(varnames.size() <= fab.nComp());


	    //	    amrex::Print() << fab.box().numPts() << "\n";


            // add ghost indicator if the fab has ghost cells
            // if(ngrow > 0)
            // {

            // }
        }
        num_domains += mf.size();
    }


}

void FabToBlueprintTopology(const amrex::Geometry& geom,
			    const amrex::FArrayBox& fab,
			    int ngrow, 
			    conduit_cpp::Node &res)
{
    int dims = BL_SPACEDIM;

    // get the details of the entire level from geom
    amrex::Box level_box = geom.Domain();

    int level_nx = level_box.size()[0];
    int level_ny = level_box.size()[1];
    int level_nz = dims > 2 ? level_box.size()[2] : 0;

    amrex::Real level_x_min = geom.ProbLo(0);
    amrex::Real level_x_max = geom.ProbHi(0);

    amrex::Real level_y_min = geom.ProbLo(1);
    amrex::Real level_y_max = geom.ProbHi(1);

    amrex::Real level_z_min = dims > 2 ? geom.ProbLo(2) : 0;
    amrex::Real level_z_max = dims > 2 ? geom.ProbHi(2) : 0;

    // geom.CellSize()[i] == (level_x_max - level_x_min) / float64(level_nx);
    amrex::Real level_dx = geom.CellSize()[0];

    // geom.CellSize()[j] == (level_y_max - level_y_min) / float64(level_ny);
    amrex::Real level_dy = geom.CellSize()[1];

    // geom.CellSize()[k] == (level_z_max - level_z_min) / float64(level_nz)
    amrex::Real level_dz = dims > 2 ? geom.CellSize()[2] : 0.0;

    // now extract the FAB details
    const amrex::Box &fab_box = fab.box();


    int i_min = fab_box.smallEnd(0);//+ngrow;
    int i_max = fab_box.bigEnd(0);//-ngrow;

    int j_min = fab_box.smallEnd(1);//+ngrow;
    int j_max = fab_box.bigEnd(1);//-ngrow;

    int k_min = dims > 2 ? fab_box.smallEnd(2) : 0;
    int k_max = dims > 2 ? fab_box.bigEnd(2) : 0;

    amrex::Print() << "x : " << i_min << " " << i_max << "\n";
    amrex::Print() << "y : " << j_min << " " << j_max << "\n";
    amrex::Print() << "z : " << k_min << " " << k_max << "\n";

    int nx = (i_max - i_min + 1);
    int ny = (j_max - j_min + 1);
    int nz = dims > 2 ? (k_max - k_min +1) : 1;

    amrex::Real x_min = level_x_min + level_dx * i_min;
    //float64 x_max = level_x_min + level_dx * i_max;

    amrex::Real y_min = level_y_min + level_dy * j_min;
    //float64 y_max = level_y_min + level_dy * j_max;

    amrex::Real z_min = dims > 2 ? level_z_min + level_dz * k_min : 0.0;
    //float64 z_max = dims > 2 ? level_z_min + level_dz * k_max : 0.0;

    // create uniform coordset
    // (which also holds all implicit details needed for the topology)
    res["coordsets/coords/type"] = "uniform";
    res["coordsets/coords/dims/i"] = nx+1;  // the +1 is because the mfs are cell centered so using 'element'
    res["coordsets/coords/dims/j"] = ny+1;

    res["coordsets/coords/spacing/dx"] = level_dx;
    res["coordsets/coords/spacing/dy"] = level_dy;

    res["coordsets/coords/origin/x"] = x_min;
    res["coordsets/coords/origin/y"] = y_min;

    if(dims > 2)
    {
      res["coordsets/coords/dims/k"] = nz+1;
      res["coordsets/coords/spacing/dz"] = level_dz;
      res["coordsets/coords/origin/z"] = z_min;
    }

    // create a rectilinear topology that refs our coordset
    res["topologies/mesh/type"] = "uniform";
    res["topologies/mesh/coordset"] = "coords";


    // add logical elements origin
    res["topologies/mesh/elements/origin/i0"] = i_min;
    res["topologies/mesh/elements/origin/j0"] = j_min;
    if( dims > 2)
    {
        res["topologies/mesh/elements/origin/k0"] = k_min;
    }

}

bool Nestsets(const int level,
              const int n_levels,
              const amrex::FArrayBox &fab,
              const amrex::Vector<const amrex::BoxArray*> box_arrays,
              const amrex::Vector<amrex::IntVect> &ref_ratio,
              const amrex::Vector<int> &domain_offsets,
              conduit_cpp::Node &nestset)
{

    nestset["association"] = "element";
    nestset["topology"] = "topo";

    const int dims = BL_SPACEDIM;
    const amrex::Box &box = fab.box();

    bool valid = false;
    if(level > 0)
    {
      // check for parents
      std::vector<std::pair<int,amrex::Box> > isects
        = box_arrays[level-1]->intersections(amrex::coarsen(box, ref_ratio[level-1]));

      for(int b = 0; b < isects.size(); ++b)
      {
        valid = true;
        // get parent box in terms of this level
        Box parent = amrex::refine(isects[b].second, ref_ratio[level-1]);
        Box overlap = box & parent;
        int parent_id = isects[b].first + domain_offsets[level-1];

        const std::string& w_name = amrex::Concatenate("window_",
                                                        parent_id,
                                                        4);
        conduit_cpp::Node window = nestset["windows/"+w_name];
        window["domain_id"] = parent_id;
        window["domain_type"] = "parent";
        // box coordinates are global to the level,
        // but the the window is local to this box so
        // subtract the current box origin
        window["origin/i"] = overlap.smallEnd()[0] - box.smallEnd()[0];
        window["origin/j"] = overlap.smallEnd()[1] - box.smallEnd()[1];
        if(dims == 3)
        {
            window["origin/k"] = overlap.smallEnd()[2] - box.smallEnd()[2];
        }
        window["dims/i"] = overlap.size()[0];
        window["dims/j"] = overlap.size()[1];
        if(dims == 3)
        {
            window["dims/k"] = overlap.size()[2];
        }
        window["ratio/i"] = ref_ratio[level-1][0];
        window["ratio/j"] = ref_ratio[level-1][1];
        if(dims == 3)
        {
            window["ratio/k"] = ref_ratio[level-1][2];
        }
      }
    }

    if(level < n_levels - 1)
    {
      // check for children
      std::vector<std::pair<int,amrex::Box> > isects
        = box_arrays[level+1]->intersections(amrex::refine(box, ref_ratio[level]));

      for(int b = 0; b < isects.size(); ++b)
      {
        valid = true;
        // get child box in terms of this level
	amrex::Box child = amrex::coarsen(isects[b].second, ref_ratio[level]);
        int child_id = isects[b].first + domain_offsets[level+1];
	amrex::Box overlap = box & child;

        const std::string& w_name = amrex::Concatenate("window_",
                                                        child_id,
                                                        4);

        conduit_cpp::Node window = nestset["windows/"+w_name];
        window["domain_id"] = child_id;
        window["domain_type"] = "child";
        // box coordinates are global to the level,
        // but the the window is local to this box so
        // subtract the current box origin
        window["origin/i"] = overlap.smallEnd()[0] - box.smallEnd()[0];
        window["origin/j"] = overlap.smallEnd()[1] - box.smallEnd()[1];
        if(dims == 3)
        {
            window["origin/k"] = overlap.smallEnd()[2] - box.smallEnd()[2];
        }
        window["dims/i"] = overlap.size()[0];
        window["dims/j"] = overlap.size()[1];
        if(dims == 3)
        {
            window["dims/k"] = overlap.size()[2];
        }
        window["ratio/i"] = ref_ratio[level][0];
        window["ratio/j"] = ref_ratio[level][1];
        if(dims == 3)
        {
            window["ratio/k"] = ref_ratio[level][2];
        }
      }
    }
    return valid;
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

