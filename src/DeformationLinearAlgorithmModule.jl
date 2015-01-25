module DeformationLinearAlgorithmModule

using JFinEALE

function linearstatics(modeldata::ModelDataDictionary)
    # Algorithm for static linear deformation (stress) analysis.
    #
    # Arguments
    # modeldata = dictionary  with fields as follows.
    #
    # model_data.fens = finite element node set (mandatory)
    #
    # For each region (connected piece of the domain made of a particular material),
    # mandatory:
    # model_data.region= array of dictionaries, each region 
    #           needs to define
    #     fes= finite element set that covers the region
    #     integration_rule =integration rule
    #     property = material property object
    #     
    #        If the material orientation matrix is not the identity and needs
    #        to be supplied,  include the attribute
    #     Rm= constant orientation matrix or a handle  to a function to compute
    #           the orientation matrix (see the class femm_base).
    #
    #
    # For essential boundary conditions (optional):
    # model_data.boundary_conditions.essential = cell array of struct,
    #           each piece of surface with essential boundary condition gets one
    #           element of the array with a struct with the attributes
    #           as recognized by the set_ebc() method of nodal_field
    #     component=which component is affected (if not supplied, 
    #           or if supplied as empty default is all components)
    #     displacement=fixed (prescribed) displacement (scalar, typically 0.0),
    #           or a function that needs to return all components prescribed at a node 
    #     fes = finite element set on the boundary to which 
    #                       the condition applies
    #               or alternatively
    #     node_list = list of nodes on the boundary to which 
    #                       the condition applies
    #           Only one of fes and node_list needs to be given.
    #    
    # For traction boundary conditions (optional):
    # model_data.boundary_conditions.traction = cell array of struct,
    #           each piece of surface with traction boundary condition gets one
    #           element of the array with a struct with the attributes
    #     traction=traction (vector), supply a zero for component in which 
    #           the boundary condition is inactive
    #     fes = finite element set on the boundary to which 
    #                       the condition applies
    #     integration_rule= integration rule
    #
    # For body loads (optional):
    # model_data.body_load = cell array of struct,
    #          each piece of the domain can have each its own body load
    #     force  = force density vector
    #     fes = finite element set to which the load applies
    #     integration_rule= integration rule
    #
    # For multi point constraints (MPC) (optional):
    # model_data.mpc= cell array of structs, each for one MPC.
    #      node_list = list of node numbers involved in the MPC,
    #      dof_list= numbers of degrees of freedom for the nodes above,
    #      umultipliers=multipliers for the nodes above, 
    #      penfact=the penalty factor to multiply  the constraint matrix,
    #          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
    #          where m_i is the multiplier.
    #
    # 
    # Control parameters:
    # The following attributes  may be supplied as fields of the model_data struct:
    #      modelreduction = What is the model reduction scheme?  Default is 3-D.
    #      renumber = true or false flag (default is true)
    #      renumbering_method = optionally choose the renumbering
    #           method  ('symrcm' or 'symamd')
    #      factorize = should the solution method be based on the
    #           Choleski factorization? true or false flag (default is true)
    # 
    # Output
    # model_data =  structure that was given on input updated with:
    # The struct model_data on output incorporates in addition to the input the fields
    #     geom = geometry field
    #     u = computed displacement field
    #
    
    # What is the model reduction scheme?  Default is 3-D.
    modelreduction = get(modeldata, "modelreduction", DeformationModelReduction3D); 
    
    renumber = get(modeldata, "renumber", true); # Should we renumber?
    #Renumbering_options =struct( [] );
    
    # # Should we renumber the nodes to minimize the cost of the solution of
    # # the coupled linear algebraic equations?
    # if (renumber) 
    #     renumbering_method = 'symamd'; # default choice
    #     if ( isfield(model_data,'renumbering_method'))
    #         renumbering_method  =model_data.renumbering_method;;
    #     end
    #     # Run the renumbering algorithm
    #     model_data =renumber_mesh(model_data, renumbering_method);;
    #     # Save  the renumbering  (permutation of the nodes)
    #     clear  Renumbering_options; Renumbering_options.node_perm  =model_data.node_perm;
    # end
    
    factorize = get(modeldata, "factorize", true); # Should we factorize?
    
    
    # Extract the nodes
    fens=get(()->error("Must get fens!"), modeldata, "fens")
    
    # Construct the geometry field
    geom = NodalField(name ="geom",data =fens.xyz)
    
    # Construct the displacement field
    u = NodalField(name ="u",data =zeros(nfens(geom),ndofn(geom)))
    
    # Apply the essential boundary conditions on the displacement field
    boundary_conditions = get(()->error("Must get boundary conditions!"), modeldata, "boundary_conditions");
    essential = get(boundary_conditions, "essential", nothing);
    if (essential!= nothing)
        for j=1:length(essential)
            fes = get(essential[j], "fes", nothing);
            if (fes!= nothing)
                fenids= connectednodes(fes);
            else
                fenids= get(()->error("Must get node list!"), essential[j], "node_list");
            end
            component=get(essential[j], "component", 1:ndofn(u)); # which components?  Default is all
            displacement=get(essential[j], "displacement", nothing);
            u_fixed=zeros(JFFlt,length(fenids),length(component));
            if (displacement!= nothing)
                if (typeof(displacement)==Function)
                    for k=1:length(fenids)
                        u_fixed[k,:] = displacement(geom.values[fenids[k],:]);
                    end
                else
                    fill! (u_fixed,displacement); # presumably displacement is a number
                end
            end
            for k=1:length(component)
                for n=1:length(fenids)
                    setebc!(u,fenids[n],true,k,u_fixed[n,k]);
                end                
            end            
            applyebc!(u);
        end
    end
    
    
    # Number the equations
    numberdofs!(u)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<
    
    # Initialize the heat loads vector
    F =zeros(JFFlt,u.nfreedofs);
    
    # Create the finite element models for the regions (if not supplied as input)
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        if (femm==nothing) #need to construct the FEMM
            # Construct the property  and material  objects
            property=get(()->error("Must get property!"), region[i], "property");
            mater=MaterialDeformationLinear(property);
            Rm=get(region[i], "Rm", MaterialOrientation());
            # This is the model object for the current region: note that we supply 
            # integration  rule and the  material orientation matrix
            fes=get(()->error("Must get finite elements!"), region[i], "fes");
            integration_rule=get(()->error("Must get integration rule!"), region[i], "integration_rule");
            femm=FEMMDeformationLinear(FEMMBase(fes,integration_rule,Rm),mater);
            setindex!(region[i], femm, "femm");
        end
    end
    setindex!(modeldata, region, "region"); # put the generated data back

    
    # Construct the system stiffness matrix 
    K=  spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        # # Give the  FEMM a chance  to precompute  geometry-related quantities
        # region.femm = associate_geometry(region.femm,geom);
        # Add up all the conductivity matrices for all the regions        
        K = K + stiffness(modelreduction, femm, geom, u);
        # Loads due to the essential boundary conditions on the displacement field
        boundary_conditions = get(modeldata, "boundary_conditions", nothing);
        essential = get(boundary_conditions, "essential", nothing);
        if (essential!= nothing) # there was at least one EBC applied
            F = F + nzebcloadsstiffness(modelreduction, femm, geom, u);
        end
    end
    
    # # Process the body load
    # body_load = get(modeldata, "body_load", nothing);
    # if (body_load  !=nothing)
    #     for j=1:length(model_data.body_load)
    #         body_load =model_data.body_load{j};
    #         femm = femm_deformation_linear (struct ('material',[],...
    #             'fes',body_load.fes,...
    #             'integration_rule',body_load.integration_rule));
    #         fi= force_intensity(struct('magn',body_load.force));
    #         F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 3);        
    #     end
    #     clear body_load fi  femm
    # end
    
    # Process the traction boundary condition
    traction = get(boundary_conditions, "traction", nothing);
    if (traction != nothing)
        for j=1:length(traction)
            traction_vector= get(()->error("Must get traction value!"), traction[j], "traction_vector");
            if (typeof(traction_vector)==Function)
                fi = ForceIntensity(JFFlt,traction_vector);
            else
                fi = ForceIntensity(traction_vector);
            end
            femm=get(()->error("Must get femm for the traction boundary condition!"), traction[j], "femm");
            # Note the sign  which reflects the formula (negative sign
            # in front of the integral)
            F = F + distribloads(femm, geom, u, fi, 2);
        end 
    end
    
    # # Process the nodal force boundary condition
    # if (isfield(model_data.boundary_conditions, 'nodal_force' ))
    #     for j=1:length(model_data.boundary_conditions.nodal_force)
    #         nodal_force =model_data.boundary_conditions.nodal_force{j};
    #         femm = femm_deformation_linear (struct ('material',[],...
    #             'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
    #             'integration_rule',point_rule));
    #         fi= force_intensity(struct('magn',nodal_force.force));
    #         F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 0);          
    #     end
    #     clear nodal_force fi femm
    # end
    
    # # Apply multi point constraints
    # if isfield(model_data,'mpc')
    #     for i=1:length(model_data.mpc)
    #         mpc =model_data.mpc{i};
    #         dofnums=0*mpc.umultipliers;# Construct an array of the degree of freedom numbers
    #         for kx=1:length(mpc.node_list)
    #             dofnums(kx)=u.dofnums(mpc.node_list(kx),mpc.dof_list(kx));
    #         end
    #         # Now call the utility function to calculate the constraint matrix
    #         [Kmpc,Fmpc]=apply_penalty_mpc(u.nfreedofs,dofnums,mpc.umultipliers,0.0,mpc.penfact);
    #         K = K + Kmpc;
    #         F = F + Fmpc;
    #     end
    #     clear Kmpc Fmpc
    # end
    
    # Solve the system of linear algebraic equations
    if (factorize)
        # Factorize the left-hand side matrix (Choleski)
        K=cholfact(K);
    end
    U=  K\F;
    scattersysvec! (u,U[:])

    
    # Update the model data
    setindex!(modeldata, geom, "geom");
    setindex!(modeldata, u, "u");
    setindex!(modeldata, dot(F,U)/2, "work");
    return modeldata            # ... And return the updated model data
end
    

function exportdeformation(modeldata::ModelDataDictionary)
    # Export the deformation for visualization in Paraview.
    # 
    # Arguments
    # model_data= model data as produced by linearstatics()
    # model_data.postprocessing = optional struct with optional fields
    #     u_scale = deflection scale, default 1.0;
    #     camera  = camera, default is [] which means use the default 
    #           orientation of the view;
    #     draw_mesh= should the mesh be rendered?  Boolean.  Default true.
    #     cmap= colormap (default: jet)
    #     boundary_only= should we plot only the boundary?  Boolean.  Default true.
    #     edgecolor= color used for the mesh edges
    # Output
    # model_data = structure on input updated with
    # model_data.postprocessing.gv=graphic viewer used to display the data
    #    

    # Defaults
    boundary_only= false;
    ffile= "deformation"
    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing); 
    if (postprocessing!=nothing)
        boundary_only=  get(postprocessing, "boundary_only", boundary_only);
          ffile=  get(postprocessing, "file", ffile); 
    end
    
    fens=get(()->error("Must get fens!"), modeldata, "fens")
    geom = get(()->error("Must get geometry field!"), modeldata, "geom");  
    u = get(()->error("Must get displacement field!"), modeldata, "u"); 
    
    # Plot the surface for each region
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        rfile = ffile * "$i" * ".vtk";
        if boundary_only
            bfes= meshboundary(femm.femmbase.fes);
            vtkexportmesh (rfile, fens, bfes;  vectors=u.values, vectors_name ="u")
        else
            vtkexportmesh (rfile, fens, femm.femmbase.fes;  vectors=u.values, vectors_name ="u")
        end
    end
    
    return modeldata
end


function exportstress(modeldata::ModelDataDictionary)
    # Export the stress/deformation for visualization in Paraview.
    # 
    # Arguments
    # model_data= model data as produced by linearstatics()
    # model_data.postprocessing = optional struct with optional fields
    #     u_scale = deflection scale, default 1.0;
    #     camera  = camera, default is [] which means use the default 
    #           orientation of the view;
    #     draw_mesh= should the mesh be rendered?  Boolean.  Default true.
    #     cmap= colormap (default: jet)
    #     boundary_only= should we plot only the boundary?  Boolean.  Default true.
    #     edgecolor= color used for the mesh edges
    # Output
    # model_data = structure on input updated with
    # model_data.postprocessing.gv=graphic viewer used to display the data
    #    

    # Defaults
    boundary_only= false;
    ffile= "stress"
    output=:Cauchy
    component=1
    outputRm=MaterialOrientation()
    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing); 
    if (postprocessing!=nothing)
        boundary_only=  get(postprocessing, "boundary_only", boundary_only);
        ffile=  get(postprocessing, "file", ffile);
        output=  get(postprocessing, "output", output);
        component=  get(postprocessing, "component", component);
        outputRm=  get(postprocessing, "outputRm", outputRm);
    end
    
    fens=get(()->error("Must get fens!"), modeldata, "fens")
    geom = get(()->error("Must get geometry field!"), modeldata, "geom");  
    u = get(()->error("Must get displacement field!"), modeldata, "u"); 
    modelreduction = get(()->error("Must get model reduction!"), modeldata, "modelreduction"); 

    
    if (typeof(component)==Symbol)
        componentnum=stresscomponentmap(modelreduction)[component]
    else
        componentnum=component
    end
    
    # Plot the surface for each region
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        rfile = ffile * string(output) * string(component) * "$i" * ".vtk";
        fld= fieldfromintegpoints(modelreduction,femm, 
                          geom, u, output, componentnum)
        if boundary_only
            bfes= meshboundary(femm.femmbase.fes);
            vtkexportmesh (rfile, fens, bfes;
                           scalars=fld.values, scalars_name=string(output)*string(component),
                           vectors=u.values, vectors_name ="u")
        else
            vtkexportmesh (rfile, fens, femm.femmbase.fes;
                           scalars=fld.values, scalars_name=string(output)*string(component),
                           vectors=u.values, vectors_name ="u")
        end
    end
    
    return modeldata
end


end
