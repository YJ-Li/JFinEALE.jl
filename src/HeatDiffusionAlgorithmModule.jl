module HeatDiffusionAlgorithmModule

using JFFoundationModule
using FESetModule
using NodalFieldModule
using ForceIntensityModule
using PropertyHeatDiffusionModule
using MaterialHeatDiffusionModule
using FEMMBaseModule
using FEMMHeatDiffusionModule

function steadystate(modeldata::ModelDataDictionary)
# Steady-state heat conduction solver.
#
# Arguments: 
#
# modeldata = struct  with fields as follows.
#
# modeldata.fens = finite element node set (mandatory)
#
# For each region (connected piece of the domain made of a particular material),
# mandatory:
# modeldata.region= cell array of struct with the attributes, each region 
#           gets a struct with attributes
#     region.conductivity = material conductivity
#     region.Q = material internal heat generation rate
#     region.fes= finite element set that covers the region
#     region.integration_rule =integration rule
#     
# For essential boundary conditions (optional):
# modeldata.boundary_conditions.essential = cell array of struct,
#           each piece of surface with essential boundary condition gets one
#           element of the array with a struct with the attributes
#     essential.temperature=fixed (prescribed) temperature (scalar),  or
#           handle to a function with signature
#               function T =f(x)
#     essential.fes = finite element set on the boundary to which 
#                       the condition applies
#               or alternatively
#     essential.node_list = list of nodes on the boundary to which 
#                       the condition applies
#               Only one of essential.fes and essential.node_list needs a given.
#            
# For convection boundary conditions (optional):
# modeldata.boundary_conditions.convection = cell array of struct,
#           each piece of surface with convection boundary condition gets one
#           element of the array with a struct with the attributes
#     convection.ambient_temperature=ambient temperature (scalar)
#     convection.surface_transfer_coefficient  = surface heat transfer coefficient
#     convection.fes = finite element set on the boundary to which 
#                       the condition applies
#     convection.integration_rule= integration rule
#            
# For flux boundary conditions (optional):
# modeldata.boundary_conditions.flux = cell array of struct,
#           each piece of surface with flux boundary condition gets one
#           element of the array with a struct with the attributes
#     flux.normal_flux=ambient temperature (scalar)
#     flux.fes = finite element set on the boundary to which 
#                       the condition applies
#     flux.integration_rule= integration rule
#    
# Control parameters:
# modeldata.renumber = true or false flag (default is true)
# modeldata.renumbering_method = optionally choose the renumbering 
#       method  (symrcm or symamd)
#
# Output:
# modeldata= the struct on input is augmented with 
#     for each region and for each boundary condition of the convection or 
#         flux type the finite element model is added to modeldata.  For 
#         instance modeldata.region{1}.fem or, if convection BC is defined, 
#         modeldata.boundary_conditions.convection{1}.fem 
#     renumbering information, if renumbering was requested      
#     geom =the nodal field that is the geometry
#     temp =the nodal field that is the computed temperature 
    
#     Control parameters
    renumber = get(modeldata, "renumber", true); # Should we renumber?
    # Renumbering_options.node_perm =struct( [] );
    
    # # Should we renumber the nodes to minimize the cost of the solution of 
    # # the coupled linear algebraic equations?
    # if (renumber)
    #     renumbering_method = get(modeldata, "renumbering_method", "symamd"); # default choice
    #     # Run the renumbering algorithm
    #     modeldata =renumbermesh(modeldata, renumbering_method);;
    #     # Save  the renumbering  (permutation of the nodes)
    #     Renumbering_options.node_perm  =get(modeldata, "node_perm", []);
    # end
    
    # Extract the nodes
    fens=get(()->error("Must get fens!"), modeldata, "fens")
    
    # Construct the geometry field
    geom = NodalField(name ="geom",data =fens.xyz)
     
    # Construct the temperature field
    temp = NodalField(name ="temp",data =zeros(size(fens.xyz,1),1))
    
    # Apply the essential boundary conditions on the temperature field
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
            temperature=get(essential[j], "temperature", nothing);
            T_fixed=zeros(JFFlt,length(fenids));
            if (temperature!= nothing)
                if (typeof(temperature)==Function)
                    for k=1:length(fenids)
                        T_fixed[k] = temperature(geom.values[fenids[k],:])[1];
                    end
                else
                    fill! (T_fixed,temperature);
                end
            end
            setebc!(temp,fenids[:],trues(length(fenids)),fenids[:]*0+1,T_fixed);
            applyebc!(temp);
        end
    end
    
    # Number the equations
    numberdofs!(temp)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<
    
    # Initialize the heat loads vector
    F =zeros(JFFlt,temp.nfreedofs);
    
    # Create the heat diffusion finite element models for the regions (if not supplied as input)
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        if (femm==nothing) #need to construct the FEMM
            # Construct the property  and material  objects
            conductivity=get(()->error("Must get conductivity!"), region[i], "conductivity");
            prop=PropertyHeatDiffusion(conductivity,0.0)
            mater=MaterialHeatDiffusion (prop);
            Rm=get(region[i], "Rm", nothing);
            # This is the model object for the current region: note that we supply 
            # integration  rule and the  material orientation matrix
            fes=get(()->error("Must get finite elements!"), region[i], "fes");
            integration_rule=get(()->error("Must get integration rule!"), region[i], "integration_rule");
            femm=FEMMHeatDiffusion(FEMMBase(fes,integration_rule,Rm),mater);
            setindex!(region[i], femm, "femm");
        end
    end
    setindex!(modeldata, region, "region"); # put the generated data back
    
    # # Create the finite element models for the boundary conditions
    # # Convection boundary conditions:
    # if (isfield(modeldata.boundary_conditions, 'convection' ))
    #     amb = 0*clone(temp,'amb'); # create the ambient temperature field
    #     for j=1:length(modeldata.boundary_conditions.convection)
    #         convection =modeldata.boundary_conditions.convection{j};
    #         # Note that we need to supply the surface heat transfer coefficient
    #         convection.fem = femm_heat_diffusion (struct ('material',[],...
    #             'fes',convection.fes,...
    #             'integration_rule',convection.integration_rule,...
    #             'surface_transfer', convection.surface_transfer_coefficient));
    #         modeldata.boundary_conditions.convection{j} =convection;    
    #     end
    #     clear convection 
    # end
    # #  Flux boundary condition:
    # if (isfield(modeldata.boundary_conditions, 'flux' ))
    #     for j=1:length(modeldata.boundary_conditions.flux)
    #         flux =modeldata.boundary_conditions.flux{j};
    #         flux.femm = femm_heat_diffusion (struct ('material',[],...
    #             'fes',flux.fes,...
    #             'integration_rule',flux.integration_rule));
    #         modeldata.boundary_conditions.flux{j}=flux;    
    #     end
    #     clear flux fi  femm
    # end
     
     
    # Construct the system conductivity matrix 
    K=  spzeros(temp.nfreedofs,temp.nfreedofs); # (all zeros, for the moment)
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        # Add up all the conductivity matrices for all the regions        
        K = K + conductivity(femm, geom, temp);
        Q=get(region[i], "Q", [0.0]);
        if (typeof(Q)==Function)
            fi = ForceIntensity(zeros(JFFlt,1,1),Q);
        else
            fi = ForceIntensity(Q);
        end
        F = F + distribloads(femm.femmbase, geom, temp, fi, 3);
        # Loads due to the essential boundary conditions on the temperature field
        boundary_conditions = get(modeldata, "boundary_conditions", nothing);
        essential = get(boundary_conditions, "essential", nothing);
        if (essential!= nothing)
            F = F + nzebcloadsconductivity(femm, geom, temp);
        end
    end
     
    # # Process the convection boundary condition
    # if (isfield(modeldata.boundary_conditions, 'convection' ))
    #     amb = 0*clone(temp,'amb'); # create the ambient temperature field
    #     for j=1:length(modeldata.boundary_conditions.convection)
    #         convection =modeldata.boundary_conditions.convection{j};
    #         # Apply the fixed ambient temperature
    #         fenids= connected_nodes(convection.fes);
    #         fixed=ones(length(fenids),1);
    #         comp=[];
    #         T_fixed =0;# default
    #         if (isfield( convection, 'ambient_temperature' ))
    #             T_fixed = convection.ambient_temperature;
    #         end
    #         val=zeros(length(fenids),1)+T_fixed;
    #         # This looks like we are setting essential boundary conditions, 
    #         # but in reality we are modifying the ambient temperature
    #         amb = set_ebc(amb, fenids, fixed, comp, val);
    #         amb = apply_ebc (amb);
    #         K = K + surface_transfer(convection.fem, sysmat_assembler_sparse, geom, temp);   
    #         F = F + surface_transfer_loads(convection.fem, sysvec_assembler, geom, temp, amb); 
    #         # Note that EBC will contribute through the surface heat transfer matrix            
    #         if (isfield(modeldata.boundary_conditions, 'essential' ))
    #             F = F + nz_ebc_loads_surface_transfer(convection.fem, sysvec_assembler, geom, temp);
    #         end
    #     end
    #     clear convection fenids fixed comp T_fixed  val 
    # end
    
    # # Process the flux boundary condition
    flux = get(boundary_conditions, "flux", nothing);
    if (flux != nothing)
        for j=1:length(flux)
            normal_flux= get(()->error("Must get normal flux value!"), flux[j], "normal_flux");
            fi = ForceIntensity(normal_flux);
            femm=get(()->error("Must get femm for the flux boundary condition!"), flux[j], "femm");
            # Note the sign  which reflects the formula (negative sign
            # in front of the integral)
            F = F - distribloads(femm, geom, temp, fi, 2);
        end 
    end

    # Solve for the temperatures
    K=cholfact(K);
    U=  K\F;
    scattersysvec! (temp,U[:])
   
    # Update the model data
    setindex!(modeldata, geom, "geom");
    setindex!(modeldata, temp, "temp");
    return modeldata            # ... And return the updated model data
end


end
