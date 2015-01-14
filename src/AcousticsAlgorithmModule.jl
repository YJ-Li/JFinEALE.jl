module AcousticsAlgorithmModule

using JFinEALE.JFFoundationModule
using JFinEALE.FESetModule
using JFinEALE.NodalFieldModule
using JFinEALE.ForceIntensityModule
using JFinEALE.PropertyAcousticFluidModule
using JFinEALE.MaterialAcousticFluidModule
using JFinEALE.FEMMBaseModule
using JFinEALE.FEMMAcousticsModule


function steadystate(modeldata::ModelDataDictionary)
    # Steady-state acoustics solver.
    #
    #
    # Arguments: 
    #
    # modeldata = dictionary with fields as follows.
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
    # dP/dn=-rho*a ... The normal derivative of the pressure in terms of the acceleration
    # modeldata.boundary_conditions.flux = array of dictionaries,
    #           each piece of surface with flux boundary condition gets one
    #           element of the array with a dictionary with the attributes
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
    #     P    =the nodal field that is the computed pressure 
    
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
    
    # Construct the acoustic pressure field
    P = NodalField(name ="P",data =zeros(JFCplxFlt,size(fens.xyz,1),1))

    # The following properties are the same for all regions and all boundaries (ABC)
    bulk_modulus=get(()->error("Must get bulk modulus!"), modeldata, "bulk_modulus");
    mass_density=get(()->error("Must get mass density!"), modeldata, "mass_density");
    omega=get(()->error("Must get angular frequency!"), modeldata, "omega");
            
    
    # Apply the essential boundary conditions on the acoustic pressure field
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
            pressure=get(essential[j], "pressure", nothing);
            Pfixed=zeros(JFFlt,length(fenids));
            if (pressure!= nothing)
                if (typeof(pressure)==Function)
                    for k=1:length(fenids)
                        Pfixed[k] = pressure(geom.values[fenids[k],:])[1];
                    end
                else
                    fill! (Pfixed,pressure);
                end
            end
            setebc!(P,fenids[:],trues(length(fenids)),fenids[:]*0+1,Pfixed);
            applyebc!(P);
        end
    end
    
    # Number the equations
    numberdofs!(P)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<
    
    # Initialize the heat loads vector
    F =zeros(JFFlt,P.nfreedofs);
    
    # Create the acoustics finite element models for the regions (if not supplied as input)
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        if (femm==nothing) #need to construct the FEMM
            # Construct the property  and material  objects
            prop=PropertyAcousticFluid(bulk_modulus,mass_density)
            mater=MaterialAcousticFluid(prop);
            # This is the model object for the current region
            fes=get(()->error("Must get finite elements!"), region[i], "fes");
            integration_rule=get(()->error("Must get integration rule!"), region[i], "integration_rule");
            femm=FEMMAcoustics(FEMMBase(fes,integration_rule), mater);
            setindex!(region[i], femm, "femm");
        end
    end
    setindex!(modeldata, region, "region"); # put the generated data back
    
    # Create the finite element models for the boundary conditions
    # ABC boundary conditions:
    ABC = get(boundary_conditions, "ABC", nothing);
    if (ABC != nothing)
        # Construct the property  and material  objects
        prop=PropertyAcousticFluid(bulk_modulus,mass_density)
        mater=MaterialAcousticFluid(prop);
        for i=1:length(ABC)
            # This is the model object for the ABC surface
            fes=get(()->error("Must get finite elements!"), ABC[i], "fes");
            integration_rule=get(()->error("Must get integration rule!"), ABC[i], "integration_rule");
            femm=FEMMAcoustics(FEMMBase(fes,integration_rule),mater);
            setindex!(ABC[i], femm, "femm"); # save the created FEMM
        end
    end
    setindex!(modeldata, ABC, "ABC"); # put the generated data back
    
    #  Flux boundary condition:
    flux = get(boundary_conditions, "flux", nothing);
    if (flux != nothing)
        for i=1:length(flux)
            # This is the model object for the Neumann BC surface
            fes=get(()->error("Must get finite elements!"), flux[i], "fes");
            integration_rule=get(()->error("Must get integration rule!"), flux[i], "integration_rule");
            femm= FEMMBase(fes,integration_rule);
            setindex!(flux[i], femm, "femm"); # save the created FEMM
        end
    end
    setindex!(modeldata, flux, "flux"); # put the generated data back
    
    
    # Construct the system acoustic mass and stiffness matrix
    # and the absorbing boundary condition (ABC) matrix
    C=  spzeros(P.nfreedofs,P.nfreedofs); # (all zeros, for the moment)
    S=  spzeros(P.nfreedofs,P.nfreedofs); # (all zeros, for the moment)
    D=  spzeros(P.nfreedofs,P.nfreedofs); # (all zeros, for the moment)
    region=get(()->error("Must get region!"), modeldata, "region")
    for i=1:length(region)
        femm=get(region[i], "femm", nothing);
        # Add up all the acoustic mass matrices for all the regions        
        C = C + acousticmass(femm, geom, P);
        S = S + acousticstiffness(femm, geom, P);
        # Loads due to the essential boundary conditions on the pressure field
        boundary_conditions = get(modeldata, "boundary_conditions", nothing);
        essential = get(boundary_conditions, "essential", nothing);
        if (essential!= nothing)
            F = F + nzebcloadsacousticmass(femm, geom, P);
            F = F + nzebcloadsacousticstiffness(femm, geom, P);
        end
    end

    # Compute the ABC matrices    
    ABC = get(boundary_conditions, "ABC", nothing);
    if (ABC != nothing)
        for j=1:length(ABC)
            femm=get(ABC[j], "femm", nothing);
            D = D + acousticABC(femm, geom, P);
        end
    end
    
    # Process the flux boundary condition:
    # dP/dn=-rho*a ... The normal derivative of the pressure in terms of the acceleration
    flux = get(boundary_conditions, "flux", nothing);
    if (flux != nothing)
        for j=1:length(flux)
            normal_flux= get(()->error("Must get normal flux value!"), flux[j], "normal_flux");
            fi = ForceIntensity(normal_flux);
            femm=get(()->error("Must get femm for the flux boundary condition!"), flux[j], "femm");
            F = F + distribloads(femm, geom, P, fi, 2);
        end 
    end

    # Solve for the pressures
    K=lufact((-omega^2*S +omega*1.0im*D + C));
    vP=  K\F;
    scattersysvec! (P,vP[:])
    
    # Update the model data
    setindex!(modeldata, geom, "geom");
    setindex!(modeldata, P, "P");
    return modeldata            # ... And return the updated model data
end

end
