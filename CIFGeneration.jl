using DataFrames
using CSV
using Printf
using Statistics
using LinearAlgebra
using Xtals
using LightGraphs
using GraphPlot

#if length(ARGS) != 1
#    error("Pass xyz file as argument")
#end
#filelist = readdir("C:\\Users\\kgeri\\Documents\\GitHub\\MolecularSynthesis\\examples")
#filelist = readdir("C:\\Users\\kgeri\\Documents\\GradSchool\\Research\\examples")
filelist = readdir("C:\\Users\\zhang\\Desktop")
for filename in filelist
   
    if endswith(filename, ".mol") == true
    #xyz_filename = ARGS[1]
        mol_filename = filename

        atoms, bonds, bondtypes = read_mol(joinpath(homedir(), "Desktop" ,mol_filename))
        
        r = maximum([distance(atoms, i, j) for i = 1:atoms.n, j = 1:atoms.n])
        box_length = 3*r
        box = Box(box_length,box_length,box_length)
        atoms=Frac(atoms, box)

        crystal = Crystal("NewTest.cif", box, atoms, Charges{Frac}(0), bonds, Xtals.SymmetryInfo())

        ##################################################################################################
        ### Function to identify carboxyl. When species is carbon, checks if it has two oxygen neighbors. 
        ### If it does, returns oxygen id's, anchor id, and oxygen's hydrogen id

        function identify_carboxyl(crystal::Crystal, a::Int64)
        
            species = crystal.atoms.species[a]    #Identify species of selected atom
            nbs = neighbors(crystal.bonds, a)     #Identify neighbors of atom
            
            if species == :C                      #Identifying carboxylate starts with carbon
                oxygen_counter = 0
                oxygen_id = Int64[]
                X_id = 0
                H_id = 0
                for nb in nbs                     #Iterating through neighbors via identifying number
                    species_nb = crystal.atoms.species[nb] #Get species of neighbor
                    if species_nb == :O
                        
                        next_nbs = neighbors(crystal.bonds, nb)  #Get neighbors of oxygen (finding hydrogen)
                        
                        for next_nb in next_nbs
                            species_next_nb = crystal.atoms.species[next_nb]  #Get species of oxygen's neighbors
                            if species_next_nb == :H     #Identify hydrogen neighbor
                                H_id = next_nb
                            else
                            end
                        end
                        
                        oxygen_counter += 1
                        push!(oxygen_id, nb)         #Enter oxygen id into array
                    else
                        X_id = nb                     #If neighbor of carboxylate carbon isn't oxygen, must be anchor
                    end
                end
                
                if oxygen_counter == 2   #Assuming carboxylate carbon is only carbon with two oxygens bonded
                
                    return true, oxygen_id, X_id, H_id
                else 
                    return false, [], 0, 0
                end
            else
                return false, [], 0, 0
            end

        end

        ##################################################################################################
        keep = [true for i = 1:crystal.atoms.n]

    X_species = Symbol[]
    X_ids = Int64[]

    for a = 1:crystal.atoms.n
       
   

        is_carboxyl, oxygen_ids, X_id, H_id = identify_carboxyl(crystal, a)
        
        if is_carboxyl == true 
        
            push!(X_species, crystal.atoms.species[X_id])
            push!(X_ids, X_id)
            crystal.atoms.species[X_id] = :X
            @assert length(oxygen_ids) == 2
            keep[oxygen_ids] .= false
            keep[a] = false
            keep[H_id] = false
            
        end


    end

     ##################################################################################################
     tobacco_crystal = crystal.atoms[keep]
     write_cif(crystal, "tobacco_crystal.cif")
     tobacco_crystal = crystal[BitArray(keep)]

     ##################################################################################################
     function center_mass(crystal::Crystal)
        xf_cm = sum(crystal.atoms.coords.xf, dims = 2)
        frctn_xf_cm = xf_cm ./ crystal.atoms.n
        new_coords = crystal.atoms.coords.xf .- frctn_xf_cm
        crystal.atoms.coords.xf .= mod.(crystal.atoms.coords.xf, 1)
        crystal.atoms.coords.xf .= new_coords
            
        return crystal
    end

    ##################################################################################################
    tobacco_crystal = center_mass(tobacco_crystal)

    ##################################################################################################
    function write_cif_Kai(crystal::Crystal, filename::AbstractString, X_species::Array{Symbol,1}, X_ids::Array{Int64,1};
         fractional_coords::Bool=true, number_atoms::Bool=true)
     if has_charges(crystal)
         if crystal.atoms.n != crystal.charges.n
           error("write_cif assumes equal numbers of Charges and Atoms (or zero charges)")
         end
         if ! isapprox(crystal.charges.coords, crystal.atoms.coords)
           error("write_cif needs coords of atoms and charges to correspond.")
         end
     end

     # TODO is this labeling necessary for the bonds, arthur?
     # create dictionary for tracking label numbers
     label_numbers = Dict{Symbol, Int}()
     for atom in crystal.atoms.species
       if !haskey(label_numbers, atom)
           label_numbers[atom] = 1
       end
     end

     # append ".cif" to filename if it doesn't already have the extension
     if ! occursin(".cif", filename)
       filename *= ".cif"
     end
     cif_file = open(filename, "w")
     # first line should be data_xtalname_PM
     if crystal.name == ""
         @printf(cif_file, "data_PM\n")
     else
       # don't include file extension!
       @printf(cif_file, "data_%s_PM\n", split(crystal.name, ".")[1])
     end

     @printf(cif_file, "_symmetry_space_group_name_H-M\t'%s'\n", crystal.symmetry.space_group)

     @printf(cif_file, "_cell_length_a\t%f\n", crystal.box.a)
     @printf(cif_file, "_cell_length_b\t%f\n", crystal.box.b)
     @printf(cif_file, "_cell_length_c\t%f\n", crystal.box.c)

     @printf(cif_file, "_cell_angle_alpha\t%f\n", crystal.box.α * 180.0 / pi)
     @printf(cif_file, "_cell_angle_beta\t%f\n", crystal.box.β * 180.0 / pi)
     @printf(cif_file, "_cell_angle_gamma\t%f\n", crystal.box.γ * 180.0 / pi)

     @printf(cif_file, "_symmetry_Int_Tables_number 1\n\n")
     @printf(cif_file, "loop_\n_symmetry_equiv_pos_as_xyz\n")
     for i in 1:size(crystal.symmetry.operations, 2)
            @printf(cif_file, "'%s,%s,%s'\n", crystal.symmetry.operations[:, i]...)
     end
     @printf(cif_file, "\n")

     @printf(cif_file, "loop_\n_atom_site_label\n_atom_site_type_symbol\n")
     if fractional_coords
       @printf(cif_file, "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
     else
       @printf(cif_file, "_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\n")
     end
     high_precision_charges = false # if, for neutrality, need high-precision charges
     if has_charges(crystal)
       @printf(cif_file, "_atom_site_charge\n")
       # if crystal will not be charge neutral to a 1e-5 tolerance when loading it
       #    into PorousMaterials.jl, then write higher-precision charges
       if abs(sum(round.(crystal.charges.q, digits=6))) > NET_CHARGE_TOL
           @info "writing high-precision charges for " * filename * ".\n"
           high_precision_charges = true
       end
     end

     idx_to_label = Array{AbstractString, 1}(undef, crystal.atoms.n)
     nb_of_X = 0
     for i = 1:crystal.atoms.n
       # print label and type symbol
       
         if crystal.atoms.species[i] == :X
           nb_of_X += 1

           element = X_species[nb_of_X]
           
           @printf(cif_file, "%s\t%s\t", string(crystal.atoms.species[i]) *
                string(i),
               element)
         else
           @printf(cif_file, "%s\t%s\t", string(crystal.atoms.species[i]) *
               string(i),
               crystal.atoms.species[i])
         end
       
        # store label for this atom idx
         idx_to_label[i] = string(crystal.atoms.species[i]) *
                   string(i)
         # increment label
         label_numbers[crystal.atoms.species[i]] += 1
         if fractional_coords
           @printf(cif_file, "%f\t%f\t%f", crystal.atoms.coords.xf[:, i]...)
         else
           @printf(cif_file, "%f\t%f\t%f", (crystal.box.f_to_c * crystal.atoms.coords.xf[:, i])...)
         end
         if has_charges(crystal)
             if high_precision_charges
               @printf(cif_file, "\t%.10f\n", crystal.charges.q[i])
             else
               @printf(cif_file, "\t%f\n", crystal.charges.q[i])
             end
         else
         @printf(cif_file, "\n")
        end
     end

     # only print bond information if it is in the crystal
     if ne(crystal.bonds) > 0
       if ! number_atoms
            error("must label atoms with numbers to write bond information.\n")
       end
       # print column names for bond information
       @printf(cif_file, "\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n_ccdc_geom_bond_type\n")

       for edge in collect(edges(crystal.bonds))
           dxf = crystal.atoms.coords.xf[:, edge.src] - crystal.atoms.coords.xf[:, edge.dst]
           nearest_image!(dxf)
           @printf(cif_file, "%s\t%s\t%0.5f\t%s\n", idx_to_label[edge.src], idx_to_label[edge.dst],
                   norm(dxf), ". A")
       end
     end
     close(cif_file)
    end

        ##################################################################################################
        infer_bonds!(tobacco_crystal, false)
        write_cif_Kai(tobacco_crystal, "tobacco_crystal_2.cif", X_species, X_ids)


        ##################################################################################################
        filename = replace(filename, "mol" => "cif")
        write_cif_Kai(tobacco_crystal, filename, X_species, X_ids)

    else
    end
end