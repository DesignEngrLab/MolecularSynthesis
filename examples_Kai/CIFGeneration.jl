using Xtals, LightGraphs, Test, Printf, LinearAlgebra

#if length(ARGS) != 1
#    error("Pass xyz file as argument")
#end
#filelist = readdir("C:\\Users\\kgeri\\Documents\\GitHub\\MolecularSynthesis\\examples")
#path_to_mol_files = joinpath(homedir(), "Documents", "GitHub", "MolecularSynthesis", "examples")

# those 2 lines are used on my own laptop in windows
#filelist = readdir("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\examples")
#path_to_mol_files = joinpath(homedir(), "source", "repos", "MolecularSynthesis", "examples")


# those 2 lines are used on my own laptop in ubuntu
# C:\Users\zhang\source\repos\1\testZeo_onUbuntu\examples
filelist = readdir("/mnt/c/Users/zhang/source/repos/1/testZeo_onUbuntu/examples")
path_to_mol_files = joinpath("/mnt","c","Users","zhang", "source", "repos", "1", "testZeo_onUbuntu","examples")

# those 2 lines are used on HPC
#/nfs/hpc/share/zhangho2/MolecularSynthesis/examples
#filelist = readdir("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples")
#path_to_mol_files = joinpath("/nfs","hpc","share","zhangho2", "MolecularSynthesis","examples")




#path_to_mol_files = readdir("C:\\Users\\zhang\\Desktop")
for filename in filelist
   
    if endswith(filename, ".mol") == true
    
     # convert mol to a crystal so we can have bond functionality.
      function mol_to_xtal(mol_file::String)
        # read in mole file
        atoms, bonds, bondtypes = read_mol(joinpath(path_to_mol_files, mol_file))
    
        # construct box that makes tobacco work
        max_r = maximum([distance(atoms, i, j) for i = 1:atoms.n, j = 1:atoms.n])
        L = 3 * max_r
        box = Box(L, L, L)
    
        # convert atoms to fractional coords w this box
        atoms = Frac(atoms, box)
    
        # construct a crystal
        crystal = Crystal(split(mol_file, ".mol")[1], 
                  box, 
                  atoms, 
                  Charges{Frac}(0), 
                  bonds, 
                  Xtals.SymmetryInfo()
        )
        return crystal
      end  
        
      #######################################################################################
      # must be a carbon
      # must be connected two two oxygens.
      function is_C_carboxyl(xtal::Crystal, a::Int64)
        species = xtal.atoms.species[a]    # species of selected atom
        nbs = neighbors(xtal.bonds, a)     # neighbors of atom

        if species != :C
        return false
        end

        # get # oxygen neighbors
        nb_O_neighbors = 0
        for a in nbs
          if xtal.atoms.species[a] == :O
              nb_O_neighbors += 1
          end
        end
        if nb_O_neighbors == 2
          return true
        else
          return false
        end
      end

      # get ids of neighbors of atom a that are species x
      function ids_x_neighbors(xtal::Crystal, a::Int64, x::Symbol)
        ids = Int[]
        nbs = neighbors(xtal.bonds, a)     # neighbors of atom
        for a in nbs
            if xtal.atoms.species[a] == x
                push!(ids, a)
            end
        end
        return ids
      end

      function ids_X_atoms(xtal::Crystal)
        ids = Int[]
        for a = 1:xtal.atoms.n
            if is_C_carboxyl(xtal, a)
                ids_C_nbs = ids_x_neighbors(xtal, a, :C)
                ids_N_nbs = ids_x_neighbors(xtal, a, :N)
                ids_CN_nbs = vcat(ids_C_nbs, ids_N_nbs)
                if length(ids_CN_nbs) != 1
                    error("wtf, carboxylate C is not bonded to another carbon or nitrogen.")
                end
                push!(ids, ids_CN_nbs[1])
            end
        end
        return ids
      end

      # get ids to remove: 
      #    (i) carboxylate groups
      #    (ii) hydrogens connected to the oxygens of the carboxylate groups
      function ids_carboxylate(xtal::Crystal)
        ids = Int[]
        for a = 1:xtal.atoms.n
            if is_C_carboxyl(xtal, a)
                push!(ids, a)
                ids_Os = ids_x_neighbors(xtal, a, :O)
                for id_O in ids_Os
                    push!(ids, id_O)
                    ids_Hs = ids_x_neighbors(xtal, id_O, :H)
                    if length(ids_Hs) > 0
                        push!(ids, ids_Hs[1])
                    elseif length(ids_Hs) > 1
                        error("wtf")
                    end
                end
            end
        end
        return ids
      end

      ##################################################################################
      function xtal_to_tobacco_xtal(xtal::Crystal)
        # get X atom ids, those that are C's connected to carboxylate C
        ids_X = ids_X_atoms(xtal)
        
        # find atoms to keep (those that aren't caboxylates)
        ids_keep = [i for i = 1:xtal.atoms.n if ! (i in ids_carboxylate(xtal))]
        
        # relabel X atoms 
        for id_X in ids_X
            X_species = xtal.atoms.species[id_X]
            xtal.atoms.species[id_X] = Symbol("X_" * String(X_species))
        end
        
        tobacco_xtal = xtal[ids_keep]
        return tobacco_xtal
      end

      ###############################################################################
      # MAIN
      # read in .mol file of a linker
      # strip off carboxylates
      # label C atom that was connected to the carboxylate as `:X`.

      xtal = mol_to_xtal(filename)

      tobacco_xtal = xtal_to_tobacco_xtal(xtal)

      #write_xyz(tobacco_xtal)
      #write_bond_information(tobacco_xtal)

      #write_cif(tobacco_xtal, xtal.name * "_tobacco.cif")

      ###############################################################################
      function center!(xtal::Crystal)
        # geometric center
        xf_center = sum(xtal.atoms.coords.xf, dims=2) / xtal.atoms.n
        
        xtal.atoms.coords.xf .= mod.(xtal.atoms.coords.xf .- xf_center, 1.0)
        return nothing
      end

      ##############################################################################
      center!(tobacco_xtal)
      #write_cif(tobacco_xtal, "tbc_center.cif")

      ##############################################################################
      """
      write_cif(crystal, filename; fractional_coords=true, number_atoms=true)
      Write a `crystal::Crystal` to a .cif file with `filename::AbstractString`. If `filename` does
      not include the .cif extension, it will automatically be added. the `fractional_coords` flag
      allows us to write either fractional or Cartesian coordinates.
      """
      function write_cif_Kai(crystal::Crystal, filename::AbstractString)
          if has_charges(crystal)
              if crystal.atoms.n != crystal.charges.n
                  error("write_cif assumes equal numbers of Charges and Atoms (or zero charges)")
              end
              if ! isapprox(crystal.charges.coords, crystal.atoms.coords)
                  error("write_cif needs coords of atoms and charges to correspond.")
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
          @printf(cif_file, "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
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

          
          nb_of_X = 0
          for i = 1:crystal.atoms.n
              # print label and type symbol
              
              if split(string(crystal.atoms.species[i]), "_")[1] == "X"   
                  @printf(cif_file, "%s\t%s\t", "X" * string(i),
                      split(string(crystal.atoms.species[i]), "_")[2])
              else
                  @printf(cif_file, "%s\t%s\t", string(crystal.atoms.species[i]) *
                      string(i),
                      crystal.atoms.species[i])
              end
              
              # increment label
              
              @printf(cif_file, "%f\t%f\t%f", crystal.atoms.coords.xf[:, i]...)
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
              # print column names for bond information
              @printf(cif_file, "\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n_ccdc_geom_bond_type\n")

              for edge in collect(edges(crystal.bonds))
                  dxf = crystal.atoms.coords.xf[:, edge.src] - crystal.atoms.coords.xf[:, edge.dst]
                  nearest_image!(dxf)
                  i = edge.src
                  j = edge.dst
                  i_label = string(crystal.atoms.species[i]) * string(i)
                  j_label = string(crystal.atoms.species[j]) * string(j)
                  if split(string(crystal.atoms.species[i]), "_")[1] == "X"
                      i_label = "X" * string(i)
                  elseif split(string(crystal.atoms.species[j]), "_")[1] == "X"
                      j_label = "X" * string(j)
                  end
                  @printf(cif_file, "%s\t%s\t%0.5f\t%s\n", i_label, j_label,
                          norm(dxf), ". A")
              end
          end
          close(cif_file)
      end

      #####################################################################
      write_cif_Kai(tobacco_xtal, xtal.name * "_fer_tobacco.cif")


    else
    end



end