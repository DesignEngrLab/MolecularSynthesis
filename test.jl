using PorousMaterials
if length(ARGS) != 1
    error("Pass xyz file as argument")
end
xyz_filename = ARGS[1]
atoms = read_xyz(joinpath(homedir(), "Documents", "Grad School", "Research", xyz_filename))

print(atoms)