### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 2ec67f30-8313-11eb-010b-e50cfd286f4f
using Xtals, LightGraphs, Test, Printf, LinearAlgebra

# ╔═╡ d7773bd0-8749-11eb-13fb-4b9bce0a6457
filelist = readdir("/mnt/c/Users/zhang/source/repos/tobacco_3.0/output_cifs")

# ╔═╡ d7321c7e-8749-11eb-22db-03adffe6ea24
path_to_mol_files = joinpath("/mnt","c","Users","zhang", "source", "repos", "1", "testZeo_onUbuntu","examples")

# ╔═╡ 62f65460-8313-11eb-02aa-4bbd2c47a324
xtal = Crystal("7_fer_tobacco.cif")

# ╔═╡ 63cd0550-8313-11eb-29dc-d1498d3a6dc9
write_cssr(xtal, "ConvertionSucceed.cssr")

# ╔═╡ Cell order:
# ╠═2ec67f30-8313-11eb-010b-e50cfd286f4f
# ╠═d7773bd0-8749-11eb-13fb-4b9bce0a6457
# ╠═d7321c7e-8749-11eb-22db-03adffe6ea24
# ╠═62f65460-8313-11eb-02aa-4bbd2c47a324
# ╠═63cd0550-8313-11eb-29dc-d1498d3a6dc9
