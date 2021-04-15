### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 2ec67f30-8313-11eb-010b-e50cfd286f4f
using Xtals, LightGraphs, Test, Printf, LinearAlgebra

# ╔═╡ 0e017870-99e5-11eb-2575-9980fa4f8e47
#print(ARGS[1])

# ╔═╡ 18856d10-99e5-11eb-019e-e7c7bec209a6


# ╔═╡ aeb964fe-906e-11eb-03f1-8f1ea6b3cbd9
#filelist=readdir("C:\\Users\\zhang\\Desktop\\data\\crystals")

# ╔═╡ 9afedfd0-9b91-11eb-0ed7-abbf5732b92a
file=ARGS[1]

# ╔═╡ cbf1f730-9baa-11eb-0160-9d181e35aa32
xtal=Crystal(file)

# ╔═╡ e7f87ee0-9baa-11eb-25d9-d9f671bda39f
name1=string(split(file,".")[1],".cssr") 

# ╔═╡ f0c43c32-9baa-11eb-0b4a-f151197107ca
write_cssr(xtal, name1)

# ╔═╡ 2ec2e6d2-9bab-11eb-312b-1901866a60d6
print("finish cif to cssr convertion")

# ╔═╡ a6862910-906c-11eb-35bd-cd348e7b8811
#begin
	
#	for file in filelist		
# 		xtal = Crystal(file)
# 		name= string(split(file,".")[1],".cssr")  
#		write_cssr(xtal, name)
# 	end
	 	
#end

# ╔═╡ 62f65460-8313-11eb-02aa-4bbd2c47a324
#xtal = Crystal("7_fer_tobacco.cif")

# ╔═╡ 63cd0550-8313-11eb-29dc-d1498d3a6dc9
#write_cssr(xtal, "ConvertionSucceed.cssr")

# ╔═╡ Cell order:
# ╠═2ec67f30-8313-11eb-010b-e50cfd286f4f
# ╠═0e017870-99e5-11eb-2575-9980fa4f8e47
# ╠═18856d10-99e5-11eb-019e-e7c7bec209a6
# ╠═aeb964fe-906e-11eb-03f1-8f1ea6b3cbd9
# ╠═9afedfd0-9b91-11eb-0ed7-abbf5732b92a
# ╠═cbf1f730-9baa-11eb-0160-9d181e35aa32
# ╠═e7f87ee0-9baa-11eb-25d9-d9f671bda39f
# ╠═f0c43c32-9baa-11eb-0b4a-f151197107ca
# ╠═2ec2e6d2-9bab-11eb-312b-1901866a60d6
# ╠═a6862910-906c-11eb-35bd-cd348e7b8811
# ╠═62f65460-8313-11eb-02aa-4bbd2c47a324
# ╠═63cd0550-8313-11eb-29dc-d1498d3a6dc9
