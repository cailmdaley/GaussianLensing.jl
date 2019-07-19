
function read_Cℓs(filename="CAMB_fiducial_cosmo_scalCls.dat", dir="datafiles")
	ℓmin = readdlm(joinpath(dir, filename))[:, 1][1] 
	Dℓs = [zeros(Int(ℓmin)); 
	        readdlm("datafiles/CAMB_fiducial_cosmo_scalCls.dat")[:, 2]]
	ℓs   = 0:length(Dℓs)-1; 
	Cℓs = @. 2π/(ℓs * (ℓs + 1)) * Dℓs; Cℓs[1] = 0
	return Cℓs
end

#https://camb.readthedocs.io/en/latest/CAMBdemo.html
