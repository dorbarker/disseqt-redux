import Pkg

function make_github_url(namespace_package::String)

	full_url = "https://github.com/" * namespace_package * ".git"

end

function github_install(namespace_package::String)
	
	url = make_github_url(namespace_package)
	spec = Pkg.PackageSpec(url=url)
	Pkg.add(spec)
end

function github_registry(namespace_registry::String)
	url = make_github_url(namespace_registry)
	spec = Pkg.RegistrySpec(url=url)
	Pkg.Registry.add(spec)
end

Pkg.add("JLD")
Pkg.add("Gadfly")
Pkg.add("DataFrames")
Pkg.add("DataStrutures")
Pkg.add("Libz")
Pkg.add("CSV")

github_registry("BioJulia/BioJuliaRegistry")
github_install("rasmushenningsson/SynapseClient.jl")
github_install("rasmushenningsson/DISSEQT.jl")

Pkg.build("DISSEQT")
