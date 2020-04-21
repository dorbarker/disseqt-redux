using ArgParse
@everywhere using DISSEQT
@everywhere using DISSEQT.AlignUtils
@everywhere using DISSEQT.BamReader
@everywhere import DISSEQT.BamReader.BamFile
@everywhere using JLD: save
using Distributed

function find_aligned(path::String)
  
    paths = readdir(path, join=true)
    mask = match.(r".bam$", paths) .!= nothing
    paths[mask]

end


@everywhere function computecodonfrequencies(samplePath::String, 
                                 outFolder::String, 
                                 strands::Symbol;
                                 mappingQualityThreshold=30, 
                                 baseQualityThreshold=30,
                                 removeAmbiguous=true,
                                 method=:Newton,
                                 newtonRegularization=1e-6,
                                 maxIter=10000,
                                 outFormat=:JLD)

    sampleName = samplePath |> basename |> splitext |> first

	println("Computing codon frequencies: $sampleName $strands")

	bamFile = BamFile(samplePath)
	freqs,positions,coverage = mlcodonfreqs(bamFile,strands=strands,
	                                        mappingQualityThreshold=mappingQualityThreshold,
	                                        baseQualityThreshold=baseQualityThreshold,
	                                        log=log,removeAmbiguous=removeAmbiguous,method=method,
	                                        newtonRegularization=newtonRegularization,
	                                        maxIter=maxIter)

	fileprefix = joinpath(outFolder,sampleName)
	d = Dict{String,Any}("codonFreqs"=>freqs,
	                          "positions"=>positions,
	                          "coverage"=>coverage,
	                          "segmentInfo"=>sequences(bamFile),
	                          "strands"=>string(strands),
	                          "mappingQualityThreshold" => mappingQualityThreshold,
	                          "baseQualityThreshold" => baseQualityThreshold,
	                          "removeAmbiguousBases" => removeAmbiguous,
	                          "algorithm" => "Maximum Likelihood",
	                          "solver" => string(method))
	method==:Newton && (d["newtonRegularization"] = newtonRegularization)


	#typeof(outFormat) <: AbstractArray || (outFormat = [outFormat])
	save("$fileprefix.jld", d, compress=true)

	#String(take!(log))
end

function arguments()

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--bams"
            help = "Directory containing .bam files"
            required = true
        "--swarms"
            help = "Directory to which the swarm files will be saved"
            required = true


    end
    return parse_args(s)
end
 
function main()

    args = arguments()
# --- Setup --------------------------------------------------------------------

    # Specify which strands to include in each run of the swarm computations. 
    strands = [:both, :forward, :reverse]

    both = joinpath(args["swarms"], "both")
    fwd = joinpath(args["swarms"], "forward")
    rev = joinpath(args["swarms"], "reverse")

    mkpath(args["swarms"]) 
    :both in strands && mkpath(both)
    :forward in strands && mkpath(fwd) # subfolder for forward strand
    :reverse in strands && mkpath(rev) # subfolder for reverse strand

# --- Compute Mutant Swarms ----------------------------------------------------
    # find samples
    samplePaths = find_aligned(args["bams"])
    
    strand2folder = Dict(:both=>both, :forward=>fwd, :reverse=>rev)
    # compute codon frequencies
    
    arg_length = 1:length(samplePaths) 
    for strand in strands
        outFolder = strand2folder[strand]

        println("--- Computing mutant swarms for $strand strand(s) ---")

        pmap(x -> computecodonfrequencies(x, outFolder, strand), samplePaths; distributed=true)

    end
end

main()
