using ArgParse
using DISSEQT
using DISSEQT.AlignUtils
using DISSEQT.SynapseTools
using DISSEQT.BamReader
using JLD: save

function listaligned(path::AbstractString, runName::AbstractString)

    paths = readdir(path, join=true)
    names = readdir(path, join=false)

    # filter by extension
    mask = match.(r".bam$",names) .!= nothing
    paths[mask], names[mask]
end

function find_aligned(path::AbstractString, runName::AbstractString)
	paths, names = listaligned(path, runName)
	names = map(x->x[1:end-4], names) # keep only matches and remove ".bam" ending
	paths, names
end

function computecodonfrequencies_derp(samplePath::String, sampleName::String, 
                                 outFolder::String; 
                                 strands=:both, mappingQualityThreshold=30, 
                                 baseQualityThreshold=30, removeAmbiguous=true,
                                 method=:Newton, newtonRegularization=1e-6,
                                 maxIter=10000,
                                 outFormat=:JLD)
	log = IOBuffer()

	println(log, "Computing codon frequencies: ", sampleName)
	startTime = time()

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

	duration = time()-startTime
	println(log, "Done in ", duration, "s.")
	String(take!(log))
end

function computecodonfrequencies(samplePaths::Vector{String},
                                 sampleNames::Vector{String}, 
                                 outFolder::String; 
                                 log=stdout, bamDir="bam",
                                 kwargs...)

	for (sample_path, sample_name) in zip(samplePaths, sampleNames)
		println("Computing codon frequencies for sample $sample_name")
		logStr = computecodonfrequencies_derp(sample_path, sample_name, outFolder;kwargs...)
		print(log, logStr); flush(log)
		println("Finished computing codon frequencies for sample $sample_name")

		
	end
	nothing
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

    runName = "664V9AAXX"
    
    # Local logFile.
    logFile = "MutantSwarms.log"

    # Specify which strands to include in each run of the swarm computations. 
    strands = [:both, :forward, :reverse]

    both = joinpath(args["swarms"], "both")
    fwd = joinpath(args["swarms"], "forward")
    rev = joinpath(args["swarms"], "reverse")

    mkpath(args["swarms"]) 
    :forward in strands && mkpath(fwd) # subfolder for forward strand
    :reverse in strands && mkpath(rev) # subfolder for reverse strand

# --- Compute Mutant Swarms ----------------------------------------------------
    # find samples
    samplePaths, sampleNames = find_aligned(args["bams"], runName)
    
    strand2folder = Dict(:both=>both, :forward=>fwd, :reverse=>rev)
    # compute codon frequencies
    log = open(logFile, "w")
    
    for strand in strands
        outFolder = strand2folder[strand]

        println(log,"--- Computing mutant swarms for $strand strand(s) ---"); flush(log)
        computecodonfrequencies(samplePaths,
                                sampleNames,
                                outFolder,
                                bamDir=args["bams"],
                                strands=strand,
                                log=log)
    end
    close(log)
end

main()
