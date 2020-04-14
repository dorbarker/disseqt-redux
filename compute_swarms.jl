using Distributed
#nprocs()==1 && addprocs()
#using SynapseClient
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
 
function main()
# --- Synapse login ------------------------------------------------------------
    # syn = nothing # set syn=nothing if working locally
    # syn = SynapseClient.login()

# --- Setup --------------------------------------------------------------------
    # Set clean=true to rerun swarm inference. For clean=false, swarm inference will only be run if the log file is missing (i.e. no previous swarm inference was done).
    clean = true 

    runName = "664V9AAXX"
    projectFolder = "664V9AAXX"
    
    alignmentFolder = joinpath(projectFolder, "alignment")

    # Set uploadPath="synapseID" to upload files after running. Should point to "MyProject/Analysis/Alignment" folder. Set upload=nothing to skip uploading.
    uploadPath = alignmentFolder

    # Where to find sample .bam files. Synapse Folder ID or local folder. Synapse ID should point to "MyProject/Analysis/Alignment/Bam".
    bamPath = joinpath(projectFolder, "bams")

    # Name of the run. (Should corrsespond to a subfolder of bamPath if in Synapse.)

    # Local logFile.
    logFile = "MutantSwarms.log"

    # Specify which strands to include in each run of the swarm computations. 
    strands = [:both, :forward, :reverse]

# --- Cleanup ------------------------------------------------------------------
    clean = clean || !isfile(logFile)

    mutant_swarms_dir = joinpath(projectFolder, "mutant_swarms")
    fwd = joinpath(mutant_swarms_dir, "forward")
    rev = joinpath(mutant_swarms_dir, "reverse")
    if clean
        makecleanfolder(mutant_swarms_dir) || return
        :forward in strands && mkdir(fwd) # subfolder for forward strand
        :reverse in strands && mkdir(rev) # subfolder for reverse strand
    else
        println("Skipping Swarm Inference [clean=false]")
    end


# --- Compute Mutant Swarms ----------------------------------------------------
    # find samples
    samplePaths, sampleNames = find_aligned(bamPath, runName)
    
    strand2folder = Dict(:both=>mutant_swarms_dir,:forward=>fwd:reverse=>rev)
    if clean
        # compute codon frequencies
        log = open(logFile, "w")
        for strand in strands
            outFolder = strand2folder[strand]

            println(log,"--- Computing mutant swarms for $strand strand(s) ---"); flush(log)
            computecodonfrequencies(samplePaths, sampleNames, outFolder, bamDir=bamPath, strands=strand, log=log)
        end
        close(log)
    end
end

main()
