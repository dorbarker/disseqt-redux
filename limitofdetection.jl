using ArgParse
using DISSEQT
using DISSEQT.AnnotatedArrays
using Statistics
using DataFrames
using CSV
using JLD

function arguments()
    s = ArgParseSettings()

    @add_arg_table! s begin

        "--bams"
            help = "Directory containing .bam files"
            required = true

        "--swarms"
            help = "Directory containing mutant swarms"
            required = true

        "--metadata"
            help = "CSV metadata file"
            required = true

        "--features"
            help = "CSV gene features file"
            required = true        

        "--output"
            help = "JLD file path"
            required = true

    end
    return parse_args(s)
end


function loadswarms(swarms::String, bams::String)
    # replaces broken function from DISSEQT
    
    
    swarm_files = readdir(swarms, join=true)
 
    sample_names = map(x -> x |> basename |> splitext |> first, swarm_files)

    consensus_paths = map(x -> joinpath(bams, x * "_consensus.fasta"), sample_names)

    # loadswarm() from DISSEQT
    loaded_swarms = map((sample, swarm, consensus) -> loadswarm(sample, swarm, consensus), sample_names, swarm_files, consensus_paths)
    cat(loaded_swarms..., dims = 5) 

end

function main()

    args = arguments()

    metadata = CSV.read(args["metadata"]) 

    features = CSV.read(args["features"]; missingstrings=["","NA"])

    positions = featurepositions(features, "ORF") # coding positions

    # load swarms
    swarms, swarmsF, swarmsR = map(x -> loadswarms(x, args["bams"]), readdir(args["swarms"], join=true)) 

    coverage = dropdims(swarms[:coverage]; dims=(1,2,3))
    
    meanCoverage = mean(coverage[positions,:]; dims=1)[:]
    println(meanCoverage)
    # filter by mean read coverage at coding sites
    sampleMask = meanCoverage .>= 1000

    println(sampleMask)

    println("Removed ", count(.~sampleMask), " samples due to low coverage.")
    swarms  = swarms[:,:,:,:,sampleMask]
    swarmsF = swarmsF[:,:,:,:,sampleMask]
    swarmsR = swarmsR[:,:,:,:,sampleMask]

    metadata = metadata[sampleMask,:]


    println("Total number of samples used: ", count(sampleMask))
    
    
    lod = limitofdetection(swarms, swarmsF, swarmsR, groupID=metadata[!, :Run])
    lodDict = Dict( string(p[1])=>p[2] for p in pairs(convert(Dict,lod)) )
    
    save(args["output"], lodDict)

    nothing
end

main();

