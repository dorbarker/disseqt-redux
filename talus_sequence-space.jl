using DISSEQT
using DISSEQT.Plots
using DISSEQT.AnnotatedArrays
using Statistics
using DataFrames
using CSV
using JLD
using Gadfly
using Cairo
using Fontconfig


function getreferencegenomes(referenceFolder::AbstractString) 

    paths = readdir(referenceFolder, join=true)
    names = map(x -> x |> basename |> splitext |> first, paths)
    referenceGenomes = map(loadfasta, paths) 
    
    if all(x->length(x)==1, referenceGenomes)
        referenceGenomes = map(x->x[1][2], referenceGenomes) # NB: Only for genomes with 1 segment, get the sequence part
    end

    names, referenceGenomes
end

function calculate_limits(swarms::Array, limitOfDetection::AnnotatedArray)

    # get rid of annotations to simplify math stuff
    # arrays are 4x4x4xPxN - nuc3 x nuc2 x nuc1 x positions x samples
    s = convert(Array, swarms)
    α = convert(Array, limitOfDetection)

    # Put a lower limit on the limit of detection
    α = max.(α,1e-3)

    # apply transformation
    X = log2.(s.+α)

    # center variables
    X = X .- mean(X;dims=5) # remove mean over samples


    # reshape to VxN - variables x samples
    X = reshape(X, 64*size(X,4), size(X,5))

    s, α, X
end

function talus(swarms::Array, limitOfDetection::AnnotatedArray, nbrDims::Int)

    s, α, X = calculate_limits(swarms, limitOfDetection)

    # Non-intrusive prefiltering of variables that will never contribute
    σThreshold = 1e-2 # same as the lower limit we use for projection score
    σ = dropdims(std(X;dims=2);dims=2)
    σ = σ/maximum(σ) # normalize
    X = X[σ.>=σThreshold,:]


    nbrDims = 40

    pl = talusplot(X, 1:nbrDims, Guide.title("Talus plot"))
    pl
end

function sequence_space(dwarms::Array, limitOfDetection::AnnotatedArray, nbrDims::Int)

    s, α, X = calculate_limits(swarms, limitOfDetection)
    

    # 0.005 because all optima are above this value
    σThresholds = 10.0.^range(log10(0.005),stop=0,length=1000) 
    U,Σ,V,projectionScores,signalDimensions = smssvd(X, nbrDims, σThresholds; nbrIter=100)
    
        
    # Low-Rank Representation, Sample names, positions and codon list
    seqSpaceRep = "plots/seqspacerepresentation.jld"
    save(seqSpaceRep, "U",U, "Σ",Σ, "V",V,
                      "SampleID", convert(Vector{String},metadata[:SampleID]),
                      "positions", swarms[:position][:],
                      "codons",    swarms[:codon][:],
                      "consensus", dropdims(swarms[:consensus]; dims=(1,2,3)) )


    # Projection Score Plot
    sigEnd = cumsum(signalDimensions)
    nbrSignals = length(signalDimensions)
    signalDesc = string.(1:nbrSignals, " (", signalDimensions, "d)")
    df = DataFrame(Sigma=repeat(σThresholds',nbrSignals,1)[:], ProjectionScore=projectionScores[sigEnd,:][:], SignalNbr=repeat(signalDesc,1,length(σThresholds))[:])
    coords = Coord.cartesian(xmin=log10(σThresholds[1]), xmax=log10(σThresholds[end]), ymin=0)
    pl = plot(df,x=:Sigma,y=:ProjectionScore,color=:SignalNbr,Geom.line,coords,
              Scale.x_log10,Scale.color_discrete(),Guide.xlabel("σ Threshold"),Guide.ylabel("Projection Score"),Guide.colorkey(title="Signal"))


    # make SMSSVD plots

    virusColors = Dict("WT"  =>colorant"black",
                       "Less"=>colorant"blue",
                       "More"=>colorant"green",
                       "Stop"=>colorant"red")
    virusColorMap = groupcolors(metadata[:Reference],virusColors)

    mutagenColors = Dict("5FU" =>RGB(1,0,0),
                         "AML" =>RGB(0,1,0),
                         "AZC" =>RGB(0,0,1),
                         "Mn"  =>RGB(0.8,0.8,0),
                         "Mock"=>RGB(1,0,1),
                         "RIBA"=>RGB(0,1,1))
    mutagenColorMap = groupcolors(metadata[:Mutagen],mutagenColors)

    P = pairwisescatterplot(V, metadata[:Reference], virusColorMap,
                               metadata[:Mutagen],   mutagenColorMap,
                               point_size=0.6mm)
 
end


function main()
    

    # projectFolder should point to MyProject folder (containing subfolders Metadata and Analysis)
    projectFolder   = "syn11639899" # VignuzziLabPublic/Projects/FitnessLandscapes
    analysisFolder  = joinpath(projectFolder, "Analysis")
    alignmentFolder = joinpath(analysisFolder, "Alignment")
    referenceFolder = joinpath(alignmentFolder, "reference_genomes") 

    
    isdir("MutantSwarms") || mkdir("MutantSwarms")
    isdir("plots") || mkdir("plots")

    metadataID = joinpath(projectFolder, "Metadata", "fitness landscapes invitro.csv")
    metadata = CSV.read(metadataID)

    # skip samples with unknown dose
    no_unknown_dose_metadata = metadata[ .~ismissing.(metadata[:Dose]), : ] 
    
    # Get Reference Genomes
    referenceNames, referenceGenomes = getreferencegenomes(referenceFolder)
	
    # get CVB3 feature list
    featureFileID = "CVB3_features.csv"
    features = CSV.read(featureFileID; missingstrings=["","NA"])
    positions = featurepositions(features, "ORF") # coding positions

    # load swarms
    swarms  = loadswarm(metadata[:SampleID],metadata[:SwarmPath],metadata[:ConsensusPath])

    # filter by positions
    swarms = swarms[:,:,:,inset(:position,positions),:]


    coverage = dropdims(swarms[:coverage]; dims=(1,2,3))
    meanCoverage = mean(coverage; dims=1)[:]

    # filter by mean read coverage at coding sites
    sampleMask = meanCoverage .>= 1000
    println("Removed ", count(.~sampleMask), " samples due to low coverage.")
    swarms  = swarms[:,:,:,:,sampleMask]
    metadata = metadata[sampleMask,:]


    # remove mixed samples
    mixes = estimatemix(swarms, referenceGenomes)
    mixMask = maximum(mixes;dims=2)[:] .>= 0.98
    println("Removed ", count(.~mixMask), " mixed samples.")
    swarms  = swarms[:,:,:,:,mixMask]
    metadata = metadata[mixMask,:]



    println("Total number of samples used: ", size(swarms,5))


    # load limit of detection
    limitOfDetectionID = joinpath(analysisFolder,"LimitOfDetection","fitness landscapes invitro","limitofdetection.jld")
    limitOfDetection = convert(AnnotatedArray,load(limitOfDetectionID, downloadLocation="MutantSwarms"))
    limitOfDetection = limitOfDetection[:,:,:,inset(:position,positions)] # filter by positions
    @assert isequal(swarms[:position][:], limitOfDetection[:position][:]) # the positions must be the same

end
