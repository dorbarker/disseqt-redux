using ArgParse
using SynapseClient
using DISSEQT
using DISSEQT.SynapseTools
using DISSEQT.Plots
using DISSEQT.AnnotatedArrays
using LinearAlgebra
using DataFrames
using CSV
using Statistics
using JLD
using Gadfly
using Colors

function getsamplefitness(fitness_path::String)

    fitness = CSV.read(fitness_path)
    rename!(fitness, :MetadataID => :FitnessID)

# put mutagen doses into low,medium,high categories
function addlevel!(metadata)
    levelDict = Dict(("5FU","50uM")=>"Low", ("5FU","100uM")=>"Medium", ("5FU","150uM")=>"Medium", ("5FU","200uM")=>"High", ("AML","50uM")=>"Low", ("AML","100uM")=>"Medium", ("AML","200uM")=>"High", ("AZC","50uM")=>"Low", ("AZC","100uM")=>"Medium", ("AZC","200uM")=>"Medium", ("AZC","300uM")=>"High", ("Mn","0.25mM")=>"Medium", ("Mn","0.33mM")=>"Medium", ("Mn","0.5mM")=>"High", ("Mock","0uM")=>"Low", ("RIBA","50uM")=>"Low", ("RIBA","100uM")=>"Medium", ("RIBA","200uM")=>"Medium", ("RIBA","300uM")=>"High")
    metadata[:Level] = map( (m,d)->get(levelDict, (m,d), missing), metadata[:Mutagen], metadata[:Dose] )
end

function arguments()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "-x"
            required = true
            help = "Metadata column name to use as the X-variable"

        "-y"
            required = true
            help = "Metadata column name to use as the Y-variable"

        "--metadata"
            required = true
            help = "CSV metadata file"
        
        "--sequence-space-pca"
            required = true
            help = "JLD file containing sequence space analysis"
        
        "--fitness"
            required = true
            help = "CSV describing predicted fitness" 

        "--output-fitness"
            required = true
            help = "CSV file of fitness predictions"

        "--outdir-plots"
            required = true
            help = "Output plot directory"

    end

    args = parse_args(s)
    
    args["x"] = Symbol(args["x"])
    args["y"] = Symbol(args["y"])
    
    
    args
end

function get_colours(colours_path::String, metadata::String, xvar::Symbol, yvar::Symbol)
    
    colours = CSV.read(colours_path)
    
    x_colours = Dict()
    y_colours = Dict()
    
    for (k, v) in eachrow(colours)
        
        colour = parse(Colorant, v)
        if k in metadata[!, xvar]  # the sample names
            x_colours[k] = colour
        
        elseif k in metadata[!, yvar]
            y_colours[k] = colour

        else
            continue
        end

    end

    x_colours, y_colours

end

function main()
   
    args = arguments()
 
    metadata = CSV.read(args["metadata"])


    seqSpaceRepID = getchildbyname(syn, analysisFolder, "SequenceSpace", uploadName, "SequenceSpaceRepresentation", "seqspacerepresentation.jld")
    seqSpaceRepDict = load(localpath(syn, seqSpaceRepID))
    U = seqSpaceRepDict["U"]
    Σ = seqSpaceRepDict["Σ"]
    V = seqSpaceRepDict["V"]
    sampleID  = seqSpaceRepDict["SampleID"]
    positions = seqSpaceRepDict["positions"]
    codons    = seqSpaceRepDict["codons"]
    consensus = permutedims(seqSpaceRepDict["consensus"],(2,1)) # transpose to get Samples x Codons

    nbrDims = length(Σ)
   
    Y = V.*Σ' # sample representation scaled by singular values



    # also load the Consensus-based PCA representation
    seqSpaceRepConsensusDict = JLD.load(args["sequence-space-pca"])

    UConsensus = seqSpaceRepConsensusDict["U"]
    ΣConsensus = seqSpaceRepConsensusDict["Σ"]
    VConsensus = seqSpaceRepConsensusDict["V"]
    @assert sampleID == seqSpaceRepConsensusDict["SampleID"] # make sure the same samples were used
    nbrDimsConsensus = length(ΣConsensus)

    YConsensus = VConsensus.*ΣConsensus' # sample representation scaled by singular values



    # filter metadata by samples included in the Sequence Space Representation
    ind = indexin(sampleID, metadata[:SampleID])
    @assert all(ind.!=0) "Unknown SampleID"
    metadata = metadata[ind,:]




    addlevel!(metadata) # categorize dosages as low, medium, high


    # get fitness
    appendfitness!(metadata, getsamplefitness(args["fitness"]))


    # setup dependencies

    x_colours, y_colours = get_colours(args["colours"], args["metadata"], args["x"], args["y"])

    # virusColors = Dict("WT"  =>colorant"black",
    #                    "Less"=>colorant"blue",
    #                    "More"=>colorant"green",
    #                    "Stop"=>colorant"red")
    # virusColorMap = groupcolors(metadata[:Reference],virusColors)

    # mutagenColors = Dict("5FU" =>RGB(1,0,0),
    #                      "AML" =>RGB(0,1,0),
    #                      "AZC" =>RGB(0,0,1),
    #                      "Mn"  =>RGB(0.8,0.8,0),
    #                      "Mock"=>RGB(1,0,1),
    #                      "RIBA"=>RGB(0,1,1))
    # mutagenColorMap = groupcolors(metadata[:Mutagen],mutagenColors)

    x_colour_map = groupcolors(metadata[args["x"]], x_colours)
    y_colour_map = groupcolors(metadata[args["y"]], y_colours)

    # Isomap representation
    K = Y*Y'
    D = sqrt.(diag(K) .+ diag(K)' .- 2*K) # sample distance matrix
    Z,stress,converged,dists = kruskalisomap(D, 2, 2, maximum(D)/20)
    converged || error("Kruskal's stress MDS did not converge.")

    # Isomap representation (Consensus PCA)
    KConsensus = YConsensus*YConsensus'
    DConsensus = sqrt.(max.(0., diag(KConsensus) .+ diag(KConsensus)' .- 2*KConsensus)) # sample distance matrix (make sure it's nonnegative to handle rounding errors for colocated points)
    ZConsensus,_,converged,_ = kruskalisomap(DConsensus, 2, 2, maximum(DConsensus)/20)
    converged || error("Kruskal's stress MDS did not converge for Consensus PCA.")

    # prediction
    fitnessMask = .~ismissing.(metadata[:Fitness])

    # make a table with all predictors
    predictionResults = DataFrame(Name=String[], Type=String[], VarianceExplained=Float64[], NbrSamples=Int[])


    # run twice, the second time excluding samples that constitute singleton groups in the "Lineage, Mutagen, Dose" group predictor
    mask = copy(fitnessMask)
    for i=1:2
        fitness = collect(skipmissing(metadata[mask,:Fitness]))

        # Isomap landscape
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(Z[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(Z[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, Z[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] Isomap (2d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Reduced space landscape
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(Y[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(Y[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, Y[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] SMSSVD ($(nbrDims)d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )


        # Isomap landscape (Consensus PCA)
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(ZConsensus[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(ZConsensus[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, ZConsensus[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] Consensus Isomap (2d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Reduced space landscape  (Consensus PCA)
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(YConsensus[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(YConsensus[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, YConsensus[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] Consensus PCA ($(nbrDimsConsensus)d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )


        # Nearest neighbor predictor: Reduced dimension
        f = leaveoneoutpredict(NearestNeighborModel, Y[mask,:], fitness)
        push!(predictionResults, ("[NN] SMSSVD ($(nbrDims)d)", "Nearest Neighbor", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Nearest neighbor predictor: Consensus
        f = leaveoneoutpredict(NearestNeighborModel, consensus[mask,:], fitness)
        push!(predictionResults, ("[NN] Consensus", "Nearest Neighbor", 1.0.-varianceunexplained(f,fitness), count(mask)) )


        # Group predictor: xvar, yvar 
        f = leaveoneoutpredict(GroupModel, metadata[mask,[args["x"], args["y"]]], fitness)
        push!(predictionResults, ("[G] Lineage/Mutagen", "Group", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Group predictor: Lineage, Level
        # f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Level]], fitness)
        # push!(predictionResults, ("[G] Lineage/Dose", "Group", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # # Group predictor: Lineage, Drug, Dose
        # #f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Mutagen,:Dose]], fitness)
        # f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Mutagen,:Level]], fitness)
        push!(predictionResults, ("[G] Lineage/Mutagen/Dose", "Group", 1.0.-varianceunexplained(f,fitness), count(mask)) )
        newMask = i==1 ? .~isnan.(f) : trues(size(f)) # For next iteration, exclude singleton groups
        i==1 && println(count(.~newMask), " singleton groups removed.")

        mask[mask] = newMask
    end

    println(predictionResults)
    CSV.write(args["output-csv"], predictionResults; missingstring="NA")

    N = div(size(predictionResults,1),2)
    xMax = maximum(predictionResults[.~isnan.(predictionResults[:VarianceExplained]), :VarianceExplained])

    # with all samples
    df = predictionResults[1:N,:]
    df = df[.~isnan.(df[:VarianceExplained]),:] # remove predictors that couldn't be used ([G] Lineage/Mutagen/Dose)
    plotArgs = (layer(df, x=:VarianceExplained, color=:Name, Geom.bar(orientation=:horizontal)),
                Coord.Cartesian(xmin=0, xmax=xMax, yflip=true),
                Guide.yticks(ticks=nothing),
                Guide.xlabel("Variance Explained"), Guide.ylabel(""))
    
    themeArgs = (:key_title_font_size=>2*11pt, 
                 :key_label_font_size=>2*8pt,
                 :major_label_font_size=>2*11pt,
                 :minor_label_font_size=>2*8pt,
                 :point_label_font_size=>2*8pt)

    pl1 = plot(plotArgs..., Theme(;themeArgs...))
    saveplot(pl1, [:png,:svg,:pdf], joinpath(args["outdir-plots"], "fitnesspredictions_allsamples"))

    pl2 = plot(plotArgs..., Theme(key_position=:none; themeArgs...)) # no legend
    saveplot(pl2, [:png,:svg,:pdf], joinpath(args["outdir-plots"], "fitnesspredictions_allsamples_nolegend"))



    # with singleton groups removed
    df = predictionResults[N+1:end,:]
    df = df[.~isnan.(df[:VarianceExplained]),:] # remove predictors that couldn't be used (none)
    plotArgs = (layer(df, x=:VarianceExplained, color=:Name, Geom.bar(orientation=:horizontal)),
                Coord.Cartesian(xmin=0, xmax=xMax, yflip=true),
                Guide.yticks(ticks=nothing),
                Guide.xlabel("Variance Explained"), Guide.ylabel(""))

    pl3 = plot(plotArgs..., Theme(;themeArgs...))
    saveplot(pl3, [:png,:svg,:pdf], joinpath(args["outdir-plots"], "fitnesspredictions"))

    pl4 = plot(plotArgs..., Theme(key_position=:none; themeArgs...)) # no legend
    saveplot(pl4, [:png,:svg,:pdf], joinpath(args["outdir-plots"], "fitnesspredictions_nolegend"))

end

end

main()
