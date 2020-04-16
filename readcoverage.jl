using ArgParse
using DataFrames
using Gadfly
using DataStructures
using Libz

import DISSEQT.BamReader
import DISSEQT.readcoverage

# get a read coverage matrix per reference sequence
function readcoverage(inputs::AbstractArray; kwargs...) 
	nbrSamples = length(inputs)

	coverage = map(x->readcoverage(x;kwargs...),inputs) # Vector{Vector{Vector{Int}}} : Samples-ReferenceSequence-Coverage

	# reorganize to create one matrix per reference sequence
	nbrSeqsPerSample = map(length, coverage)
	nbrSeqs = nbrSeqsPerSample[1]

	@assert all(nbrSeqsPerSample .== nbrSeqs) "All samples must have the same number of reference sequences."

	covMatrices = Vector{Matrix{Int}}(undef,nbrSeqs)
	for i=1:nbrSeqs
		seqLengthPerSample = map( x->length(x[i]), coverage )
		seqLength = maximum(seqLengthPerSample) # To avoid failing if some samples have indels

		# gather read coverage in matrix filled with zeros at the end
		M = zeros(nbrSamples, seqLength)
		for j=1:nbrSamples
			s = seqLengthPerSample[j]
			M[j,1:s] = coverage[j][i]
		end
		covMatrices[i] = M
	end
	covMatrices
end

function listaligned(path::AbstractString)

    paths = readdir(path, join=true)
    names = readdir(path, join=false)

    # filter by extension
    mask = match.(r".bam$",names) .!= nothing
    paths[mask], names[mask]
end

function find_aligned(path::AbstractString)
	paths, names = listaligned(path)
	names = map(x->x[1:end-4], names) # keep only matches and remove ".bam" ending
	paths, names
end

function referencefromlog(logFile::AbstractString)
    logFileText = readlines(logFile) # load entire text file
    map!(rstrip, logFileText, logFileText)

    # find lines of type
    # Assigned reference "RefName" to sample "SampleName".

    pattern = r"^Assigned reference \"([^\"]+)\" to sample \"([^\"]+)\".$"
    matches = map( x->match(pattern,x), logFileText )
    matches = matches[map!(x->x!=nothing, falses(length(matches)), matches)] # remove non-matches

    logRefs = [splitext(m.captures[1])[1] for m in matches] # reference name without file ending (.fasta)
    logSamples = [m.captures[2] for m in matches]           # sample names

    logSamples, logRefs
end

function referencefromlog(logFile::AbstractString, sampleNames::AbstractArray)
    logSamples, logRefs = referencefromlog(logFile)

    references = Vector{String}(undef,length(sampleNames))
    for (i,name) in enumerate(sampleNames)
        ind = findfirst(isequal(name), logSamples)
        references[i] = ind!=0 ? logRefs[ind] : ""
    end
    references
end


function referencefromlog(logFile::AbstractString, sampleNames::AbstractString) 
    referencefromlog(logFile, [sampleNames])[1]
end


# Name of the run. (Should corrsespond to a subfolder of bamPath if in Synapse.)
runName = "664V9AAXX"

# Should point to MyProject/Analysis/Alignment
#alignmentFolder = "syn18694207" # VignuzziLabPublic/Projects/FitnessLandscapes/Analysis/Alignment

function _savecoverageplot(layer,format,outPath,args...)
    if format==:png
        args = (Theme(background_color=colorant"white"), args...)
    end

    pl = plot(layer,args...)

    filename = string(outPath,'.',format)
    if format==:png
        draw(PNG(filename, 29.7cm, 21cm), pl)
    elseif format==:pdf
        draw(PDF(filename, 29.7cm, 21cm), pl)
    elseif format==:svg
        draw(SVG(filename, 29.7cm, 21cm), pl)
    else
        error("Unknown image format $format")
    end
    filename
end


function coverageplots(bamPath, outputFolder; log=stdout)
# --- Setup --------------------------------------------------------------------

    outFormats=[:png, :pdf]
    alignLogFile = "AlignUtils.log"
# --- Cleanup ------------------------------------------------------------------
    mkpath(outputFolder) 

# --- Make plots ---------------------------------------------------------------
    # find samples
    samplePaths, sampleNames = find_aligned(bamPath)
    references = referencefromlog(alignLogFile, sampleNames)

    # split samples by reference sequence
    for reference in unique(references)
        mask = references.==reference
        nbrSamples = count(mask)
        println(log, "Processing $nbrSamples samples for reference \"$reference\".")

        filesCreated, segments = _coverageplots(samplePaths[mask], 
						sampleNames[mask],
						outputFolder,
						"readcoverage_$(basename(reference))",
						outFormats)
    end
end

function _coverageplots(samplePaths, sampleNames, outputFolder, outName, outFormats )
    println("Computing read coverage...")
    coverage = readcoverage(samplePaths)
    println("Done.")

    # check if these samples have one or multiple segments
    nbrSegments = length(coverage)

    segmentNames = [""]
    if nbrSegments>1
        seqs = sequences(BamFile(samplePaths[1])) # array of (segmentName,length)
        segmentNames = String[x[1] for x in seqs]
    end

    outFiles = Vector{String}()
    outSegment = Vector{String}()

    println("Plotting...")
    for (i,segmentName) in enumerate(segmentNames)
        # Create DataFrame with all the data for the plot
        cov = coverage[i]
        maximum(cov)==0 && (cov[1]=1) # Gadfly fix for logscale when no sample has data.
        nbrSamples = size(cov,1)
        nbrPos     = size(cov,2)

        pos = repeat( (1:nbrPos)', nbrSamples, 1)
        name = repeat(sampleNames, 1, nbrPos)
        df = DataFrame(position=pos[:], coverage=cov[:], name=name[:])

        xTickStep = Int( 10^ceil(log10(nbrPos/2))/10 )
        xTicks = xTickStep:xTickStep:nbrPos


        covLayer = layer(df,x=:position,y=:coverage,color=:name,Geom.line)
        covPlotParams = (Scale.y_log10, Coord.cartesian(xmin=1,xmax=nbrPos,ymin=0), Guide.xticks(ticks=collect(xTicks)))

        filename = outName
        if !isempty(segmentName)
            filename = "$(filename)_$segmentName"
        end

        for f in outFormats
            push!(outFiles, _savecoverageplot(covLayer,f,joinpath(outputFolder,filename),covPlotParams...))
            push!(outSegment, segmentName)
        end
    end
    println("Done.")

    outFiles, outSegment
end

function arguments()

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--bams"
            help = "Directory containing .bam files"
            required = true
        "--output"
            help = "Output directory"
            required = true          
    end

    return parse_args(s)

end

function main()
    args = arguments()

    coverageplots(args["bams"], args["output"])

end

main()
