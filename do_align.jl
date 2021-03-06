using ArgParse

using DISSEQT
using DISSEQT.AlignUtils

# I fully realize this is madness
import DISSEQT.AlignUtils.singlefastq
import DISSEQT.AlignUtils.concatfastq
import DISSEQT.AlignUtils.removetempfiles
import DISSEQT.AlignUtils.trimmedfastqname
import DISSEQT.AlignUtils.bwaalign
import DISSEQT.AlignUtils.trimfastq
import DISSEQT.AlignUtils.bwaindexfilescreated
import DISSEQT.AlignUtils.bamsort
import DISSEQT.AlignUtils.bamindex
import DISSEQT.AlignUtils.bwaindex
import DISSEQT.AlignUtils.bwaindexfilescreated

using JLD
using Dates


mutable struct DiskBuffer
	filename::String
	stream#::IOStream
end

DiskBuffer() = DiskBuffer(tempname(),nothing)
openbuf(db::DiskBuffer) = db.stream = open(db.filename, "w")
function closebuf(db::DiskBuffer)::String
	db.stream == nothing && return ""
	close(db.stream)
	str = read(db.filename, String)
	rm(db.filename) # delete file
	str
end

function printiferror(log,str)
	str = strip(str)
	isempty(str) || println("ERROR: ", str)
end

# only prints if str is nonempty
function printifwarning(log,str)
	str = strip(str)
	isempty(str) || println("WARNING: ", str)
end

function printifinfo(log,str)
	str = strip(str)
	isempty(str) || println(str)
end

function getreferenceinfo(refs::Vector, reference_genome::String)
    #[getreferenceinfo(r, referenceFolder) for r in refs]
    
    # Weird redundancy for backward compatibility
    [(r[1], r[2], reference_genome, reference_genome) for r in refs]
end

function align_sample(sample::Sample, adapters::String,
	                   outFolder::String, tempFolder::String; 
	                   log::IO=devnull, globalLog::IO=devnull, 
	                   maxAlignIterations::Int=5, consensusMinSupport::Int=100,
	                   consensusIndelMinSupport::Int=consensusMinSupport,
	                   keepUnmapped::Bool=true,
	                   threads::Int=1)
	tempFiles = String[]
	# make sure fastq-files are concatenated
	ret, singleFastq = singlefastq(sample.fastq, sample.fastq2, joinpath(tempFolder, sample.name), tempFiles, log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Fastq concatenation failed.")
		
	end
	
	trimmedFastq = trimmedfastqname(singleFastq, joinpath(tempFolder,sample.name), tempFiles)
	ret = trimfastq(singleFastq, trimmedFastq, adapters, log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Fastq trimming failed.")

	end


	unsortedBam = joinpath(tempFolder, "$(sample.name)_unsorted_temp.bam")
	push!(tempFiles,unsortedBam) # we will overwrite in successive iterations

	currReference = sample.referencePathLocal
	for i=1:maxAlignIterations

		# align
		ret = bwaalign(trimmedFastq, unsortedBam, currReference; nbrThreads=threads, log=log, keepUnmapped=keepUnmapped)
		if ret != 0
			removetempfiles(tempFiles, log=log)
			printiferror(globalLog, "$(sample.name): Alignment failed.")

		end


		# check if consensus equals reference
		println(log, "--- Computing consensus (minimum support for substitutions=$consensusMinSupport and indels=$consensusIndelMinSupport) ---")
		cons = consensus(unsortedBam, currReference, minSupport=consensusMinSupport, indelMinSupport=consensusIndelMinSupport)
		#ref  = collect(FastaReader(currReference))
		ref  = loadfasta(currReference)
		if isequal(cons, ref) 
			println(log, "Consensus equal to reference."); flush(log)
			break # done if they are equal
		end
		println(log, "Consensus not equal to reference."); flush(log)

		# otherwise save the new consensus as reference
		currReference = joinpath(tempFolder, "$(sample.name)_consensus_$i.fasta")
		push!(tempFiles,currReference)
		savefasta(currReference, cons)

		# and generate a bwa index for it...
		append!(tempFiles,bwaindexfilescreated(currReference))
		ret = bwaindex(currReference; log=log)
		if ret != 0
			removetempfiles(tempFiles, log=log)
			printiferror(globalLog, "$(sample.name): BWA indexing of consensus failed.")
			return ret	
		end



		if i==maxAlignIterations
			printiferror(log, "Consensus did not converge within $maxAlignIterations iterations.")
			printiferror(globalLog, "$(sample.name): Consensus failed to converge.")
			flush(log)
			return 1
		end
	end

	# copy consensus
	sample.consensus = joinpath(outFolder, sample.name * "_consensus.fasta")
	cp(currReference,sample.consensus,force=true)


	sample.bam = joinpath(outFolder, "$(sample.name).bam")
	ret = bamsort(unsortedBam, sample.bam; nbrThreads=threads, log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Bam sorting failed.")

	end

	ret = bamindex(sample.bam; log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Bam indexing failed.")

	end

	removetempfiles(tempFiles, log=log)

	println(globalLog, "$(sample.name) aligned successfully.")
	return sample 
end




# helper function made to be called by pmap
function align_single(sample::Sample, adapters::String,
                       outFolder::String, tempFolder::String,
                       maxAlignIterations::Int, consensusMinSupport::Int, 
	                   consensusIndelMinSupport::Int,
	                   keepUnmapped::Bool,
                       threads::Int)
	log = open(joinpath(outFolder, "$(sample.name).log"), "w")
	globalLog = IOBuffer()

	try 
		sample = align_sample(sample, adapters, outFolder, tempFolder, 
			          log=log, globalLog=globalLog, maxAlignIterations=maxAlignIterations,
			          consensusMinSupport=consensusMinSupport,
			          consensusIndelMinSupport=consensusIndelMinSupport,
		              keepUnmapped=keepUnmapped,
			          threads=threads)
	catch err
		printiferror(globalLog, string(err))
	end

	close(log)
	sample, String(take!(globalLog))
end

function find_samples(fastqPath::String)
    
    sample_directories = readdir(fastqPath)

    samples = Vector{Sample}(undef, length(sample_directories))

    for (i, sample) in enumerate(sample_directories)
        
        sample_fastqs = readdir(joinpath(fastqPath, sample), join=true)

        samples[i] = Sample(sample, sample_fastqs)

    end
    
    samples
end

function align_samples(samples::Vector{Sample}, adapters::String,
                        outFolder::String, tempFolder::String,
                        log::IO;
                        maxAlignIterations::Int=5, consensusMinSupport::Int=1, 
                        consensusIndelMinSupport::Int=consensusMinSupport,
                        keepUnmapped::Bool=true,
                        threads::Int=4)

	
	# Make sure all the references have bwa index files - before we start threading!!!
	refsPathsLocal = unique([s.referencePathLocal for s in samples])
	
	for r in refsPathsLocal
		if bwaindex(r) != 0
			println(log, "Indexing of reference \"$r\" failed.")
			return
		end
	end

	aligned_samples = Vector{Sample}(undef, length(samples))
	for (i, s) in enumerate(samples)
		# println("Aligning sample $(s.name)")
		aligned_sample, logStr = align_single(s, adapters, outFolder, tempFolder, maxAlignIterations, consensusMinSupport, consensusIndelMinSupport, keepUnmapped, threads)
		aligned_samples[i] = aligned_sample
		print(log, logStr); flush(log)
		# println("Finished aligning sample $(s.name)")
	end
	# --------------------------------------------------------------------------
	aligned_samples
end

function arguments()
	
    s = ArgParseSettings()
    
    @add_arg_table! s begin

        "--fastqs"
            help = "Directory containing .fastq.gz files"
            required = true

        "--bams"
            help = "Directory to which .bam files will be written"    
            required = true

        "--adapters"
            help = "FASTA-formatted file containing adater sequences"
            required = true

        "--reference-genome"
            help = "FASTA-formattedreference genome"
            required = true

        "--consensus-min-support"
            help = "Minimum consensus support during alignment"
            arg_type = Int  # this might need to be a float
            default = 100

        "--threads"
            help = "Number of cpu threads to use [1]"
            arg_type = Int
            default = 1
    end

    return parse_args(s)
end
 
function main()

    args = arguments()
    
# --- Setup --------------------------------------------------------------------

    refs = [(r".*", args["reference-genome"])]
    refs = getreferenceinfo(refs, args["reference-genome"]) 

    # Local logFile.
    logFile = "AlignUtils.log"

    # Local file with sample info. Required for Uploading of a previously run alignment.
    sampleInfoFile = "sampleInfo.jld"

    mkpath(args["bams"]) 
    mkpath("temp") 

# --- Alignment ----------------------------------------------------------------

    samples = Sample[]
    log = open(logFile,"w")
    isfile(sampleInfoFile) && rm(sampleInfoFile)

    start = time()
    println(log,"Starting alignment batch run at $(now())"); flush(log)
    samples = find_samples(args["fastqs"])

    assign_reference!(samples, refs, log=log); flush(log)

    aligned_samples = align_samples(samples, 
                                    args["adapters"],
                                    args["bams"],
                                    "temp",
                                    log,
                                    threads=args["threads"],
                                    consensusMinSupport=args["consensus-min-support"]); flush(log)

    reference_sanity_check(aligned_samples, refs, log=log)

    duration = time()-start
    println(log,"Finished alignment batch run in $(duration)s."); flush(log)

    close(log)

    save(sampleInfoFile, "samples", samples, compress=true)

end

main()
