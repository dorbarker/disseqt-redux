import Cairo
export Cairo

import Compose
import Fontconfig

using Distributed
nprocs()==1 && addprocs()
using DISSEQT
using DISSEQT.AlignUtils
#using SynapseClient
#using DISSEQT.SynapseTools
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
	isempty(str) || println(log, "ERROR: ", str)
end

# only prints if str is nonempty
function printifwarning(log,str)
	str = strip(str)
	isempty(str) || println(log, "WARNING: ", str)
end

function printifinfo(log,str)
	str = strip(str)
	isempty(str) || println(log, str)
end


function getreferenceinfo(ref::Tuple{T,String}, referenceFolder::String) where {T}
    path = joinpath(referenceFolder, ref[2])
    ref[1], ref[2], path, path
end

function getreferenceinfo(refs::Vector, referenceFolder::String)
    [getreferenceinfo(r, referenceFolder) for r in refs]
end

function bwaversion()
	"foo"
end

function samtoolsversion()
	"bar"
end

function fastqmcfversion()
	"baz"
end


function align_sample!(sample::Sample, adapters::String,
	                   outFolder::String, tempFolder::String; 
	                   log::IO=devnull, globalLog::IO=devnull, 
	                   maxAlignIterations::Int=5, consensusMinSupport::Int=100,
	                   consensusIndelMinSupport::Int=consensusMinSupport,
	                   keepUnmapped::Bool=true,
	                   nbrThreads::Int=4)
	tempFiles = String[]

	# make sure fastq-files are concatenated
	ret,singleFastq = singlefastq(sample.fastqLocal,sample.fastqLocal2,joinpath(tempFolder,sample.name),tempFiles,log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Fastq concatenation failed.")
		return ret
	end
	
	trimmedFastq = trimmedfastqname(singleFastq, joinpath(tempFolder,sample.name), tempFiles)
	ret = trimfastq(singleFastq, trimmedFastq, adapters, log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Fastq trimming failed.")
		return ret	
	end


	unsortedBam = joinpath(tempFolder, "$(sample.name)_unsorted_temp.bam")
	push!(tempFiles,unsortedBam) # we will overwrite in successive iterations

	currReference = sample.referencePathLocal
	for i=1:maxAlignIterations

		# align
		ret = bwaalign(trimmedFastq, unsortedBam, currReference; nbrThreads=nbrThreads, log=log, keepUnmapped=keepUnmapped)
		if ret != 0
			removetempfiles(tempFiles, log=log)
			printiferror(globalLog, "$(sample.name): Alignment failed.")
			return ret	
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
	ret = bamsort(unsortedBam, sample.bam; nbrThreads=nbrThreads, log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Bam sorting failed.")
		return ret
	end

	ret = bamindex(sample.bam; log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Bam indexing failed.")
		return ret
	end

	removetempfiles(tempFiles, log=log)

	println(globalLog, "$(sample.name) aligned successfully.")
	return 0
end




# helper function made to be called by pmap
function align_single!(sample::Sample, adapters::String,
                       outFolder::String, tempFolder::String,
                       maxAlignIterations::Int, consensusMinSupport::Int, 
	                   consensusIndelMinSupport::Int,
	                   keepUnmapped::Bool,
                       nbrThreads::Int)
	log = open(joinpath(outFolder, "$(sample.name).log"), "w")
	globalLog = IOBuffer()

	try 
		align_sample!(sample, adapters, outFolder, tempFolder, 
			          log=log, globalLog=globalLog, maxAlignIterations=maxAlignIterations,
			          consensusMinSupport=consensusMinSupport,
			          consensusIndelMinSupport=consensusIndelMinSupport,
		              keepUnmapped=keepUnmapped,
			          nbrThreads=nbrThreads)
	catch err
		printiferror(globalLog, string(err))
	end

	close(log)
	sample, String(take!(globalLog))
end

 
# utility function that tries to identify sample names in folder were each sample might have multiple fastq files
function find_samples(fastqPath::String, namePrefix::String="", pattern::Regex=r".+(?=_L\d+_R1_\d+.fastq.gz$)"; log=devnull)
	
	if namePrefix != "" && namePrefix[1] != '_'
		namePrefix = namePrefix * "_"
	end

	if isempty(fastqPath)
		error("Could not find \"$runName\" in \"$fastqPath\".")
	end

	fileNames = readdir(fastqPath, join=false)
	filePaths = readdir(fastqPath, join=true)

	# find files matching pattern
	# matches = map( x->match(pattern,x), fileNames )
	# mask = falses(matches)
	# map!( x->x!=nothing, mask, matches );
	matches = match.(pattern,fileNames)
	mask = matches.!=nothing

	# log warning for non-matching files
	for f in fileNames[.~mask]
		printifwarning(log, "File \"$f\" was ignored.")
	end


	# remove files that do not match pattern
	fileNames = fileNames[mask]
	filePaths = filePaths[mask]
	matches = convert(Vector{RegexMatch}, matches[mask])

	# extract the matching part of the regex as strings
	matchingNames = Vector{String}(undef,size(matches))
	map!( x->x.match, matchingNames, matches )

	# put all files with the same match together (and add the run name to the sample name)
	samples = [ Sample("$namePrefix$u", filePaths[matchingNames.==u]) for u in unique(matchingNames) ]


	for s in samples
		n = length(s.fastq)
		println(log, "Found sample \"$(s.name)\" with $n fastq files.")
	end
	samples
end

function bwaindex(reference::String; log=devnull)
	println(log, "--- Indexing with bwa index ---"); flush(log)
	out = DiskBuffer()
	out2 = DiskBuffer() # bwa index writes info to stderr...

	ret = 0
	try
		run(pipeline(`bwa index $reference`,stdout=openbuf(out),stderr=openbuf(out2)))
	catch
		ret = 1
	end

	printifinfo(log,closebuf(out))
	printifinfo(log,closebuf(out2))
	ret != 0 && printiferror(log, "BWA indexing failed.")
	flush(log)
	ret
end

function align_samples!(samples::Vector{Sample}, adapters::String,
                        outFolder::String, tempFolder::String,
                        log::IO;
                        maxAlignIterations::Int=5, consensusMinSupport::Int=1, 
                        consensusIndelMinSupport::Int=consensusMinSupport,
                        keepUnmapped::Bool=true,
                        nbrThreads::Int=4)

	println(log, "[bwa] ", bwaversion())
	println(log, "[samtools] ", samtoolsversion())
	println(log, "[fastq-mcf] ", fastqmcfversion())
	
	# Make sure all the references have bwa index files - before we start threading!!!
	refsPathsLocal = unique([s.referencePathLocal for s in samples])
	
	for r in refsPathsLocal
		if bwaindex(r) != 0
			println(log, "Indexing of reference \"$r\" failed.")
			return
		end
	end


	# Threading modelled after the pmap implementation in http://docs.julialang.org/en/release-0.4/manual/parallel-computing/ (NB: not the same as actual pmap())
	# --------------------------------------------------------------------------
	procList = procs() # find the available processes
	n = length(samples)
	i = 1
	# function to produce the next work item from the queue.
	# in this case it's just an index.
	nextidx() = (idx=i; i+=1; idx)
	@sync begin
		for p in procList
			if p != myid() || length(procList)==1
				@async begin
					while true
						idx = nextidx()
						if idx > n
							break
						end

						s = samples[idx]

						println("Aligning sample $(s.name)")

						# align in worker thread
						samples[idx], logStr = fetch(@spawnat p align_single!(s, adapters, outFolder, tempFolder, maxAlignIterations, consensusMinSupport, consensusIndelMinSupport, keepUnmapped, nbrThreads))
						print(log, logStr); flush(log)
						println("Finished aligning sample $(s.name)")
					end
				end
			end
		end
	end
	# --------------------------------------------------------------------------
end
 
function main()
# --- Synapse login ------------------------------------------------------------

    # syn = SynapseClient.login()
    syn=nothing
    isdir("synapsecache") || mkdir("synapsecache")

# --- Setup --------------------------------------------------------------------
    # Set clean=true to rerun alignment. For clean=false, alignment will only be run if the log file is missing (i.e. no previous alignment was done).
    clean = true 

    # Name of the run. (Should corrsespond to a subfolder of fastqPath.)
    runName = "664V9AAXX"
    projectFolder = runName # VignuzziLabPublic/Projects/FitnessLandscapes
    alignmentFolder = joinpath(projectFolder, "Analysis", "Alignment")

    # Set uploadPath="synapseID" to upload files after running. Should point to "MyProject/Analysis/Alignment" folder. Set upload=nothing to skip uploading.
    uploadPath = alignmentFolder

    # Where to find sample .fastq files. Synapse Folder ID or local folder. Synapse ID should point to "MyProject/Raw Data/Sequencing".
    fastqPath = joinpath(projectFolder, "fastqs")


    # This prefix will be added to all sample names. Normally same as runName. Set to "" if the .fastq files already have this prefix.
    namePrefix = runName

    # Where to find reference genomes. Synapse Folder ID or local folder. [Reference genomes will be uploaded if local.]
    # Synapse ID should point to "MyProject/Analysis/Alignment/ReferenceGenomes".
    referenceFolder = joinpath(projectFolder, "reference_genomes")

    # Rules for matching sample IDs to references.
    # The rule can either be a Regex or a function taking the sample ID and returning true/false.
    refs = [(r"_p\d+_WT_", "WT.fasta"),
            (r"_p\d+_299_", "CVB3_299.fasta"), 
            (r"_p\d+_372_", "CVB3_372.fasta")]
            #(r"Undetermined", "phiX174.fasta")] # Always keep this unless PhiX wasn't used in the sequencing.
    refs = getreferenceinfo(refs, referenceFolder) # Get path/synapse id and local path info.

    # File with adapters. Synapse File ID or local file. [Uploaded if local.]
    adapters = joinpath(projectFolder, "adapters", "adapters.fa")

    # Local logFile.
    logFile = "AlignUtils.log"

    # Local file with sample info. Required for Uploading of a previously run alignment.
    sampleInfoFile = "sampleInfo.jld"

    # Minimum support for updating consensus during alignment. (Otherwise fallback to reference sequence.)
    consensusMinSupport = 100

# --- Cleanup ------------------------------------------------------------------
    clean = clean || !isfile(logFile)

    if clean
        makecleanfolder("bam") || return
        makecleanfolder("temp") || return
    else
        println("Skipping Alignment [clean=false]")
    end


# --- Alignment ----------------------------------------------------------------

    samples = Sample[]
    if clean
        log = open(logFile,"w")
        isfile(sampleInfoFile) && rm(sampleInfoFile)

        start = time()
        println(log,"Starting alignment batch run at $(now())"); flush(log)
        samples = find_samples(fastqPath, namePrefix, log=log); flush(log)

        assign_reference!(samples, refs,log=log); flush(log)

        align_samples!(samples, adapters, joinpath(projectFolder, "bam"), "temp",log,consensusMinSupport=consensusMinSupport); flush(log)
        reference_sanity_check(samples, refs, log=log)

        duration = time()-start
        println(log,"Finished alignment batch run in $(duration)s."); flush(log)

        close(log)

        save(sampleInfoFile, "samples", samples, compress=true)
    end

# --- Upload -------------------------------------------------------------------
#     if uploadPath != nothing
#         if !clean # reload sample info if necessary
#             try
#                 samples = load(sampleInfoFile, "samples")
#             catch
#                 println("Sample info file is missing or corrupt. Please rerun alignment.")
#             end
#         end
# 
#         refPaths = [r[3] for r in refs]
#         uploadaligned(syn, uploadPath, runName, @__FILE__, logFile, adapters, refPaths, samples)
#     end


end

main()
