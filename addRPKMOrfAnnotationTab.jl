#!/usr/bin/env julia

using ArgParse
using DataFrames

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "ann"
        help = "rap orf counts annotation file"
        required = true
        arg_type = String
         "bed"
        help = "bed file of gene positions"
        required = true
        arg_type = String
    end
    
    return parse_args(s)
end

function loadGeneLengths(bed)
    println(STDERR, "Loading gene lengths...")
    geneLengths = Dict{String,Int}()
    fp = open(bed)
    for line in eachline(fp)
        contig, start, stop, exon, score, strand = split(chomp(line), '\t')
        start = int(start)
        stop = int(stop)
        gene, rest = split(exon, ':')
        if !haskey(geneLengths, gene)
            geneLengths[gene] = 0
        end
        geneLengths[gene] += (stop - start) + 1
    end
    close(fp)
    return geneLengths
end

function rpkm(num, totalMapped, len)
    if num != nothing && len != nothing
        return int(1e11*float(num)/(totalMapped*len))/1e2
    else
        return 0.0
    end
end

function addRPKM(ann, samples, lens)
    println(STDERR, "Calculating RPKMs...")
    withRPKM = copy(ann)
    genes = withRPKM[1]
    for sample in samples
        tot = reduce(+, ann[sample])
        rpkms = Float64[]
        for gene in genes
            counts = first(withRPKM[withRPKM[1].==gene, sample])
            rpkmVal = rpkm(counts, tot, lens[gene])
            push!(rpkms, rpkmVal)
        end
        rpkmName = symbol(replace(string(sample),"_Count","")*"_RPKM")
        withRPKM[rpkmName] = rpkms
    end
    return sort(withRPKM, cols=:Total_Counts, rev=true)
end

function addTotalSample(ann, samples)
    genes = ann[1]
    totalCounts = Int[]
    for gene in genes
        rowtot = 0
        for sample in samples
            rowtot += first(ann[ann[1].==gene, sample])
        end
        push!(totalCounts, rowtot)
    end
    ann[symbol("Total_Counts")] = totalCounts
end

function main()
    parsed_args = parse_commandline()
    println(STDERR, "Loading data...")
    ann = readtable(parsed_args["ann"], separator = '\t')
    samples = filter(x->contains(string(x), "_Count"), names(ann))
    if isempty(samples)
        samples = names(ann)[34:end] # from John's "rap"
    end
    addTotalSample(ann, samples)
    lens = loadGeneLengths(parsed_args["bed"])
    withRPKM = addRPKM(ann, samples, lens)
    outFile = replace(replace(parsed_args["ann"],".tab",""),".tsv","")*"_rpkm.tsv"
    writetable(outFile, withRPKM, separator = '\t')
end

main()
