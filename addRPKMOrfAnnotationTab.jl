#!/usr/bin/env julia

using ArgParse

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

function totalMapped(ann, samples)
    println(STDERR, "Calculating total mapped reads...")
    total = 0
    for sample in samples
        total += reduce(+, ann[sample])
    end
    return total
end

function rpkm(num, totalMapped, len)
    if num != nothing && len != nothing
        return int(1e11*float(num)/(totalMapped*len))/1e2
    else
        return 0.0
    end
end

function addRPKM(ann, samples, lens, total)
    println(STDERR, "Calculating RPKMs...")
    withRPKM = copy(ann)
    genes = withRPKM[:orf_id]
    for sample in samples
        rpkms = Float64[]
        for gene in genes
            counts = first(withRPKM[withRPKM[:orf_id].==gene, sample])
            rpkmVal = rpkm(counts, total, lens[gene])
            push!(rpkms, rpkmVal)
        end
        rpkmName = symbol(string(sample)*"_RPKM")
        withRPKM[rpkmName] = rpkms
    end
    return sort(withRPKM, cols=:Total_Counts, rev=true)
end

function addTotalSample(ann, samples)
    genes = ann[:orf_id]
    totalCounts = Int[]
    for gene in genes
        rowtot = 0
        for sample in samples
            rowtot += first(ann[ann[:orf_id].==gene, sample])
        end
        push!(totalCounts, rowtot)
    end
    ann[symbol("Total_Counts")] = totalCounts
end

using DataFrames

function main()
    parsed_args = parse_commandline()
    println(STDERR, "Loading data...")
    ann = readtable(parsed_args["ann"], separator = '\t')
    samples = names(ann)[34:length(names(ann))]
    addTotalSample(ann, samples)
    lens = loadGeneLengths(parsed_args["bed"])
    total = totalMapped(ann, samples)
    withRPKM = addRPKM(ann, samples, lens, total)
    outFile = first(split(parsed_args["ann"],".tab"))*"_rpkm.tab"
    writetable(outFile, withRPKM, separator = '\t')
end

main()
