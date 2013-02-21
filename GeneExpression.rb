# routines to process gene expression data

class GeneExpression
  # generate a poisson random number
  def GeneExpression.poisson(lambda)
    return lambda if (lambda > 700)
    l = Math.exp(-lambda)
    k = 0
    p = 1
    while p >= l
      k += 1
      u = rand
      p = p *u
    end
    k -1
  end
  
  # calculate Strekel (2000) R value
  def GeneExpression.rvalue(total, counts)
    r = 0
    sum = counts.reduce(:+)
    freq = sum/(1.0*total.reduce(:+))
    total.size.times do |i|
      if (!counts[i].nil? && counts[i] > 0)
        r += (counts[i]*Math.log((counts[i])/(total[i]*freq)))
      end
    end
    (100*r).to_i/100.0 # round nicely
  end
end

# return list of samples, descs from a given experiment in the sample_desc table
def samplesFromExperiment(db, experiment)
  STDERR << "Loading samples from database...\n"
  sql = "SELECT DISTINCT sample, desc FROM sample_desc WHERE experiment=?" 
  sdesc = Hash.new
  db.execute(sql, experiment) do |row|
    sdesc[row.first] = row.last.to_s
  end
  return sdesc.keys.sort, sdesc
end

# return list of counts for experiment by sample and transcript
def countsFromExperiment(db, experiment)
  STDERR << "Loading counts from database...\n"
  sql = "SELECT sample, transcript, count FROM counts WHERE experiment=? GROUP BY sample, transcript"
  counts = Hash.new
  totals = Hash.new
  db.execute(sql, experiment) do |row|
    sample, transcript, count = row
    counts[sample] = Hash.new if !counts[sample]
    totals[sample] = 0 if !totals[sample]
    counts[sample][transcript] = count.to_i
    totals[sample] += count.to_i
  end
  return counts, totals
end 

# return annotation info for organism
def annFromOrg(db, org)
  STDERR << "Loading annotation from database...\n"
  sql = "SELECT transcript, jgi,ncbi,annotation,swissprot,swissprot_annotation,"
  sql += "kegg,kegg_def,gos_cluster,phytax_cluster,phytax,ko,ko_def,pfam,pfam_def,"
  sql += "pfam2go,tfam,tfam_def,tfam2go,cog,cog_def,plant,animal FROM annotation "
  sql += "WHERE org=?"
  ann = Hash.new
  db.execute(sql, org).each do |row|
    transcript = row.shift
    ann[transcript] = row
  end
  ann
end

#return set of fields in annotation table
def ann_headers
  ["jgi","ncbi","annotation","swissprot","swissprot_annotation","kegg","kegg_def",
    "gos_cluster","phytax_cluster","phytax","ko","ko_def","pfam","pfam_def",
    "pfam2go","tfam","tfam_def","tfam2go","cog","cog_def","plant","animal"]
end


#rank hash by values descending
def rankHash(standings)
  ranks = Hash.new
  rank = 1
  standings.keys.sort{|x,y| standings[y]<=>standings[x]}.each do |key|
    ranks[key] = rank
    rank += 1
  end
  ranks
end


# load counts file into db
def loadCounts(file, db, org)
  experiment, sample = file.split("-")
  STDERR << "Loading " << file << "...\n"
  genes = Hash.new
  db.execute("SELECT gene_id, transcript FROM gtf WHERE org = ?", org) do |row|
    gene_id, transcript = row
    genes[gene_id] = transcript
  end
  File.new(file).each do |line|
    next if line=~/^#/
    gene, chr, strand, start, stop, count, length, rpkm, expressed_exons, transcripts = line.chomp.split("\t")
    if (genes[gene].nil?)
      STDERR << "Error loading: #{$_}. No gene found\n"
      exit(1)
    else
      db.execute("INSERT INTO counts VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", experiment, 
          sample, genes[gene], chr, strand, start, stop, count, length, rpkm, 
            expressed_exons, transcripts)
    end
  end
end

# return length of genes for organism
def gene_lengths(db, org)
  lengths = Hash.new
  db.execute("SELECT transcript, sum(abs(1+start-stop)) FROM gtf WHERE org = ? GROUP BY transcript", org).each do |row|
    transcript, len = row
    lengths[transcript] = len
  end
  lengths
end

# compute gene-based rpkm from counts, total_mapped_reads, gene_exons_length
def gene_rpkm(counts, total_mapped_reads, gene_length)
  return sprintf("%.2f",(10e8*counts.to_f)/(total_mapped_reads*gene_length)).to_f
end

# compute array of gene-based rpkms (and maximums) from counts, total, gene lengths
def rpkmsFromCounts(counts, totals, lengths)
  STDERR << "Calculating RPKMs...\n"
  rpkms = Hash.new
  maxrpkms = Hash.new
  counts.keys.each do |sample|
    rpkms[sample] = Hash.new if !rpkms[sample]
    counts[sample].keys.each do |transcript|
      maxrpkms[transcript] = 0 if !maxrpkms[transcript]
      rpkms[sample][transcript] = gene_rpkm(counts[sample][transcript].to_i, totals[sample], lengths[transcript])
      maxrpkms[transcript] = rpkms[sample][transcript] if rpkms[sample][transcript] > maxrpkms[transcript]
    end
  end
  return rpkms, maxrpkms
end
