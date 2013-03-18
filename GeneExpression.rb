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
  sql = "SELECT transcript, kegg_hit, kegg_desc, kegg_pathway, ko, ko_desc, "
  sql += "ko_pathway, ec, uniprot, organelle, organelle_id, organelle_species, "
  sql += "organelle_evalue, best_hit, best_hit_gos_core_cluster, best_hit_species, "
  sql += "best_hit_taxon_id, best_hit_group, pfams, pfams_desc, tigrfams, "
  sql += "tigrfams_desc, kog, kog_desc, kog_class, kog_group, transmembrane_domains "
  sql += "FROM annotation WHERE org=?"
  ann = Hash.new
  db.execute(sql, org).each do |row|
    transcript = row.shift
    ann[transcript] = row
  end
  ann
end

#return set of fields in annotation table
def ann_headers
  ["kegg hit", "kegg desc", "kegg pathway", 
    "ko", "ko desc", "ko pathway", "ec", "uniprot", "organelle",
    "organelle id", "organelle_species", "organelle evalue",
    "best hit", "best hit gos core cluster", "best hit species",
    "best hit taxon id", "best hit group", "pfams", "pfams desc",
    "tigrfams", "tigrfams desc", "kog", "kog desc", "kog class",
    "kog group", "transmembrane domains"]
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
  db.execute("SELECT transcript, length FROM gen_length WHERE org = ?", org).each do |row|
    transcript, len = row
    lengths[transcript] = len.to_i
  end
  lengths
end

# compute gene-based rpkm from counts, total_mapped_reads, gene_exons_length
def gene_rpkm(counts, total_mapped_reads, gene_length)
  return sprintf("%.2f",(1e9*counts.to_f)/(total_mapped_reads*gene_length)).to_f
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
