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
  sql = "SELECT transcript, contig_ssurna_taxonomy,contig_ssurna_evalue,"
  sql += "kegg_hit,kegg_desc,kegg_pathway,ko,ko_desc,ko_pathway,ec,"
  sql += "uniprot,kog_id,kog_desc,kog_class,kog_group,organelle,"
  sql += "organelle_id,organelle_species,organelle_evalue,"
  sql += "best_hit,best_hit_gos_core_cluster,best_hit_species,"
  sql += "best_hit_taxon_id,best_hit_group,pfams,pfams_desc,"
  sql += "tigrfams,tigrfams_desc,transmembrane_domains "
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
  ["contig_ssuRNA_taxonomy", "contig_ssuRNA_e-value", "kegg_hit", "kegg_desc", 
    "kegg_pathway", "KO", "KO_desc", "KO_pathway", "EC", "uniprot", 
    "KOG_id", "KOG_desc", "KOG_class", "KOG_group", "organelle", 
    "organelle_id", "organelle_species", "organelle_e-value", "best_hit", 
    "best_hit_GOS_core_cluster", "best_hit_species", "best_hit_taxon_id", 
    "best_hit_group", "PFams", "PFams_desc", "TIGRFams", "TIGRFams_desc", 
    "transmembrane_domains"]
end


#rank hash by values descending
def rankHash(standings)
  ranks = Hash.new
  rank = 1
  standings.keys.each do |key|
    standings[key] = 0 if !standings[key]
  end
  standings.keys.sort{|x,y| standings[y]<=>standings[x]}.each do |key|
    ranks[key] = rank
    rank += 1
  end
  ranks
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
def gene_rpkm(transcript, counts, total_mapped_reads, gene_length)
  begin
    sprintf("%.2f",(1e9*counts.to_f)/(total_mapped_reads*gene_length)).to_f
  rescue
    STDERR << "Problem with " << transcript << " counts " << counts << " total " << total_mapped_reads << " len " << gene_length << "\n"
  end
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
      rpkms[sample][transcript] = gene_rpkm(transcript, counts[sample][transcript].to_i, totals[sample], lengths[transcript])
      maxrpkms[transcript] = rpkms[sample][transcript] if rpkms[sample][transcript] > maxrpkms[transcript]
    end
  end
  return rpkms, maxrpkms
end
