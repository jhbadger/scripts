# check to see if needed programs installed
def required_programs(array)
  array.each do |item|
      prog = `which #{item}`
      if prog == ""
        STDERR << "You don't have #{item} on your path, but it is required\n"
        exit(1)
      end
  end
end


# compute gene-based rpkm from counts, total_mapped_reads, gene_exons_length
def rpkm(counts, total_mapped_reads, gene_length)
    if counts && gene_length
      sprintf("%.2f",(1e9*counts.to_f)/(total_mapped_reads*gene_length)).to_f
    else
      0.0
    end
end

# fill out gene_counts with counts, rpkm
def computeRPKM(db, experiment, sample, gene_lengths)
	counts = Hash.new
	total = 0
	db.query("SELECT gene, SUM(counts) FROM exon_counts WHERE experiment=? AND sample=? GROUP BY gene", experiment, sample).each do |row|
		gene, count = row
		counts[gene] = count
		total += count
	end
	counts.keys.each do |gene|
		grpkm = rpkm(counts[gene], total, gene_lengths[gene])
		db.query("INSERT INTO gene_counts VALUES(?,?,?,?,?)", experiment, sample, gene, counts[gene], grpkm)
	end
end

# write bed file for bedtools counting
def writeBed(db, org, fname)
	out = File.new(fname, "w")
	db.query("SELECT contig, start, stop, exon, 0, strand FROM exons where org=? ORDER by contig", org).each do |row|
		out.print row.join("\t") + "\n"
	end
	out.close
	fname
end

# create hash of transcript lengths by summing exon lengths
def gene_lengths(db, org)
	lengths = Hash.new
	db.query("SELECT gene, SUM (1 + stop - start) FROM exons WHERE org = ? GROUP BY gene", org).each do |row|
		gene, len = row
		lengths[gene] = len
	end
	lengths
end

class SQLite3::Database
	# add count method
	def count(sql, *others)
		sql = "SELECT COUNT(*) FROM "+sql if !sql.index(/^SELECT/i)
		get_first_value(sql, *others)
	end
end