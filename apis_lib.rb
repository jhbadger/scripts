# helper functions

# returns true if file likely to be DNA, false otherwise
def isDNA?(fasta)
  seq = File.read(fasta, 10000).split("\n").grep(/^[^>]/).join
  seq.count("AGTCN").to_f / seq.length > 0.90
end

# returns species from phylodb-formattted string
def headerSpecies(header)
  header.split("::")[1].split("||")[0].tr("()","")
end

# returns pure seguid from phylodb-formatted string
def headerSeguid(header)
  header.gsub(">","").gsub("lcl|","").split(" ")[0]
end

#returns first word from header
def headerName(header)
  header.split(" ")[0]
end

class Array # additional methods for Array class
  # return majority consensus for counts array
  def majority
    consensus = []
    size.times do |i|
      total = 0.0
      self[i].values.each{|val|total+=val}
      name = self[i].keys.sort {|x,y| self[i][y] <=> self[i][x]}.first
      
      if (self[i][name]/total > 0.5)
        consensus[i] = name
      else
        consensus[i] = "Mixed"
      end
    end
    return consensus
  end

  # return absolute consensus for counts array
  def absolute
    consensus = []
    size.times do |i|
      if (self[i].size == 1)
        consensus.push(self[i].keys.first)
      else
        consensus.push("Mixed")
      end
    end
    return consensus
  end
end

class NewickTree # Additional methods for NewickTree class
  # returns array of consensus taxonomy at each relative level of tree
  def consensusTax(pep, taxonomy, ruleMaj)
    consensus = []
    return  [] if relatives(pep).nil?
    relatives(pep).each do |list|
      counts = []
      list.each do |relative|
        seqid, sp = relative.split("__")
        groups = taxonomy[sp]
        next if groups.nil?
        groups.size.times do |i|
          counts[i] = Hash.new if counts[i].nil?
          counts[i][groups[i]] = 0 if counts[i][groups[i]].nil?
          counts[i][groups[i]] += 1
        end
      end
      if (ruleMaj)
        consensus.push(counts.majority)
      else
        consensus.push(counts.absolute)
      end
    end
    return consensus
  end

  # returns classification of node based on taxonomy
  def createClassification(pep, exclude, taxonomy, ruleMaj)
    cons = consensusTax(pep, taxonomy, ruleMaj)
    lines = []
    cons.each do |line|
      excluded = false
      excluded = exclude.to_a.collect{|x| line.grep(/#{x}/)}.flatten
      lines.push(line) if excluded.empty?
    end
    first = lines[0]
    first=[nil,nil,nil,nil,nil,nil,nil] if first.nil?
    if (lines[1].nil?)
      second = nil
    else
      second = lines[1]
    end
    mixed = false
    classification = []
    7.times do |level|
      mixed = true if first[level] == "Mixed"
      first[level] = "Mixed" if mixed
      if (first[level] == "Mixed" || second.nil? || first[level] == second[level])
        outgroup = 0
      else
        outgroup = 1
      end
      first[level] = "Undefined" if first[level].nil?
      classification.push(first[level][0..45])
      classification.push(outgroup)
    end
    return classification
  end
end

# function to generate seguid link to metarep for drawing tree
def segLink(entry)
  metalink = "http://www.jcvi.org/phylo-metarep/phylodb/seguid/"
  return metalink + entry.gsub(":","<>").gsub("/", "[]").gsub("+","()")
end

# returns mga called proteins from DNA
def asProt(fasta, verbose)
	header = nil
  orfs = Hash.new
  STDERR << "Running mga to find ORFS...\n" if verbose
  `mga #{fasta}`.split("\n").each do |line|
  	if (line =~/^#/ && (line !~ /gc =/ && line !~ /self:/))
    	header = line.chomp.split("# ")[1].split(" ").first
    elsif (line =~/^gene/)
    	n, s, e, strand, frame = line.chomp.split(" ")
      orfs[header] = [] if (orfs[header].nil?)
      orfs[header].push("#{s} #{e} #{strand} #{frame}")
    end
  end
  STDERR << "Writing peptides...\n" if verbose
  pep = fasta + ".pep"
  out = File.new(pep, "w")
  Bio::FlatFile.new(Bio::FastaFormat, ZFile.new(fasta)).each do |seq|
  	if (orfs[seq.full_id])
    	id = seq.full_id
      seq = Bio::Sequence::NA.new(seq.seq)
      orfs[id].each do |orf|
      	s, e, strand, frame = orf.split(" ")
        s = s.to_i
        e = e.to_i
        frame = frame.to_i + 1
        subseq = seq.subseq(s, e)
        next if (subseq.length < 3*minOrf)
        id.gsub!("/","!")
        subseq = subseq.complement if strand == "-"
       	trans = subseq.translate(frame, 11)
        out.print trans.to_fasta("#{id}_#{s}_#{e}_#{frame}_#{strand}", 60)
      end
    end
    out.close
   end
   pep
 end

# load taxonomy and return species-based hash
def loadTaxonomy(taxf, verbose)
  # recusively walk up tax tree
  def recurseTaxonomy(tax, current)
   name, parent, rank = tax[current]
    if (name.nil? || parent.to_i == 1)
      []
    else
      recurseTaxonomy(tax, parent).to_a + [name] 
    end
  end
  STDERR << "Loading taxonomy...\n" if verbose
  tax = Hash.new
  sp = Hash.new
  if taxf =~/.gz/ # compressed file
    intax = IO.popen("gunzip -c #{taxf}")
  else
    inTax = File.new(taxf)
  end
  inTax.each do |line|
    current, name, parent, rank = line.chomp.split("\t")
    name = name.tr("()","")
    tax[current.to_i] = [name, parent.to_i, rank]
    sp[name] = current.to_i if rank == "species"
  end
  inTax.close
  taxonomy = Hash.new
  sp.keys.each do |s|
    taxonomy[s.gsub(" ","_")] = recurseTaxonomy(tax, sp[s])
  end
  taxonomy
end

# load peptides (if needed) into sqlite db, returning db handle
def loadPeps(fasta, verbose)
	db = SQLite3::Database.new(File.basename(fasta) + "_pep.db")
	begin
		db.execute("CREATE TABLE peptides (name, seq, processed)")
		db.execute("CREATE UNIQUE INDEX name_idx ON peptides(name)")
		db.execute("CREATE INDEX processed_idx ON peptides(processed)")
	rescue
	end
	if db.get_first_value("SELECT count(*) FROM peptides") == 0
		STDERR << "Loading peptides...\n" if verbose
		Bio::FlatFile.new(Bio::FastaFormat, File.new(fasta)).each do |seq|
			name = headerName(seq.definition)
      begin
			 db.execute("INSERT INTO peptides VALUES(?,?,?)", name, seq.seq, 0)
      rescue
        STDERR << "Error #{$!} loading peptide #{pep}...\n" if verbose
      end
		end
	end
	db
end

# load BLAST output into sqlite db, returning db handle
def loadBlast(fasta, blasts, ethresh, maxhits, verbose)
  db = SQLite3::Database.new(File.basename(fasta) + "_blast.db")
  begin
    db.execute("CREATE TABLE blast (name, hit, evalue, score)")
    db.execute("CREATE INDEX name_idx ON blast(name)")
  rescue
  end
  if db.get_first_value("SELECT count(*) FROM blast") == 0
    STDERR << "Loading blast results...\n" if verbose
    counts = Hash.new
    sortCmd = "sort -t $'\t' -k1 -k12 -r -n -u"
    `cat #{blasts.join(" ")} | #{sortCmd}`.split("\n").each do |row|
      if row !~ /^QUERY/
        row = row.split("\t")
        name, hit, evalue, score = row[0], row[1], row[10].to_f, row[11].to_f
        counts[name] = 0 if !counts[name]
        if evalue <= ethresh && counts[name] < maxhits
          begin
            db.execute("INSERT INTO blast VALUES(?,?,?,?)", name, hit, evalue, score)
            counts[name] += 1
          rescue
            STDERR << "Error #{$!} loading blast for #{pep}...\n" if verbose
          end
        end
      end
    end
  end
  db
end

# creates alignment db and returns handle to it
def createAlignmentDB(fasta)
  db = SQLite3::Database.new(File.basename(fasta) + "_alignment.db")
  begin
    db.execute("CREATE TABLE alignment (name, alignment)")
    db.execute("CREATE UNIQUE INDEX name_idx ON alignment(name)")
  rescue
  end
  db
end

# creates trees db and returns handle to it
def createTreesDB(fasta)
  db = SQLite3::Database.new(File.basename(fasta) + "_tree.db")
  begin
    db.execute("CREATE TABLE trees (name, tree)")
    db.execute("CREATE UNIQUE INDEX name_idx ON trees(name)")
  rescue
  end
  db
end

# return seqs from fastacmd formatted blast database
def fetchSeqs(blastids, database)
  seqs = []
  `fastacmd -d #{database} -s "#{blastids.join(',')}"`.split(/^>/).each do |seq|
    lines = seq.split("\n")
    if (!lines.empty?)
      header = lines.shift
      out = ">" + headerSeguid(header) + "__" + headerSpecies(header) + "\n"
      out += lines.join("\n")
      seqs.push(out)
    end
  end
  seqs
end

# runs muscle to align sequences, populates db returns aligment
def align(db, pep, peps, blastids, database, verbose)
  STDERR << "Making alignment for " << pep << "...\n" if verbose
  seq = peps.get_first_value("SELECT seq FROM peptides WHERE name = ?", pep)
  homologs = fetchSeqs(blastids, database)
  hom = pep + ".hom"
  afa = pep + ".afa"
  out = File.new(hom, "w")
  out.print ">" + pep + "\n" + seq  + "\n" + homologs.join("\n")
  out.close
  begin
    system("muscle -in '#{hom}' -out '#{afa}' -quiet")
    db.execute("INSERT INTO alignment VALUES(?,?)", pep, File.read(afa))
    File.unlink(hom)
  rescue
    STDERR << "Error #{$!} aligning " << pep << "...\n" if verbose
  end
  afa
end

# produces stock format needed for quicktree from fastaFormat
def fasta2Stockholm(alignFile)
  stock = alignFile + ".stock"
  stockf = File.new(stock, "w")
  stockf.printf("# STOCKHOLM 1.0\n")
  align = Hash.new 
  aSize = 0
  nSize = 0
  Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each do |seq|
    name = headerName(seq.definition)
    align[name] = seq.seq
    aSize = seq.seq.length
    nSize = name.size if (nSize < name.size)
  end
  0.step(aSize, 50) do |i|
    stockf.printf("\n")
    align.keys.sort.each do |key|
      stockf.printf("%-#{nSize}s %s\n", key, align[key][i..i+49]) 
    end
  end
  stockf.printf("//\n")
  stockf.close
  stock
end

# infers tree by desired method, populates db, returns tree
def infer(db, afa, pep, method, verbose)
  STDERR << "Making tree for " << pep << "...\n" if verbose
  tree = nil
  if (method == "nj")
    begin
      stock = fasta2Stockholm(afa)
      tree = NewickTree.new(`quickTree -boot 100 '#{stock}'`.tr("\n",""))
      tree.midpointRoot
    rescue
      STDERR << "Error #{$!} inferring nj tree for " << pep << "...\n" if verbose
    end
  end
  begin
    db.execute("INSERT INTO trees VALUES(?,?)", pep, tree.to_s)
    File.unlink(afa, stock)
  rescue
    STDERR << "Error #{$!} writing tree for " << pep << "...\n" if verbose
  end
  tree
end

# creates classification db and returns handle to it
def createClassificationDB(fasta)
  ranks = ["kingdom", "phylum", "class", "ord", "family", "genus", "species"]
  cols = ranks.collect{|x| [x, x + "_outgroup"]}.flatten
  begin
    db = SQLite3::Database.new(File.basename(fasta) + "_classification.db")
    db.execute("CREATE TABLE classification (name, #{cols.join(",")})")
    db.execute("CREATE UNIQUE INDEX name_idx ON classification(name)")
    cols.each do |col|
      db.execute("CREATE INDEX #{col}_idx ON classification(#{col})")
    end
  rescue
  end
  db
end

# classifies tree and populates db
def classify(db, pep, tree, ruleMaj, exclude, taxonomy, verbose)
  STDERR << "Making classification for " << pep << "...\n" if verbose
  classification = tree.createClassification(pep, exclude, taxonomy, ruleMaj)
  begin
    line = ([pep] + classification).collect{|x| "\"#{x}\""}.join(",")
    db.execute("INSERT INTO classification VALUES(#{line})")
  rescue
    STDERR << "Error #{$!} writing classification for " << pep << "...\n" if verbose
  end
end
