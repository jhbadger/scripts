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

# get taxonomy array (or string) for taxon
def getTaxonomy(taxon, taxonomy, joined=false)
  seqid, sp = taxon.split("__")
  tx=taxonomy[sp]
  if tx.nil?
    nil
  elsif joined
    return tx.join("; ")
  else
    tx
  end
end

class NewickNode # additional methods for NewickNode class
  # return array of arrays of taxa representing relatives at each level
  def relatives
    relatives = []
    node = self
    while(!node.nil? && !node.parent.nil?)
      relatives.push(node.parent.taxa - node.taxa)
      node = node.parent
    end
    return relatives
  end

  # returns array of consensus taxonomy at each relative level of tree
  def consensusTax(taxonomy, ruleMaj = false)
    consensus = []
    rels = relatives
    return  [] if rels.nil?
    rels.each do |list|
      counts = []
      list.each do |relative|
        groups = getTaxonomy(relative, taxonomy)
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
end

class NewickTree # Additional methods for NewickTree class
  # returns classification of node based on taxonomy
  def createClassification(pep, exclude, taxonomy, ruleMaj)
    node = findNode(pep)
    return nil if node.nil?
    cons = node.consensusTax(taxonomy, ruleMaj)
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
    sp[name] = current.to_i if rank == "species" || rank == "no rank"
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

def getBlastHits(blastIdx, pep)
  pid = pep.full_id
  matches = []
  if f=blastIdx.seek(pep.full_id)
    f.each do |line|
      bid, match = line.chomp.split("\t")
      if (bid != pid)
        break
      else
        matches.push(match)
      end
    end
  end
  matches
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

# runs muscle to align sequences, returns alignment as row
def align(pep, blastids, database, verbose)
  STDERR << "Making alignment for " << pep.full_id << "...\n" if verbose
  homologs = fetchSeqs(blastids, database)
  hom = pep.full_id + ".hom"
  afa = pep.full_id + ".afa"
  out = File.new(hom, "w")
  out.print pep.to_s + homologs.join("\n")
  out.close
  begin
    afa = `muscle -in '#{hom}' -quiet`.gsub("\n","%%")
  rescue
    STDERR << "Error #{$!} aligning " << pep.full_id << "...\n" if verbose
  end
  File.unlink(hom)
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

# converts pointless stockholm alignment to a useful fasta one
def stockholm2Fasta(alignFile)
  afa = alignFile + ".afa"
  afaf = File.new(afa, "w")
  seqs = Hash.new
  start = false
  File.new(alignFile).each do |line|
    if (line =~ /^#|\/\//)
      start = true
      next
    end
    next if !start
    name, ali = line.split(" ")
    next if (!start || ali.nil?)
    ali.gsub!(".","-")
    seqs[name] = "" if (seqs[name].nil?)
    seqs[name] += ali += "\n"
  end
  seqs.keys.each do |name|
   afaf.printf(">%s\n%s",name,seqs[name])
 end
 afaf.close
 afa
end

# infers tree by desired method, populates db, returns tree
def infer(db, afa, pep, method, verbose)
  STDERR << "Making tree for " << pep << "...\n" if verbose
  tree = nil
  if (method == "nj")
    begin
      stock = fasta2Stockholm(afa)
      tree = NewickTree.new(`quicktree -boot 100 '#{stock}'`.tr("\n",""))
      tree.midpointRoot
    rescue
      STDERR << "Error #{$!} inferring nj tree for " << pep << "...\n" if verbose
    end
  end
  begin
    db.execute("REPLACE INTO trees VALUES(?,?)", pep, tree.to_s)
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
    db.execute("REPLACE INTO classification VALUES(#{line})")
  rescue
    STDERR << "Error #{$!} writing classification for " << pep << "...\n" if verbose
  end
end

def trimAlignment(trimFile, alignFile, maxGapFract = 0.5, exclude = nil)
  if (File.exist?(alignFile) && !File.exist?(trimFile))
    seqs = []
    badCols = []
    len = 0
    Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each do |seq|
      seq.data.tr!("\n","")
      seqs.push(seq)
    end
    seqs[0].data.length.times do |i|
      gapNum = 0
      count = 0
      seqs.each do |seq|
        next if (exclude && seq.full_id =~/#{exclude}/)
        gapNum += 1 if (seq.data[i].chr == "-" || seq.data[i].chr == "?" || seq.data[i].chr == ".")
        count += 1
      end
      badCols.push(i) if (gapNum > maxGapFract*count)
    end
    out = File.new(trimFile, "w")
    seqs.each do |seq|
      badCols.each do |col|
        seq.data[col] = "!"
      end
      seq.data.tr!("!","")
      len = seq.data.length
      out.print Bio::Sequence.auto(seq.data).to_fasta(seq.definition, 60)
    end
    out.close
    return len if len > 0
  end
  return nil
end

# back aligns dna to pep alignment and puts it in dnaAlign
def backAlign(dna, pepAlign, dnaAlign)
  pep = Hash.new
  Bio::FlatFile.new(Bio::FastaFormat, File.new(pepAlign)).each do |seq|
    pep[seq.full_id] = seq.seq
  end
  dnaAlign = File.new(dnaAlign, "w")
  Bio::FlatFile.new(Bio::FastaFormat, File.new(dna)).each do |seq|
    if (!pep[seq.full_id])
      raise "No #{seq.full_id} in #{pepAlign}\n"
    end
    dseq = ""
    j = 0
    pep[seq.full_id].length.times do |i|
      c = pep[seq.full_id][i].chr
      if (c == "-")
        dseq += "---"
      else
        dseq += seq.seq[j..j+2]
        j += 3
      end
    end
    dnaAlign.print Bio::Sequence::NA.new(dseq).to_fasta(seq.definition, 60)
  end
  dnaAlign.close
end

# calc average percent identity of an fasta alignment
def calcPercentIdent(fasta)
  pos = nil
  idents = []
  len = nil
  counts = 0
  Bio::FlatFile.new(Bio::FastaFormat, File.new(fasta)).each do |seq1|
    len = seq1.length if len.nil?
    Bio::FlatFile.new(Bio::FastaFormat, File.new(fasta)).each do |seq2|
      next if seq2.full_id == seq1.full_id
      idents.push(0)
      seq1.length.times do |i|
        idents[idents.size - 1] += 1 if (seq1.seq[i] == seq2.seq[i])
      end
    end
  end
  tot = 0
  idents.each {|ident| tot+=ident}
  avIdent = (tot * 100 / idents.size) / len
 return avIdent
end

# given a NewickTree and an alignment add ML distances
def estimateMLBranchLengths(tree, alignFile, tmpdir)
  outgroup = tree.taxa.sort.last
  tree.reroot(tree.findNode(outgroup))
  bClades = tree.clades(true)
  fasta2Phylip(alignFile, "#{tmpdir}/infile")
  tree.write("#{tmpdir}/intree")
  treepuzzle = "puzzle infile intree"
  system("cd #{tmpdir};echo  \"y\" | #{treepuzzle} > /dev/null")
  tree = NewickTree.fromFile("#{tmpdir}/intree.tree")
  tree.reroot(tree.findNode(outgroup))
  tree.addBootStrap(bClades)
  File.unlink(tmpdir+"/intree", tmpdir+"/intree.tree", tmpdir+"/infile")
  File.unlink(tmpdir+"/intree.dist", tmpdir+"/intree.puzzle")
  return tree
end

def fasta2Phylip(alignFile, phyFile) 
  seqs = Hash.new
  name = nil
  inFile = File.new(alignFile)
  inFile.each do |line|
    line.chomp!
    line.tr!("*","")
    if (line =~ /^>/)
      name = line[1..line.length].split(";").pop
      seqs[name] = ""
    else
      seqs[name] += line.gsub(".","-")
    end
  end
  inFile.close
  phy = File.new(phyFile, "w")
  lineLen = 60
  phy.printf("%d %d\n", seqs.size, seqs[name].length)
  pos = 0
  while (pos < seqs[name].length)
    seqs.keys.sort.each do |name|
      if (pos == 0)
        phy.printf("%-10s ", name)
      end
      phy.printf("%s\n", seqs[name][pos..pos + lineLen - 1])
    end
    pos += lineLen
    phy.printf("\n")
  end
  phy.close
end

def removeAA(trimFile, alignFile, aaList)
  if (File.exist?(alignFile) && !File.exist?(trimFile))
    seqs = []
    badCols = []
    Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each do |seq|
      seq.data.tr!("\n","")
      seqs.push(seq)
    end
    seqs[0].data.length.times do |i|
      bad = 0
      seqs.each do |seq|
  bad += 1 if aaList.include?(seq.data[i].chr)
      end
      badCols.push(i) if (bad == seqs.size)
    end
    out = File.new(trimFile, "w")
    seqs.each do |seq|
      badCols.each do |col|
  seq.data[col] = "!"
      end
      seq.data.tr!("!","")
      out.print Bio::Sequence.auto(seq.data).to_fasta(seq.definition, 60)
    end
    out.close
  end
end

def fasta2Nexus(alignFile, dna, nexFile = nil)
  seqs = Hash.new
  name = nil
  seqFile = File.new(alignFile)
  Bio::FlatFile.new(Bio::FastaFormat, seqFile).each do |seq|
    seqs[seq.full_id] = seq.seq.gsub("?","-").gsub(".","-")
  end
  seqFile.close
  if (dna)
    type = "NUC"
  else
    type = "PROT"
  end
  if (nexFile.nil?)
    out = STDOUT
  else
    out = File.new(nexFile, "w")
  end 
  lineLen = 40
  aLen = seqs[seqs.keys.first].size
  out.print "#NEXUS\nBEGIN DATA;\n"
  out.print "DIMENSIONS NTAX=#{seqs.size} NCHAR=#{aLen};\n"
  out.print "FORMAT DATATYPE=#{type} INTERLEAVE MISSING=-;\n"
  out.print "MATRIX\n"
  pos = 0
  while (pos < aLen)
    seqs.keys.sort.each do |name|
      out.printf("%35s ", name)
      out.printf("%s\n", seqs[name][pos..pos + lineLen - 1])
    end
    pos += lineLen 
    out.printf("\n")
  end
  out.print ";\nEND;\n"
  out.close if nexFile
end

def aliasFasta(fasta, ali, out, outgroup = nil, trim = false)
  outFile = File.new(out, "w")
  aliFile = File.new(ali, "w") if (!ali.nil?)
  aliHash = Hash.new
  orfNum = "SID0000001"
  Bio::FlatFile.new(Bio::FastaFormat, File.new(fasta)).each do |seq|
    newName = orfNum
    name = seq.definition.split(" ").first
    newName = "SID0000000" if (outgroup == name)
    aliFile.printf("%s\t%s\n", newName, seq.definition) if (ali)
    aliHash[newName] = seq.definition.tr("(),:","_")
    seq.definition = newName
    outFile.print seq
    orfNum = orfNum.succ if (outgroup != name)
  end
  outFile.close
  aliFile.close if (ali)
  if (@trim)
    trimAlignment(out+"_trim", out)
    File.unlink(out)
    File.rename(out+"_trim", out)
  end
  return aliHash
end

# do this to avoid splitting on "|"
class Bio::FastaFormat
  def full_id
    return definition.split(" ").first
  end
end

# run command on grid with project number
def qsystem(cmd, project)
  qsub = "qsub -P #{project} -e stderr -cwd -o stdout "
  system("#{qsub} \"#{cmd}\"")
end

# generate command line from trollop opts, minus unwanted options
def cmdLine(prog, opts, exclude)
    cmd = prog
    keys = opts.keys - exclude - [:help]
    keys.each do |key|
      k = key.to_s
      cmd += " --#{key} #{opts[key]}" if !k.index("_given") && opts[key]
    end
    cmd
end

# gets rid of a directory
def cleanup(dir)
  system("rm -rf #{dir}") 
end

class FlatFileDB
  def initialize(fileName, col = 0, sep = "\t", verbose = false)
    @fileName = fileName
    @file = File.new(@fileName)
    @idxfile = fileName + ".idx"
    @verbose = verbose
    buildIndex(col, sep) if !File.exists?(@idxfile) || File.mtime(@idxfile) < File.mtime(@fileName)
    readIndex(col, sep)
  end
  def buildIndex(col, sep)
    begin
      ifile = File.new(@idxfile, "w")
    rescue
      STDERR << "Permissions issue creating index file " << @idxfile << "\n" 
      exit(1)
    end
    STDERR << "Building index for " << @fileName << "...\n" if @verbose
    @idx = Hash.new
    @file.rewind
    oldfield = nil
    pos = @file.tell
    @file.each do |line| 
      fields = line.chomp.split(sep)
      ifield = fields[col]
      ifile.printf("%s\t%d\n", ifield, pos) if ifield != oldfield
      oldfield = ifield
      pos = @file.tell
    end
    ifile.close
  end
  def readIndex(col, sep)
    STDERR << "Loading index for " << @fileName << " into memory...\n" if @verbose
    @idx = Hash.new
    File.new(@idxfile).each do |line|
      field, pos = line.chomp.split(sep)
      @idx[field] = pos.to_i
    end
  end
  def seek(key)
    if (@idx[key])
      @file.seek(@idx[key])
      @file
    else
      nil
    end
  end
end