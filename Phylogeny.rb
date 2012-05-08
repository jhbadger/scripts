require 'Newick'
require 'bio'

$VERBOSE = false

def align(alignFile, homologFile)
  if (File.exist?(homologFile) && !File.exist?(alignFile))
    system("muscle -quiet -in #{homologFile} -out #{alignFile}")
  end
end

def gblocks(trimFile, alignFile)
  ali = aliasFasta(alignFile, nil, alignFile + ".out")
  system("Gblocks #{alignFile}.out -t p >/dev/null")
  out = File.new(trimFile, "w")
  gb = alignFile + ".out-gb"
  len = nil
  Bio::FlatFile.new(Bio::FastaFormat, File.new(gb)).each {|seq|
    nosp = seq.seq.gsub(" ","")
    out.print Bio::Sequence.auto(nosp).to_fasta(ali[seq.entry_id], 60)
    len = nosp.length if len.nil?
  }
  out.close
  File.unlink(alignFile + ".out", alignFile + ".out-gb", alignFile + 
              ".out-gb.htm")
  return len
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
        next if (exclude && seq.entry_id =~/#{exclude}/)
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
    pep[seq.entry_id] = seq.seq
  end
  dnaAlign = File.new(dnaAlign, "w")
  Bio::FlatFile.new(Bio::FastaFormat, File.new(dna)).each do |seq|
    if (!pep[seq.entry_id])
      raise "No #{seq.entry_id} in #{pepAlign}\n"
    end
    dseq = ""
    j = 0
    pep[seq.entry_id].length.times do |i|
      c = pep[seq.entry_id][i].chr
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
      next if seq2.entry_id == seq1.entry_id
      idents.push(0)
      seq1.length.times {|i|
        idents[idents.size - 1] += 1 if (seq1.seq[i] == seq2.seq[i])
      }
    end
  end
  tot = 0
  idents.each {|ident| tot+=ident}
  avIdent = (tot * 100 / idents.size) / len
 return avIdent
end

def removeAA(trimFile, alignFile, aaList)
  if (File.exist?(alignFile) && !File.exist?(trimFile))
    seqs = []
    badCols = []
    Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each {|seq|
      seq.data.tr!("\n","")
      seqs.push(seq)
    }
    seqs[0].data.length.times {|i|
      bad = 0
      seqs.each {|seq|
	bad += 1 if aaList.include?(seq.data[i].chr)
      }
      badCols.push(i) if (bad == seqs.size)
    }
    out = File.new(trimFile, "w")
    seqs.each {|seq|
      badCols.each {|col|
	seq.data[col] = "!"
      }
      seq.data.tr!("!","")
      out.print Bio::Sequence.auto(seq.data).to_fasta(seq.definition, 60)
    }
    out.close
  end
end

def fasta2Nexus(alignFile, dna, nexFile = nil)
  seqs = Hash.new
  name = nil
  seqFile = File.new(alignFile)
  Bio::FlatFile.new(Bio::FastaFormat, seqFile).each {|seq|
    seqs[seq.entry_id] = seq.seq.gsub("?","-").gsub(".","-")
  }
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
    seqs.keys.sort.each {|name|
      out.printf("%35s ", name)
      out.printf("%s\n", seqs[name][pos..pos + lineLen - 1])
    }
    pos += lineLen 
    out.printf("\n")
  end
  out.print ";\nEND;\n"
  out.close if nexFile
end

def fasta2Phylip(alignFile, phyFile) 
  seqs = Hash.new
  name = nil
  inFile = File.new(alignFile)
  inFile.each {|line|
    line.chomp!
    line.tr!("*","")
    if (line =~ /^>/)
      name = line[1..line.length].split(";").pop
      seqs[name] = ""
    else
      seqs[name] += line.gsub(".","-")
    end
  }
  inFile.close
  phy = File.new(phyFile, "w")
  lineLen = 60
  phy.printf("%d %d\n", seqs.size, seqs[name].length)
  pos = 0
  while (pos < seqs[name].length)
    seqs.keys.sort.each {|name|
      if (pos == 0)
        phy.printf("%-10s ", name)
      end
      phy.printf("%s\n", seqs[name][pos..pos + lineLen - 1])
    }
    pos += lineLen
    phy.printf("\n")
  end
  phy.close
end

# converts fasta alignment into absurd PAML format
def fasta2PAML(alignFile, pamlFile) 
  seqs = []
  length = nil
  Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each {|seq|
    length = seq.seq.length
    seqs.push(Bio::Sequence::AA.new(seq.seq).to_fasta(seq.entry_id,60))
  }
  paml = File.new(pamlFile, "w")
  paml.printf(" %d %d\n", seqs.size, length)
  seqs.each {|seq|
    paml.print seq
  }
  paml.close
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

def drawTree(treeFile)
  if (File.exist?(treeFile) && !File.exist?(treeFile + ".pdf"))
    STDERR.printf("Drawing #{treeFile}...\n")
    tree = NewickTree.fromFile(treeFile)
    tree.draw(treeFile + ".pdf")
  end
end

def makeProMLTree(treeFile, alignFile, name)
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    phyFile = fasta2Phylip(alignFile, "infile")
    File.unlink("outfile") if (File.exist?("outfile"))
    File.unlink("outtree") if (File.exist?("outtree"))
    system('echo -e "Y\n" | proml')
    tree = NewickTree.fromFile("outtree")
    tree.midpointRoot
    tree.write(treeFile)
    File.unlink("outfile", "infile")
  end
end

def makeProtFastTree(treeFile, alignFile, outGroup = nil)
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    system("FastTree #{alignFile} > #{treeFile}")
    tree = NewickTree.fromFile(treeFile)
    if (outGroup.nil?)
      tree.midpointRoot
    else
      tree.reroot(tree.findNode(outGroup))
    end
    tree.write(treeFile)
  end
end

def makeNucFastTree(treeFile, alignFile, outGroup = nil)
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    system("FastTree -gtr -nt #{alignFile} > #{treeFile}")
    tree = NewickTree.fromFile(treeFile)
    if (outGroup.nil?)
      tree.midpointRoot
    else
      tree.reroot(tree.findNode(outGroup))
    end
    tree.write(treeFile)
  end
end

def makeWeighborTree(treeFile, alignFile, name)
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    phyFile = fasta2Phylip(alignFile, "infile")
    File.unlink("outfile") if (File.exist?("outfile"))
    system('echo -e "P\nY\n" | protdist')
    system("weighbor -b 18 -i outfile > outtree")
    tree = NewickTree.fromFile("outtree")
    tree.midpointRoot
    tree.write(treeFile)
    File.unlink("outtree", "infile", "outfile")
    drawTree(treeFile)
  end
end


def makePuzzleTree(treeFile, alignFile, name)
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    phyFile = fasta2Phylip(alignFile, "infile")
    File.unlink("outfile") if (File.exist?("outfile"))
    File.unlink("outtree") if (File.exist?("outtree"))
    system('echo -e "Y\n" | treepuzzle infile')
    tree = NewickTree.fromFile("infile.tree")
    tree.midpointRoot
    tree.write(treeFile)
    File.unlink("infile.tree", "infile", "infile.puzzle", "infile.dist")
    drawTree(treeFile)
  end
end

def fasta2Stockholm(stockFile, alignFile)
  stock = File.new(stockFile, "w")
  stock.printf("# STOCKHOLM 1.0\n")
  align = Hash.new 
  aSize = 0
  nSize = 0
  Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each {|seq|
    align[seq.entry_id] = seq.seq
    aSize = seq.seq.length
    nSize = seq.entry_id.size if (nSize < seq.entry_id.size)
  }
  0.step(aSize, 50) {|i|
    stock.printf("\n")
    align.keys.sort.each {|key|
      stock.printf("%-#{nSize}s %s\n", key, align[key][i..i+49]) 
    }
  }
  stock.printf("//\n")
  stock.close
  return stockFile
end

def stockholm2Fasta(fastaFile, stockFile)
  seqs = Hash.new
  File.new(stockFile).each do |line|
    if (line =~/^[A-Z|a-z|0-9]/)
      name, data = line.chomp.split(" ")
      data = data.gsub(".","-")
      if (seqs[name].nil?)
        seqs[name] = data
      else
        seqs[name] += data
      end
    end
  end
  fasta = File.new(fastaFile, "w")
  seqs.keys.each do |key|
    fasta.print ">#{key}\n#{seqs[key].gsub("*","").gsub(Regexp.new(".{1,60}"), "\\0\n")}"
  end
  fasta.close
  return fastaFile
end



def makeQuickNJTree(treeFile, alignFile, name = nil, boot = false, 
		    process = "")
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    STDERR.printf("making NJ tree for #{name}...\n") if (!name.nil?)
    stock = fasta2Stockholm(alignFile + ".stock", alignFile)
    if (boot)
      quickTree = "quicktree -boot 100"
    else
      quickTree = "quicktree"
    end
    tree = NewickTree.new(`#{quickTree} '#{stock}'`.tr("\n",""))
    tree.midpointRoot
    tree.addECnums(alignFile) if (process.index("ecdb"))
    tree.write(treeFile)
    File.unlink(stock)
  end
end

# makes NJ tree in absurd PAML format
def makeQuickPAMLTree(treeFile, alignFile)
  if (File.exist?(alignFile) && !File.exist?(treeFile))
    stock = fasta2Stockholm(alignFile + ".stock", alignFile)
    tree = NewickTree.new(`quicktree #{stock}`.tr("\n",""))
    tree.midpointRoot
    out = File.new(treeFile, "w")
    out.printf(" %d %d\n%s\n", tree.taxa.size, 1, tree.to_s)
    out.close
    File.unlink(stock)
  end
end

def aliasFasta(fasta, ali, out, outgroup = nil, trim = false)
  outFile = File.new(out, "w")
  aliFile = File.new(ali, "w") if (!ali.nil?)
  aliHash = Hash.new
  orfNum = "SID0000001"
  Bio::FlatFile.new(Bio::FastaFormat, File.new(fasta)).each {|seq|
    newName = orfNum
    name = seq.definition.split(" ").first
    newName = "SID0000000" if (outgroup == name)
    aliFile.printf("%s\t%s\n", newName, seq.definition) if (ali)
    aliHash[newName] = seq.definition.tr("(),:","_")
    seq.definition = newName
    outFile.print seq
    orfNum = orfNum.succ if (outgroup != name)
  }
  outFile.close
  aliFile.close if (ali)
  if (@trim)
    trimAlignment(out+"_trim", out)
    File.unlink(out)
    File.rename(out+"_trim", out)
  end
  return aliHash
end

def makeQuickDistTree(treeFile, distFile, aliasHash)
  if (File.exist?(distFile) && !File.exist?(treeFile))
    STDERR.printf("making NJ tree for #{distFile}...\n")
    quickTree = "quicktree"
    tree = NewickTree.new(`#{quickTree} -in m #{distFile}`.tr("\n",""))
    tree.midpointRoot
    tree.unAlias(aliasHash)
    tree.write(treeFile)
  end
end

def writeCloseEC(treeFile, eclist, neighborFile, name)
  top = nil
  tree = NewickTree.fromFile(treeFile)
  neighbor = File.new(neighborFile, "w")
  relatives = tree.relatives(name)
  if (!relatives.nil?)
    relatives.each {|list|
      allECs = []
      list.each {|relative|
        ecs = relative.split("_").last.split(" ")
        ecs.each {|ec| allECs.push(ec)}
      }
      bestEC = allECs.mostCommon
      neighbor.printf("%s\n", bestEC)
    }
    neighbor.close
  end
  return top
end

def classifyNeighbors(neighborFile, combine = nil, exclude = nil)
  groupClass = ["kingdom", "phylum", "class", "order", "family", "genus"]
  neighbors = []
  consensus = Hash.new
  File.new(neighborFile).each {|line|
    neighbors.push(line.chomp.split("\t"))
  }
  level1 = nil
  level2 = nil
  neighbors.size.times {|i|
    if (exclude.nil? || neighbors[i].grep(/#{exclude}/).empty?)
      level1, level2 = neighbors[i], neighbors[i + 1]
      break
    else
      gr = neighbors[i].grep(/#{exclude}/).first
      STDERR.printf("Excluding #{gr} as requested\n")
    end
  }
  if (neighbors[0].nil?)
    return consensus
  end
  [groupClass.size, neighbors[0].size].min.times {|i|
    group = groupClass[i]
    next if (level1.nil?)
    if (combine)
      if (level1[i] == "Mixed")
	consensus[group] = "Closest relative is unresolved at #{group} level"
      else
	consensus[group] = level1[i]
      end
    else
      if (level1[i] == "Mixed")
	consensus[group] = "Closest relative is unresolved at #{group} level"
      elsif ((!level2.nil? && !level2[i].nil?) && level1[i] != level2[i])
	consensus[group] = "Outgroup of #{level1[i]}"
      elsif (!level1[i].nil?)
	consensus[group] = "Contained within #{level1[i]}"
      end
    end
  }
  return consensus
end 

def addSpecies(qseq, treeFile, spHash)
  species = Hash.new
  tree = NewickTree.fromFile(treeFile)
  spHash.keys.each do |key|
    species[key] =  key + "__#{spHash[key].tr("():,","")}"
  end
  tree.unAlias(species)
  tree.write(treeFile)
end

def addDefinitions(qseq, treeFile, alignFile)
  definitions = Hash.new
  tree = NewickTree.fromFile(treeFile)
  Bio::FlatFile.new(Bio::FastaFormat, File.new(alignFile)).each {|seq|
    next if qseq.definition == seq.definition
    name = seq.definition.split(" ").first
    definitions[name] = seq.definition.tr('()[] ",', "_")
  }
  tree.unAlias(definitions)
  tree.write(treeFile)
end

# returns true if file is DNA/RNA 
def isDNA?(file)
  string = ""
  count = 0
  File.new(file).each {|line|
    next if (line =~ /^>/)
    string += line.chomp.upcase.gsub("-","").gsub(".","")
    count += 1
    break if (count == 100)
  }
  agtc = 1.0 * string.count("AGTCU")
  if (agtc / string.size > 0.90)
    return true
  else
    return false
  end
end
