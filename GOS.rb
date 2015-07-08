# routines useful for processing GOS data

require 'Stats'

$VERBOSE = nil

#build csv file of pep info for a directory
def buildInfoTable(out, hmm, dir)
  pepAnn = Hash.new
  pepTax = Hash.new
  pepHmm = Hash.new
  pepHmmAnn = Hash.new
  STDERR.printf("Processing %s...\n", dir)
  STDERR.printf("\tLoading HMM Data...\n")
  `grep #{dir} #{hmm}`.each {|line|
    hmmId, id, sample, ann = line.chomp.split("\t")
    pepHmm[id] = hmmId.gsub(",","")
    pepHmmAnn[id] = ann.gsub(",","")
  }
  mainPep = dir + "/" + dir + ".pep"
  if !File.file?(mainPep)
    STDERR.printf("Error: no %s\n", mainPep)
    exit(1)
  end
  STDERR.printf("\tLoading main peptide file...\n")
  `grep ">" #{mainPep}`.each {|line|
    id = line.chomp.tr(">","").split(" ").first
    pepAnn[id] = ""
    pepTax[id] = "NO_TREE"
  }
  annPep = dir + "/" + "ann.pep"
  if !File.file?(annPep)
    STDERR.printf("Error: no %s\n", annPep)
    exit(1)
  end
  STDERR.printf("\tLoading annotated peptide file...\n")
  `grep ">" #{annPep}`.each {|line|
    id, ann = line.chomp.tr(">","").split(" ", 2)
    pepAnn[id] = ann.gsub(",","") if (!ann.nil?)
    pepTax[id] = ""
  }
  STDERR.printf("\tLoading taxonomy...\n")
  ranks =  ["kingdom","phylum","class","order","family","genus","species"]
  ranks.each {|level|
    Dir.glob(dir + "/*" + level + ".html").each {|html|
      tax = nil
      if (html =~/is_unresolved/)
        tax = "Mixed"
      elsif (html =~/(Contained_within[_]*|Outgroup_of[_]*)(.*)_#{level}/)
        tax = $2.gsub(",","")
      else
        STDERR.printf("Error: %s cannot be parsed at %s\n", html, level)
        exit(1)
      end
      File.new(html).each {|line|
        if (line =~/\"blast\/(.*).blastp/)
          id = $1
          pepTax[id] += tax + "; " if (!pepTax[id].include?("Mixed"))
        end
      }
    }
  }
  sample = $1 if (dir =~/(GS[^-]*)/)
  sample.gsub!(/IOLG|IOSM|IOVIR/,"")
  size = $1 if (dir =~/-([1-9|-]*(kb|KB))/)
  sample = $1 if sample == "GS" && (dir =~/(GS-[0-9]*)/)
  STDERR.printf("\tWriting data...\n")
  pepAnn.keys.sort.each {|key|
    out.printf("%s,%s,%s,%s,%s,%s,%s,%s\n",key, sample, size, 
               pepAnn[key], pepTax[key], pepHmm[key], pepHmmAnn[key], dir) 
  }
end

# return filter size based on sample name
def classifySample(sample)
  if sample =~/VIR|viral/
    return "VIR"
  elsif sample =~/LG/ || sample =~/01a/
    return 3.0
  elsif sample =~/SM/ || sample =~/01b/ || sample =~/GS-25-/ || sample =~/GOS108XLRVAL-4F-1-400_FG5HGJH01/ || sample =~/GOS108XLRVAL-4F-1-400_FG5HGJH02/ || sample =~/GOS108XLRVAL-4F-1-400_FJGGSX101/
    return 0.8
  elsif sample =~/00[b|c|d]/
    return 0.22
  else
    return 0.1
  end
end

# return site name based on sample
def siteName(sample)
  s = sample
  s = s.gsub("GS0","GS")
  s = s.gsub("GS-","GS")
  s = s.gsub("GOS", "GS")
  s = s.gsub(/^0[0-9]+-/,"")
  s = s.gsub("IOSM","")
  s = s.gsub("IOLG","") 
  s = s.gsub("IOVIR","")
  s = s.gsub("XLRVAL","")
  s = s.gsub("viral","")
  s = s.gsub("LG","") 
  return s.split("-").first
end


# load metadata
def loadMetaData(csv)
  headers = nil
  values = Hash.new
  File.new(csv).each {|line|
    fields = line.chomp.split(",")
    if (headers.nil?)
      headers = fields
    else
      values[fields.first] = Hash.new
      headers.size.times {|i|
	values[fields.first][headers[i]] = fields[i]
      }
    end
  }
  return values
end

# compute correlation of counts with metadata
def corr(hash, metadata, property)
  counts = []
  meta = []
  hash.keys.each {|key|
    if (!metadata[key].nil? && !metadata[key][property].nil?)
      if (metadata[key][property] != "NA")
	meta.push(metadata[key][property].to_f)
	counts.push(hash[key])
      end
    end
  }
  return counts.correlationWith(meta)
end

# compute taxonomy breakdown for given db and group
def clusterTaxonDist(db, kingdom)
  dist = Hash.new
  ranks = Hash.new
  counts = Hash.new
  STDERR << "Computing taxononomy distro for " << kingdom << "...\n"
  db.execute("SELECT phylum, class, ord, family, genus, cluster_num FROM classification " +
    "JOIN cluster ON cluster.seq_name = classification.seq_name WHERE kingdom = '#{kingdom}'").each do |row|
    phylum, cl, ord, family, genus, cluster_num = row
    group = phylum
    counts[cluster_num] = 0 if (!counts[cluster_num])
    counts[cluster_num] += 1
    if (kingdom == "Bacteria")
      if (genus =~/Pelagibacter/)
        group = "Pelagibacter"
      elsif (genus =~/Prochlorococcus/)
        group = "Prochlorococcus"
      elsif (genus =~/Synechococcus/)
        group = "Synechococcus"
      elsif (cl == "Alphaproteobacteria")
        group = "Other alphaproteobacteria"
      elsif (group == "Cyanobacteria")
        group = "Other cyanobacteria"
      elsif (group =~/Bacteroidetes|Chlorobi/)
        group = "Bacteroidetes/Chlorobi"
      elsif (cl =~/Betaproteobacteria/)
        group = "Betaproteobacteria"
      elsif (cl =~/Gammaproteobacteria/)
        group = "Gammaproteobacteria"
      elsif (cl =~/Deltaproteobacteria/)
        group = "Deltaproteobacteria"
      elsif (cl =~/Epsilonproteobacteria/)
        group = "Epsilonproteobacteria"
      elsif (group =~/Proteobacteria/)
        group = "Other proteobacteria"
      elsif (group =~/Rhodobacterales/)
        group = "Rhodobacterales"
      elsif (group =~/Actinobacteria/)
        group = "Actinobacteria"  
      elsif (group =~/Firmicutes/)
        group = "Firmicutes" 
      elsif (group =~/Chlamydiae|Verrucomicrobia/)
        group = "Chlamydiae/Verrucomicrobia"
      elsif (group =~/Spirochaetes/)
        group = "Spirochaetes"
      elsif (group =~/Thermotogae/)
        group = "Thermotogae"
      elsif (group =~/Planctomycetes/)
        group = "Planctomycetes"  
      elsif (group =~/Unknown/)
        group = "Unknown" 
      else
        group = "Other bacteria"
      end
    elsif (kingdom == "Eukaryota")
      group = phylum
    elsif (kingdom == "Viruses")
      group = phylum
      group = "Other viruses" if (group =~/\(phylum\)/)
    end
    next if group == "Mixed"
    group = kingdom + ": " + group
    dist[cluster_num] = Hash.new if (!dist[cluster_num])
    dist[cluster_num][group] = 0 if (!dist[cluster_num][group])
    dist[cluster_num][group] += 1
    ranks[group] = 0 if (!ranks[group])
    ranks[group] += 1
  end
  dist.keys.each do |cluster|
    sum = dist[cluster].values.reduce(:+)
    dist[cluster].keys.each do |key|
      dist[cluster][key] = ((dist[cluster][key]*1000) / sum) / 10.0
    end
  end
  return dist, ranks.keys.sort {|x, y| ranks[y] <=> ranks[x]}, counts
end


