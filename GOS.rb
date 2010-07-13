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
  if sample =~/VIR/
    return "VIR"
  elsif sample =~/LG/ || sample =~/01a/
    return 3.0
  elsif sample =~/SM/ || sample =~/01b/ || sample =~/GS-25-/
    return 0.8
  elsif sample =~/00[b|c|d]/
    return 0.22
  else
    return 0.1
  end
end



# compute counts of peptides in different samples and filters
def printSampleCounts(table)
  counts = Hash.new
  File.new(table).each {|line|
    fields = line.chomp.split(",") 
    site = fields[1] 
    tax = fields[4] 
    dir = fields[7]
    filter = classifySample(dir)
    counts[site] = Hash.new if counts[site].nil?
    counts[site][filter] = 0 if counts[site][filter].nil?
    counts[site][filter] += 1
  }
  counts.keys.sort.each {|sample|
    printf("sampleCounts[\"%s\"] = Hash.new\n", sample)
    counts[sample].keys.sort {|x,y| y.to_s <=> x.to_s}.each {|filter|
      printf("sampleCounts[\"%s\"][%s] = %d\n", sample, filter, counts[sample][filter])
    }
  }
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

$sampleCounts = Hash.new
$sampleCounts["GS108"] = Hash.new
$sampleCounts["GS108"]["VIR"] = 1642
$sampleCounts["GS108"][3.0] = 51974
$sampleCounts["GS108"][0.8] = 57484
$sampleCounts["GS108"][0.1] = 76521
$sampleCounts["GS109"] = Hash.new
$sampleCounts["GS109"][0.1] = 82440
$sampleCounts["GS110"] = Hash.new
$sampleCounts["GS110"]["VIR"] = 28165
$sampleCounts["GS110"][3.0] = 37007
$sampleCounts["GS110"][0.8] = 53653
$sampleCounts["GS110"][0.1] = 133104
$sampleCounts["GS111"] = Hash.new
$sampleCounts["GS111"][0.1] = 80207
$sampleCounts["GS112"] = Hash.new
$sampleCounts["GS112"]["VIR"] = 62095
$sampleCounts["GS112"][3.0] = 37243
$sampleCounts["GS112"][0.8] = 54682
$sampleCounts["GS112"][0.1] = 133753
$sampleCounts["GS113"] = Hash.new
$sampleCounts["GS113"][0.1] = 149292
$sampleCounts["GS114"] = Hash.new
$sampleCounts["GS114"][0.1] = 464694
$sampleCounts["GS115"] = Hash.new
$sampleCounts["GS115"][0.1] = 84015
$sampleCounts["GS116"] = Hash.new
$sampleCounts["GS116"][0.1] = 82174
$sampleCounts["GS117"] = Hash.new
$sampleCounts["GS117"]["VIR"] = 59320
$sampleCounts["GS117"][3.0] = 55582
$sampleCounts["GS117"][0.8] = 59839
$sampleCounts["GS117"][0.1] = 449479
$sampleCounts["GS119"] = Hash.new
$sampleCounts["GS119"][0.1] = 84430
$sampleCounts["GS120"] = Hash.new
$sampleCounts["GS120"][0.1] = 71993
$sampleCounts["GS121"] = Hash.new
$sampleCounts["GS121"][0.1] = 153121
$sampleCounts["GS122"] = Hash.new
$sampleCounts["GS122"]["VIR"] = 61580
$sampleCounts["GS122"][3.0] = 29827
$sampleCounts["GS122"][0.8] = 47329
$sampleCounts["GS122"][0.1] = 141516
$sampleCounts["GS123"] = Hash.new
$sampleCounts["GS123"][0.1] = 144005
$sampleCounts["GS148"] = Hash.new
$sampleCounts["GS148"][0.1] = 134498
$sampleCounts["GS149"] = Hash.new
$sampleCounts["GS149"][0.1] = 144505
$sampleCounts["GS48"] = Hash.new
$sampleCounts["GS48"][3.0] = 9751
$sampleCounts["GS48"][0.8] = 54196
$sampleCounts["GS48"][0.1] = 115843
$sampleCounts["GS49"] = Hash.new
$sampleCounts["GS49"][0.1] = 121429

