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


$dataset = Hash.new
$dataset["GS001a"] = Hash.new
$dataset["GS001b"] = Hash.new
$dataset["GS001c"] = Hash.new
$dataset["GS00a"] = Hash.new
$dataset["GS00b"] = Hash.new
$dataset["GS00c"] = Hash.new
$dataset["GS00d"] = Hash.new
$dataset["GS002"] = Hash.new
$dataset["GS003"] = Hash.new
$dataset["GS004"] = Hash.new
$dataset["GS005"] = Hash.new
$dataset["GS006"] = Hash.new
$dataset["GS007"] = Hash.new
$dataset["GS008"] = Hash.new
$dataset["GS009"] = Hash.new
$dataset["GS010"] = Hash.new
$dataset["GS011"] = Hash.new
$dataset["GS012"] = Hash.new
$dataset["GS013"] = Hash.new
$dataset["GS014"] = Hash.new
$dataset["GS015"] = Hash.new
$dataset["GS016"] = Hash.new
$dataset["GS017"] = Hash.new
$dataset["GS018"] = Hash.new
$dataset["GS019"] = Hash.new
$dataset["GS020"] = Hash.new
$dataset["GS021"] = Hash.new
$dataset["GS022"] = Hash.new
$dataset["GS023"] = Hash.new
$dataset["GS025"] = Hash.new
$dataset["GS026"] = Hash.new
$dataset["GS027"] = Hash.new
$dataset["GS028"] = Hash.new
$dataset["GS029"] = Hash.new
$dataset["GS030"] = Hash.new
$dataset["GS031"] = Hash.new
$dataset["GS032"] = Hash.new
$dataset["GS033"] = Hash.new
$dataset["GS034"] = Hash.new
$dataset["GS035"] = Hash.new
$dataset["GS036"] = Hash.new
$dataset["GS037"] = Hash.new
$dataset["GS038"] = Hash.new
$dataset["GS039"] = Hash.new
$dataset["GS040"] = Hash.new
$dataset["GS041"] = Hash.new
$dataset["GS042"] = Hash.new
$dataset["GS043"] = Hash.new
$dataset["GS044"] = Hash.new
$dataset["GS045"] = Hash.new
$dataset["GS046"] = Hash.new
$dataset["GS047"] = Hash.new
$dataset["GS048"] = Hash.new
$dataset["GS049"] = Hash.new
$dataset["GS050"] = Hash.new
$dataset["GS051"] = Hash.new
$dataset["GS052"] = Hash.new
$dataset["GS053"] = Hash.new
$dataset["GS054"] = Hash.new
$dataset["GS055"] = Hash.new
$dataset["GS057"] = Hash.new
$dataset["GS058"] = Hash.new
$dataset["GS059"] = Hash.new
$dataset["GS060"] = Hash.new
$dataset["GS061"] = Hash.new
$dataset["GS062"] = Hash.new
$dataset["GS063"] = Hash.new
$dataset["GS064"] = Hash.new
$dataset["GS065"] = Hash.new
$dataset["GS066"] = Hash.new
$dataset["GS067"] = Hash.new
$dataset["GS068"] = Hash.new
$dataset["GS069"] = Hash.new
$dataset["GS070"] = Hash.new
$dataset["GS071"] = Hash.new
$dataset["GS072"] = Hash.new
$dataset["GS073"] = Hash.new
$dataset["GS074"] = Hash.new
$dataset["GS075"] = Hash.new
$dataset["GS076"] = Hash.new
$dataset["GS077"] = Hash.new
$dataset["GS078"] = Hash.new
$dataset["GS079"] = Hash.new
$dataset["GS080"] = Hash.new
$dataset["GS082"] = Hash.new
$dataset["GS083"] = Hash.new
$dataset["GS084"] = Hash.new
$dataset["GS086"] = Hash.new
$dataset["GS088"] = Hash.new
$dataset["GS089"] = Hash.new
$dataset["GS090"] = Hash.new
$dataset["GS091"] = Hash.new
$dataset["GS093"] = Hash.new
$dataset["GS094"] = Hash.new
$dataset["GS098"] = Hash.new
$dataset["GS099"] = Hash.new
$dataset["GS100"] = Hash.new
$dataset["GS102"] = Hash.new
$dataset["GS103"] = Hash.new
$dataset["GS108"] = Hash.new
$dataset["GS109"] = Hash.new
$dataset["GS110"] = Hash.new
$dataset["GS111"] = Hash.new
$dataset["GS112"] = Hash.new
$dataset["GS113"] = Hash.new
$dataset["GS114"] = Hash.new
$dataset["GS115"] = Hash.new
$dataset["GS116"] = Hash.new
$dataset["GS117"] = Hash.new
$dataset["GS119"] = Hash.new
$dataset["GS120"] = Hash.new
$dataset["GS121"] = Hash.new
$dataset["GS122"] = Hash.new
$dataset["GS123"] = Hash.new
$dataset["GS124"] = Hash.new
$dataset["GS125"] = Hash.new
$dataset["GS126"] = Hash.new
$dataset["GS127"] = Hash.new
$dataset["GS128"] = Hash.new
$dataset["GS130"] = Hash.new
$dataset["GS131"] = Hash.new
$dataset["GS132"] = Hash.new
$dataset["GS133"] = Hash.new
$dataset["GS134"] = Hash.new
$dataset["GS135"] = Hash.new
$dataset["GS136"] = Hash.new
$dataset["GS137"] = Hash.new
$dataset["GS138"] = Hash.new
$dataset["GS139"] = Hash.new
$dataset["GS140"] = Hash.new
$dataset["GS141"] = Hash.new
$dataset["GS142"] = Hash.new
$dataset["GS143"] = Hash.new
$dataset["GS144"] = Hash.new
$dataset["GS145"] = Hash.new
$dataset["GS146"] = Hash.new
$dataset["GS147"] = Hash.new
$dataset["GS148"] = Hash.new
$dataset["GS149"] = Hash.new
$dataset["GS200"] = Hash.new
$dataset["GS201"] = Hash.new
$dataset["GS202"] = Hash.new
$dataset["GS203"] = Hash.new
$dataset["GS204"] = Hash.new
$dataset["GS205"] = Hash.new
$dataset["GS215"] = Hash.new
$dataset["GS216"] = Hash.new
$dataset["GS217"] = Hash.new
$dataset["GS218"] = Hash.new
$dataset["GS219"] = Hash.new
$dataset["GS220"] = Hash.new
$dataset["GS221"] = Hash.new
$dataset["GS222"] = Hash.new
$dataset["GS223"] = Hash.new
$dataset["GS224"] = Hash.new
$dataset["GS226"] = Hash.new
$dataset["GS237"] = Hash.new
$dataset["GS238"] = Hash.new
$dataset["GS239"] = Hash.new
$dataset["GS240"] = Hash.new
$dataset["GS241"] = Hash.new
$dataset["GS242"] = Hash.new
$dataset["GS243"] = Hash.new
$dataset["GS244"] = Hash.new
$dataset["GS245"] = Hash.new
$dataset["GS246"] = Hash.new
$dataset["GS247"] = Hash.new
$dataset["GS249"] = Hash.new
$dataset["GS250"] = Hash.new
$dataset["GS251"] = Hash.new
$dataset["GS252"] = Hash.new
$dataset["GS253"] = Hash.new
$dataset["GS254"] = Hash.new
$dataset["GS257"] = Hash.new
$dataset["GS258"] = Hash.new
$dataset["GS259"] = Hash.new
$dataset["GS260"] = Hash.new
$dataset["GS262"] = Hash.new
$dataset["GS263"] = Hash.new
$dataset["GS264"] = Hash.new
$dataset["GS265"] = Hash.new
$dataset["GS266"] = Hash.new
$dataset["GS267"] = Hash.new
$dataset["GS268"] = Hash.new
$dataset["GS269"] = Hash.new
$dataset["GS270"] = Hash.new
$dataset["GS271"] = Hash.new
$dataset["GS272"] = Hash.new
$dataset["GS277"] = Hash.new
$dataset["GS278"] = Hash.new
$dataset["GS279"] = Hash.new
$dataset["GS280"] = Hash.new
$dataset["GS281"] = Hash.new
$dataset["GS282"] = Hash.new
$dataset["GS283"] = Hash.new
$dataset["GS284"] = Hash.new
$dataset["GS285"] = Hash.new
$dataset["GS286"] = Hash.new
$dataset["GS299"] = Hash.new
$dataset["GS300"] = Hash.new
$dataset["GS301"] = Hash.new
$dataset["GS302"] = Hash.new
$dataset["GS303"] = Hash.new
$dataset["GS304"] = Hash.new
$dataset["GS305"] = Hash.new
$dataset["GS306"] = Hash.new
$dataset["GS307"] = Hash.new
$dataset["GS308"] = Hash.new
$dataset["GS309"] = Hash.new
$dataset["GS310"] = Hash.new
$dataset["GS311"] = Hash.new
$dataset["GS312"] = Hash.new
$dataset["GS313"] = Hash.new
$dataset["GS320"] = Hash.new
$dataset["GS321"] = Hash.new
$dataset["GS322"] = Hash.new
$dataset["GS323"] = Hash.new
$dataset["GS324"] = Hash.new
$dataset["GS325"] = Hash.new
$dataset["GS326"] = Hash.new
$dataset["GS327"] = Hash.new
$dataset["GS355"] = Hash.new
$dataset["GS362"] = Hash.new
$dataset["GS364"] = Hash.new
$dataset["GS367"] = Hash.new
$dataset["GS69"] = Hash.new
$dataset["GS355"][0.1] = "355gs01um"
$dataset["GS362"][0.1] = "362gs01um"
$dataset["GS364"][0.1] = "364gs01um"
$dataset["GS367"][0.1] = "367gs01um"
$dataset["GS108"][0.1] = "GOS108XLRVAL-4F-1-400_FG5HGJH01"
$dataset["GS108"][0.1] = "GOS108XLRVAL-4F-1-400_FG5HGJH02"
$dataset["GS108"][0.1] = "GOS108XLRVAL-4F-1-400_FJGGSX101"
$dataset["GS112"][0.1] = "GOS112XLRVAL-4F-1-400_FG67BMZ01"
$dataset["GS112"][0.1] = "GOS112XLRVAL-4F-1-400_FG67BMZ02"
$dataset["GS117"][0.1] = "GOS117XLRVAL-4F-1-400_FG87RRN01"
$dataset["GS117"][0.1] = "GOS117XLRVAL-4F-1-400_FG87RRN02"
$dataset["GS117"][0.1] = "GOS117XLRVAL-4F-1-400_FJGGSX102"
$dataset["GS108"][0.1] = "GS108-4F-01-695_FWXVZTC02"
$dataset["GS110"][0.1] = "GS110-4F-01-732_FWZPJ9401"
$dataset["GS110"][3.0] = "GS110LG-4F-01-901_FWKLEGN01"
$dataset["GS117"][3.0] = "GS117LG-4F-01-693_FWZPJ9402"
$dataset["GS127"][0.1] = "GS127-0p1um-DNA-fragment-tenth-F4WYUEW02"
$dataset["GS131"][0.1] = "GS131-0p1um-DNA-fragment-tenth-F4WYUEW02"
$dataset["GS133"][0.1] = "GS133-0p1um-DNA-pairend-half-FYUA4HS01"
$dataset["GS135"][0.1] = "GS135-0p1um-DNA-fragment-tenth-F4JX21H01"
$dataset["GS141"][0.1] = "GS141-0p1um-DNA-fragment-tenth-F4JX21H01"
$dataset["GS143"][0.1] = "GS143-0p1um-DNA-fragment-tenth-F4JX21H01"
$dataset["GS148"][0.1] = "GS148-4F-01-697_FV9FT7F01"
$dataset["GS149"][0.1] = "GS149-4F-01-699_FV9FT7F02"
$dataset["GS200"][0.1] = "GS200-0p1um-DNA-fragment-tenth-F4JX21H01"
$dataset["GS299"][0.8] = "GS299-0p8um-DNA-fragment-quarter-F46EHO302"
$dataset["GS299"][3.0] = "GS299-3p0um-DNA-fragment-quarter-F46EHO302"
$dataset["GS300"][0.8] = "GS300-0p8um-DNA-fragment-quarter-F4CFQKI02"
$dataset["GS300"][3.0] = "GS300-3p0um-DNA-fragment-quarter-F4CFQKI02"
$dataset["GS301"][0.8] = "GS301-0p8um-DNA-fragment-quarter-F5HANXW01"
$dataset["GS302"][0.8] = "GS302-0p8um-DNA-fragment-quarter-F5OB7B201"
$dataset["GS302"][3.0] = "GS302-3p0um-DNA-fragment-quarter-F5OB7B201"
$dataset["GS303"][0.8] = "GS303-0p8um-DNA-fragment-quarter-F566QMO01"
$dataset["GS304"][3.0] = "GS304-3p0um-DNA-fragment-quarter-F566QMO01"
$dataset["GS305"][0.8] = "GS305-0p8um-DNA-fragment-quarter-F5HANXW01"
$dataset["GS306"][0.8] = "GS306-0p8um-DNA-fragment-quarter-F5KYWCV02"
$dataset["GS306"][3.0] = "GS306-3p0um-DNA-fragment-quarter-F5KYWCV02"
$dataset["GS307"][0.8] = "GS307-0p8um-DNA-fragment-quarter-F5BNKFK01"
$dataset["GS307"][3.0] = "GS307-3p0um-DNA-fragment-quarter-F5BNKFK01"
$dataset["GS308"][0.8] = "GS308-0p8um-DNA-fragment-quarter-F51BN4L02"
$dataset["GS309"][0.8] = "GS309-0p8um-DNA-fragment-half-F5KNHNV02"
$dataset["GS320"][0.8] = "GS320-0p8um-DNA-fragment-quarter-MIDpool-FRDLTJ302"
$dataset["GS320"][3.0] = "GS320-3p0um-DNA-fragment-quarter-MIDpool-FRDLTJ302"
$dataset["GS321"][0.8] = "GS321-0p8um-DNA-fragment-quarter-MIDpool-FRHLGZN02"
$dataset["GS321"][3.0] = "GS321-3p0um-DNA-fragment-quarter-MIDpool-FRHLGZN02"
$dataset["GS322"][0.8] = "GS322-0p8um-DNA-fragment-quarter-MIDpool-FZTCAPU02"
$dataset["GS322"][3.0] = "GS322-3p0um-DNA-fragment-quarter-MIDpool-FZTCAPU02"
$dataset["GS323"][0.8] = "GS323-0p8um-DNA-fragment-quarter-F5BNKFK02"
$dataset["GS323"][3.0] = "GS323-3p0um-DNA-fragment-quarter-F5BNKFK02"
$dataset["GS324"][0.8] = "GS324-0p8um-DNA-fragment-quarter-F3TO1RL02"
$dataset["GS324"][3.0] = "GS324-3p0um-DNA-fragment-quarter-F3TO1RL02"
$dataset["GS326"][3.0] = "GS326-3p0um-DNA-fragment-quarter-F51BN4L02"
$dataset["GS327"][0.8] = "GS327-0p8um-DNA-fragment-quarter-MIDpool-FZWZ24A02"
$dataset["GS327"][3.0] = "GS327-3p0um-DNA-fragment-quarter-MIDpool-FZWZ24A02"
$dataset["GS059"][0.1] = "GS59-0p1um-DNA-fragment-tenth-F4WYUEW02"
$dataset["GS064"][0.1] = "GS64-0p1um-DNA-fragment-tenth-F4WYUEW02"
$dataset["GS068"][0.1] = "GS68-0p1um-DNA-fragment-tenth-F4WYUEW02"
$dataset["GS100"][0.1] = "01-GS100-G-1p6-2kb_mod"
$dataset["GS108"][0.8] = "01-GS108-G-4-6kb"
$dataset["GS112"][0.1] = "01-GS112-G-4-6kb"
$dataset["GS114"][0.1] = "01-GS114-G-4-6kb"
$dataset["GS117"][0.1] = "01-GS117-G-3-5kb"
$dataset["GS122"][0.1] = "01-GS122-G-4-6kb"
$dataset["GS069"][0.1] = "01-GS69-G-3-5kb_mod"
$dataset["GS100"][0.1] = "02-GS100-G-4-6kb_mod"
$dataset["GS114"][0.1] = "02-GS114-G-2-3kb"
$dataset["GS69"][0.1] = "02-GS69-G-5-6kb_mod"
$dataset["GS114"][0.1] = "03-GS114-G-1p6-2kb"
$dataset["GS00a"][0.1] = "GS-00a-01-01-2P5KB"
$dataset["GS00b"][0.22] = "GS-00b-01-01-10P0KB"
$dataset["GS00c"][0.22] = "GS-00c-01-01-3P5KB"
$dataset["GS00d"][0.22] = "GS-00d-01-01-2P0KB"
$dataset["GS001a"][3.0] = "GS-01a-01-01-2P0KB"
$dataset["GS001b"][0.8] = "GS-01b-01-01-2P0-4P0KB"
$dataset["GS001c"][0.1] = "GS-01c-01-01-3P0KB"
$dataset["GS002"][0.1] = "GS-02-01-01-1P6KB"
$dataset["GS003"][0.1] = "GS-03-01-01-2P2KB"
$dataset["GS004"][0.1] = "GS-04-01-01-1P4-1P8KB"
$dataset["GS005"][0.1] = "GS-05-01-01-1P6KB"
$dataset["GS006"][0.1] = "GS-06-01-01-1P8KB"
$dataset["GS007"][0.1] = "GS-07-01-01-1P6KB"
$dataset["GS008"][0.1] = "GS-08-01-01-1P6KB"
$dataset["GS009"][0.1] = "GS-09-01-01-1P2-1P4KB"
$dataset["GS010"][0.1] = "GS-10-01-01-1P5KB"
$dataset["GS102"][0.1] = "GS102-G-01-3-4kb_mod"
$dataset["GS103"][0.1] = "GS103-G-01-4-6kb_mod"
$dataset["GS108"][0.1] = "GS108-G-02-4-6KB"
$dataset["GS109"][0.1] = "GS109-G-01-3-4KB"
$dataset["GS011"][0.1] = "GS-11-01-01-1P3-1P8KB"
$dataset["GS110"][0.1] = "GS110-G-01-4-6kb"
$dataset["GS110"][0.1] = "GS110-G-02-4-6KB"
$dataset["GS111"][0.1] = "GS111-G-02-4-6KB"
$dataset["GS112"][0.1] = "GS112-G-02-4-6KB"
$dataset["GS113"][0.1] = "GS113-G-01-4-6KB"
$dataset["GS113"][0.1] = "GS113-G-02-6-8KB"
$dataset["GS114"][0.1] = "GS114-G-02-4-6KB"
$dataset["GS114"][0.1] = "GS114-G-03-2-3KB"
$dataset["GS114"][0.1] = "GS114-G-04-1P6-2KB"
$dataset["GS115"][0.1] = "GS115-G-01-4-6KB"
$dataset["GS115"][0.1] = "GS115-G-02-4-6KB"
$dataset["GS116"][0.1] = "GS116-G-01-4-6KB"
$dataset["GS116"][0.1] = "GS116-G-02-4-6KB"
$dataset["GS117"][0.1] = "GS117-G-02-3-5KB"
$dataset["GS119"][0.1] = "GS119-G-01-3-4KB"
$dataset["GS012"][0.1] = "GS-12-01-01-1P3-1P6KB"
$dataset["GS012"][0.1] = "GS-12-02-01-1P6-2P0KB"
$dataset["GS120"][0.1] = "GS120-G-02-4-6KB"
$dataset["GS121"][0.1] = "GS121-G-01-4-6KB"
$dataset["GS122"][0.1] = "GS122-G-02-4-6KB"
$dataset["GS123"][0.1] = "GS123-G-01-4-6KB"
$dataset["GS124"][0.1] = "GS124-G-01-4-6kb_mod"
$dataset["GS124"][0.1] = "GS124-G-02-6-8kb_mod"
$dataset["GS125"][0.1] = "GS125-G-01-4-6kb_mod"
$dataset["GS126"][0.1] = "GS126-G-01-4-6kb_mod"
$dataset["GS126"][0.1] = "GS126-G-02-6-8kb_mod"
$dataset["GS128"][0.1] = "GS128-G-01-6-8kb_mod"
$dataset["GS013"][0.1] = "GS-13-02-01-1P6-2P0KB"
$dataset["GS130"][0.1] = "GS130-G-01-6-8kb_mod"
$dataset["GS132"][0.1] = "GS132-G-01-6-8kb_mod"
$dataset["GS134"][0.1] = "GS134-G-01-6-8kb_mod"
$dataset["GS136"][0.1] = "GS136-G-01-6-8kb_mod"
$dataset["GS137"][0.1] = "GS137-G-01-6-8KB_mod"
$dataset["GS138"][0.1] = "GS138-G-01-6-8kb_mod"
$dataset["GS139"][0.1] = "GS139-G-01-4-6KB_mod"
$dataset["GS014"][0.1] = "GS-14-01-01-1P6-2P0KB"
$dataset["GS140"][0.1] = "GS140-G-01-6-8kb_mod"
$dataset["GS142"][0.1] = "GS142-G-01-6-8kb_mod"
$dataset["GS142"][0.1] = "GS142-G-02-8-10kb_mod"
$dataset["GS144"][0.1] = "GS144-G-01-4-6kb_mod"
$dataset["GS145"][0.1] = "GS145-G-01-4-6KB_mod"
$dataset["GS146"][0.1] = "GS146-G-01-8-10kb_mod"
$dataset["GS147"][0.1] = "GS147-G-01-8-10kb_mod"
$dataset["GS148"][0.1] = "GS148-G-01-4-6kb"
$dataset["GS148"][0.1] = "GS148-G-02-6-8kb"
$dataset["GS149"][0.1] = "GS149-G-01-4-6kb"
$dataset["GS149"][0.1] = "GS149-G-02-6-8kb"
$dataset["GS015"][0.1] = "GS-15-01-01-1P3-1P6KB"
$dataset["GS016"][0.1] = "GS-16-01-01-1P3-1P6KB"
$dataset["GS017"][0.1] = "GS-17-01-01-1P8-2P2KB"
$dataset["GS018"][0.1] = "GS-18-01-01-1P6-2P0KB"
$dataset["GS019"][0.1] = "GS-19-01-01-1P3-1P6KB"
$dataset["GS020"][0.1] = "GS-20-01-01-1P3-1P8KB"
$dataset["GS201"][0.1] = "GS201-G-01-4-6KB_mod"
$dataset["GS202"][0.1] = "GS202-G-01-4-6KB_mod"
$dataset["GS203"][0.1] = "GS203-G-01-4-6KB"
$dataset["GS204"][0.1] = "GS204-G-01-4-6KB_mod"
$dataset["GS205"][0.1] = "GS205-G-01-3-4KB_mod"
$dataset["GS021"][0.1] = "GS-21-01-01-1P3-1P8KB"
$dataset["GS215"][0.1] = "GS215-G-01-3-4KB_mod"
$dataset["GS215"][0.1] = "GS215VA-G-02-3-4KB_mod"
$dataset["GS215"][0.1] = "GS215VU-G-01-3-4KB_mod"
$dataset["GS216"][0.1] = "GS216-G-01-4-6KB_mod"
$dataset["GS217"][0.1] = "GS217-G-01-4-6KB"
$dataset["GS218"][0.1] = "GS218-G-01-4-6KB_mod"
$dataset["GS219"][0.1] = "GS219-G-01-3-4KB_mod"
$dataset["GS022"][0.1] = "GS-22-01-01-1P3-1P6KB"
$dataset["GS220"][0.1] = "GS220-G-01-3-4KB"
$dataset["GS221"][0.1] = "GS221-G-01-3-4KB"
$dataset["GS221"][0.1] = "GS221-G-01-4-6KB"
$dataset["GS222"][0.1] = "GS222-G-01-4-6KB_mod"
$dataset["GS223"][0.1] = "GS223-G-01-4-6KB_mod"
$dataset["GS224"][0.1] = "GS224-G-01-4-6KB"
$dataset["GS226"][0.1] = "GS226-G-01-4-6KB"
$dataset["GS023"][0.1] = "GS-23-01-01-1P8-2P2KB"
$dataset["GS237"][0.1] = "GS237-G-01-3-4KB_mod"
$dataset["GS238"][0.1] = "GS238-G-01-4-6KB_mod"
$dataset["GS239"][0.1] = "GS239-G-01-4-6KB_mod"
$dataset["GS240"][0.1] = "GS240-G-01-3-4KB_mod"
$dataset["GS241"][0.1] = "GS241-G-01-4-6KB_mod"
$dataset["GS242"][0.1] = "GS242-G-01-3-4KB_mod"
$dataset["GS243"][0.1] = "GS243-G-01-4-6KB_mod"
$dataset["GS244"][0.1] = "GS244-G-01-4-6KB_mod"
$dataset["GS245"][0.8] = "GS245-G-01-4-6KB_mod"
$dataset["GS246"][0.1] = "GS246-G-01-4-6KB_mod"
$dataset["GS247"][0.1] = "GS247-G-01-4-6KB_mod"
$dataset["GS249"][0.1] = "GS249-G-01-2-3KB_mod"
$dataset["GS025"][0.8] = "GS-25-01-01-1P8-2P2KB"
$dataset["GS250"][0.1] = "GS250-G-01-4-6KB_mod"
$dataset["GS251"][0.1] = "GS251-G-01-4-6KB"
$dataset["GS252"][0.1] = "GS252-G-01-3-4KB"
$dataset["GS253"][0.1] = "GS253-G-01-4-6KB"
$dataset["GS254"][0.1] = "GS254-G-01-4-6KB"
$dataset["GS257"][0.1] = "GS257-G-01-4-6KB"
$dataset["GS258"][0.1] = "GS258-G-01-4-6KB"
$dataset["GS259"][0.1] = "GS259-G-01-4-6KB"
$dataset["GS026"][0.1] = "GS-26-01-01-1P3-1P6KB"
$dataset["GS026"][0.1] = "GS-26-02-01-1P0-1P3KB"
$dataset["GS260"][0.1] = "GS260-G-01-4-6KB"
$dataset["GS262"][0.1] = "GS262-G-01-4-6KB_mod"
$dataset["GS263"][0.1] = "GS263-G-01-4-6KB"
$dataset["GS264"][0.1] = "GS264-G-02-4-6KB_mod"
$dataset["GS265"][0.1] = "GS265-G-01-4-6KB"
$dataset["GS266"][0.1] = "GS266-G-01-4-6KB"
$dataset["GS267"][0.1] = "GS267-G-01-4-6KB"
$dataset["GS268"][0.1] = "GS268-G-01-4-6KB"
$dataset["GS269"][0.1] = "GS269-G-01-3-4KB"
$dataset["GS027"][0.1] = "GS-27-01-01-1P8-2P0KB"
$dataset["GS270"][0.1] = "GS270-G-01-4-6KB"
$dataset["GS271"][0.1] = "GS271-G-01-4-6KB"
$dataset["GS272"][0.1] = "GS272-G-01-4-6KB"
$dataset["GS277"][0.1] = "GS277-G-01-4-6KB"
$dataset["GS278"][0.1] = "GS278-G-01-4-6KB"
$dataset["GS279"][0.1] = "GS279-G-01-4-6KB"
$dataset["GS028"][0.1] = "GS-28-01-01-1P6-2P0KB"
$dataset["GS280"][0.1] = "GS280-G-01-4-6KB"
$dataset["GS281"][0.1] = "GS281-G-01-4-6KB"
$dataset["GS282"][0.1] = "GS282-G-01-4-6KB"
$dataset["GS283"][0.1] = "GS283-G-01-4-6KB"
$dataset["GS284"][0.1] = "GS284-G-01-4-6KB"
$dataset["GS285"][0.1] = "GS285-G-01-4-6KB"
$dataset["GS286"][0.1] = "GS286-G-01-4-6KB"
$dataset["GS029"][0.1] = "GS-29-01-01-1P0-1P3KB"
$dataset["GS299"][0.1] = "GS299-G-01-4-6KB_mod"
$dataset["GS030"][0.1] = "GS-30-01-01-1P3-1P6KB"
$dataset["GS030"][0.1] = "GS-30-02-01-1P0-1P3KB"
$dataset["GS300"][0.1] = "GS300-G-01-3-4KB_mod"
$dataset["GS301"][0.1] = "GS301-G-01-4-6KB_mod"
$dataset["GS302"][0.1] = "GS302-G-01-4-6KB_mod"
$dataset["GS305"][0.1] = "GS305-G-01-6-8KB_mod"
$dataset["GS306"][0.1] = "GS306-G-01-3-4KB_mod"
$dataset["GS307"][0.1] = "GS307-G-01-3-4KB_mod"
$dataset["GS308"][0.1] = "GS308-G-01-4-6KB_mod"
$dataset["GS309"][0.1] = "GS309-G-01-4-6KB_mod"
$dataset["GS031"][0.1] = "GS-31-01-01-1P3-1P8KB"
$dataset["GS310"][0.1] = "GS310-G-01-4-6KB"
$dataset["GS311"][0.1] = "GS311-G-01-4-6KB"
$dataset["GS312"][0.1] = "GS312-G-01-4-6KB"
$dataset["GS313"][0.1] = "GS313-G-01-4-6KB"
$dataset["GS032"][0.1] = "GS-32-01-01-1P3-1P6KB"
$dataset["GS032"][0.1] = "GS-32-02-01-1P0-1P3KB"
$dataset["GS320"][0.1] = "GS320-G-01-4-6KB_mod"
$dataset["GS321"][0.1] = "GS321-G-01-4-6KB_mod"
$dataset["GS322"][0.1] = "GS322-G-01-4-6KB_mod"
$dataset["GS323"][0.1] = "GS323-G-01-4-6KB_mod"
$dataset["GS324"][0.1] = "GS324-G-01-4-6KB"
$dataset["GS325"][0.1] = "GS325-G-01-4-6KB_mod"
$dataset["GS326"][0.1] = "GS326-G-01-4-6KB_mod"
$dataset["GS327"][0.1] = "GS327-G-01-4-6KB_mod"
$dataset["GS033"][0.1] = "GS-33-01-01-1P3-1P8KB"
$dataset["GS033"][0.1] = "GS33-F-01-40kb"
$dataset["GS034"][0.1] = "GS-34-01-01-1P5-1P8KB"
$dataset["GS035"][0.1] = "GS-35-01-01-1P5KB"
$dataset["GS036"][0.1] = "GS-36-01-01-2P2KB"
$dataset["GS037"][0.1] = "GS-37-01-01-1P5KB"
$dataset["GS038"][0.1] = "GS-38-01-01-1P5KB_mod"
$dataset["GS039"][0.1] = "GS-39-01-01-2P0KB_mod"
$dataset["GS040"][0.1] = "GS-40-01-01-1P6-2P0KB_mod"
$dataset["GS041"][0.1] = "GS-41-01-01-1P6-2P0KB_mod"
$dataset["GS042"][0.1] = "GS-42-01-01-1P6-2P0KB_mod"
$dataset["GS043"][0.1] = "GS-43-01-01-1P6-2P0KB_mod"
$dataset["GS044"][0.1] = "GS-44-01-01-1P6-2P0KB_mod"
$dataset["GS045"][0.1] = "GS-45-01-01-1P6-2P0KB_mod"
$dataset["GS046"][0.1] = "GS-46-01-01-1P0-1P3KB_mod"
$dataset["GS047"][0.1] = "GS-47-01-01-1P6-2P0KB"
$dataset["GS048"][0.1] = "GS-48-01-01-1P6-2P0KB"
$dataset["GS048"][0.8] = "GSIOSM048-G-01-8-10KB"
$dataset["GS048"][3.0] = "GSIOLG48-G-01-8-10KB"
$dataset["GS049"][0.1] = "GS-49-01-01-1P6-2P0KB"
$dataset["GS050"][0.1] = "GS-50-01-01-1P3-1P6KB_mod"
$dataset["GS051"][0.1] = "GS-51-01-01-1P6-2P0KB"
$dataset["GS052"][0.1] = "GS52-G-01-2-3KB_mod"
$dataset["GS053"][0.1] = "GS53-G-01-3-4KB_mod"
$dataset["GS054"][0.1] = "GS54-G-01-4-6KB"
$dataset["GS055"][0.1] = "GS55-G-01-4-6KB_mod"
$dataset["GS057"][0.1] = "GS57-G-02-4-6KB"
$dataset["GS058"][0.1] = "GS58-F-01-40kb_mod"
$dataset["GS058"][0.1] = "GS58-F-02-40kb_mod"
$dataset["GS058"][0.1] = "GS58FPL-D-01-0-40KB_mod"
$dataset["GS058"][0.1] = "GS58-G-01-8-10kb_mod"
$dataset["GS060"][0.1] = "GS60-G-01-4-6KB"
$dataset["GS061"][0.1] = "GS61-G-01-6-8KB_mod"
$dataset["GS062"][0.1] = "GS62-G-01-5kb_mod"
$dataset["GS063"][0.1] = "GS63-G-01-4-6KB"
$dataset["GS065"][0.8] = "GS65-G-01-10-12KB_mod"
$dataset["GS066"][0.1] = "GS66-G-01-10-12KB_mod"
$dataset["GS067"][0.1] = "GS67-G-01-6-8KB_mod"
$dataset["GS069"][0.1] = "GS69-F-01-40kb_mod"
$dataset["GS069"][0.1] = "GS69-F-02-40kb_mod"
$dataset["GS070"][0.1] = "GS70-G-01-8-10kb_mod"
$dataset["GS071"][0.1] = "GS71-G-01-4-6KB_mod"
$dataset["GS072"][0.1] = "GS72-G-01-10-12kb_mod"
$dataset["GS073"][0.1] = "GS73-G-01-6-8KB_mod"
$dataset["GS074"][0.1] = "GS74-G-01-4-6KB_mod"
$dataset["GS075"][0.1] = "GS75-G-01-3-4KB_mod"
$dataset["GS076"][0.1] = "GS76-G-01-10-12KB_mod"
$dataset["GS077"][0.1] = "GS77-G-01-4-6KB_mod"
$dataset["GS078"][0.1] = "GS78-G-01-3-5KB_mod"
$dataset["GS079"][0.1] = "GS79-G-01-4-6KB_mod"
$dataset["GS080"][0.1] = "GS80-G-01-10-12KB_mod"
$dataset["GS082"][0.1] = "GS82-G-01-4-6KB_mod"
$dataset["GS083"][0.1] = "GS83-G-01-4-6KB_mod"
$dataset["GS084"][0.1] = "GS83-G-02-3-4KB"
$dataset["GS084"][0.1] = "GS84-02-2-3kb_mod"
$dataset["GS084"][0.1] = "GS84-F-01-40kb_mod"
$dataset["GS084"][0.1] = "GS84-F-02-40kb_mod"
$dataset["GS084"][0.1] = "GS84-G-01-1-2kb_mod"
$dataset["GS084"][0.1] = "GS84-G-03-1kb_mod"
$dataset["GS084"][0.1] = "GS84-G-04-2kb_mod"
$dataset["GS086"][0.1] = "GS86-G-01-10-12KB_mod"
$dataset["GS088"][0.1] = "GS88-F-01-40kb_mod"
$dataset["GS088"][0.1] = "GS88-F-02-40kb_mod"
$dataset["GS088"][0.1] = "GS88-G-01-3kb_mod"
$dataset["GS088"][0.1] = "GS88-G-02-3KB"
$dataset["GS088"][0.1] = "GS88-G-03-3KB_mod"
$dataset["GS089"][0.1] = "GS89-F-01-40kb_mod"
$dataset["GS089"][0.1] = "GS89-F-02-40kb_mod"
$dataset["GS089"][0.1] = "GS89-G-01-3kb_mod"
$dataset["GS090"][0.1] = "GS90-G-01-10-12KB_mod"
$dataset["GS091"][0.1] = "GS91-G-01-6-8KB_mod"
$dataset["GS093"][0.1] = "GS93-G-01-8-10KB_mod"
$dataset["GS094"][0.1] = "GS94-G-01-4-6KB_mod"
$dataset["GS098"][0.1] = "GS98-G-01-3-4KB_mod"
$dataset["GS099"][0.1] = "GS99-G-01-4-6kb_mod"
$dataset["GS099"][0.1] = "GS99-G-02-2-3KB_mod"
$dataset["GS099"][0.1] = "GS99-G-03-1P6-2KB_mod"
$dataset["GS108"][3.0] = "GSIOLG108-G-01-2-4KB"
$dataset["GS110"][3.0] = "GSIOLG110-G-01-3-4KB"
$dataset["GS112"][3.0] = "GSIOLG112-G-01-3-4KB"
$dataset["GS114"][3.0] = "GSIOLG114-G-01-4-6KB_mod"
$dataset["GS117"][3.0] = "GSIOLG117-G-01-3-4KB"
$dataset["GS122"][3.0] = "GSIOLG122-G-01-3-4KB"
$dataset["GS108"][0.8] = "GSIOSM108-G-01-3-4KB"
$dataset["GS110"][0.8] = "GSIOSM110-G-01-3-4KB"
$dataset["GS112"][0.8] = "GSIOSM112-G-01-3-4KB"
$dataset["GS114"][0.8] = "GSIOSM114-G-01-4-6KB"
$dataset["GS117"][0.8] = "GSIOSM117-G-01-3-4KB"
$dataset["GS122"][0.8] = "GSIOSM122-G-01-3-4KB"

$genomeEquiv=Hash.new
$genomeEquiv["GS00b"]=Hash.new
$genomeEquiv["GS00c"]=Hash.new
$genomeEquiv["GS00d"]=Hash.new
$genomeEquiv["GS001b"]=Hash.new
$genomeEquiv["GS001c"]=Hash.new
$genomeEquiv["GS002"]=Hash.new
$genomeEquiv["GS003"]=Hash.new
$genomeEquiv["GS004"]=Hash.new
$genomeEquiv["GS005"]=Hash.new
$genomeEquiv["GS006"]=Hash.new
$genomeEquiv["GS007"]=Hash.new
$genomeEquiv["GS008"]=Hash.new
$genomeEquiv["GS009"]=Hash.new
$genomeEquiv["GS010"]=Hash.new
$genomeEquiv["GS011"]=Hash.new
$genomeEquiv["GS012"]=Hash.new
$genomeEquiv["GS013"]=Hash.new
$genomeEquiv["GS014"]=Hash.new
$genomeEquiv["GS015"]=Hash.new
$genomeEquiv["GS016"]=Hash.new
$genomeEquiv["GS017"]=Hash.new
$genomeEquiv["GS018"]=Hash.new
$genomeEquiv["GS019"]=Hash.new
$genomeEquiv["GS021"]=Hash.new
$genomeEquiv["GS022"]=Hash.new
$genomeEquiv["GS023"]=Hash.new
$genomeEquiv["GS026"]=Hash.new
$genomeEquiv["GS027"]=Hash.new
$genomeEquiv["GS028"]=Hash.new
$genomeEquiv["GS029"]=Hash.new
$genomeEquiv["GS030"]=Hash.new
$genomeEquiv["GS031"]=Hash.new
$genomeEquiv["GS034"]=Hash.new
$genomeEquiv["GS035"]=Hash.new
$genomeEquiv["GS036"]=Hash.new
$genomeEquiv["GS037"]=Hash.new
$genomeEquiv["GS047"]=Hash.new
$genomeEquiv["GS049"]=Hash.new
$genomeEquiv["GS051"]=Hash.new
$genomeEquiv["GS109"]=Hash.new
$genomeEquiv["GS111"]=Hash.new
$genomeEquiv["GS113"]=Hash.new
$genomeEquiv["GS114"]=Hash.new
$genomeEquiv["GS115"]=Hash.new
$genomeEquiv["GS116"]=Hash.new
$genomeEquiv["GS119"]=Hash.new
$genomeEquiv["GS120"]=Hash.new
$genomeEquiv["GS123"]=Hash.new
$genomeEquiv["GS148"]=Hash.new
$genomeEquiv["GS149"]=Hash.new
$genomeEquiv["GS048"]=Hash.new
$genomeEquiv["GS108"]=Hash.new
$genomeEquiv["GS110"]=Hash.new
$genomeEquiv["GS112"]=Hash.new
$genomeEquiv["GS117"]=Hash.new
$genomeEquiv["GS122"]=Hash.new
$genomeEquiv["GS025"]=Hash.new
$genomeEquiv["GS001a"]=Hash.new
$genomeEquiv["GS00b"][0.22]=159.210985
$genomeEquiv["GS00c"][0.22]=167.431308
$genomeEquiv["GS00d"][0.22]=169.785216
$genomeEquiv["GS001b"][0.8]=19.037984
$genomeEquiv["GS001c"][0.1]=55.568896
$genomeEquiv["GS002"][0.1]=47.949944
$genomeEquiv["GS003"][0.1]=22.834564
$genomeEquiv["GS004"][0.1]=24.534753
$genomeEquiv["GS005"][0.1]=19.816728
$genomeEquiv["GS006"][0.1]=28.225878
$genomeEquiv["GS007"][0.1]=18.922631
$genomeEquiv["GS008"][0.1]=65.235233
$genomeEquiv["GS009"][0.1]=34.566327
$genomeEquiv["GS010"][0.1]=40.672082
$genomeEquiv["GS011"][0.1]=65.230282
$genomeEquiv["GS012"][0.1]=60.909534
$genomeEquiv["GS013"][0.1]=39.725499
$genomeEquiv["GS014"][0.1]=61.577814
$genomeEquiv["GS015"][0.1]=67.520127
$genomeEquiv["GS016"][0.1]=59.127882
$genomeEquiv["GS017"][0.1]=119.019226
$genomeEquiv["GS018"][0.1]=65.744026
$genomeEquiv["GS019"][0.1]=68.327888
$genomeEquiv["GS021"][0.1]=63.397924
$genomeEquiv["GS022"][0.1]=62.817992
$genomeEquiv["GS023"][0.1]=67.519106
$genomeEquiv["GS026"][0.1]=54.221367
$genomeEquiv["GS027"][0.1]=114.755242
$genomeEquiv["GS028"][0.1]=97.852425
$genomeEquiv["GS029"][0.1]=70.157047
$genomeEquiv["GS030"][0.1]=216.357411
$genomeEquiv["GS031"][0.1]=242.024775
$genomeEquiv["GS034"][0.1]=68.711675
$genomeEquiv["GS035"][0.1]=61.616293
$genomeEquiv["GS036"][0.1]=32.463854
$genomeEquiv["GS037"][0.1]=36.256087
$genomeEquiv["GS047"][0.1]=32.839701
$genomeEquiv["GS049"][0.1]=55.115866
$genomeEquiv["GS051"][0.1]=56.845883
$genomeEquiv["GS109"][0.1]=33.110194
$genomeEquiv["GS111"][0.1]=31.836716
$genomeEquiv["GS113"][0.1]=58.945258
$genomeEquiv["GS114"][0.1]=180.896572
$genomeEquiv["GS115"][0.1]=30.208242
$genomeEquiv["GS116"][0.1]=33.675949
$genomeEquiv["GS119"][0.1]=35.167868
$genomeEquiv["GS120"][0.1]=30.085377
$genomeEquiv["GS123"][0.1]=63.477831
$genomeEquiv["GS148"][0.1]=54.126628
$genomeEquiv["GS149"][0.1]=66.03469
$genomeEquiv["GS048"][0.1]=42.609171
$genomeEquiv["GS108"][0.1]=29.25399
$genomeEquiv["GS110"][0.1]=52.918064
$genomeEquiv["GS112"][0.1]=52.798134
$genomeEquiv["GS117"][0.1]=171.489805
$genomeEquiv["GS122"][0.1]=56.210333
$genomeEquiv["GS025"][0.8]=18.441473
$genomeEquiv["GS048"][0.8]=13.933269
$genomeEquiv["GS108"][0.8]=16.417576
$genomeEquiv["GS110"][0.8]=12.233203
$genomeEquiv["GS112"][0.8]=10.071848
$genomeEquiv["GS117"][0.8]=13.931343
$genomeEquiv["GS122"][0.8]=7.613603
$genomeEquiv["GS001a"][3.0]=14.863094
$genomeEquiv["GS048"][3.0]=1.812964
$genomeEquiv["GS108"][3.0]=8.088258
$genomeEquiv["GS110"][3.0]=3.03309
$genomeEquiv["GS112"][3.0]=2.876536
$genomeEquiv["GS117"][3.0]=11.022205
$genomeEquiv["GS122"][3.0]=1.904963

