#db wrapper
require 'ZFile'

# Exception raised when dataset not found
class DatasetNotFoundError < RuntimeError
end

class DBwrapper 
  def close
    @connection.close if @connection
    @connected = false
  end
  def query(sql)
    errors = 0
    begin
      connect if (!@connected)
      @query = @connection.query(sql)
    rescue
      10.times {|errors|
        if (errors < 10)
          STDERR.printf("%s\n", sql)
          STDERR.printf("Error %s: Reconnecting...\n", $!)
          errors += 1
          connect
        else
          STDERR.printf("Giving up after %s\n", $!)
        end
      }
    end
    return self
  end
  def each
    @query.each {|row|
      yield row
    }
  end

  # return hash of apis db objects (mysql or sqlite) depending on source string
  def  self.dbsFromSources(sources, user, password)
    dbs = Hash.new
    sources.each {|source|
      if (source =~/.db/) # sqlite
        dbs[File.basename(source, ".db")] = SQLite.new(db)
      elsif (source =~/::/)
        host, dbname = source.split("::")
        if (!dbs[dbname])
          dbs[dbname] = MyDB.new(host, dbname, user, password)
          dbs[dbname].close
        end
      else
        host = source
        db = MyDB.new(host, "", user, password)
        db.query("SHOW DATABASES").each {|row|
          dbname = row[0]
          if (dbname =~/_apis/ && !dbs[dbname])
            dbs[dbname] = MyDB.new(host, dbname, user, password)
            dbs[dbname].close
          end
        }
        db.close
      end
    }
    return dbs
  end

  # return hash of human readable dataset names from metadata table where it exists
  def self.populateMetaName(dbs)
    metaName = Hash.new
    dbs.keys.each do |dbname|
      metaName[dbname] = Hash.new
      table = dbs[dbname].get("SHOW TABLES WHERE tables_in_#{dbname} = 'metadata'").to_s
      if (table != "")
        dbs[dbname].query("SELECT dataset, value FROM metadata WHERE prop='name'").each do |row|
          metaName[dbname][row.first] = row.last
        end
        dbs[dbname].query("SELECT dataset, value FROM metadata WHERE prop='location'").each do |row|
          metaName[dbname][row.first] += " " + row.last
        end
      end
    end
    return metaName
  end
  
  def count(condition)
    connect if (!@connected)
    query = "SELECT count(*) FROM #{condition}"
    return query(query).each {|row|
      return row[0].to_i
     }
  end
  
  # returns true if exists one or more records matching condition
  def exists?(condition)
    num = count(condition)
    if (num > 0)
      return true
    else
      return false
    end
  end

  # Insert one or more records into table
  def insert(table, records)
    return if records.empty? || records.first.empty?
    connect if (!@connected)
    val = ""
    records.each {|record|
      val << "("
      record.each {|value|
        val << "'" << value.to_s.tr("'","") << "'"
        val << ","
      }
      val.chop!
      val << "),"
    }
    val.chop!
    query("INSERT INTO #{table} VALUES #{val}")
  end
  
  def selectFirst(query)
    data = []
    query(query + " limit 1").each {|row|
      data = row
    }
    return data
  end
  
  def blastdb
    if (!@blastdb)
        query("select * from #{@combodb}.apisdbs").each {|row|
        @blastdb = row.first
      }
      if (@blastdb.nil?)
        STDERR.printf("Error: looks like there isn't a current blastdb for %s\n", db)
        exit(1)
        end
    end
    return @blastdb
  end
  
  def tax
    if (!@tax)
      @tax = Hash.new
      query("select name, species, strain, taxonomy, supergroup, form from #{@combodb}.contigs").each {|row|
        name, species, strain, taxonomy, supergroup, form = row
        @tax[name] = Hash.new
        @tax[name]["species"] = species
        @tax[name]["strain"] = strain
        if (form == "Mitochondria")
          taxonomy = "Bacteria; Proteobacteria; Alphaproteobacteria; Rickettsiales; Rickettsiaceae; Rickettsieae; Mitochondrion;"
        elsif (form == "Plastid")  
          taxonomy = "Bacteria; Cyanobacteria; Prochlorophytes; Prochlorococcaceae; Chloroplast;  Chloroplast;  Chloroplast;"
        end
        @tax[name]["taxonomy"] = taxonomy
        @tax[name]["format"] = form
        if (!supergroup.nil? && supergroup.length > 3)
          @tax[name]["supergroup"] = supergroup 
        end
      }
    end
    return @tax
  end
  def fetchProtID(id)
    fetched = false
    while(!fetched)
      begin
        query("select proteins.name, annotation, species, proteins.seq from #{@combodb}.contigs, #{@combodb}.proteins where contig_name = contigs.name and proteins.name = '#{id}'").each {|row|
          name, annotation, species, seq = row
          fetched = true
          seq.gsub!("*","")
          header = "#{name} #{annotation} {#{species}}" 
          return ">#{header}\n#{seq.gsub(Regexp.new(".{1,60}"), "\\0\n")}", seq.length  
        }
        return nil
      rescue
        connect
      end
    end
  end
  def fetchTranID(id)
    fetched = false
    while(!fetched)
      begin
        query("select transcripts.name, proteins.annotation, species, transcripts.seq from #{@combodb}.contigs, #{@combodb}.transcripts, #{@combodb}.proteins where transcripts.contig_name = contigs.name and proteins.name = transcripts.name and transcripts.name = '#{id}'").each {|row|
          name, annotation, species, seq = row
          fetched = true
          seq.gsub!("*","")
          header = "#{name} #{annotation} {#{species}}" 
          return ">#{header}\n#{seq.gsub(Regexp.new(".{1,60}"), "\\0\n")}", seq.length  
        }
        return nil
      rescue
        connect
      end
    end
  end
  def fetchContigSeq(id)
    fetched = false
    while(!fetched)
      begin
        query("select name, species, form, seq from #{@combodb}.contigs where name = '#{id}'").each {|row|
          name, species, form, seq = row
          fetched = true
          seq.gsub!("*","")
          header = "#{name} [#{form}] {#{species}}" 
          return ">#{header}\n#{seq.gsub(Regexp.new(".{1,60}"), "\\0\n")}"
        }
        return nil
      rescue
        connect
      end
    end
  end
  def fetchrRNA(id)
    fetched = false
    while(!fetched)
      begin
        q = "SELECT rrnas.name, annotation, species, rrnas.seq FROM "
        q += "#{@combodb}.contigs, #{@combodb}.rrnas WHERE " 
        q += "contig_name = contigs.name AND rrnas.name = '#{id}'"
        p q
        query(q).each {|row|
          name, species, form, seq = row
          fetched = true
          seq.gsub!("*","")
          header = "#{name} #{annotation} {#{species}}" 
          return ">#{header}\n#{seq.gsub(Regexp.new(".{1,60}"), "\\0\n")}"
        }
        return nil
      rescue
        connect
      end
    end
  end
  def fetchFunction(id)
    get("select annotation from #{@combodb}.proteins where proteins.name = '#{id}'")[0]
  end
  # returns best BLAST hit for a given protein in a dataset
  def fetchTopBlastHit(name, dataset)
    query = "SELECT subject_name, score FROM blast "
    query += "WHERE dataset = \"#{dataset}\" AND seq_name = \"#{name}\""
    query += " ORDER BY score DESC LIMIT 1"
    sname, score = get(query)
    return sname
  end
  # returns best BLAST hit for a given protein in a dataset
  def fetchTree(name, dataset)
    query = "SELECT tree FROM tree "
    query += "WHERE dataset = \"#{dataset}\" AND seq_name = \"#{name}\""
    get(query)
  end
  # returns APIS annotation for a given protein in a dataset
  def fetchAPISAnnotation(name, dataset)
    query = "SELECT annotation FROM annotation "
    query += "WHERE dataset = \"#{dataset}\" AND seq_name = \"#{name}\""
    query += " AND source=\"APIS\""
    get(query)
  end
  # returns taxon id, species, taxonomy
  def fetchSpeciesInfo(name)
    get("SELECT t2.taxon_id, t2.species, t2.taxonomy FROM proteins t1 INNER JOIN " +
    "contigs t2 ON t1.contig_name = t2.name WHERE t1.name = '#{name}'")
  end
  # fetches swiss prot name and definition
  def fetchSwiss(name)
    get("SELECT t1.hit, t2.definition FROM blast t1 INNER JOIN swiss t2 ON t1.hit = t2.swiss " +
    "WHERE t1.source = 'phylodb' AND t1.target = 'swiss' AND t1.blast = 'blastp' AND " +
    "t1.query = '#{name}'")
  end
  # fetches Kegg ortholog name and definition
  def fetchKegg(name)
    pathway = nil
    ko, rest = get("SELECT t2.ko FROM blast t1 INNER JOIN kegg t2 ON t1.hit = t2.kegg_seq_id WHERE t1.source = 'phylodb' 
    AND t1.target = 'kegg' AND t1.blast = 'blastp' AND t2.ko IS NOT NULL AND t1.query = '#{name}'")
    if (ko)
      pathway, rest = get("SELECT pathway FROM pathways WHERE ko = '#{ko}'")
    end
    return [ko, pathway]
  end
  # fetches Cluster info
  def fetchClusters(name)
    get("SELECT gos_cluster_id, phytax_cluster_id, plant, animal, pfam, tigrfams FROM map WHERE " + 
    "phylodb_seq_id = '#{name}'")
  end
  # fetches go from pfam
  def fetchPfam2Go(pfam)
    return "" if pfam.nil?
    pfam = $1 if (pfam =~/(.*)\.\d+/) # no decimals
    go = ""
    query("SELECT go FROM pfam2go WHERE pfam = '#{pfam}'").each do |row|
      go += ";" if (go != "")
      go += row[0]
    end
    return go
  end
  # fetches ann for pfam
  def fetchPfamAnn(pfam)
    get("SELECT definition FROM pfam WHERE pfam='#{pfam}'")
  end
  # fetches ann for tigrfam
  def fetchTigrfamAnn(tigrfam)
    get("SELECT definition FROM tigrfams WHERE tigrfams='#{tigrfam}'")
  end
  # fetches go from tigrfam
  def fetchTigrFam2Go(tigrfam)
    return "" if tigrfam.nil?
    go = ""
    query("SELECT go FROM tigrfams2go WHERE tigrfams = '#{tigrfam}'").each do |row|
      go += ";" if (go != "")
      go += row[0]
    end
    return go
  end
  def fetchAnnotation(name)
    get("SELECT annotation FROM proteins WHERE name = \"#{name}\"")
  end
  # returns taxonomy for phylodb (combodb, etc.) protein
  def fetchTaxonomy(name)
    id, contig = name.split("-")
    get("SELECT taxonomy FROM contigs WHERE name = \"#{contig}\"")
  end

  # return the annotation of the closest sequence to id on tree 
  def findClosestFunction(tree, id)
    begin
      tree.relatives(id).each {|match|
        acc, contig = match.first.split("-")
        contig, rest = contig.split("__")
        match_id = acc + "-" + contig
        function = fetchFunction(match_id)
        if (!function.nil? && function.split(" ").size > 1 && 
            function !~/unnamed/ && function !~/unknown/ && 
            function !~/numExons/ && function !~/^\{/)
          return function
        end
      }
      return false
    rescue
      return false
    end
  end
  
  # delete contig(s) matching condition and all associated records
  def deleteContig(condition)
    query("select name, species from #{@combodb}.contigs where #{condition}").each {|row|
      name, species  = row
      STDERR.printf("Deleting %s (%s)...\n", name, species)
      query("delete from proteins where contig_name = '#{name}'")
      query("delete from transcripts where contig_name = '#{name}'")
      query("delete from geneorders where contig_name = '#{name}'")
      query("delete from rrnas where contig_name = '#{name}'")
      query("delete from contigs where name = '#{name}'")
    }
  end
  # build up taxonomy line from taxonomy table given taxid 
  def buildTaxFromTaxId(taxid, string = false, verbose = false)
    levels = ["kingdom", "phylum", "class", "order", "family", 
              "genus", "species"]
    name = ""
    tax = [""]*7
    while (name != "root")
      query = "select parent_id, name, rank from taxonomy WHERE tax_id = #{taxid}"
      pid, name, rank = selectFirst(query)
      STDERR.printf("%d\t%d\t%s\t%s\n", taxid, pid, name, rank) #if verbose
      return nil if pid.nil?
      pos = levels.index(rank)
      if (pos.nil?)
        pos = 0 if name == "Viruses" || name == "Viroids"
        pos = 1 if name =~ /viruses/
      end
      tax[pos] = name.tr(",()[]'\"/","") if (pos)
      taxid = pid
    end
    6.step(0, -1) {|i|
      if (tax[i] == "")
        tax[i] = tax[i + 1].split(" (").first + " (" + levels[i] + ")"
      end
    }
    if (string)
      tline = ""
      tax.each {|lev|
        tline += lev
        tline += "; " if lev != tax.last
      }
      return tline
    else
      return tax
    end
  end
  
  # create dataset if it doesn't exist
  def createDataset(dataset, owner, date, database, comments = "", group = "",
                    username = "", password = "")
    if (!exists?("dataset where dataset='#{dataset}'"))
      insert("dataset", [[dataset, owner, date, database, comments, group,
                          username, password]])
    end
  end
  
  # delete dataset(s) matching condition and all associated records
  def deleteDataset(condition, keepBlast = false)
    query("select dataset from dataset where #{condition}").each {|row|
      dataset, rest  = row
      ["alignment", "annotation", "tree", "classification", "blast",
       "sequence", "dataset"].each {|tbl|
        next if (tbl == "blast" && keepBlast)
        STDERR.printf("Deleting data in %s for %s...\n", tbl, dataset)
        query("delete from #{tbl} where dataset = '#{dataset}'")
      }
    }
  end
  
  # load peptides from file
  def loadPeptides(prot, dataset, include = false)
    if (!exists?("sequence where dataset='#{dataset}'") || include)
      STDERR.printf("Loading Peptides from %s...\n", prot)
    else
      return
    end
    seqs = []
    inputs = []
    count = 0
    Bio::FlatFile.new(Bio::FastaFormat, ZFile.new(prot)).each {|seq|
      id = seq.entry_id.gsub("(","").gsub(")","")
      seqs.push([id, dataset, seq.seq, 0])
      count += 1
      name, rest = seq.definition.split(" ", 2)
      if (!rest.nil? && rest.length > 3)
        inputs.push([id, dataset, rest.strip, 'input'])
      end
      if (count % 10000 == 0)
        insert("sequence", seqs)
        insert("annotation", inputs)
        seqs = []
        inputs = []
        STDERR.printf("Loaded sequence %d...\n", count)
      end
    }
    if (seqs.size > 0)
      insert("sequence", seqs)
      insert("annotation", inputs)
    end
  end
  
  # return processed state of peptide
  def processed?(seq_name, dataset)
    query("select processed from sequence where seq_name = '#{seq_name}' AND dataset='#{dataset}'").each {|row|
      return true if row[0] == "1"
    }
    return false
  end
  
  # set processed state for peptide
  def setProcessed(seq_name, dataset)
    query("update sequence set processed=1 where seq_name = '#{seq_name}' AND dataset='#{dataset}'") 
  end
  
  # unset processed state for peptide
  def setUnProcessed(seq_name, dataset)
    query("update sequence set processed=0 where seq_name = '#{seq_name}' AND dataset='#{dataset}'") 
  end
  
  # return blast info for given seq_name + dataset and evalue
  def fetchBlast(seq_name, dataset, evalue, maxTree, tax)
    homologs = []
    query("select subject_name, subject_length, evalue from blast where seq_name='#{seq_name}' and dataset='#{dataset}' and evalue <= #{evalue} order by evalue").each {|row|
      homologs.push(row[0]) if (homologs.size < maxTree && !homologs.include?(row[0]) && row[1].to_i < 2000)
    }
    return homologs
  end
  
  # delete blast records for seq with evalue > threshold
  def deleteBlast(seq_name, dataset, evalue)
    query("delete from blast where seq_name='#{seq_name}' and dataset='#{dataset}' and evalue > #{evalue}")
  end
  # stores alignment, deleting previous if there
  def createAlignment(seq_name, dataset, afa)
    query("delete from alignment where seq_name = '#{seq_name}' and dataset = '#{dataset}'")
    lines = []
    Bio::FlatFile.new(Bio::FastaFormat, File.new(afa)).each {|aseq|
      lines.push([seq_name, dataset, aseq.entry_id, 
                  aseq.definition.split(" ", 2).last, aseq.seq])
    }
    insert("alignment", lines)
  end
  
  
  # returns dataset for peptide, with optional list of possible datasets
  def getDataset(seq_name, datasets = nil)
    if (datasets.nil?)
      query("select dataset from dataset").each {|dataset|
        datasets.push(dataset[0])
      }
    end
    datasets.each {|dataset|
      if (exists?("sequence where dataset='#{dataset}' and seq_name='#{seq_name}'"))
        return dataset
      end
    }
    raise DatasetNotFoundError, "Cannot find dataset for #{seq_name}"
  end
  
  # stores tree
  def createTree(seq_name, dataset, tree)
    if (exists?("tree where seq_name='#{seq_name}' AND dataset = '#{dataset}'"))
        query("UPDATE tree SET tree = '#{tree}' WHERE seq_name='#{seq_name}' AND dataset = '#{dataset}'")
      else
        insert("tree", [[seq_name, dataset, tree]])
      end
  end
  
  # interprets tree, creating classification
  def createClassification(tree, name, dataset, exclude, ruleMaj)
    cons = consensusTax(tree, name, ruleMaj)
    lines = []
    cons.each {|line|
      lines.push(line) if (line.grep(/#{exclude}/).empty? || exclude.nil?)
    }
    first = lines[0]
    first=[nil,nil,nil,nil,nil,nil,nil] if first.nil?
    if (lines[1].nil?)
      second = nil
    else
      second = lines[1]
    end
    mixed = false
    classification = [name, dataset]
    7.times {|level|
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
    }
    return classification
  end
  
  # creates phylogenomic annotation
  def createAnnotation(tree, seq_name, dataset)
    function  = findClosestFunction(tree, seq_name)
    if (function)
      return [seq_name, dataset, function.strip, "APIS"]
    else
      return nil
    end
  end

  # returns array of consensus taxonomy at each relative level of tree
  def consensusTax(tree, taxon, ruleMaj)
    consensus = []
    return  [] if (tree.relatives(taxon).nil?)
    tree.relatives(taxon).each {|list|
      counts = []
      list.each {|relative|
        acc, contig = relative.split("-")
	contig, rest = contig.split("__")
	next if (tax[contig].nil? || tax[contig]["species"].nil?)
	groups = tax[contig]["taxonomy"].split(/; |;/)
        groups.size.times {|i|
          counts[i] = Hash.new if counts[i].nil?
          counts[i][groups[i]] = 0 if counts[i][groups[i]].nil?
          counts[i][groups[i]] += 1
        }
      }
      if (ruleMaj)
        consensus.push(counts.majority)
      else
        consensus.push(counts.absolute)
      end
    }
    return consensus
  end
end

class MyDB < DBwrapper
  def initialize(host, dbname, user, password, combodb = nil)
    require "mysql"
    @host = host
    @db = dbname
    @user = user
    @password = password
    @connected = false
    if (combodb.nil?)
      @combodb = dbname
    else
      @combodb = combodb
      connect
    end
  end
  def connect
    @connected = false
    while(!@connected)
      begin 
        @connection = Mysql::new(@host, @user, @password, @db) 
        @connected = true
      rescue
        #STDERR.printf("%s: waiting...\n", $!)
        sleep 0.2
      end
    end
  end
  def get(query)
    connect if !@connected
    return @connection.query(query).fetch_row
  end
  # creates and inserts classification row
  def createClassification(tree, name, dataset, exclude, ruleMaj)
    classification = super(tree, name, dataset, exclude, ruleMaj)
    where = "WHERE seq_name='#{name}' AND dataset = '#{dataset}'"
    old = []
    query("SELECT * FROM classification #{where}").each {|row|
      old = row
    }
    if (!old.empty?)
      cl = classification.dup
      name, dataset = cl.shift, cl.shift
      name, dataset = old.shift, old.shift
      ranks = ["kingdom", "phylum", "class", "ord", "family", "genus", 
               "species"]
      ranks.each {|rank|
        val = cl.shift
        out = cl.shift
        oval = old.shift
        oout = old.shift.to_i
        if (val != oval)
          query("UPDATE classification SET #{rank} = '#{val}' #{where}")
        end
        if (out != oout)
          query("UPDATE classification SET #{rank}_outgroup = #{out} #{where}")
        end
      }
    else
      insert("classification", [classification])
    end
  end  
  # creates and inserts annotation string
  def createAnnotation(tree, seq_name, dataset)
    function = super(tree, seq_name, dataset)
    if (function)
      where = "WHERE seq_name='#{seq_name}' AND dataset = '#{dataset}' AND source = 'APIS'"
      if (exists?("annotation #{where}"))
        query("UPDATE annotation SET annotation = '#{function[-2]}' #{where}")
      else
        insert("annotation", [function]) 
      end
    end
  end
end

class SQLite < DBwrapper
  class NoSQLiteDB < RuntimeError
  end
  def initialize(file, mode = "r")
    require 'sqlite3'
    if (File.exist?(file))
      @db = file
    elsif (mode == "r")
      raise(NoSQLiteDB, "No such database as #{file}!")
    end
    connect
  end
  def connect
    @connection = SQLite3::Database.new(@db)
    @connected = true
  end
  # Insert one or more records into table, overridden because SQLite doesn't
  # allow multiple inserts
  def insert(table, records)
    records.each {|record|
      super(table, [record])
    }
  end
  def count(condition)
    return @connection.get_first_value("select count(*) from " + condition).to_i
  end
  # return first row of result immediately
  def get(query)
     return @connection.get_first_row(query)
  end
  # returns array of consensus taxonomy at each relative level of tree
  def consensusTax(tree, taxon, ruleMaj)
    consensus = []
    return  [] if (tree.relatives(taxon).nil?)
    tree.relatives(taxon).each {|list|
      counts = []
      list.each {|relative|
        query =  "SELECT taxonomy FROM apis_proteins "
        query += "WHERE protein_id=#{relative}"
        groups = @connection.get_first_value(query).split(/; |;/) 
        groups.size.times {|i|
          counts[i] = Hash.new if counts[i].nil?
          counts[i][groups[i]] = 0 if counts[i][groups[i]].nil?
          counts[i][groups[i]] += 1
        }
      }
      if (ruleMaj)
        consensus.push(counts.majority)
      else
        consensus.push(counts.absolute)
      end
    }
    return consensus
  end
  # returns function of protein based on either protein id or name
  def fetchFunction(id)
    if (id =~/[^0-9]/)
      where = "protein_name = '#{id}'"
    else
      where = "protein_id = #{id}"
    end
    query = "SELECT annotation FROM apis_proteins WHERE #{where}"
    return @connection.get_first_value(query)
  end
  def findClosestFunction(tree, id)
    begin
      tree.relatives(id).each {|match|
        function = fetchFunction(match.first)
        if (!function.nil? && function.split(" ").size > 1 && 
            function !~/unnamed/ && function !~/unknown/ && 
            function !~/numExons/ && function !~/^\{/)
          return function
        end
      }
      return false
    rescue
      return false
    end
  end
  # creates a version of the tree with protein names & species
  def SpeciesTree(tree)
    ali = Hash.new
    tree.taxa.each {|taxon|
      if (taxon !~/[^0-9]/)
        query = "SELECT protein_name, taxonomy FROM apis_proteins WHERE "
        query += "protein_id = #{taxon}"
        name, taxonomy = @connection.get_first_row(query)
        if (name)
          species = taxonomy.split(/; |;/).last.gsub(" ","_")
          ali[taxon] = name + "__" + species
        end
      end
    }
    return tree.unAlias(ali)
  end
end

class Array
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
  def mostCommon
    count = Hash.new
    self.each do |el|
      if (count[el].nil?)
	      count[el] = 1
      else
	      count[el] += 1
      end
    end
    sorted = count.keys.sort {|a, b| count[b] <=> count[a]}
    return sorted[0]
  end
end
