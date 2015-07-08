# Represents result for a single match in Btab format

require 'ZFile'
class BtabMatch
  # name of match
  attr_reader :name
  # range of query match
  attr_reader :qrange
  # range of subject match
  attr_reader :srange
  # percent identical
  attr_reader :percentId
  # percent similar
  attr_reader :percentSim
  # score
  attr_reader :score
  # subject description
  attr_reader :desc
  # reading frame
  attr_reader :frame
  # strand +/-
  attr_reader :strand
  # query match length
  attr_reader :qlen
  # subject match length
  attr_reader :slen
  # evalue
  attr_reader :evalue
  # species of match
  attr_reader :species
  # annotation of match
  attr_reader :annotation

  def initialize(name, qleft, qright, sleft, sright,
                 percentId, percentSim, score, desc, frame,
                 strand, evalue)
    @name = name
    @qrange = [qleft.to_i, qright.to_i]
    @srange = [sleft.to_i, sright.to_i]
    @percentId = percentId.to_f
    @percentSim = percentSim.to_f
    @score = score.to_i
    @desc = desc
    @frame = frame
    @strand = strand
    @slen = (srange.last - srange.first).abs + 1
    @qlen = (qrange.last - qrange.first).abs + 1
    $VERBOSE = nil
    @evalue = evalue.to_f
    $VERBOSE = true
    if (@desc =~/\{([^\}]*)/)
      @species = $1.strip
    else
      @species = nil
    end
    if (@species && @desc.index("{") > 0)
      @annotation = @desc[0..@desc.index("{") - 2]
    end
  end
end

# Represents result for a single query in Btab format
class BtabQuery 
  # name of query
  attr_reader :name
  # length of query
  attr_reader :length
  # array of matches
  attr_reader :matches
  # array of lines
  attr_reader :lines

  def initialize(name, length)
    @name = name
    @length = length.to_i
    @matches = []
    @lines = []
  end

  # add a match to the query
  def addMatch(subject, qleft, qright, sleft, sright, percentId, 
               percentSim, score, desc, frame, strand, evalue, line)
    match = BtabMatch.new(subject, qleft, qright, sleft, sright,
                          percentId, percentSim, score, desc, frame,
                          strand, evalue)
    @matches.push(match)
    @lines.push(line)
  end

  # list full length put. orthologs w/ evalue better or equal as given
  def listOrths(evalue = 1e-10, fraction = 0.5)
    names = []
    @matches.each {|match|
      next if (match.nil? || match.percentId > 98)
      if (match.evalue <= evalue)
        mlen = (match.srange.first - match.srange.last).abs
        next if ((mlen * 1.0) / @length < fraction)
	names.push(match.name)
      end
    }
    return names
  end

  # list  top full length ortholog
  def bestOrth(fraction = 0.5, percentThreshold = 98)
    @matches.each {|match|
      next if (match.percentId > percentThreshold)
      mlen = (match.srange.first - match.srange.last).abs
      next if ((mlen * 1.0) / @length < fraction)
      return match
    }
    return nil
  end

  # list each match
  def each
    @matches.each {|match|
      yield match
    }
  end
end

# Represents BLAST output in Btab format
class Btab 
  # date of BLAST run
  attr_reader :date
  # program used
  attr_reader :program
  # database used
  attr_reader :database
  # queries if all loaded
  attr_reader :queries

  def initialize(file)
    @date = nil
    @program = nil
    @database = nil 
    @query = nil
    @file = ZFile.new(file)
  end

  # load all query info into hash
  def loadAll
    @queries = Hash.new
    self.each {|query|
      @queries[query.name] = query if (!query.nil?)
    }
  end

  # step through each query
  def each
    oldQ = nil
    @file.each do |line|
      query, date, qlen, program, database, subject, qleft, qright, 
      sleft, sright, percentId, percentSim, dummy, score, dummy2,
      desc, frame, strand, slen, evalue = line.chomp.split("\t")
      desc = dummy2 if (desc.nil? || desc.length < 10) # fix timelogic btab
      @date = date if (!@date)
      @program = program if (!@program)
      @database = database if (!@database)
      if (query != oldQ)
        yield @query if !oldQ.nil?
        @query = BtabQuery.new(query, qlen)
      end
      @query.addMatch(subject, qleft, qright, sleft, sright,
                      percentId, percentSim, score, desc,
                      frame, strand, evalue, line)
      oldQ = query
    end
    yield @query if (!@query.nil?)
  end 
    
end

