class Hit
  attr_reader :name
  attr_reader :pvalue
  attr_reader :qmatch
  attr_reader :smatch
  attr_reader :description
  attr_reader :sstart
  attr_reader :send
  attr_reader :qstart
  attr_reader :qend
  attr_reader :identity
  def initialize(fp, name, desc = nil)
    @fp = fp
    @name = name
    @description = desc
    @qmatch = ""
    @smatch = ""
    @sstart = 0
    @send = 0
    @qstart = 0
    @qend = 0
    @pvalue = nil
    @ident = 0
    @fp.each {|line|
      break if (line =~/^BLAST|^>/)
      if (line =~ /Expect = ([^ ]*)/)
        if (@pvalue.nil?)
          @pvalue = $1.to_f
        else
          break
        end
      end
      if (line =~ /^Query:/)	
	fields = line.split(" ")
        @qstart = fields[1].to_i if (@qstart == 0)
	@qmatch += fields[2]
        @qend = fields[3]
      elsif (line =~ /^Sbjct:/)
	fields = line.split(" ")
        @sstart = fields[1].to_i if (@sstart == 0)
	@smatch += fields[2]
        @send = fields[3]
      elsif (line =~ /Identities = ([0-9]*)\/([0-9]*)/)
        @identity = 100*$1.to_i/$2.to_i
      end
    }
  end
  def strand
    if (@send.to_i > @sstart.to_i)
      return "+"
    else
      return "-"
    end
  end
  def coverage(length)
    return (100*@smatch.length)/length
  end
end

class Query
  attr_reader :name
  attr_reader :description
  attr_reader :length
  def initialize(fp, name, desc = nil)
    @fp = fp
    @name = name
    @description = desc
    @firstHit = nil
    @fp.each {|line|
      if line =~/([0-9]*) letters/
        @length = $1.to_i
        break
      end
    }
  end
  def each
    seen = false
    @fp.each {|line|
      break if (line =~/^BLAST/)
      if (line =~/>([^ ]*)(.*)/)
	yield Hit.new(@fp, $1.chomp, $2.chomp)
	seen = true 
      end
    }
    return nil if (!seen)
  end
  def firstHit
    if (@firstHit.nil?)
      self.each {|hit|
	@firstHit = hit
	break
      }
    end
    return @firstHit
  end
end

class Blast
  def initialize(file)
    @fp = File.new(file)
  end
  def each
    seen = false
    @fp.each {|line|
      if (line =~ /^Query=[\ ]*([^ ]*)(.*)/)
	yield Query.new(@fp, $1.chomp, $2.chomp)
	seen = true
      end
    }
    return nil if (!seen)
  end
end

class String
  def to_fasta(header, len)
    self.tr!("\n","")
    seq = ">#{header}\n" 
    seq += self.gsub(Regexp.new(".{1,#{len}}"), "\\0\n")
    return seq
  end
end
