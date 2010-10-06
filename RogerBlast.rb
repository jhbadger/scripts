

require 'ZFile'

class RogerBlast
  attr_reader :query, :score, :evalue, :ident, :strand, :target, :length
  def initialize(file)
    @file = ZFile.new(file)
  end
  def each
    inQuery = false
    inTarget = false
    @query = ""
    @target = ""
    oldLine = ""
    @file.each do |line|
      if (line =~/^Query=/)
        @query = line.chomp.split("Query= ", 2).last
        inQuery = true
      elsif (inQuery && line == "\n")
        inQuery = false
        if (@query =~/([0-9]+) letters/)
          @length = $1.to_i
        end
      elsif (inQuery)
        @query += " " + line.strip.chomp
      elsif (line =~/Score =  ([0-9]+) bits \([0-9]+\), Expect = (.*)/)
        @score = $1.to_i
        @evalue = $2
        @evalue = "1" + @evalue if (@evalue =~/^e/)
        @evalue = @evalue.to_f
      elsif (line =~/Identities = ([0-9]+)\/([0-9]+)/)
        @ident = ($1.to_i*1000/$2.to_i)/10.0
      elsif (line =~/Strand = (.*)/)
        @strand = $1.split("/ ").last
        yield self
      elsif (line =~/^>/)
        @target = line.chomp.split(">", 2).last
        inTarget = true
      elsif (inTarget && line == "\n")
        inTarget = false
       elsif (inTarget)
         if (oldLine.strip.size < 57) 
           @target += " " 
         end
         @target += line.strip.chomp 
      end
      oldLine = line
    end
  end
end