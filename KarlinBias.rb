require 'bio'
include Bio

class CodonVector
  attr :codonCounts
  attr :aaCounts
  attr :definition
  def initialize(seq = nil)
    @aaCounts = Hash.new
    @codonCounts = Hash.new
    @definition = ""
    totAA = 0
    CodonTable[11].keys.each {|codon|
      aa = CodonTable[11][codon]
      next if (aa == "*")
      @aaCounts[aa] = 0 if (@aaCounts[aa].nil?)
      @codonCounts[aa] = Hash.new if (@codonCounts[aa].nil?)
      @codonCounts[aa][codon] = 0
    }
    if (!seq.nil?)
      @definition = seq.definition
      0.step(seq.length - 3, 3) {|i|
	codon = seq.seq[i, 3].to_s.downcase
	aa = CodonTable[11][codon]
	next if (aa == "*" || aa == "X" || aa.nil?)
	@aaCounts[aa] += 1
	totAA += 1
	@codonCounts[aa][codon] += 1
      }
      @aaCounts.keys.each {|aa|
	next if (@aaCounts[aa] == 0)
	@codonCounts[aa].keys.each {|codon|
	  @codonCounts[aa][codon] /= (1.0*@aaCounts[aa])
	}
	@aaCounts[aa] /= (1.0*totAA)
      }
    end
  end
  def +(vector)
    sum = CodonVector.new
    sum.codonCounts.keys.each {|aa|
      sum.aaCounts[aa] = vector.aaCounts[aa] + self.aaCounts[aa]
      sum.codonCounts[aa].keys.each {|codon|
	sum.codonCounts[aa][codon] += (vector.codonCounts[aa][codon] +
				       self.codonCounts[aa][codon])
      }
    }
    return sum
  end
  def -(vector)
    difference = CodonVector.new
    difference.codonCounts.keys.each {|aa|
      difference.aaCounts[aa] = (vector.aaCounts[aa] - self.aaCounts[aa]).abs
      difference.codonCounts[aa].keys.each {|codon|
	difference.codonCounts[aa][codon] += (vector.codonCounts[aa][codon] -
				       self.codonCounts[aa][codon]).abs
      }
    }
    return difference
  end
  def *(factor)
    product = CodonVector.new
    product.codonCounts.keys.each {|aa|
      product.aaCounts[aa] = self.aaCounts[aa] * factor
      product.codonCounts[aa].keys.each {|codon|
	product.codonCounts[aa][codon] = self.codonCounts[aa][codon] * factor
      }
    }
    return product
  end
  def CodonVector.average(vectors)
    sum = CodonVector.new
    vectors.each {|vector|
      sum = sum + vector
    }
    return sum * (1.0 / vectors.size)
  end
  def definition=(value)
    @definition = value
  end
  def CodonVector.medianBias(vectors)
    averageVector = CodonVector.average(vectors)
    bias = []
    vectors.each {|vector|
      bias.push(vector.karlinBias(averageVector))
    }
    return bias.sort[bias.size / 2]
  end
  def karlinBias(group)
    difference = self - group
    bias = 0
    @aaCounts.keys.each {|aa|
      subBias = 0
      difference.codonCounts[aa].keys.each {|codon|
	subBias += difference.codonCounts[aa][codon]
      }
      bias += (@aaCounts[aa]*subBias)
    }
    return bias
  end
  def CodonVector.findRibosomal(vectors)
    set = []
    vectors.each {|vector|
      if (vector.definition =~ /rpl|rps|rpm|ribosomal protein/i)
	set.push(vector)
      end
    }
    return CodonVector.average(set)
  end
  def CodonVector.findChaperone(vectors)
    set = []
    vectors.each {|vector|
      if (vector.definition =~ /groEL|dnaK|Ion1|Ion2|HSP|groES|clpX|ppiC/i ||
	  vector.definition =~ /grpE|msrA|ftsH-3|trxA|chaperon/i)
	set.push(vector)
      end
    }
    return CodonVector.average(set)
  end
  def CodonVector.findTrans(vectors)
    set = []
    vectors.each {|vector|
      if (vector.definition =~ /EF-Ts|EF-G|EF-TU|rpoB|rpoC|rpo-C|typA|bipA/i ||
	  vector.definition =~ /infB|EF-P|greA|nusG|nusA|rpoA|infC|frr/i)
	set.push(vector)
      end
    }
    return CodonVector.average(set)
  end
  def expression(average, ribosomal, chaperone, trans)
    return karlinBias(average) / (0.5*karlinBias(ribosomal) +
				  0.25*karlinBias(chaperone) +
				  0.25*karlinBias(trans))
  end
  def alien(medBias, ribosomal, chaperone, trans, level = 0.25)
    thresh = medBias + level
    rb = karlinBias(ribosomal)
    cb = karlinBias(chaperone)
    tb = karlinBias(trans)
    if (rb >= thresh && cb >= thresh && tb >= thresh)
      return (rb + cb + tb) / 3
    else
      return false
    end
  end
end
