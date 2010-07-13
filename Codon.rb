# class to calculate codon bias, dinucleotide frequencies, etc.
class Codon

def initialize
  @seqs = []
end  

# adds a sequence to the set
def add(seq)
  begin
    @seqs.push(seq.seq)
  rescue
    @seqs.push(seq) # if seq only String not Sequence
  end
end

# calculates 64-vector for sequence
def calc(seq, vector = nil)
  vector = Hash.new if (vector.nil?)
  0.step(seq.length - 3, 3) {|i|
    codon = seq[i, 3]
    vector[codon] = 0 if (vector[codon].nil?)
    vector[codon] += 1
  }
  return vector
end

# calculates 64-vector for all seqs in set
def calcAll
  vector = Hash.new
  @seqs.each {|seq|
    calc(seq, vector)
  }
  return vector
end

end



