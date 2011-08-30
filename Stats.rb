class Array
  def average
    tot = 0
    self.size.times {|i|
      tot += self[i]
    }
    return tot/(self.size*1.0)
  end
  def median
    sorted = self.sort
    mid = sorted.size / 2
    if (self.size % 2 == 1)
      return sorted[mid]
    else
      return (sorted[mid] + sorted[mid - 1]) / 2.0
    end
  end
  def sumit
    tot = 0
    self.size.times {|i|
      tot += self[i]
    }
    return tot
  end
  def product
    tot = 1
    self.size.times {|i|
      tot *= self[i]
    }
    return tot
  end
  def stdErrMean
    mean = self.average
    tot = 0
    self.size.times {|i|
      tot += (self[i] - mean)**2
    }
    variance = tot / (1.0*self.size)
    return Math.sqrt(variance)
  end
  def stdDev
    mean = self.average
    tot = 0
   self.size.times {|i|
      tot += (self[i] - mean)**2
    }
    variance = tot / (self.size - 1.0)
    return Math.sqrt(variance)
  end
  def zScores
    mean = self.average
    dist = self.stdErrMean
    zScores = Array.new
   self.size.times {|i|
      zScores[i] = ((self[i] - mean) / dist)
    }   
    return zScores
  end
  def correlationWith(scores)
    xZ = self.zScores
    yZ = scores.zScores
    r = 0
   self.size.times {|i|
      r += (xZ[i]*yZ[i])
    }
    r /= self.size
    return r
  end
  def euclidDist(array)
    d = 0
    self.size.times {|i|
      d += (self[i] - array[i])**2
    }
    return Math.sqrt(d)
  end
end

class Fixnum
  def factorial
    product = 1
    1.upto(self) {|i|
      product = product * i ;
    }
    return product
  end
  def choose(k)
    return self.factorial / (k.factorial * (self - k).factorial)
  end
  def confidenceInterval(p)
    stdDev = Math.sqrt(p*(1-p))/Math.sqrt(self)   
  end
end

def normalP(x) 
  # Abramowitz & Stegun 26.2.19
  d1 = 0.0498673470
  d2 = 0.0211410061
  d3 = 0.0032776263
  d4 = 0.0000380036
  d5 = 0.0000488906
  d6 = 0.0000053830
  a = x.abs
  t = 1.0 + a*(d1+a*(d2+a*(d3+a*(d4+a*(d5+a*d6)))))
  # to 16th power
  t *= t;  t *= t;  t *= t;  t *= t
  t = 1.0 / (t+t)  # the MINUS 16th
  return t
end

def computeZScore(observed, expected, n)
  p = 1.0*observed / n
  pe = 1.0*expected / n
  z = (p - pe) / (Math.sqrt(pe*(1-pe))/Math.sqrt(n))
  return z
end

def approxBinom(x, y, p)
  expected = (x + y) * p
  return 2*normalP(computeZScore(x, expected, x + y))
end



