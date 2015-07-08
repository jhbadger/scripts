require 'Stats'

class Histogram

  def initialize(data, binNum = 10)
    delta = (data.max - data.min) / (binNum * 1.0)
    @bins = []
    binMax = data.min
    dataMax = data.max
    bin = nil
    @total = 0
    data.sort.each {|datum|
      if (datum >= binMax && datum != dataMax)
	bin = []
	@bins.push(bin)
	binMax += delta
      end
      bin.push(datum)
      @total += 1.0
    }
  end

  def to_s(title = "", marked = [], legend = "")
    s = ""
    s += "Histogram of #{title}\n\n" if (title != "")
    @bins.each {|bin|
      barLen = (bin.size / @total) * 40
      markLen = 0
      binMin = bin.min
      binMax = bin.max
      if (!marked.empty?)
	markNum = 0
	marked.each {|mark|
	  if (mark >= binMin && mark < binMax)
	    markNum += 1.0
	  end
	}
	markLen = barLen * (markNum / barLen)
	markLen = 1 if (markLen < 1 && markNum > 0)
	barLen -= markLen
      end
      s += sprintf("%-40s [%-8.3f, %8.3f] %6d (%3.1f%%)\n", "#" * markLen +
		   "*" * barLen, binMin, binMax, bin.size, 
		   bin.size * 100 / @total)
    }
    s += "\n# represents #{legend}\n" if (legend != "")
    return s
  end
end
