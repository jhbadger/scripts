class Chooser
  attr_reader :items, :k

  def initialize(items, k)
    @items = items
    @k = k
    @result = []    
  end
  
  def reset
    @k.times {|i|
      @result[i] = i
    }
    @result[@k] = @items.size
    (@k + 1).upto(@items.size) {|i|
      @result[i] = 0
    }
  end
  private :reset

  def size
    return @items.size
  end

  def nextSet
    good = false
    index = 0
    (@k - 1).downto(0) {|i|
      if (@result[i]+1 < @result[i+1]) 
	good = true
	index = i
	break
      end
    }
    if (good)
      @result[index] += 1
      (index + 1).upto(@k - 1) {|i|
	@result[i] = @result[i - 1] + 1
      }
      return true
    else
      return false
    end
  end
  private :nextSet

  def each
    reset
    loop do
      pick = []
      @k.times {|i|
	pick.push(@items[@result[i]])
      }
      yield pick.sort
      break if (!nextSet)
    end
  end

  def all
    reset
    perms = []
    self.each {|pick|
      perms.push(pick)
    }
    return perms
  end

end
