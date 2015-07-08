class Array
  def hex
    color = "#"
    self.each do |value|
      value = value.to_s(base = 16).upcase
      value = "0" + value if (value.length == 1)
      color += value
    end
    return color
  end
  def dist(col)
    return Math.sqrt(((col[0] - self[0])**2 + (col[1] - self[1])**2 + (col[2] - self[2])**2).to_f)
  end
  def averageDist(cols)
    dist = 0
    cols.each do |col|
      dist += self.dist(col)
    end
    return dist/cols.size
  end
end

def pickDistantColor(colors)
  maxD = 0
  bestC = nil
  100.times do |i|
    candidate = generateRandomColor([128,128,128]) 
    d = candidate.averageDist(colors)
    if (d > maxD)
      bestC = candidate
      maxD = d
    end
  end
  return bestC
end

def generateRandomColor(mix) 
  red = rand(256)
  green = rand(256)
  blue = rand(256)
  red = (red + mix[0]) / 2
  green = (green + mix[1]) / 2
  blue = (blue + mix[2]) / 2
  return [red, green, blue]
end
