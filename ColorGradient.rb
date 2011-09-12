class ColorGradient
  # initialize with two arrays of r,g,b values
  def initialize(c1, c2)
    @c1 = c1
    @c2 = c2
  end
  # ratio is fraction between 0 and 1
  def getColor(ratio)
    red = (@c1[0] * (1 - ratio) + @c2[0] * ratio).round
    green = (@c1[1] * (1 - ratio) + @c2[1] * ratio).round
    blue = (@c1[2] * (1 - ratio) + @c2[2] * ratio).round
    return [red, green, blue]
  end
end