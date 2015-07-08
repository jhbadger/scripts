class ColorGradient
  # initialize with three arrays of r,g,b values
  def initialize(c1, c2, c3)
    @c1 = c1
    @c2 = c2
    @c3 = c3
  end
  # ratio is fraction between 0 and 1
  def getColor(ratio)
    if (ratio < 0.5)
      ratio *= 2
      red = (@c1[0] * (1 - ratio) + @c2[0] * ratio).round
      green = (@c1[1] * (1 - ratio) + @c2[1] * ratio).round
      blue = (@c1[2] * (1 - ratio) + @c2[2] * ratio).round
    else
      ratio = (ratio - 0.5)*2
      red = (@c2[0] * (1 - ratio) + @c3[0] * ratio).round
      green = (@c2[1] * (1 - ratio) + @c3[1] * ratio).round
      blue = (@c2[2] * (1 - ratio) + @c3[2] * ratio).round
    end
    return [red, green, blue]
  end
end
