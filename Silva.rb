class Silva
  def initialize(tax, alignment)
    @tax = tax
    @alignment = alignment
  end
  def findTaxonomy(string)
    line = `grep #{string} #{@tax}| head -1`.chomp
    id, tax = line.split(" ", 2)
    return tax
  end
end
