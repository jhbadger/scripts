require 'rubygems'
require 'bio'
include Bio

# class to encapsulate NCBI taxonomy database
class Taxonomy
  def initialize(dir)
    STDERR.printf("Loading NCBI Taxonomy...")
    @parent = Hash.new
    @names = Hash.new
    @rank = Hash.new

    File.new(dir + "/" + "names.dmp").each {|line|
      num, name, blah, type = line.split("\t|\t")
      next if (type.index("scientific name").nil?)
      @names[num.to_i] = name
    }
    File.new(dir + "/" + "nodes.dmp").each {|line|
      num, parent, rk = line.split("\t|\t")
      name = @names[num.to_i]
      if (name.nil?)
        STDERR.printf("Error! #{num} not found!\n")
      else
          if (name == "Bacteria" || name == "Archaea" || name == "Eukaryota")
            rk = "kingdom"
          elsif rk == "kingdom"
            rk = "subkingdom"
          end
          @parent[num.to_i] = parent.to_i
          @rank[num.to_i] = rk
      end
    }
    STDERR.printf("done\n")
  end
  
  # write taxonomy.txt file to be loaded into database
  def writeAll(filename)
    out = File.new(filename, "w")
    @names.keys.sort.each {|num|
      out.printf("%d\t%d\t%s\t%s\n", num, @parent[num], @names[num],
      @rank[num]) 
    }
    out.close
  end
end

tax = Taxonomy.new("/Users/jbadger/Desktop")
tax.writeAll("/Users/jbadger/Desktop/taxonomy.txt")
