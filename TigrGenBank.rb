require 'bio'
include Bio

# class for manipulating category annotated .gb files created by TigrGenome2Genbank
class TigrGenBank
  attr_reader :product
  attr_reader :role
  def initialize(gbfile = nil)
    @product = Hash.new
    @role = Hash.new
    add(gbfile) if (!gbfile.nil?)
  end

  # populate product and role hashes
  def add(gbfile)
    compress = false
    if (gbfile.index(".gz"))
      system("PATH=/usr/bin;zcat #{gbfile} > #{gbfile}.tmp")
      gbfile = gbfile + ".tmp"
      compress = true
    end
    FlatFile.new(GenBank, File.new(gbfile)).each {|seq|
      seq.each_cds {|cds|
        @product[cds.assoc["protein_id"]] = cds.assoc["product"]
        @role[cds.assoc["protein_id"]] = cds.assoc["note"].split(";")
      }
    }
    File.unlink(gbfile) if compress
  end

  # summarize main categories, optional array of selected keys
  def summarizeMain(subKeys = nil)
    counts = Hash.new
    total = 0
    @role.keys.each {|name|
      if (subKeys.nil? || subKeys.include?(name))
        category = @role[name].first
        if (counts[category].nil?)
          counts[category] = [name]
        else
          counts[category].push(name)
        end
        total += 1.0
      end
    }
    return counts, total
  end
end
