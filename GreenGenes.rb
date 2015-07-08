require 'ZFile'
require 'rubygems'
require 'bio'

class GreenGenes
  def initialize(file)
    @file = ZFile.new(file)
  end
  def each
    sep = $/
    $/ = "END"
    @file.each {|rec|
      $/ = sep
      record = Hash.new
      rec.split("\n").each {|line|
        if line.index("=")
          key, value = line.split("=")
          record[key] = value
        end
      }
      yield record
    }
  end
end


#GreenGenes.new("/Users/jbadger/Downloads/greengenes16SrRNAgenes.txt.gz").each {|record|
#  begin
#    seq = Bio::Sequence::NA.new(record["aligned_seq"].tr(".-",""))
#    header = record["prokMSA_id"] + " " + record["ncbi_tax_string"] + "; " +
#      record["source"]
#    if (header =~/Bacteroidetes|Firmicutes|Proteobacteria|Actinobacteria/)
#      print seq.to_fasta(header, 60)
#    end
#  rescue
#  end
#}
