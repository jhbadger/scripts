require 'rubygems'
require 'DBwrapper'
require 'ZFile'
require 'bio'

class FastaDB
  def initialize(fasta)
    @name = File.basename(fasta)[0..File.basename(fasta).rindex(".") - 1]
    if (!File.exists?(@name + ".db"))
      system("touch #{@name}.db")
      @db = SQLite.new(@name + ".db")
      @db.query("CREATE TABLE seq (name text NOT NULL , annotation text NOT NULL, seq text NOT NULL)")
      @db.query("CREATE INDEX name_idx ON seq (name ASC)")
      Bio::FlatFile.new(Bio::FastaFormat, ZFile.new(fasta)).each {|seq|
        @db.insert("seq", [[seq.entry_id, seq.definition[2 + seq.entry_id.size..seq.definition.length], seq.seq]])
      }
    else
      @db = SQLite.new(@name + ".db")
    end
  end
  def getFasta(name, contigOnly = false)
    seq = @db.query("SELECT name, annotation, seq FROM seq WHERE name = '#{name}'").each {|row|
      name, ann, seq = row
      name, contig = name.split("-")
      if (contigOnly)
        header = contig
      else
        header = name + " " + ann
      end
      return ">" + header + "\n" + seq.gsub(Regexp.new(".{1,60}"), "\\0\n")
    }
  end
end
