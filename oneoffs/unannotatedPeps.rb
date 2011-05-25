require 'bio'

Bio::FlatFile.new(Bio::FastaFormat, File.new(ARGV.first)).each do |seq|
  if (seq.definition.index(seq.entry_id, 1))
    print seq
  end
end
