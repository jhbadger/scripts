
# wrapper around scripts to load right version of csv depending on ruby version
require 'csv'
if CSV.const_defined? :Reader
	# Ruby 1.8 compatible
 	require 'fastercsv'
    Object.send(:remove_const, :CSV)
    CSV = FasterCSV
  else
    # CSV is now FasterCSV in ruby 1.9
  end

