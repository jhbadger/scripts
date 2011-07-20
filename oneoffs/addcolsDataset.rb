Dir.glob("*/apis/dataset.tsv").each do |dataset|
  out = File.new(dataset, "w")
  File.new(dataset + ".orig").each do |line|
    name, owner, date_added, database_used, comments, group, username, password = line.chomp.split("\t")
    out.print [name, owner, date_added, database_used, comments, group, username, password].join("\t") + "\n"
  end
end