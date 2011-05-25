require 'csv'

csv1 = CSV.read("/Users/jbadger/Desktop/go.csv", options = {:headers => true})
csv2 = CSV.read("/Users/jbadger/Desktop/clusters.csv", options = {:headers => true})

link = Hash.new
csv2.each do |row|
  goses = row["cluster_GO"].to_s.split(",")
  goses.each do |gos|
    link[gos] = [] if link[gos].nil?
    link[gos].push(row)
  end
end

print (csv1.headers + ["Cluster","Larger","Smaller"]*20).to_csv

csv1.each do |row|
  go, rest = row["annotation"].to_s.split(" ")
  if (link[go])
    link[go].each do |clust|
      row.push(clust["annotation"])
      row.push(clust["enriched larger"])
      row.push(clust["enriched smaller"])
    end
  end
  print row.to_csv
end
