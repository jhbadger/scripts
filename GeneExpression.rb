# routines to process gene expression data

class GeneExpression
  # generate a poisson random number
  def GeneExpression.poisson(lambda)
    return lambda if (lambda > 700)
    l = Math.exp(-lambda)
    k = 0
    p = 1
    while p >= l
      k += 1
      u = rand
      p = p *u
    end
    return k -1
  end
  
  # calculate Strekel (2000) R value
  def GeneExpression.rvalue(total, counts)
    r = 0
    sum = counts.reduce(:+)
    freq = sum/(1.0*total.reduce(:+))
    total.size.times do |i|
      if (!counts[i].nil? && counts[i] > 0)
        r += (counts[i]*Math.log((counts[i])/(total[i]*freq)))
      end
    end
    return (100*r).to_i/100.0 # round nicely
  end
end

# return list of samples, descs from a given experiment in the counts table
def samplesFromExperiment(db, experiment)
  sql = "SELECT DISTINCT sample, desc FROM sample_desc " 
  sql += "WHERE experiment='#{experiment}' " 
  sdesc = Hash.new
  db.execute(sql) do |row|
    sdesc[row.first] = row.last.to_s
  end
  return sdesc.keys.sort, sdesc
end

# return set of select sum statements for given field and set of samples
def samplePivot(table, samples, field, prefix = "")
  sql = ""
  samples.each do |sample|
    if (prefix == "")
      ali = sample
    else
      ali = prefix + "_" + sample
    end
    sql += "sum(case when #{table}.sample='#{sample}' then #{field} end) as '#{ali}', "
  end
  sql.chop.chop
end

#return set of fields in annotation table
def ann_headers
  ["jgi","ncbi","annotation","swissprot","swissprot_annotation","kegg","kegg_def",
    "gos_cluster","phytax_cluster","phytax","ko","ko_def","pfam","pfam_def",
    "pfam2go","tfam","tfam_def","tfam2go","cog","cog_def","plant","animal"]
end


#rank by size list pf values
def rankValues(values)
  ranks = Hash[values.sort.uniq.reverse.each_with_index.to_a]
  values.collect{ |v| ranks[v] + 1 }
end