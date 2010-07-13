# routines to process gene expression data


# generate a poisson random number
def poisson(lambda)
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
def rvalue(total, counts)
  r = 0
  sum = total.keys.reduce(0) {|s, lib| s += counts[lib]}
  freq = sum/(1.0*total.values.reduce(:+))
  total.keys.each do |lib|
    if (!counts[lib].nil? && counts[lib] > 0)
      r += (counts[lib]*Math.log((counts[lib])/(total[lib]*freq)))
    end
  end
  return r
end
