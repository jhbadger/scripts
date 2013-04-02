# routines to aid dealing with McCrow's Tree stuff

require 'Colors'
require 'Newick'

def loadTax(taxFile)
  tax = Hash.new
  sizes = Hash.new
  File.new(taxFile).each do |line|
    id, tx = line.chomp.split(" ", 2)
    tx.gsub!(/\([0-9]*\)/, "")
    if (tx =~/;/)
      tx = tx.split(";")
      tx = tx[0] + ";" + tx[1]
      tax[id] = tx
      sizes[tx] = 0 if sizes[tx].nil?
      sizes[tx] += 1      
    end
  end
  sizes.keys.each do |key|
    if (sizes[key] < 10)
      tax.keys.each do |id|
        if (tax[id] == key)
          tax[id] = "Misc phylum"
        end
      end
    end
  end
  return tax
end

def makeColors(tax)
  colors = {"Bacteria;Firmicutes"=>[176, 23, 31], "Bacteria;Actinobacteria"=>[255, 193, 37], 
    "Bacteria;Verrucomicrobia"=>[0, 250, 154], "Bacteria;Tenericutes"=>[32, 178, 170], 
    "Bacteria;Proteobacteria"=>[184, 184, 184], "Bacteria;Deinococcus-Thermus"=>[255, 181, 197], 
    "Bacteria;Bacteroidetes"=>[180, 82, 205], "Bacteria;Acidobacteria"=>[58, 95, 205], 
    "Bacteria;Aquificae"=>[142, 142, 56]}
  col = colors.values
  tax.values.sort.uniq.each do |taxon|
    if (!colors[taxon])
      colors[taxon] = pickDistantColor(col)
      col.push(colors[taxon])
    end
  end
  colors.keys.each do |taxon|
    colors[taxon] = colors[taxon].hex
  end
  return colors
end

def makeInfo(tree, tax, colors, highlights)
  info = tree + ".info"
  infotmp = tree + ".info.tmp"
  system("treeinfo.pl #{tree} > #{infotmp}")
  infoF = File.new(info, "w")
  File.new(infotmp).each do |line|
    fields = line.split(" ", 4)
    num, dist, type, sp = line.chomp.split("\t")
    if (type == "r")
      tx = tax[sp]
      if (tx.nil?)
        s = sp.split("__")[1]
        tx = tax[s]
        if (tx.nil?)
          g, s = s.to_s.split(" ")
          tax.keys.each do |key|
            if (key=~/#{g}/)
              tx = tax[key]
              break
            end
          end
        end
      end
      if (tx.nil?)
        STDERR.printf("No taxonomy for %s\n", sp)
      end
      if (highlights[sp])
        type = "q"
      end
      line = num.to_s + "\t" + dist.to_s + "\t" + type.to_s +  "\t" + sp.to_s + "\t" + tx.to_s + "\n"
    end
    infoF.print line
  end
  infoF.printf("ttclr\tq\torange\tblue\n")
  colors.keys.each do |key|
    infoF.printf("htax\t%s\t%s\n", key, colors[key])
  end
  infoF.close
  File.unlink(infotmp)
  return info
end


def addLegend(info, colors, title)
  svg = info + ".svg"
  tmp = info + ".svgtmp"
  system("treeinfo2svg.pl #{info} > #{tmp}")
  svgF = File.new(svg, "w")
  y = -150
  File.new(tmp).each {|line|
    if line =~/<\/svg>/
      colors.keys.each {|key|
        y += 210
        svgF.printf("<rect x=\"100\" y=\"#{y}\" width=\"6000\"  height=\"300\" fill=\"#{colors[key]}\"/>
        <text x=\"100\" y=\"#{y+200}\" font-size=\"290\" font-weight=\"bold\" font-family=\"Verdana\" fill=\"black\">#{key}</text>\n")
      }
    else
      svgF.print line
    end
  }
  svgF.printf("<text x=\"7300\" y=\"700\" font-size=\"500\" font-weight=\"bold\" font-family=\"Verdana\" 
  fill=\"black\">#{title}</text>")
  svgF.printf("</svg>\n")
  svgF.close
  File.unlink(tmp)
  return svg
end
