class Mothur
  def Mothur.cleanup(dir)
    system("rm -rf #{dir}")	
  end
  
	def Mothur.classify(seqs, template = ENV["HOME"] + "/lib/mothur/silva/silva.bacteria.fasta", 
	  taxonomy=ENV["HOME"] + "/lib/mothur/silva/silva.rdp.taxonomy")
	  tmpdir = ENV["HOME"] + "/.tmp"
    Dir.mkdir(tmpdir) if !File.exists?(tmpdir)
    tmpdir += "/" + name + "_" + Time.now.to_f.to_s
    Dir.mkdir(tmpdir) if !File.exists?(tmpdir)
    out = File.new(tmpdir + "/seq.fasta", "w")
    num = 1
    seqs.each do |seq|
      out.printf(">seq%d\n%s\n", num, seq)
      num += 1
    end
    out.close
    batch = File.new(tmpdir + "/batch", "w")
    batch.printf("classify.seqs(fasta=seq.fasta, template=#{template}, taxonomy=#{taxonomy}, processors=2)")
    batch.close
    system("cd #{tmpdir};mothur batch >/dev/null")
    tax = Dir.glob("#{tmpdir}/*.taxonomy").first
    classes = []
    File.new(tax).each do |line|
      name, tx = line.chomp.split(" ", 2)
      classes.push(tx)
    end
    cleanup(tmpdir)
    return classes
	end
end
