def runGrid(cmd, project, queue = "default", memory = 2)
  ecmd = cmd.tr(" ()&;", "_")
  qsub = "qsub -P #{project} -o #{ecmd}.out -e #{ecmd}.err -cwd "
  if (queue != "default")
    qsub += "-l \"#{queue},memory=#{memory}G\" "
  else
    qsub += "-l \"memory=#{memory}G\" "
  end
  qsub += " \"#{cmd}\""
  system(qsub)
end

def blast(type, query, db, evalue, project, queue)
  blastFile = File.basename(query) + "_vs_" + File.basename(db) + ".#{type}"
  cmd = "pblastNCBI -b #{type} -p #{project} -e #{evalue} -q #{queue} #{db} #{query} -z 1000 &"
  if (!File.exists?(blastFile))
    system(cmd)
  end
end
