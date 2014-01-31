# class to encapsulate running jobs on a Grid Engine (or locally on multiple CPUs/cores)
class Grid
  attr_accessor :name
  
  def initialize(command, project = nil, memory = nil, queue = nil, dir = Dir.pwd)
    @command = command
    @name = @command.gsub("/","_").split(" ").first
    @project = project
    @memory = memory
    @queue = queue
    @dir = dir
    Dir.mkdir(@dir) if !Dir.exists?(@dir)
    @count = 0
    @files = []
  end
  
  # return name for next data file
  def next
    @count += 1
    filename = @dir + "/" + @name + "_input." + @count.to_s
    @files.push(filename)
    filename
  end
  
  # write job script
  def writeJob
    out = File.new(@name + ".com", "w")
    out.printf("\#!/bin/sh\n")
    out.printf("\#$ -t 1-%d\n", @count)
    out.printf("cd %s\n", Dir.pwd)
    filename = @dir + "/" + @name + "_input."
    out.printf("%s %s$SGE_TASK_ID\n", @command, filename)
    out.close
    File.chmod(0777, @name + ".com")
  end
  
  # submit array to grid
  def submit(sync = false, local=false, verbose=false, maxLocal=4)
    writeJob
    if local
      1.upto(@count).each do |i|
        com = @name+".com"
        cmd = "export SGE_TASK_ID=#{i};bash #{@name}.com > #{@dir}/#{com}.o#{i} 2>#{@dir}/#{com}.e#{i}"
        STDERR << "Submitting #{cmd}...\n" if verbose
        fork do
          `#{cmd}` 
        end
        Process.waitall if i%maxLocal == 0
      end  
      Process.waitall
    else
      qsub = "qsub -t 1:" + @count.to_s
      qsub += " -sync yes" if (sync)
      qsub += " -P #{@project}" if @project
      printL = false
      if (@memory)
        qsub += " -l " if !printL
        qsub += "," if printL
        qsub +=  "memory=#{@memory}"
        printL = true
      end
      if (@queue)
        qsub += " -l " if !printL
        qsub += "," if printL
        qsub +=  @queue
        printL = true
      end
      qsub += " -o #{@dir} "
      qsub += " -e #{@dir} "
      qsub += " #{Dir.pwd}/#{@name}.com"
      system(qsub)
    end
  end

  # concatenate output, error, removing input files
  def cleanup
    Dir.glob(@dir + "/" + @name + "_input*").each do |file|
      File.unlink(file)
    end
    out = File.new(@name + ".out", "w")
    Dir.glob(@dir + "/" + @name + ".com.o*").each do |file|
      File.new(file).each do |line|
        out.print line
      end
      File.unlink(file)
    end
    out.close
    err = File.new(@name + ".err", "w")
    Dir.glob(@dir + "/" + @name + ".com.e*").each do |file|
      File.new(file).each do |line|
        err.print line
      end
      File.unlink(file)
    end
    err.close
    File.unlink(@name + ".com")
  end

  def join(fileendings, dataset)
    ofiles = Hash.new
    fileendings.each do |ending|
      ofiles[ending] = File.new(dataset + "_" + ending, "a")
    end
    @files.each do |file|
      fileendings.each do |ending|
        if File.exist?(file + "_" + ending)
          File.new(file + "_" + ending).each do |line|
            ofiles[ending].print line
          end
          File.unlink(file + "_" + ending)
        end
      end
    end
    ofiles.keys.each do |key|
      ofiles[key].close
    end
    File.unlink(File.basename(@name) + ".com*") if File.exists?(File.basename(@name) + ".com*")
  end
end
