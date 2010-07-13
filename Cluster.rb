# Class to encapsulate clustering methods like cd-hit
class Cluster
  def initialize(seq, method, fraction, project = nil)
    @seq = seq
    @method = method
    @fraction = fraction
    @project = project
    return cluster
  end
  
  def cluster
    if (@method == "cd-hit" || @method == "cd-hit-est")
      @out = "#{@seq}_#{@fraction}_#{@method}"
      cmd = "#{@method} -i #{@seq} -d 0 -c #{@fraction} -o #{@out}"
      qsub = `which qsub`
      if (project.nil? || qsub == "")
        system(cmd)
      else
        
      end
      return @out
    end
  end
end
