# wrapper to handle compressed files transparently
class ZFile < File
  def self.new(fn,options="r")
    if (options == "r")
      if fn == "-"
        result = f = $stdin
        if block_given?
          result = yield f
        end
      elsif fn =~ /\.gz$/
        result = f = IO.popen("gunzip -c #{fn}")
        if block_given?
          result = yield f
        end
      elsif fn =~ /\.bz2$/
        result = f = IO.popen("bunzip2 -c #{fn}")
        if block_given?
          result = yield f
        end
      else
        result = f = File.new(fn,options)
        if block_given?
          result = yield f
        end
      end
    else
      result = f = File.new(fn,options)
      if block_given?
        result = yield f
      end
      f.close
    end
    return result
  end  
end
