require 'Stats'

class GA

  def initialize(genomeSize = 10, alphabet = "AGCT", popSize = 100, 
                 maxGen = 100, mutRate = 0.1)
    @genomeSize = genomeSize
    @alphabet = alphabet
    @popSize = popSize
    @fitness = Hash.new
    @mutRate = mutRate
    @maxGen = maxGen
    @perfect = false
    @popSize.times {
      genome = createGenome
      @fitness[genome] = fitness(genome)
    }
    @genomes = @fitness.keys.sort{|x,y| @fitness[x]<=>@fitness[y]}
  end
  
  def fitness(genome)
    target = "AGCTTAGCTT"
    @perfect = true if (genome == target)
    fitness = 0
    @genomeSize.times {|i|
      fitness += 1 if (genome[i] == target[i])
    }
    return fitness
  end

  def createGenome
    genome = ""
    @genomeSize.times {
      genome += @alphabet[rand(@alphabet.length)].chr
    }
    return genome
  end

  def asexualRep(genome)
    newGenome = " " * @genomeSize
    @genomeSize.times {|i|
      if (rand() > @mutRate)
        newGenome[i] = genome[i].chr
      else
        newGenome[i] = @alphabet[rand(@alphabet.length)].chr
      end
    }
    return newGenome
  end

  def mate(mama, dada)
    newGenome = " " * @genomeSize
    @genomeSize.times {|i|
      mallele = 
      if (rand() > @mutRate)
        if (rand() > 0.5)
          newGenome[i] = mama[i, 1]
        else
          newGenome[i] = dada[i, 1]
        end
      else
        newGenome[i] = @alphabet[rand(@alphabet.length)].chr
      end
    }
    return newGenome
  end

  def pickParent
    spinner = rand * @totFit
    seenFit = 0
    @genomes.each {|genome|
      seenFit += @fitness[genome]
      return genome if (spinner < seenFit)
    }
    return pickParent
  end

  def nextGen    
    @totFit = 0
    @genomes = @fitness.keys.sort{|x,y| @fitness[x]<=>@fitness[y]}
    @genomes.each {|genome|
      @totFit += @fitness[genome]
    }
    newPop = Hash.new
    (@popSize - 1).times {
      child = mate(pickParent(), pickParent())
      newPop[child] = fitness(child)
    }
    newPop[@bestGenome] = @fitness[@bestGenome]
    @fitness = newPop
  end

  def show(genome)
    return genome
  end

  def insert(genome)
    p @genomes[0]
    p @fitness[@genomes[0]]
  end

  def run
    gen = 0
    while(gen < @maxGen)
      gen +=1 
      @bestGenome = @fitness.keys.sort{|x, y| @fitness[x]<=>@fitness[y]}.last
      printf("%s %8.3f\n", show(@bestGenome), @fitness[@bestGenome])
        printf("%4d: av Fitness = %8.3f best Fitness = %8.3f\n", 
               gen, @fitness.values.average, @fitness.values.max)
      break if (@perfect)
      nextGen
    end
  end
end
