class SuffixArray
  def initialize(the_string)
    @the_string = the_string
    @suffix_array = Array.new
    #build the suffixes
    last_index = the_string.length-1
    (0..last_index).each do |i|
      the_suffix = the_string[i..last_index]
      the_position = i
      # << is the append (or push) operator for arrays in Ruby
      @suffix_array << { :suffix=>the_suffix, :position=>the_position }
    end
      
    #sort the suffix array
    @suffix_array.sort! { |a,b| a[:suffix] <=> b[:suffix] }
  end
  
  def find_first_substring_index(the_substring, n_mismatches = 0)
    #first one found, not necessarily first on in the string
    finder(the_substring, n_mismatches, true)
  end
  
  def find_all_substring_indices(the_substring, n_mismatches = 0)
    finder(the_substring, n_mismatches, false)
  end
  
  private
  def finder(the_substring, n_mismatches, quit_after_first_result)
    results = []
    #uses typical binary search
    high = @suffix_array.length - 1
    low = 0
    while(low <= high)
      mid = (high + low) / 2
      this_suffix = @suffix_array[mid][:suffix]
      compare_len = the_substring.length-1
      comparison = this_suffix[0..compare_len]
      
      if n_mismatches == 0
        within_n_mismatches = comparison == the_substring
      else
        within_n_mismatches = hamming_distance(the_substring, comparison) <= n_mismatches
      end
      
      if within_n_mismatches
        results << @suffix_array[mid][:position]
        return results[0] if quit_after_first_result
      end
      
      if comparison > the_substring
        high = mid - 1
      else
        low = mid + 1
      end
      
      # the comparisons order in the original version:
      #if comparison > the_substring
      # high = mid - 1
      #elsif comparison < the_substring
      # low = mid + 1
      #else
      # return @suffix_array[mid][:position]
      #end
    end
  
    if quit_after_first_result
      return nil
    else
      return results
    end
  end
  
  def hamming_distance(a, b)
    # from Mladen Jablanović's answer at http://stackoverflow.com/questions/5322428/finding-a-substring-while-allowing-for-mismatches-with-ruby
    a.chars.zip(b.chars).count{|ca, cb| ca != cb}
  end
  
  
end