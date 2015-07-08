class CDP1802
  def initialize
    @opcodes = {:INC => 0x1, :DEC => 0x2, :IRX => 0x60, :GLO => 0x8, :PLO => 0xA,
      :GHI => 0x9, :PHI=> 0xB, :LDN => 0x0, :LDA => 0x4, :LDX => 0xF0, :LDXA => 0x72,
      :LDI => 0xF8, :STR => 0x5, :STXD => 0x73, :OR => 0xF1, :ORI => 0xF9, :XOR => 0xF3,
      :XRI => 0xFB, :AND => 0xF2, :ANI => 0xFA, :SHR => 0xF6, :SHRC => 0x76, :RSHR => 0x76,
      :SHL => 0xFE, :SHLC => 0x7E, :RSHL => 0x7E, :ADD => 0xF4, :ADI => 0xFC, :ADC => 0x74,
      :ADCI => 0x7C, :SD => 0xF5, :SDI => 0xFD, :SDB => 0x75, :SDBI => 0x7D, :SM => 0xF7,
      :SMI => 0xFF, :SMB => 0x77, :SMBI => 0x7F, :BR => 0x30, :NBR => 0x38, :SKP => 0x38,
      :BZ => 0x32, :BNZ => 0x3A, :BDF => 0x33, :BPZ => 0x33, :BGE => 0x33, :BNF => 0x3B,
      :BM => 0x3B, :BL => 0x3B, :BQ => 0x31, :BNQ => 0x39, :B1 => 0x34, :BN1 => 0x3C,
      :B2 => 0x35, :BN2 => 0x3D, :B3 => 0x36, :BN3 => 0x3E, :B4 => 0x37, :BN4 => 0x3F,
      :LBR => 0xC0, :NLBR => 0xC8, :LSKP => 0xC8, :LBZ => 0xC2, :LBNZ => 0xCA,
      :LBDF => 0xC3, :LBNF => 0xCB, :LBQ => 0xC1, :LBNQ => 0xC9, :LSKP => 0xC8,
      :NLBR => 0xC8, :LSZ => 0xCE, :LSNZ => 0xC6, :LSDF => 0xCF, :LSNF => 0xC7,
      :LSQ => 0xCD, :LSNQ => 0xC5, :LSIE => 0xCC, :IDL => 0x00, :NOP => 0xC4,
      :SEP => 0xD, :SEX => 0xE, :SEQ => 0x7B, :REQ => 0x7A, :SAV => 0x78,
      :MARK => 0x79, :RET =>0x70, :DIS => 0x71, :OUT => 0x6, :INP => 0x6}
    @opdesc = {:INC => "increment reg N", :DEC => "decrement reg N",  
      :IRX => "increment reg X", :GLO => "get low reg N", :PLO => "put low reg N", 
      :GHI => "get high reg N", :PHI=> "put high reg N", :LDN => "load via N", 
      :LDA => "load advance", :LDX => "load via X", :LDXA => "load via X and advance",
      :LDI => "load immediate", :STR => "store via N", :STXD => "store via X and decrement",
      :OR => "or ", :ORI => "or immediate", :XOR => "exclusive-or", 
      :XRI => "exclusive-or immediate", :AND => "and", :ANI => "and immediate", 
      :SHR => "shift right", :SHRC => "shift right with carry", :RSHR => "ring shift right",
      :SHL => "shift left", :SHLC => "shift left with carry", :RSHL => "ring shift left",
      :ADD => "add", :ADI => "add immediate", :ADC => "add with carry", 
      :ADCI => "add with carry immediate", :SD => "subtract D", :SDI => "subtract D immediate",
      :SDB => "subtract D with borrow", :SDBI => "subtract D with borrow immediate",
      :SM => "subtract memory", :SMI => "subtract memory immediate", 
      :SMB => "subtract memory with borrow", :SMBI => "subtract memory with borrow immediate",
      :BR => "unconditional short branch", :NBR => "no short branch", :SKP => "skip",
      :BZ => "short branch if D=0", :BNZ => "short branch if D not 0", 
      :BDF => "short branch if DF=0", :BPZ => "short branch if pos or zero",
      :BGE => "short branch if equal or greater", :BNF => "short branch if DF=0",
      :BM => "short branch if minus", :BL => "short branch if less",
      :BQ => "short branch if Q=1", :BNQ => "short branch if Q=0", 
      :B1 => "short branch if EF1=1", :BN1 => "short branch if EF1=0",
      :B2 => "short branch if EF2=1", :BN2 => "short branch if EF2=0",
      :B3 => "short branch if EF3=1", :BN3 => "short branch if EF3=0",
      :B4 => "short branch if EF4=1", :BN4 => "short branch if EF4=0",
      :LBR => "long branch", :NLBR => "no long branch", :LSKP => "long skip",
      :LBZ => "long branch if D=0", :LBNZ => "long branch if D not 0",
      :LBDF => "long branch if DF=1", :LBNF => "long branch if DF=0",
      :LBQ => "long branch if Q=1", :LBNQ => "long branch if Q=0", :LSKP => "long skip",
      :NLBR => "no long branch", :LSZ => "long skip if D=0", :LSNZ => "long skip if D NOT 0",
      :LSDF => "long skip if DF=1", :LSNF => "long skip id DF=0", 
      :LSQ => "long skip if Q=1", :LSNQ => "long skip if Q=0", :LSIE => "long skip if IE=1",
      :IDL => "idle", :NOP => "no operation", :SEP => "set P", :SEX => "set X",
      :SEQ => "set Q", :REQ => "reset Q", :SAV => "save", :MARK => "push X, P to stack",
      :RET => "return", :DIS => "disable", :OUT => "output", :INP => "input"}
  end
  def formatByte(byte, base)
    if (base == 16 || base == 10)
      return byte.to_s(base).upcase
    elsif (base == 2)
      bin = byte.to_s(2)
      if (bin.size < 8)
        bin = "0"*(8-bin.size) + bin
      end
      return bin
    else
      STDERR << "Unsupported base #{base}"
      exit(1)
    end
  end
  def assemble(code, base)
    labels = Hash.new
    asm = ""
    code.split("\n").each do |line|
      data, comments = line.split("..")
      address, rest = data.split(" ", 2)
      next if (address !~/[0-9]+/)
      if (rest =~/:/)
        lab, rest = rest.split(":")
        labels[lab] = address
      end
      opcode, operand = rest.split(" ")
      if (@opcodes[opcode.to_sym])
        byte = @opcodes[opcode.to_sym]
        if (byte < 0x10)
          operand = labels[operand] if (labels[operand])
          byte = 16*byte + operand.to_i if (operand)
          asm <<  formatByte(byte, base)<< "\n"
        else
          asm << formatByte(byte, base) << "\n"
          asm << formatByte(operand.to_i, base) << "\n" if (operand)
        end
      else
        STDERR << "Error at #{address}\n"
        exit(1)
      end
    end
    return asm
  end
  def disassemble(code, base)
     code.split("\n").each do |line|
       
     end
  end
end
