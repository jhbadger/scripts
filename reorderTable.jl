#!/usr/bin/env julia

using ArgParse
using DataFrames

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "table"
        required = true
        arg_type = String
        "collist"
        required = true
        arg_type = String
        "output"
        required = true
        arg_type = String
    end
    
    return parse_args(s)
end

function main()
    opts = parse_commandline()
    cols = split(readall(opts["collist"]),"\n")
    table = readtable(opts["table"])
    df = DataFrame()
    for col in cols
        df[symbol(col)] = table[symbol(col)]
    end
    writetable(opts["output"], df, separator = '\t')
end

main()
