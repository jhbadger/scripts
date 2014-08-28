#!/usr/bin/env julia

using ArgParse
using DataFrames

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "--kind"
        help = "kind of merge to perform: left, right, inner, outer"
        default = "left"
        "table1"
        required = true
        arg_type = String
        "table2"
        required = true
        arg_type = String
        "output_name"
        required = true
        arg_type = String
        "common_field_name"
        required = true
        arg_type = String
    end
    
    return parse_args(s)
end

function cleanFrame(df)
    rows,cols = size(df)
    for i in 1:rows
        for j in 1:cols
            if isna(df[i,j])
                if contains(string(names(df)[j]),"Total")
                    df[i,j] = "0"
                else
                    df[i,j] = ""
                end
            end
        end
    end
    df
end

function main()
    opts = parse_commandline()
    types = fill(UTF8String, 2000)
    table1 = readtable(opts["table1"], nastrings=["NA"], eltypes=types)
    table2 = readtable(opts["table2"], nastrings=["NA"], eltypes=types)
    merged = join(table1, table2, on=symbol(opts["common_field_name"]),
                  kind=symbol(opts["kind"]))
    writetable(opts["output_name"], cleanFrame(merged))
end

main()
