# ExtractEHF.jl

function main()
    if length(ARGS) != 2
        println("Usage: julia ExtractEHF.jl input.txt output.csv")
        exit(1)
    end

    infile  = ARGS[1]
    outfile = ARGS[2]

    fields = ["t", "U", "F1", "F2", "KE1", "KE2", "PE"]

    num = raw"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

    kv_pattern = Regex(raw"\b([A-Za-z][A-Za-z0-9]*)\s*:\s*(" * num * raw")")

    rows = Vector{Vector{Float64}}()

    nlines = 0
    nskip = 0

    open(infile, "r") do io
        for line in eachline(io)
            nlines += 1

            values = Dict{String, Float64}()

            for m in eachmatch(kv_pattern, line)
                key = m.captures[1]
                val = parse(Float64, m.captures[2])
                values[key] = val
            end

            if all(f -> haskey(values, f), fields)
                push!(rows, [values[f] for f in fields])
            else
                nskip += 1
            end
        end
    end

    open(outfile, "w") do io
        println(io, join(fields, ","))

        for row in rows
            println(io, join(row, ","))
        end
    end

    println("Read $nlines lines")
    println("Extracted $(length(rows)) rows")
    println("Skipped $nskip lines")
end

main()
