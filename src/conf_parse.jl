"""
    conf_enc(input::AbstractString; maxn = 7, spins = 2,
    notation = Dict("s" => 0, "p" => 1, "d" => 2, "f" => 3))

Parse an input string as the configuration of an atom. An electronic 
configuration is string like "[He] 2s2 2p1", where the notation [X] indicates 
that all subshells associated with the noble gas X are fully occupied. 
The function encodes the configuration provided "[He] 2s2 2p1" and returns 
`((2.0, 2.0), (1.0,))`, where subtuples `(2.0, 2.0)` and `(1.0,)` are 
populations of all the l = 0 and l = 1 shells, i.e. `((1s², 2s²), (2p¹,))`. 
Within each subtuple populations of levels with the same azimuthal quantum 
number l and subsequent radial quantum numbers nᵣ are listed.

Optional keyword arguments are:

* `maxn`  - maximal permitted principal quantum number
* `spins` - number of spin up/spin down electrons per level (1 or 2)
* `notation` - dictionary encoding the correspondence of quantum numbers l = 0, 1... to letters: s, p, d... See [spectroscopic notation](https://en.wikipedia.org/wiki/Spectroscopic_notation)

# Examples
```jldoctest
julia> conf = conf_enc("[Ar] 3d1 4s2") # Scandium, Z = 21
((2.0, 2.0, 2.0, 2.0), (6.0, 6.0), (1.0,))

julia> sum(collect(Iterators.flatten(conf)))
21.0
```
"""  
function conf_enc(input::AbstractString; 
    maxn=7, 
    spins=2,
    notation=Dict("s" => 0, "p" => 1, "d" => 2, "f" => 3))    
    # mnemonic for the s,p,d,f,g,h,i,k notation:
    # Sober Physicists Don't Find Giraffes Hiding In Kitchens 
    
    ls = join(keys(notation))

    p = raw"(^\s*(\[(He|Ne|Ar|Kr|Xe|Rn)\]\s*$)|" * 
        raw"^\s*(\[(He|Ne|Ar|Kr|Xe|Rn)\])?(\s+[1-9]+\d*" * 
        "[" * ls * "]" * 
        raw"{1}\d+\.?\d*)+\s*$|^[1-9]+\d*" * 
        "[" * ls * "]" * 
        raw"{1}\d+\.?\d*(\s+[1-9]+\d*" * 
        "[" * ls * "]" * 
        raw"{1}\d+\.?\d*)*\s*$)+"
    
    configur_pattern = Regex(p, "i")
    subshell_pattern = Regex(raw"(\d+)([" * ls * raw"]{1})(\d+\.?\d*)")
    
    # check syntax of the input string
    if !occursin(configur_pattern, input)
        throw(DomainError(input, 
        "input electron configuration doesn't match common pattern," * 
        "such as: \"[He] 2s2 2p1\" ... or \"1s2 2s2 2p1 ...\""))
    end
      
    # include all fully occupied subshells associated with the noble gas
    input = input |> strip |> lowercase |> 
    s -> replace(s, "[rn]" => "[xe] 4f14 5d10 6s2 6p6") |>
    s -> replace(s, "[xe]" => "[kr] 4d10 5s2 5p6") |>
    s -> replace(s, "[kr]" => "[ar] 3d10 4s2 4p6") |>
    s -> replace(s, "[ar]" => "[ne] 3s2 3p6") |>
    s -> replace(s, "[ne]" => "[he] 2s2 2p6") |>
    s -> replace(s, "[he]" => "1s2")
                                       
    configuration = []
    buff = zeros(length(notation), maxn)
    uniq = Set{Tuple{Int64,Int64}}()
    lmax = 0
    
    # parse subshells, check for physical sense 
    for subshell in split(input)
        _n, _l, _occ = match(subshell_pattern, subshell).captures
        n, l, occ = parse(Int, _n), notation[_l], parse(Float64, _occ)
        
        if !in(n, 1:maxn)
            throw(DomainError(subshell,
            "Principal quantum number is out of range 1:$maxn"))
            elseif l ≥ n
            throw(DomainError(subshell,
            "Azimuthal quantum number must be l < $n"))            
            elseif occ > (2l + 1) * spins
            throw(DomainError(subshell,
            "$_l-subshell cannot be occupied by more then" * 
            "$((2l + 1) * spins) electrons"))
            elseif in((n, l), uniq)
            throw(DomainError(input, 
            "The $n$_l subshell has already been encountered before."))
        end
    
        push!(uniq, (n, l))
        if l > lmax  
            lmax = l 
        end        
        buff[l + 1,n] = occ 
    end

    # construct electronic configuration tuple
    for l in 0:lmax
        subshell = []; flag = true
        for n in maxn:-1:(l + 1)
            if buff[l + 1, n] == 0.0 && flag 
                continue
            end
            flag = false
            push!(subshell, buff[l + 1, n])
        end
        push!(configuration, tuple(reverse(subshell)...))
    end        
    return tuple(configuration...)
end

Nₑ(conf) = sum(Iterators.flatten(conf))

#=
function conf_decode(conf)
    ls = ('s', 'p', 'd', 'f')
    s = []
    for l in 1:length(conf), n in 1:length(conf[l])
            push!(s, tuple(n + l - 1, l, conf[l][n]))
    end
    sort!(s, by = x -> x[1])
    return mapreduce(x -> "$(x[1])$(ls[x[2]])$(x[3]) ", *, s, init="") |> 
    s -> replace(s, "1s2.0 2s2.0 2p6.0 3s2.0 3p6.0 3d10.0 "  * 
                    "4s2.0 4p6.0 4d10.0 4f14.0 5s2.0 5p6.0 " * 
                    "5d10.0 6s2.0 6p6.0"            => "[Rn]") |> 
    s -> replace(s, "1s2.0 2s2.0 2p6.0 3s2.0 3p6.0 3d10.0 4s2.0 " * 
                    "4p6.0 4d10.0 5s2.0 5p6.0"      => "[Xe]") |>
    s -> replace(s, "1s2.0 2s2.0 2p6.0 3s2.0 3p6.0 3d10.0 4s2.0 " * 
                    "4p6.0"                         => "[Kr]") |>
    s -> replace(s, "1s2.0 2s2.0 2p6.0 3s2.0 3p6.0" => "[Ar]") |>
    s -> replace(s, "1s2.0 2s2.0 2p6.0"             => "[Ne]") |>
    s -> replace(s, "1s2.0"                         => "[He]") |>
    strip
end
=#