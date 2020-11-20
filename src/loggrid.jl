# r = (exp(x) for x in -30:0.1:20)
# r[1]: MethodError
# usage: 
# for rᵢ in LogRange(-30:0.1:20, scale = 1.2) ... end
# dx = lgrid.dx
# TODO: rangemap

struct LogRange
    xgrid::AbstractRange{T} where T
    xscale::Real 
    dx::Real
    function LogRange(r::AbstractRange{T}; xscale = 1.0) where T
        new(r, xscale, step(r))
    end
end

Base.eltype(lg::LogRange) = Float64
Base.length(lg::LogRange) = length(lg.xgrid)
Base.first(lg::LogRange) = exp(first(lg.xgrid))
Base.last(lg::LogRange) = exp(last(lg.xgrid))

function Base.iterate(lg::LogRange, state = 1)
    state > length(lg.xgrid) && return nothing
    x, state = iterate(lg.xgrid, state)
    return exp(x), state
end

Base.getindex(lg::LogRange, i::Integer) = exp(getindex(lg.xgrid, i))

# getindex(::LogRange, ::UnitRange{Int64})
# getindex(::LogRange, ::StepRange{Int64,Int64})



