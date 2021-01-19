# rectangle rule for integration is accurate enough
∫(h, f) = h*sum(f)
#=
function simpson(f, x)
    n = length(x)
    h = [x[i] - x[i-1] for i in 2:n]
    
    result = 0.0
    for i=2:2:n-1
        hph = h[i] + h[i-1]
        result += f[i] * (   h[i]^3 + h[i-1]^3 + 3h[i]   * h[i-1] * hph) / (6h[i]   * h[i-1]) +
                f[i-1] * (2h[i-1]^3 -   h[i]^3 + 3h[i]   * h[i-1]^2)     / (6h[i-1] * hph)    +
                f[i+1] * (  2h[i]^3 - h[i-1]^3 + 3h[i-1] * h[i]^2)       / (6h[i]   * hph)
    end

    if n % 2 == 0
        result += f[n]   * (2h[n-1]^2 + 3h[n-2] * h[n-1]) / (6 * (h[n-2] + h[n-1])) +
                  f[n-1] * ( h[n-1]^2 + 3h[n-1] * h[n-2]) /  6h[n-2] -
                  f[n-2] *   h[n-1]^3 / (6h[n-2] * (h[n-2] + h[n-1]))
    end

    return result    
end
=#