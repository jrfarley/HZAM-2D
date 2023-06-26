using SpecialFunctions

sigma_comp = 0.01

function linear_density(w)
    return sqrt(pi / 2) * sigma_comp * (erf((min(w + 0.03, 1) - w) / (sqrt(2) * sigma_comp)) - erf((max(w - 0.03, 0) - w) / (sqrt(2) * sigma_comp)))
end

println(linear_density(0.5))
