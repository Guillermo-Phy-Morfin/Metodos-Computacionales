__precompile__() # Este comando es para que julia precompile el paquete

module JacobiP

export jacobi_poli

"""
jacobi_poly(n, α, β, z)

Polinomios de Jacobi en la variable z;
de orden n y parametros α, β
"""
function jacobi_poli(n, α, β, z)
    if n == 0
        return 1.0
    elseif n == 1
        return 0.5 * ((α - β) + (α + β + 2) * z)
    end

    Pnm2 = 1.0
    Pnm1 = 0.5 * ((α - β) + (α + β + 2) * z)
    Pn = 0.0

    for k in 2:n
        a_k = 2.0 * k * (k + α + β) * (2k + α + β - 2)
        b_k = (2k + α + β - 1) * ((α^2 - β^2) + (2k + α + β) * (2k + α + β - 2) * z)
        c_k = 2.0 * (k + α - 1) * (k + β - 1) * (2k + α + β)
        d_k = (2k) * (k + α + β)

        Pn = (b_k * Pnm1 - c_k * Pnm2) / a_k
        Pnm2 = Pnm1
        Pnm1 = Pn
    end

    return Pn
end

end