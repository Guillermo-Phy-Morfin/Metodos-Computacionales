__precompile__() # Este comando es para que julia precompile el paquete

module AlgebraRepSU11

export J_0_rep, J_1_rep, J_2_rep

"""
Representación del elemento J₀ para la serie principal discreta.
- `k = -2j`, `k = 2, 3, 4, ...` (se excluye `j = -1/2` por ser un punto singular para la representación)
- `d+1` es la dimensión en la que se trunca la representación.
- `sign` determina si es el eigenvalor máximo (`-1`) o el mínimo (`1`).

# Ejemplo
```
julia> J_0_rep(2, 5, -1)
```
Esto produce una matriz 6x6 con entradas complejas.
"""
function J_0_rep(k, d , sign) #El signo determina si es la representación D^+ o D^-
    A = zeros(ComplexF64, d + 1, d + 1)

    for i in 0:d
        A[i+1,i+1] = (k/2 + i)
    end

    return sign * A
end

"""
Representación del elemento J₁ para la serie principal discreta.
- `k = -2j`, `k = 2, 3, 4, ...` (se excluye `j = -1/2` por ser un punto singular para la representación)
- `d+1` es la dimensión en la que se trunca la representación.
- `sign` determina si es el eigenvalor máximo (`-1`) o el mínimo (`1`).

# Ejemplo
```
julia> J_1_rep(2, 5, -1)
```
Esto produce una matriz 6x6 con entradas complejas.
"""
function J_1_rep(k, d, sign)  #El sign determina si es D_j^+ o D_j^-
    #j = -k/2 con k=1,2,3,...
    #m = sign(k/2+n) con n=0,1,2,3,...,d
    M = zeros(ComplexF64, d + 1, d + 1) # Matriz de (d+1 x d+1)
    
    for n in 0:d
        if n < d
            M[n+2, n+1] = sqrt((k + n) * (n + 1)) 
        end
        if n > 0
            M[n, n+1] = -sqrt(n * (k + n - 1)) 
        end
    end
    
    return sign * im / 2 * M
end

"""
Representación del elemento J₂ para la serie principal discreta.
- `k = -2j`, `k = 2, 3, 4, ...` (se excluye `j = -1/2` por ser un punto singular para la representación)
- `d+1` es la dimensión en la que se trunca la representación.
- En este caso no hay signo ya que la representación es la misma para el eigenvalor máximo o el mínimo.

# Ejemplo
```
julia> J_2_rep(2, 5)
```
Esto produce una matriz 6x6 con entradas complejas.
"""
function J_2_rep(k, d)  #en este caso la representación es la misma para D^+ y D_-
    #j = -k/2 con k=1,2,3,...
    #m = sign(k/2+n) con n=0,1,2,3,...,d
    M = zeros(ComplexF64, d + 1, d + 1) # Matriz de (d+1 x d+1)
    
    for n in 0:d
        if n < d
            M[n+2, n+1] = sqrt((k + n) * (n + 1)) 
        end
        if n > 0
            M[n, n+1] = sqrt(n * (k + n - 1)) 
        end
    end
    
    return 1 / 2 * M
end

end
     