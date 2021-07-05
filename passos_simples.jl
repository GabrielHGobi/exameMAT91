module passos_simples

export euler, euler_melhorado, runge_kutta_4_ordem. runge_kutta_fehlberg

# =====================================================================
#                         Método de Euler
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# u     real, condição inicial
# f     função lado direito do PVI
# ---------------------------------------------------------------------
# Saída:
# U     solução aproximada de u' = f(t,u)

function euler(f,a::Float64,b::Float64,N::Int64,u::Float64)
    k = (b-a)/N
    U = zeros(N+1)
    t = zeros(N+1)
    U[1] = u
    t[1] = a
    for i = 1:N
        U[i+1] = U[i] + k*f(t[i],U[i])
        t[i+1] = t[i] + k
    end
    return t,U
end

# =====================================================================
#                       Método de Euler Melhorado
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# u     real, condição inicial
# f     função lado direito do PVI
# ---------------------------------------------------------------------

# Saída:
# U     solução aproximada de u' = f(t,u)

function euler_melhorado(f,a::Float64,b::Float64,N::Int64,u::Float64)
    k = (b-a)/N
    U = zeros(N+1)
    t = zeros(N+1)
    U[1] = u
    t[1] = a
    for i = 1:N
        U[i+1] = U[i] + k*f(t[i]+k/2,U[i]+k/2*f(t[i], U[i]))
        t[i+1] = t[i] + k
    end
    return t,U
end

# =====================================================================
#                  Método de Runge-Kutta de 4ª ordem
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# u     real, condição inicial
# f     função lado direito do PVI
# ---------------------------------------------------------------------

# Saída:
# U     solução aproximada de u' = f(t,u)

function runge_kutta_4_ordem(f,a::Float64,b::Float64,N::Int64,u::Float64)
    k = (b-a)/N
    U = zeros(N+1)
    t = zeros(N+1)
    U[1] = u
    t[1] = a
    for i = 1:N
        F0 = f(t[i], U[i])
        F1 = f(t[i] + k/2, U[i] + 1/2*F0)
        F2 = f(t[i] + k/2, U[i] + 1/2*F1)
        F3 = f(t[i] + k, U[i] + F2)
        U[i+1] = U[i] + k/6*(F0+2*F1+2*F2+F3)
        t[i+1] = t[i] + k
    end
    return t,U
end

# =====================================================================
#                  Método de Runge-Kutta-Fehlberg
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# Kmax  real, tamanho do passo máximo
# Kmin  real, tamanho do passo mínimo
# tol   real, tolerância máxima per
# u     real, condição inicial
# f     função lado direito do PVI
# ---------------------------------------------------------------------

# Saída:
# U     solução aproximada de u' = f(t,u)

function runge_kutta_fehlberg(f,a::Float64,b::Float64,tol::Float64, Kmax::Float64, Kmin::Float64, u::Float64)
    N = trunc(Int, (b-a)/Kmin)
    U = zeros(N+1)
    t = zeros(N+1)
    flag::Bool = true
    i = 1
    delta = 0.0
    k = Kmax
    U[1] = u
    t[1] = a
    while flag
        F0 = k*f(t[i], U[i])
        F1 = k*f(t[i] + k/4, U[i] + 1/4*F0)
        F2 = k*f(t[i] + 3*k/8, U[i] + 3/32*F0 + 9/32*F1)
        F3 = k*f(t[i] + 12*k/13, U[i] + 1932/2197*F0 - 7200/2197*F1 + 7296/2197*F2)
        F4 = k*f(t[i] + k, U[i] + 439/216*F0 - 8*F1 + 3680/513*F2 - 845/4104*F3)
        F5 = k*f(t[i] + k/2, U[i] - 8/27*F0 + 2*F1 - 3544/2565*F2 + 1859/4104*F3 - 11/40*F4)

        R = 1/k*abs(1/360*F0 - 128/4275*F2 - 2197/75240*F3 + 1/50*F4 + 2/55*F5)

        if R <= tol
            t[i+1] = t[i] + k
            U[i+1] = U[i] + 25/216*F0 + 1408/2565*F2 + 2197/4104*F3 - 1/5*F4
            i = i+1
        end

        delta = 0.84*(tol/R)^(1/4)
        if delta <= 0.1
            k = 0.1*k
        end
        if delta >= 4
            k = 4*k
        else
            k = delta*k
        end

        if k > Kmax
            k = Kmax
        end

        if t[i] >= b
            flag = false
        elseif t[i]+k > b
            k = b - t[i]
        end

        if k < Kmin
            error("Falha: k minimo")
        end
    end
    resize!(t, i)
    resize!(U, i)
    return t, U
end



end