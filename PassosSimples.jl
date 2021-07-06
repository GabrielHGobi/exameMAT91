module PassosSimples

export euler, euler_melhorado, runge_kutta_4_ordem, runge_kutta_fehlberg

# =====================================================================
#                         Método de Euler
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# x0    real, condição inicial da primeira variavel
# y0    real, condição inicial da segunda variavel
# fx    função lado direito do PVI para a primeira variavel
# fy    função lado direito do PVI para a segunda variavel
# ---------------------------------------------------------------------
# Saída:
# t, X, Y     solução aproximada do sistema linear

function euler(fx, fy,a::Float64,b::Float64,N::Int64,x0::Float64, y0::Float64)
    k = (b-a)/N
    X = zeros(N+1)
    Y = zeros(N+1)
    t = zeros(N+1)
    X[1] = x0
    Y[1] = y0
    t[1] = a
    for i = 1:N
        X[i+1] = X[i] + k*fx(t[i],X[i], Y[i])
        Y[i+1] = Y[i] + k*fy(t[i],X[i], Y[i])
        t[i+1] = t[i] + k
    end
    return t,X,Y
end

# =====================================================================
#                       Método de Euler Melhorado
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# x0     real, condição inicial da primeira variavel
# y0     real, condição inicial da segunda variavel
# fx     função lado direito do PVI para a primeira variavel
# fy     função lado direito do PVI para a segunda variavel
# ---------------------------------------------------------------------

# Saída:
# t, X, Y     solução aproximada do sistema linear

function euler_melhorado(fx, fy,a::Float64,b::Float64,N::Int64,x0::Float64, y0::Float64)
    k = (b-a)/N
    X = zeros(N+1)
    Y = zeros(N+1)
    t = zeros(N+1)
    X[1] = x0
    Y[1] = y0
    t[1] = a
    for i = 1:N
        X[i+1] = X[i] + k*fx(t[i]+k/2, X[i]+k/2*fx(t[i], X[i], Y[i]), Y[i]+k/2*fy(t[i], X[i], Y[i]))
        Y[i+1] = Y[i] + k*fy(t[i]+k/2, X[i]+k/2*fx(t[i], X[i], Y[i]), Y[i]+k/2*fy(t[i], X[i], Y[i]))
        t[i+1] = t[i] + k
    end
    return t,X,Y
end

# =====================================================================
#                  Método de Runge-Kutta de 4ª ordem
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# x0     real, condição inicial da primeira variavel
# y0     real, condição inicial da segunda variavel
# fx     função lado direito do PVI para a primeira variavel
# fy     função lado direito do PVI para a segunda variavel
# ---------------------------------------------------------------------

# Saída:
# t, X, Y     solução aproximada do sistema linear

function runge_kutta_4_ordem(fx, fy,a::Float64,b::Float64,N::Int64,x0::Float64, y0::Float64)
    k = (b-a)/N
    X = zeros(N+1)
	Y = zeros(N+1)
    t = zeros(N+1)
    X[1] = x0
	Y[1] = y0
    t[1] = a
    for i = 1:N
        F0x = fx(t[i], X[i], Y[i])
		F0y = fy(t[i], X[i], Y[i])
        F1x = fx(t[i] + k/2, X[i] + 1/2*F0x, Y[i] + 1/2*F0y)
		F1y = fy(t[i] + k/2, X[i] + 1/2*F0x, Y[i] + 1/2*F0y)
        F2x = fx(t[i] + k/2, X[i] + 1/2*F1x, Y[i] + 1/2*F1y)
		F2y = fy(t[i] + k/2, X[i] + 1/2*F1x, Y[i] + 1/2*F1y)
        F3x = fx(t[i] + k, X[i] + F2x, Y[i] + F2y)
		F3y = fy(t[i] + k, X[i] + F2x, Y[i] + F2y)
        X[i+1] = X[i] + k/6*(F0x+2*F1x+2*F2x+F3x)
		Y[i+1] = Y[i] + k/6*(F0y+2*F1y+2*F2y+F3y)
		t[i+1] = t[i] + k	
    end
    return t,X,Y
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
# x0     real, condição inicial da primeira variavel
# y0     real, condição inicial da segunda variavel
# fx     função lado direito do PVI para a primeira variavel
# fy     função lado direito do PVI para a segunda variavel
# ---------------------------------------------------------------------

# Saída:
# t, X, Y     solução aproximada do sistema linear

function runge_kutta_fehlberg(fx, fy, a::Float64,b::Float64,tol::Float64, Kmax::Float64, Kmin::Float64, x0::Float64, y0::Float64)
    N = trunc(Int, (b-a)/Kmin)
    X = zeros(N+1)
	Y = zeros(N+1)
    t = zeros(N+1)
    X[1] = x0
	Y[1] = y0
    t[1] = a
    flag::Bool = true
    i = 1
    delta = 0.0
    k = Kmax
    while flag
        F0x = k*fx(t[i], X[i], Y[i])
        F0y = k*fy(t[i], X[i], Y[i])
        F1x = k*fx(t[i] + k/4, X[i] + 1/4*F0x, Y[i] + 1/4*F0y)
        F1y = k*fy(t[i] + k/4, X[i] + 1/4*F0x, Y[i] + 1/4*F0y)
        F2x = k*fx(t[i] + 3*k/8, X[i] + 3/32*F0x + 9/32*F1x, Y[i] + 3/32*F0y + 9/32*F1y)
        F2y = k*fy(t[i] + 3*k/8, X[i] + 3/32*F0x + 9/32*F1x, Y[i] + 3/32*F0y + 9/32*F1y)
        F3x = k*fx(t[i] + 12*k/13, X[i] + 1932/2197*F0x - 7200/2197*F1x + 7296/2197*F2x, Y[i] + 1932/2197*F0y - 7200/2197*F1y + 7296/2197*F2y)
        F3y = k*fy(t[i] + 12*k/13, X[i] + 1932/2197*F0x - 7200/2197*F1x + 7296/2197*F2x, Y[i] + 1932/2197*F0y - 7200/2197*F1y + 7296/2197*F2y)
        F4x = k*fx(t[i] + k, X[i] + 439/216*F0x - 8*F1x + 3680/513*F2x - 845/4104*F3x, Y[i] + 439/216*F0y - 8*F1y + 3680/513*F2y - 845/4104*F3y)
        F4y = k*fy(t[i] + k, X[i] + 439/216*F0x - 8*F1x + 3680/513*F2x - 845/4104*F3x, Y[i] + 439/216*F0y - 8*F1y + 3680/513*F2y - 845/4104*F3y)
        F5x = k*fx(t[i] + k/2, X[i] - 8/27*F0x + 2*F1x - 3544/2565*F2x + 1859/4104*F3x - 11/40*F4x, Y[i] - 8/27*F0y + 2*F1y - 3544/2565*F2y + 1859/4104*F3y - 11/40*F4y)
        F5y = k*fy(t[i] + k/2, X[i] - 8/27*F0x + 2*F1x - 3544/2565*F2x + 1859/4104*F3x - 11/40*F4x, Y[i] - 8/27*F0y + 2*F1y - 3544/2565*F2y + 1859/4104*F3y - 11/40*F4y)

        F0 = [F0x, F0y]
        F1 = [F1x, F1y]
        F2 = [F2x, F2y]
        F3 = [F3x, F3y]
        F4 = [F4x, F4y]
        F5 = [F5x, F5y]

        R = 1/k*norm(1/360*F0 - 128/4275*F2 - 2197/75240*F3 + 1/50*F4 + 2/55*F5)

        if R <= tol
            t[i+1] = t[i] + k
            X[i+1] = X[i] + 25/216*F0x + 1408/2565*F2x + 2197/4104*F3x - 1/5*F4x
            Y[i+1] = Y[i] + 25/216*F0y + 1408/2565*F2y + 2197/4104*F3y - 1/5*F4y
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
    resize!(X, i)
    resize!(Y, i)
    return t, X, Y
end



end