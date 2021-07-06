module PassosMultiplos

export adam_bashfort, adam_moulton, preditor_corretor

# =====================================================================
#                     Métodos de Adams-Bashforth
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# r     inteiro, número de estágios (r = 2, 3, 4)
# x0    real, condição inicial da primeira variavel
# y0    real, condição inicial da segunda variavel
# fx    função lado direito do PVI para a primeira variavel
# fy    função lado direito do PVI para a segunda variavel
# ---------------------------------------------------------------------

# Saída:
# t, X, Y     solução aproximada do sistema linear

function adam_bashfort(fx,fy,a::Float64,b::Float64,N::Int64, x0::Float64, y0::Float64, r::Int64)
    k = (b-a)/N
    X = zeros(N+5)
	Y = zeros(N+5)
    t = zeros(N+5)
    X[1] = x0
	Y[1] = y0
    t[1] = a
	for i = 2:4 # Valores iniciais via RK 4ª ordem
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
	for i = 1:N	
		if r == 2
			X[i+2] = X[i+1] + k/2*( 3*fx(t[i+1], X[i+1], Y[i+1]) - fx(t[i+1], X[i+1], Y[i+1]))
			Y[i+2] = Y[i+1] + k/2*( 3*fy(t[i+1], X[i+1], Y[i+1]) - fy(t[i+1], X[i+1], Y[i+1]))
			t[i+2] = t[i+1] + k
		elseif r == 3
			X[i+3] = X[i+2] + k/12*(23*fx(t[i+2], X[i+2], Y[i+2]) - 16*fx(t[i+1], X[i+1], Y[i+1]) + 5*fx(t[i], X[i], Y[i]))
			Y[i+3] = Y[i+2] + k/12*(23*fy(t[i+2], X[i+2], Y[i+2]) - 16*fy(t[i+1], X[i+1], Y[i+1]) + 5*fy(t[i], X[i], Y[i]))
			t[i+3] = t[i+2] + k
		elseif r == 4
			X[i+4] = X[i+3] + k/4*(55*fx(t[i+3], X[i+3], Y[i+3]) - 59*fx(t[i+2], X[i+2], Y[i+3]) + 37*fx(t[i+1], X[i+1], Y[i+1]) - 9*fx(t[i], X[i], Y[i]))
			Y[i+4] = Y[i+3] + k/4*(55*fy(t[i+3], X[i+3], Y[i+3]) - 59*fy(t[i+2], X[i+2], Y[i+3]) + 37*fy(t[i+1], X[i+1], Y[i+1]) - 9*fy(t[i], X[i], Y[i]))
			t[i+4] = t[i+3] + k
		else
			error("r invalido para o Método de Adams-Bashforth")
		end
    end
    return t,X,Y
end

# =====================================================================
#                     Métodos de Adams-Moulton
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# r     inteiro, número de estágios (r = 2, 3)
# x0    real, condição inicial da primeira variavel
# y0    real, condição inicial da segunda variavel
# fx    função lado direito do PVI para a primeira variavel
# fy    função lado direito do PVI para a segunda variavel
# ---------------------------------------------------------------------

# Saída:
# t, X, Y     solução aproximada do sistema linear

function adam_moulton(fx,fy,a::Float64,b::Float64,N::Int64, x0::Float64, y0::Float64, r::Int64)
    k = (b-a)/N
    X = zeros(N+4)
	Y = zeros(N+4)
    t = zeros(N+4)
    X[1] = x0
	Y[1] = y0
    t[1] = a
	for i = 2:3 # Valores iniciais via RK 4ª ordem
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
	for i = 1:N	
		if r == 2
			X[i+2] = X[i+1] + k/12*(5*fx(t[i+2], X[i+2], Y[i+2]) + 8*fx(t[i+1], X[i+1], Y[i+1]) - fx(t[i], X[i], Y[i]));
			Y[i+2] = Y[i+1] + k/12*(5*fy(t[i+2], X[i+2], Y[i+2]) + 8*fy(t[i+1], X[i+1], Y[i+1]) - fy(t[i], X[i], Y[i]));
			t[i+2] = t[i+1] + k
		elseif r == 3
			X[i+3] = X[i+2] + k/24*(9*fx(t[i+3], X[i+3], Y[i+3]) + 19*fx(t[i+2], X[i+2], Y[i+2]) - 5*fx(t[i+1], X[i+1], Y[i+1]) + f(t[i], X[i], Y[i]))
			Y[i+3] = Y[i+2] + k/24*(9*fy(t[i+3], X[i+3], Y[i+3]) + 19*fy(t[i+2], X[i+2], Y[i+2]) - 5*fy(t[i+1], X[i+1], Y[i+1]) + f(t[i], X[i], Y[i]))
			t[i+3] = t[i+2] + k
		else
			error("r invalido para o Método de Adams-Moulton")
		end
    end
    return t,X,Y
end

# =====================================================================
#                     Método Preditor-Corretor
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

function preditor_corretor(fx,fy,a::Float64,b::Float64,N::Int64, x0::Float64, y0::Float64)
    k = (b-a)/N
    X = zeros(N+2)
	Y = zeros(N+2)
    t = zeros(N+2)
    X[1] = x0
	Y[1] = y0
    t[1] = a
	for i = 2:4 # Valores iniciais via RK 4ª ordem
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
	for i = 4:N+1
		t[i+1] = t[i] + k
		wx = X[i] + k/24*(55*fx(t[i], X[i], Y[i]) - 59*fx(t[i-1], X[i-1], Y[i-1]) + 37*fx(t[i-2], X[i-2], Y[i-2]) - 9*fx(t[i=3], X[i-3], Y[i-3]))
		wy = Y[i] + k/24*(55*fy(t[i], X[i], Y[i]) - 59*fy(t[i-1], X[i-1], Y[i-1]) + 37*fy(t[i-2], X[i-2], Y[i-2]) - 9*fy(t[i=3], X[i-3], Y[i-3]))
		X[i+1] = X[i] + k/24*(9*fx(t[i+1], wx, wy) + 19*fx(t[i], X[i], Y[i]) - 5*fx(t[i-1], Y[i-1]) + fx(t[i-2], X[i-2], Y[i-2]))
		Y[i+1] = Y[i] + k/24*(9*fy(t[i+1], wx, wy) + 19*fy(t[i], X[i], Y[i]) - 5*fy(t[i-1], Y[i-1]) + fy(t[i-2], X[i-2], Y[i-2]))
    end
    return t,X,Y
end


end