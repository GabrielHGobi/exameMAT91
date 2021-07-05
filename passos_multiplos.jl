module passos_multiplos

export adam_bashfort. adam_moulton, preditor_corretor

# =====================================================================
#                     Métodos de Adams-Bashforth
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# u     real, condição inicial
# r     inteiro, número de estágios (r = 2, 3, 4)
# f     função lado direito do PVI
# ---------------------------------------------------------------------

# Saída:
# U     solução aproximada de u' = f(t,u)

function adam_bashfort(f,a::Float64,b::Float64,N::Int64, u::Float64, r::Int64)
    k = (b-a)/N
    U = zeros(N+5)
    t = zeros(N+5)
	U[1] = u
	t[1] = a
	for i = 2:4 # Valores iniciais via RK 4ª ordem
		F0 = f(t[i-1], U[i-1])
        F1 = f(t[i-1] + k/2, U[i-1] + 1/2*F0)
        F2 = f(t[i-1] + k/2, U[i-1] + 1/2*F1)
        F3 = f(t[i-1] + k, U[i-1] + F2)
        U[i] = U[i-1] + k/6*(F0+2*F1+2*F2+F3)
        t[i] = t[i-1] + k
	end
	for i = 1:N	
		if r == 2
			U[i+2] = U[i+1] + k/2*( 3*f(t[i+1], U[i+1]) - f(t[i], U[i]) )
			t[i+2] = t[i+1] + k
		elseif r == 3
			U[i+3] = U[i+2] + k/12*(23*f(t[i+2], U[i+2]) - 16*f(t[i+1], U[i+1]) + 5*f(t[i], U[i]))
			t[i+3] = t[i+2] + k
		elseif r == 4
			U[i+4] = U[i+3] + k/4*(55*f(t[i+3], U[i+3]) - 59*f(t[i+2], U[i+2]) + 37*f(t[i+1], U[i+1]) - 9*f(t[i], U[i]))
			t[i+4] = t[i+3] + k
		else
			error("r invalido para o Método de Adams-Bashforth")
		end
    end
    return t,U
end

# =====================================================================
#                     Métodos de Adams-Moulton
# =====================================================================
# ---------------------------------------------------------------------
# Dados de entrada:
# a     real, tempo inicial
# b     real, tempo final
# N     inteiro, número de passos
# u     real, condição inicial
# r     inteiro, número de estágios (r = 2, 3)
# f     função lado direito do PVI
# ---------------------------------------------------------------------

# Saída:
# U     solução aproximada de u' = f(t,u)

function adam_moulton(f,a::Float64,b::Float64,N::Int64, u::Float64 , r::Int64)
    k = (b-a)/N
    U = zeros(N+4)
    t = zeros(N+4)
	U[1] = u
	t[1] = a
	for i = 2:3 # Valores iniciais via RK 4ª ordem
		F0 = f(t[i-1], U[i-1])
        F1 = f(t[i-1] + k/2, U[i-1] + 1/2*F0)
        F2 = f(t[i-1] + k/2, U[i-1] + 1/2*F1)
        F3 = f(t[i-1] + k, U[i-1] + F2)
        U[i] = U[i-1] + k/6*(F0+2*F1+2*F2+F3)
        t[i] = t[i-1] + k
	end
	for i = 1:N	
		if r == 2
			U[i+2] = U[i+1] + k/12*(5*f(t[i+2], U[i+2]) + 8*f(t[i+1], U[i+1]) - f(t[i], U[i]));
			t[i+2] = t[i+1] + k
		elseif r == 3
			U[i+3] = U[i+2] + k/24*(9*f(t[i+3], U[i+3]) + 19*f(t[i+2], U[i+2]) - 5*f(t[i+1], U[i+1]) + f(t[i], U[i]))
			t[i+3] = t[i+2] + k
		else
			error("r invalido para o Método de Adams-Moulton")
		end
    end
    return t,U
end

# =====================================================================
#                     Método Preditor-Corretor
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

function preditor_corretor(f,a::Float64,b::Float64,N::Int64, u::Float64)
    k = (b-a)/N
    U = zeros(N+2)
    t = zeros(N+2)
	U[1] = u
	t[1] = a
	for i = 2:4 # Valores iniciais via RK 4ª ordem
		F0 = f(t[i-1], U[i-1])
        F1 = f(t[i-1] + k/2, U[i-1] + 1/2*F0)
        F2 = f(t[i-1] + k/2, U[i-1] + 1/2*F1)
        F3 = f(t[i-1] + k, U[i-1] + F2)
        U[i] = U[i-1] + k/6*(F0+2*F1+2*F2+F3)
        t[i] = t[i-1] + k
	end
	for i = 4:N+1
		t[i+1] = t[i] + k
		w = U[i] + k/24*(55*f(t[i], U[i]) - 59*f(t[i-1], U[i-1]) + 37*f(t[i-2], U[i-2]) - 9*f(t[i=3], U[i-3]))
		U[i+1] = U[i] + k/24*(9*f(t[i+1], w) + 19*f(t[i], U[i]) - 5*f(t[i-1], U[i-1]) + f(t[i-2], U[i-2]))
    end
    return t,U
end


end