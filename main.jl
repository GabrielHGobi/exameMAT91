using Plots: display
# =========================================================================================
# *****************************************************************************************
#             			               MAT-91 2021 - Exame
#                          Aluno: Gabriel Henrique Gobi      COMP23
#                          Aluno: Jian Lucas Brito Veras      ELE23
# *****************************************************************************************
# =========================================================================================

# Tema: Problema de oscilações elétricas
# Modelar e resolver o problema de oscilações elétricas, usando os métodos vistos no curso.
# Referência: Seção 4.6 do livro ODE Architect Companion - The ultimate ODE power tool.


# Pacotes de Julia
using Printf
using Plots
using TypedTables

# Métodos para resolução de EDOs
push!(LOAD_PATH, pwd())
using PassosSimples
using PassosMultiplos

# Parâmetros do sistema
R = 1.5 # Valor da resistência em ohm
L = 12.0 # Valor da indutância em henry
C = 1.0 # Valor da capacitância em faraday
V = 1 # Valor da tensão, dependende do tempo em s, em volt

# Tempo de simulação em s
a, b = [0.0, 80.0]

# EDO do problema de oscilações elétricas (PVI)
dQ(t, Q, I) = I
dI(t, Q, I) = V/L - R/L * I - 1/(L*C) * Q

# Condições iniciais do PVI
Q0 = 0.0
I0 = 0.0

# Solução analítica (para V = cte)
alpha = R/(2*L) 
w = sqrt(1/(L*C) - (R/(2*L))^2 ) 
phi = atan(-alpha/w)
A = -V*C/cos(phi)
sol(t) = A*exp(-alpha*t)*cos(w*t+phi) + V*C
t_analitica = a:0.01:b
Q_analitica = @. sol(t_analitica)




# Resolução da EDO pelos diferentes métodos e análise



# ============================================================================
# Por métodos de passo simples

N = 50 # qtde de passos

plt_passos_simples = plot(t_analitica, Q_analitica,
		title = "Comparação de Métodos de Passo Simples (N = $N)",
        xlabel = "Tempo (s)",
        ylabel = "Carga Q no capacitor (C)",
        lcolor =:black ,
        ls =:dot,
        lw =:1.5,
        legend = :topright,
		label = "Analítica")

t, Q_euler, I_euler = PassosSimples.euler(dQ, dI, a, b, N, Q0, I0)
scatter!(t, Q_euler, 
    markershape = :diamond, 
    markersize =:3,
    mcolor = :blue,
    markerstrokewidth =:0.1,
    label= "Euler")
Q_analitica = @. sol(t)

t, Q_euler_melhorado, I_euler_melhorado = PassosSimples.euler_melhorado(dQ, dI, a, b, N, Q0, I0)
scatter!(t, Q_euler_melhorado, markershape = :circle, 
    markersize =:3,
    mcolor = :red,
    markerstrokewidth =:0.1, 
    label = "Euler Melhorado")

t, Q_rk4, I_rk4 = PassosSimples.runge_kutta_4_ordem(dQ, dI, a, b, N, Q0, I0)
scatter!(t, Q_rk4,
    markershape = :utriangle, 
    markersize =:4, 
    mcolor = :purple,
    markerstrokewidth =:0.1, 
    label= "RK 4ª ordem")

t_rkf, Q_rkf, I_rkf = PassosSimples.runge_kutta_fehlberg(dQ, dI, a, b, 0.000001, (b-a)/N, 0.01, Q0, I0)
scatter!(t_rkf, Q_rkf, 
    markershape = :star,
    markersize =:3, 
    mcolor = :yellow,
    markerstrokewidth =:0.1, 
    label = "RK Fehlberg")

savefig(plt_passos_simples, "./figures/PassoSimplesN$N")

dif_euler = broadcast(abs,Q_analitica - Q_euler)
dif_euler_melhorado = broadcast(abs,Q_analitica - Q_euler_melhorado)
dif_rk4 = broadcast(abs,Q_analitica - Q_rk4)

erros = plot(t, dif_euler, label = "Euler",
    title = "Erro (N = $N)",
    xlabel = "Tempo (s)",
    ylabel = "Erro")
plot!(t, dif_euler_melhorado, label = "Euler Melhorado")
plot!(t, dif_rk4, label = "RK 4ª ordem")
savefig(erros, "./figures/ErrosPassoSimplesN$N")

# Criando tabela para mostrar resultados
tabelaPassoSimples = Table(Tempo = t,
                    Analitica = Q_analitica,
                    Euler = Q_euler,
                    ErroEuler = dif_euler,
                    EulerMelhorado = Q_euler_melhorado,
                    ErroEulerMelhorado = dif_euler_melhorado,
                    RK4 = Q_rk4,
                    ErroRK4 = dif_rk4)
display(tabelaPassoSimples)

# ============================================================================
# Por métodos de passos múltiplos

Q_analitica = @. sol(t_analitica)

N = 50 # qtde de passos

plt_adam_b = plot(t_analitica, Q_analitica,
		title = "Comparação de estágios de Adam-Bashforth (N = $N)",
        xlabel = "Tempo (s)",
        ylabel = "Carga Q no capacitor (C)",
        lcolor =:black ,
        ls =:dot,
        lw =:1.5,
        legend = :topright,
		label = "Analítica")

t, Q_adam_b_2, I_adam_b_2 = PassosMultiplos.adam_bashfort(dQ, dI, a, b, N, Q0, I0, 2)
scatter!(t, Q_adam_b_2, 
    markershape = :diamond, 
    markersize =:3,
    mcolor = :blue,
    markerstrokewidth =:0.1,
    label= "r = 2")
Q_analitica = @. sol(t)

t, Q_adam_b_3, I_adam_b_3 = PassosMultiplos.adam_bashfort(dQ, dI, a, b, N, Q0, I0, 3)
scatter!(t, Q_adam_b_3, 
    markershape = :circle, 
    markersize =:3,
    mcolor = :red,
    markerstrokewidth =:0.1,
    label= "r = 3")

t, Q_adam_b_4, I_adam_b_4 = PassosMultiplos.adam_bashfort(dQ, dI, a, b, N, Q0, I0, 4)
scatter!(t, Q_adam_b_4, 
    markershape = :utriangle, 
    markersize =:3,
    mcolor = :purple,
    markerstrokewidth =:0.1,
    label= "r = 4")    
    
savefig(plt_adam_b, "./figures/PassoMultiploAdamBashforthN$N")



Q_analitica = @. sol(t_analitica)
plt_adam_m = plot(t_analitica, Q_analitica,
		title = "Comparação de estágios de Adam-Moulton (N = $N)",
        xlabel = "Tempo (s)",
        ylabel = "Carga Q no capacitor (C)",
        lcolor =:black ,
        ls =:dot,
        lw =:1.5,
        legend = :topright,
		label = "Analítica")

t, Q_adam_m_2, I_adam_m_2 = PassosMultiplos.adam_moulton(dQ, dI, a, b, N, Q0, I0, 2)
scatter!(t, Q_adam_m_2, 
    markershape = :diamond, 
    markersize =:3,
    mcolor = :blue,
    markerstrokewidth =:0.1,
    label= "r = 2")
Q_analitica = @. sol(t)

t, Q_adam_m_3, I_adam_m_3 = PassosMultiplos.adam_moulton(dQ, dI, a, b, N, Q0, I0, 3)
scatter!(t, Q_adam_m_3, 
    markershape = :circle, 
    markersize =:3,
    mcolor = :red,
    markerstrokewidth =:0.1,
    label= "r = 3")

savefig(plt_adam_m, "./figures/PassoMultiploAdamMoultonN$N")


# ==============================================

Q_analitica = @. sol(t_analitica)
plt_passos_multiplos = plot(t_analitica, Q_analitica,
		title = "Comparação de Métodos de Passos Múltiplos (N = $N)",
        xlabel = "Tempo (s)",
        ylabel = "Carga Q no capacitor (C)",
        lcolor =:black ,
        ls =:dot,
        lw =:1.5,
        legend = :topright,
		label = "Analítica")

t, Q_adam_b_4, I_adam_b_4 = PassosMultiplos.adam_bashfort(dQ, dI, a, b, N, Q0, I0, 4)
scatter!(t, Q_adam_b_4, 
    markershape = :diamond, 
    markersize =:3,
    mcolor = :blue,
    markerstrokewidth =:0.1,
    label= "Adam-Bashforth (r=4)")
Q_analitica = @. sol(t)

t, Q_adam_m_3, I_adam_m_3 = PassosMultiplos.adam_moulton(dQ, dI, a, b, N, Q0, I0, 3)
scatter!(t, Q_adam_m_3, 
    markershape = :circle, 
    markersize =:3,
    mcolor = :red,
    markerstrokewidth =:0.1, 
    label = "Adam-Moulton (r=3)")

t, Q_preditor_corretor, I_preditor_corretor = PassosMultiplos.preditor_corretor(dQ, dI, a, b, N, Q0, I0)
scatter!(t, Q_preditor_corretor,
    markershape = :utriangle, 
    markersize =:4, 
    mcolor = :purple,
    markerstrokewidth =:0.1, 
    label= "Preditor-Corretor")
savefig(plt_passos_multiplos, "./figures/PassoMultiploN$N")


dif_adam_b = broadcast(abs,Q_analitica - Q_adam_b_4)
dif_adam_m = broadcast(abs,Q_analitica - Q_adam_m_3)
dif_preditor_corretor = broadcast(abs,Q_analitica - Q_preditor_corretor)

erros = plot(t, dif_adam_b, label = "Adam-Bashforth",
    title = "Erro (N = $N)",
    xlabel = "Tempo (s)",
    ylabel = "Erro")
plot!(t, dif_adam_m, label = "Adam-Moulton")
plot!(t, dif_preditor_corretor, label = "Preditor-Corretor")
savefig(erros, "./figures/ErrosPassoMultiploN$N")

# Criando tabela para mostrar resultados
tabelaPassoMultiplo = Table(Tempo = t,
                        Analitica = Q_analitica,
                        AdamBashforth = Q_adam_b_4,
                        ErroAdamB = dif_adam_b,
                        AdamMoulton = Q_adam_m_3,
                        ErroAdamM = dif_adam_m,
                        PreditorCorretor = Q_preditor_corretor,
                        ErroPreditorCorretor = dif_preditor_corretor)
display(tabelaPassoMultiplo)



# ===============================================================



Q_analitica = @. sol(t_analitica)

N = 25 # qtde de passos

plt_final = plot(t_analitica, Q_analitica,
		title = "Comparação final (N = $N)",
        xlabel = "Tempo (s)",
        ylabel = "Carga Q no capacitor (C)",
        lcolor =:black ,
        ls =:dot,
        lw =:1.5,
        legend = :topright,
		label = "Analítica")


t_rkf, Q_rkf, I_rkf = PassosSimples.runge_kutta_fehlberg(dQ, dI, a, b, 100.0, (b-a)/N, 0.00001, Q0, I0)
scatter!(t_rkf, Q_rkf, 
    markershape = :diamond,
    markersize =:3, 
    mcolor = :blue,
    markerstrokewidth =:0.1, 
    label = "RK Fehlberg")
Q_analitica = @. sol(t)

println(length(Q_rkf))

t, Q_preditor_corretor, I_preditor_corretor = PassosMultiplos.preditor_corretor(dQ, dI, a, b, N, Q0, I0)
scatter!(t, Q_preditor_corretor,
    markershape = :circle, 
    markersize =:4, 
    mcolor = :red,
    markerstrokewidth =:0.1, 
    label= "Preditor-Corretor")

savefig(plt_final, "./figures/FinalN$N")

