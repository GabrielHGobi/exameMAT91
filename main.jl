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

# Funções para resolução de EDOs
include("passos_simples.jl")
include("passos_multiplos.jl")

# Parâmetros do sistema
R = 1.5 # Valor da resistência em ohm
L = 12.0 # Valor da indutância em henry
C = 1.0 # Valor da capacitância em faraday
V(t) = 10.0*sin(t) # Valor da tensão, dependende do tempo em s, em volt

# Tempo de simulação em s
a, b = [0.0, 100.0]

# EDO do problema de oscilações elétricas (PVI)
dQ(t, Q, I) = I
dI(t, Q, I) = V(t)/L - R/L * I - 1/(L*C) * Q

# Condições iniciais do PVI
Q0 = 0.0
I0 = 0.0

# Resolução da EDO pelos diferentes métodos e análise

N = 200 # qtde de passos
t, Q, I = passos_simples.euler_melhorado(dQ, dI, a, b, N, Q0, I0)
plot(t, Q)



