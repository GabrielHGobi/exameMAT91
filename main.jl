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

module main

# Pacotes de Julia
using Printf
using Plots

# Funções para resolução de EDOs
include("passos_simples.jl")
include("passos_multiplos.jl")
# include("filename.jl")

# EDO do problema de oscilações elétricas (PVI)

R = 820 # Valor da resistência em ohm
L = 0.0047 # Valor da indutância em henry
C = 3.3*10^(-9) # Valor da capacitância em faraday
V(t) = 12 # Valor da tensão, dependende do tempo em s, em volt


# Resolução da EDO pelos diferentes métodos e análise

dQ(Q, I, t) = I
dI(Q, I, t) = V(t)/L - R/L * I + 1/(L*C) * Q

Q0 = 0
I0 = 0

end