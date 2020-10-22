# -*- coding: utf-8 -*-

from funcoes_butt import *
from funcoes_cheb1 import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# Boas Vindas
print("\nBem Vindo(a)!\n")
print("--------------------------------------------------------------------")
print("Tipos de Filtros:")
print("(BW) Butterworth ou (C1) Chebyshev Tipo 1\n")
print("(PB) - Passa-Baixa ou (PA) - Passa-Alta")
print("--------------------------------------------------------------------")

# Pega o tipo de filtro e os pontos de projeto
tipo1 = input("Digite a SIGLA do tipo que deseja (BW) ou (C1): ")
tipo2 = input("Digite a SIGLA do tipo que deseja (PB) ou (PA): ")
print("\nAgora vamos aos Pontos de Projeto ... ")
Wp = float(input("Digite a Frequência de Passagem (Wp): "))
Ws = float(input("Digite a Frequência de Rejeição (Ws): "))
Ap = float(input("Digite a Atenuação de Passagem (Ap): "))
As = float(input("Digite a Atenuação de Rejeição (As): "))

ApAlt = -1
AsAlt = int(As - abs(Ap) - 1)

# Inicializa um Objeto da Classe Butterworth
if tipo1 == "BW":
    filtro = butterworth(tipo2, Wp, Ws, ApAlt, AsAlt)
else:
    filtro = chebyshev(tipo2, Wp, Ws, ApAlt, AsAlt)

# Define e exibe Função de Transferência
H = filtro.func_tranf()

G = 10**((Ap+1)/20)
H.num = H.num * G

print("--------------------------------------------------------------------")
print("Sua Função de Transferência H(s) é:")
print(H)
print("--------------------------------------------------------------------")

# Plotar Gráficos
print("Deseja Plotar os gráficos de Bode?")
resposta = input("'s' ou 'n' (sem aspas): ")
if resposta == 's':
        w, y, phase = H.bode(w = np.arange(0, max(Wp, Ws)+1000, step = 1))
        plt.figure(1)
        plt.grid(True)
        plt.title("Butterworth PB Ativo")
        plt.scatter(Wp, Ap)
        plt.scatter(Ws, As)
        plt.xlim(0, max(Wp, Ws)+1000)
        plt.ylim(min(Ap, As)-5, max(Ap, As)+5)
        plt.plot(w, y)
    
        # Plotagem da Fase
        plt.figure(2)
        plt.grid(True)
        plt.title("Butterworth PB Ativo")
        plt.xlim(0, 11000)
        plt.ylim(0, 100)
        plt.plot(w, phase, 'r')
elif resposta == 'n':
    print("Gráfico não iniciado!")
else:
    print("Resposta não identificada.")
print("--------------------------------------------------------------------\n")
