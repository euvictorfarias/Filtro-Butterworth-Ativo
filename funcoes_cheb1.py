# -*- coding: utf-8 -*-

#Importa Bibliotecas Necessárias
from scipy import signal
import numpy as np


class chebyshev:
    
    # Iniciação do Objeto e seus parâmetros
    def __init__(self, tipo, Wp1, Ws1, Ap, As, Wp2 = 0, Ws2 = 0):
        self.tipo = tipo
        self.Wp = Wp1
        self.Ws = Ws1
        self.Ap = Ap
        self.As = As
        self.e = chebyshev.constProp(self)
        self.N = chebyshev.ordem(self)
        
    
    # Constante de Proporcionalidade
    def constProp(self):
        e = np.sqrt(pow(10, (-0.1*self.Ap)) - 1)
        self.e = e
        return e
    
    
    # Essa função define e retorna a ordem do filtro
    def ordem(self):
        if self.tipo == "PB":
            n = np.arccosh(np.sqrt(pow(10, (-0.1*self.As)) - 1) / self.e
                           ) / np.arccosh(self.Ws/self.Wp)
        elif self.tipo == "PA":
            n = np.arccosh(np.sqrt(pow(10, (-0.1*self.As)) - 1) / self.e
                           ) / np.arccosh(self.Wp/self.Ws)
        N = int(np.ceil(n))
        self.N = N
        return N
    
    
    # Essa função define e retorna as raízes do denominador da FT
    def raizes_unit(self):
        Sk = list()
        Ok = list()
        Wk = list()
        for k in range(1, self.N+1):
            Ok.append(-np.sinh((1/self.N) * np.arcsinh(1/self.e) 
                               ) * np.sin(np.pi/(2*self.N)*(2*k - 1)))
            Wk.append(np.cosh((1/self.N) * np.arcsinh(1/self.e) 
                              ) * np.cos(np.pi/(2*self.N)*(2*k - 1)))
            Sk.append(complex(Ok[k-1], Wk[k-1]))
        self.Sk = Sk
        return Sk
    
    
    # Essa função define e retorna a FT do filtro
    def func_tranf(self):
        chebyshev.raizes_unit(self)
        poli = np.poly(self.Sk)
        coef = poli.real
        D = list()
        aux = 0
        for i in range(-self.N, 1):
            D.append(coef[aux]*pow(self.Wp, i))
            aux = aux + 1
            
        if self.N % 2 != 0:
            if self.tipo == "PB":
                num, den = signal.lp2lp(coef[-1], coef, self.Wp)
                H = signal.TransferFunction(num[-1], den)
            elif self.tipo == "PA":
                num, den = signal.lp2hp(coef[-1], coef, self.Wp)
                H = signal.TransferFunction(num, den)
        else:
            aux = 1 / np.sqrt(1 + self.e**2)
            if self.tipo == "PB":
                num, den = signal.lp2lp(coef[-1], coef, self.Wp)
                H = signal.TransferFunction(num * aux, den)
            elif self.tipo == "PA":
                num, den = signal.lp2hp(coef[-1], coef, self.Wp)
                H = signal.TransferFunction(num * aux, den)
        self.H = H
        return H

    
    
    
    
    
    
    