# -*- coding: utf-8 -*-

#Importa Bibliotecas Necessárias
from math import log10, pi, sin, cos, ceil
from scipy import signal
import numpy as np


class butterworth:
    
    # Essa função inicializa os valores do filtro
    def __init__(self, tipo, Wp1, Ws1, Ap, As):
        self.tipo = tipo
        self.Wp = Wp1
        self.Ws = Ws1
        self.Ap = Ap
        self.As = As
        self.N = butterworth.ordem(self)
    
    
    
    # Essa função define e retorna a ordem do filtro
    def ordem(self):
        if self.tipo == "PB":
            n = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Ws/self.Wp))
        elif self.tipo == "PA":
            n = log10( (pow(10, (-self.As/10)) - 1) / (pow(10,(
                -self.Ap/10)) - 1)) / (2*log10(self.Wp/self.Ws))
        N = ceil(n)
        self.N = int(N)
        return int(N)
        
    
    # Essa função define e retorna a frequência de corte do filtro
    def freq_corte(self):
        if self.tipo == "PB":
            Wc = self.Wp / (-1 + 10**(-self.Ap/10))**(1/(2*self.N))
            self.Wc = Wc
            return Wc
        elif self.tipo == "PA":
            Wc = self.Wp * pow((pow(10, (-self.Ap/10))-1), (1/(2*self.N)))
            self.Wc = Wc
            return Wc

    
    # Essa função define e retorna as raízes do denominador da FT
    def raizes_unit(self):
        Sk = list()
        for k in range(1, self.N + 1):
            Sk.append( complex(-sin( (pi*(2*k-1)) / (2*self.N) ), cos( 
                (pi*(2*k-1)) / (2*self.N) )))
        self.Sk = Sk
        return Sk
    
    
    # Essa função define e retorna a FT do filtro
    def func_tranf(self):
        butterworth.raizes_unit(self)
        poli = list()
        poli = np.poly(self.Sk)
        coefReal = poli.real
        self.Wc = butterworth.freq_corte(self)
        
        if self.tipo == "PB":
            num, den = signal.lp2lp(coefReal[-1], coefReal, self.Wc)
            H = signal.TransferFunction(num, den)
        elif self.tipo == "PA":
            num, den = signal.lp2hp(coefReal[-1], coefReal, self.Wc)
            H = signal.TransferFunction(num, den)
        self.H = H
        return H
    

    
    
    
    
    
    
    
    
    
    
    
    