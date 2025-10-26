import numpy as np
from inspect import signature

class V_ff:
    def __init__(self):
        self.output = 0.0
        # self.w_2 = [1, 1]
        # self.x_2 = [-1/np.sqrt(3), 1/np.sqrt(3)]
        # self.w_3 = [5/9, 8/9, 5/9]
        # self.x_3 = [-3/5, 0, 3/5]

        self.w = [None, None, [1, 1], [5/9, 8/9, 5/9]]
        self.x = [None, None, [-1/np.sqrt(3), 1/np.sqrt(3)], [-np.sqrt(3/5), 0, np.sqrt(3/5)]]
        self.mode = 2
    
    # funkcja do wywolania, sama okresla czy 1d czy 2d
    def quadrature(self, ff, a, b, a_2=None, b_2=None, n=100, mode=2):
        self.output = 0.0
        params_count = len(signature(ff).parameters)
        if mode != 2 and mode != 3:
            print("mozna calkowac tylko w trybie 2 lub 3 wagi; automatycznie ustawiono 2")
        else:
            self.mode = mode

        if params_count==1:
            self.quadrature_1d(ff,a,b,n)  
        else:
            if (a_2, b_2) == (None, None):
                self.quadrature_2d(ff,a,b,a,b,n)
            else:
                self.quadrature_2d(ff,a,b,a_2,b_2,n)  
    
    def quadrature_1d(self, ff, a, b, n):
        step = (b-a)/n
        for _ in range(n):
            b = a + step
            self.output += self.quadrature_1d_subdiv(ff, self.w[self.mode],self.x[self.mode], a, b)
            a = a + step
        
    def quadrature_1d_subdiv(self, ff, w, x, a, b):
        sum = 0
        for w_, x_ in zip(w, x):
            sum += w_ * ff(((b - a) * x_ + (b + a)) * 0.5)
        
        return (b-a) * sum / 2
    
    def quadrature_2d(self, ff, a, b, a_2, b_2, n):
        step = (b-a)/n
        step_2 = (b_2-a_2)/n
        a_2_save = a_2
        for _ in range(n):
            b = a + step
            a_2 = a_2_save
            for _ in range(n):
                b_2 = a_2 + step_2
                self.output += self.quadrature_2d_subdiv(ff, self.w[self.mode],self.x[self.mode], a, b, a_2, b_2)
                a_2 = a_2 + step_2
            a = a + step

    def quadrature_2d_subdiv(self, ff, w, x, a, b, a_2, b_2):
        sum = 0
        for wi, xi in zip(w, x):
            for wj, xj in zip(w, x):
                arg_y = ((b_2 - a_2) * xi + (b_2 + a_2)) * 0.5
                arg_x = ((b - a) * xj + (b + a)) * 0.5
                sum += wi * wj * ff(arg_y, arg_x)
        
        # to jest ten jacaobian
        # w sumie logiczne ze trzeba dwa razy skalowac po dwoch wymiarach
        jacobian = ((b_2 - a_2) / 2.0) * ((b - a) / 2.0)
        return sum * jacobian
        