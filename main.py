import math
from scipy.interpolate import interp1d


class Performance:
    def __init__(self, io, scr, s1r, tl, kd, kd_performance, hn, kt, x, r, omega_r, cd, factor_delta_u, last_period=10,
                 interval=0.01):
        self.factor_delta_u = factor_delta_u
        self.interval = interval
        self.last_period = last_period
        self.cd = cd
        self.omega_r = omega_r
        self.r = r
        self.x = x
        self.kt = kt
        self.hn = hn
        self.kd = kd
        self.kd_performance = kd_performance
        self.tl = tl
        self.s1r = s1r
        self.scr = scr
        self.io = io

        fa = 1
        fv = 1
        na = 1
        nv = 1
        self.scs = self.scr * fa * na
        self.s1s = self.s1r * fv * nv
        self.ts = self.s1s / self.scs
        self.t0 = self.ts * 0.2

        self.t1 = kt * (hn ** x)
        if qh := (2 / 3) * (self.ts / self.t1) < 1:
            if qh > 0.5:
                self.qh = qh
            else:
                self.qh = 0.5
        else:
            self.qh = 1

        self.beta_i = 0.05

        self.sa_ta = self.spectrum_point(self.t1, self.kd)[1]
        self.cs = self.sa_ta / self.r
        self.delta_e = 100 * self.hn * self.factor_delta_u / self.cd
        self.delta_u = 100 * self.hn * self.factor_delta_u

        self.sa_e = self.cs * self.omega_r
        self.sd_e = self.delta_e
        self.sa_u = self.sa_e

        self.mu = None
        self.beta_h = None
        self.beta_effective = None
        self.beta_d = None
        self.w_effective = None
        self.t_effective = None
        self.sa_t_effective = None
        self.sd_t_effective = None

        self.sd_u = self.performance_point()

        self.performance_curve = [[0, 0],
                                  [self.sa_e, self.sd_e],
                                  [self.sa_u, self.sd_u]]

        print(self.performance_curve)

        # self.design_spectrum = list(self.complete_spectrum())

    # def complete_spectrum(self):
    #     t = 0
    #     while t < self.last_period:
    #         yield self.spectrum_point(t)

    def spectrum_point(self, t, kd):
        scd = self.scs * kd
        s1d = self.s1s * kd

        if t < self.t0:
            sa = scd * (0.4 + (0.6 * (t / self.t0)))
        elif t < self.ts:
            sa = scd
        elif t < self.tl:
            sa = min(s1d / t, scd)
        else:
            sa = (s1d / (t ** 2)) * self.tl
        t += self.interval
        return [t, sa]

    def beta_d_finder(self, beta_effective):
        # sourcery skip: inline-immediately-returned-variable
        beta_effective = beta_effective * 100
        effective_damping_table = [2, 5, 10, 20, 30, 40, 50]
        damping_coefficient_table = [0.8000, 1.0000, 1.2000, 1.5000, 1.7000, 1.9000, 2.0000]

        answer1 = None

        if beta_effective < 5:
            answer1 = ((damping_coefficient_table[1] - damping_coefficient_table[0]) /
                       (effective_damping_table[1] - effective_damping_table[0]) *
                       (beta_effective - effective_damping_table[0])) + damping_coefficient_table[0]
        elif beta_effective < 10:
            answer1 = ((damping_coefficient_table[2] - damping_coefficient_table[1]) /
                       (effective_damping_table[2] - effective_damping_table[1]) *
                       (beta_effective - effective_damping_table[1])) + damping_coefficient_table[1]
        elif beta_effective < 20:
            answer1 = ((damping_coefficient_table[3] - damping_coefficient_table[2]) /
                       (effective_damping_table[3] - effective_damping_table[2]) *
                       (beta_effective - effective_damping_table[2])) + damping_coefficient_table[2]
        elif beta_effective < 30:
            answer1 = ((damping_coefficient_table[4] - damping_coefficient_table[3]) /
                       (effective_damping_table[4] - effective_damping_table[3]) *
                       (beta_effective - effective_damping_table[3])) + damping_coefficient_table[3]
        elif beta_effective < 40:
            answer1 = ((damping_coefficient_table[5] - damping_coefficient_table[4]) /
                       (effective_damping_table[5] - effective_damping_table[4]) *
                       (beta_effective - effective_damping_table[4])) + damping_coefficient_table[4]
        else:
            answer1 = ((damping_coefficient_table[6] - damping_coefficient_table[5]) /
                       (effective_damping_table[6] - effective_damping_table[5]) *
                       (beta_effective - effective_damping_table[5])) + damping_coefficient_table[5]

        y_interp = interp1d(effective_damping_table, damping_coefficient_table)
        return y_interp(beta_effective)

    def performance_point(self):
        sd_u = self.sd_e
        while sd_u < self.delta_u:
            mu = sd_u / self.sd_e
            beta_h = self.qh * ((2 / math.pi) - self.beta_i) * (1 - (1 / mu))
            beta_effective = self.beta_i + beta_h
            beta_d = self.beta_d_finder(beta_effective)
            w_effective = (self.sa_u * 980 / sd_u) ** 0.5
            t_effective = 2 * math.pi / w_effective
            sa_t_effective = self.spectrum_point(t_effective, self.kd_performance)[1] / beta_d
            sd_t_effective = sa_t_effective * 980 / (w_effective ** 2)

            # equal = math.isclose(sd_u, sd_t_effective, rel_tol=0.0001), sd_u, beta_effective, sd_t_effective
            if math.isclose(sd_u, sd_t_effective, rel_tol=0.0001):
                self.mu = mu
                self.beta_h = beta_h
                self.beta_effective = beta_effective
                self.beta_d = beta_d
                self.w_effective = w_effective
                self.t_effective = t_effective
                self.sa_t_effective = sa_t_effective
                self.sd_t_effective = sd_t_effective
                break
            sd_u += 0.0001

        return sd_u


if __name__ == '__main__':
    # Chinautla - Suelo tipo D
    # Estructura Ordinaria, Concreto reforzado DL
    espectro_test = Performance(4.1, 1.43, 0.89, 3.48, 0.66, 0.8, 12.8, 0.049, 0.75, 4, 3, 2.5, 0.01, last_period=6)
    print(espectro_test.spectrum_point(6, 0.66))
    print(espectro_test.spectrum_point(6, 0.8))
