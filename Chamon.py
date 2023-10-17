import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from scipy.integrate import quad
from scipy.integrate import simps
from scipy import interpolate
import math
rc("text", usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['figure.figsize'] = (8, 5)
plt.rcParams['font.size'] = 20

class Chamon:
    def __init__(self, z1, a1, z2, a2, E, sigma_data, pot_data):
        """ Input parameters are atomic and mass number of the projectile and target, E_CM.
            sigma_data must contain experimental values in 3 columns: theta_CM, sigma, delta_sigma.
            pot_data must contain potential data in 3 columns: R, V_C, V_SPP2.

            This class is used to convert from theta_CM to a (distance of closest approach).
            Default values for b (impact parameter) are from 0 to 30 fm in 0.1 fm steps. """
        self.z1         = z1
        self.z2         = z2
        self.a1         = a1
        self.a2         = a2
        self.E          = E
        self.mu         = a1*a2/(a1+a2) * 939
        self.v          = np.sqrt(2 * self.E / self.mu)
        self.theta      = sigma_data[:,0] * np.pi/180
        self.sigma      = sigma_data[:,1]
        self.sigma_err  = sigma_data[:,2]
        self.R          = pot_data[:,0]
        self.Vc         = pot_data[:,1]
        self.Vspp2      = pot_data[:,2]
        self.pot        = self.Vc + self.Vspp2

    def a_coul(self):
        """ Analytical formula to convert from theta_CM to a -- Coulomb only"""
        return 0.5*self.z1*self.z2*1.44/self.E*(1+1/np.sin(self.theta/2))

    def s_coul(self, MSR=0):
        """ a_coul but considering the size of the nucleus
            -> if MSR (mean square radius) = 0, calculates nuclear radius using heavy nuclei systematics:
               R = 1.31A^(1/3) - 0.84
            -> otherwise, calculates the radius using:
               R = sqrt(5/3)*MSR """
        a  = self.a_coul()
        if MSR==0:
            R1 = 1.31*self.a1**(1/3)-0.84
            R2 = 1.31*self.a2**(1/3)-0.84
        else:
            R1 = np.sqrt(5/3)*MSR
            R2 = np.sqrt(5/3)*3.96
        s = a - (R1 + R2)
        return s

    def potential(self, b = 30):
        """ Plots the nuclear + Coulomb potential """
        E_arr   = np.zeros_like(self.pot) + self.E
        L       = self.mu*self.v*b

        plt.figure()
        plt.plot(self.R, E_arr, color="red", linestyle="--")
        plt.plot(self.R, self.pot + L ** 2 / (2 * self.mu * self.R ** 2))
        plt.xlabel("r (fm)")
        plt.ylabel("V(r) (MeV)")
        plt.ylim(-200,200)
        plt.tight_layout()
        plt.show()

    def calculate_a(self):
        """ Calculates a for each value of b """
        b      = np.arange(0,30,0.1)
        E_arr  = np.zeros_like(self.pot) + self.E
        f_defl = []
        a      = []
        for i in b:
            L     = self.mu*self.v*i
            idx   = np.argwhere(np.diff(np.sign(self.pot + L ** 2 / (2 * self.mu * self.R ** 2) - E_arr))).flatten()
            a_aux = self.R[idx][-1]
            a.append(a_aux)
        return a

    def deflection(self, plot=True):
        """ Calculates and plots the deflection function """
        E_arr = np.zeros_like(self.pot) + self.E
        a = self.calculate_a()
        b = np.arange(0,30,0.1)
        f_deflec       = []
        for i in b:
            L = self.mu*self.v*i
            idx = np.argwhere(np.diff(np.sign(self.pot + L ** 2 / (2 * self.mu * self.R ** 2) - E_arr))).flatten()
            a_aux = self.R[idx][-1]
            f_to_integrate = []
            for k in range(len(self.R)):
                if self.R[k] <= a_aux + 0.5:
                    pass
                else:
                    fun = L / (self.R[k] ** 2 * np.sqrt(2 * self.mu * (self.E - self.pot[k] - L ** 2 / (2 * self.mu * self.R[k] ** 2))))
                    f_to_integrate.append(fun)
            I1 = simps(f_to_integrate, dx=self.R[1]-self.R[0])
            f = lambda r: L/(r**2*np.sqrt(2*self.mu*(self.E-self.z1*self.z2*1.44/r-L**2/(2*self.mu*r**2))))
            I2 = quad(f, self.R[-1], np.inf)
            I = I1 + I2[0]
            f_deflec.append((np.pi - 2*I)*180/np.pi)

        if plot:
            plt.figure()
            plt.plot(b, f_deflec)
            plt.xlabel("b (fm)")
            plt.ylabel(r"$\theta$ (deg)")
            plt.tight_layout()
            plt.show()
        return f_deflec

    def inv_deflection(self, plot=True):
        """ Calculates and plots the inverse of the deflection function """
        a              = self.calculate_a()
        f_deflec       = self.deflection(plot=False)
        a_nuc          = []
        a_coul         = []
        f_deflec_nuc   = []
        f_deflec_coul  = []
        for k in range(len(a)):
            if a[k]<5:
                a_nuc.append(a[k])
                f_deflec_nuc.append(f_deflec[k])
            else:
                a_coul.append(a[k])
                f_deflec_coul.append(f_deflec[k])
        a_inter        = np.linspace(a_nuc[-1], a_coul[0], 100)
        a_inter        = a_inter[1:-1]
        f_deflec_inter = np.linspace(f_deflec_nuc[-1], f_deflec_coul[0], 100)
        f_deflec_inter = f_deflec_inter[1:-1]
        a_nuc          = np.array(a_nuc)
        f_deflec_nuc   = np.array(f_deflec_nuc)

        if plot:
            plt.figure()
            plt.plot(f_deflec_nuc, a_nuc, color="red", label="Nuclear component")
            plt.plot(f_deflec_coul, a_coul, color="blue", label="Coulomb component")
            plt.plot(f_deflec_inter, a_inter, color="green", label="Interference component")
            plt.legend(fontsize=14, loc="upper left")
            plt.ylabel("a (fm)")
            plt.xlabel(r"$\theta$ (deg)")
            plt.tight_layout()
            plt.show()
        return a_inter, f_deflec_inter, a_nuc, f_deflec_nuc

    def a_inter(self):
        a, f, _,_ = self.inv_deflection(plot=False)
        coef      = np.polyfit(f, a, 1)
        poly1d_fn = np.poly1d(coef)
        return poly1d_fn(self.theta*180/np.pi)

    def a_nuc(self):
        _,_, a, f = self.inv_deflection(plot=False)
        coef      = np.polyfit(f, a, 1)
        poly1d_fn = np.poly1d(coef)
        return poly1d_fn(self.theta*180/np.pi)

    def plot_sigma_coul(self):
        """ Plots sigma vs a (coulomb component) """
        x = self.a_coul()
        plt.errorbar(x, self.sigma, yerr=self.sigma_err, label="$E_{CM}$ = " + str(self.E) + " MeV", linestyle="", capsize=3, marker="x", color="tab:red")
        plt.axhline(y=1, color="black", linestyle="-.", linewidth=1)
        plt.grid(alpha=0.4, axis="x")
        plt.xlabel("a (fm)")
        plt.ylabel("$\sigma/\sigma_R$")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_sigma_inter(self):
        """ Plots sigma vs a (interference component) """
        x = self.a_inter()
        plt.errorbar(x, self.sigma, yerr=self.sigma_err, label="$E_{CM}$ = " + str(self.E) + " MeV", linestyle="",capsize=3, marker="x", color="tab:red")
        plt.axhline(y=1, color="black", linestyle="-.", linewidth=1)
        plt.grid(alpha=0.4, axis="x")
        plt.xlabel("a (fm)")
        plt.ylabel("$\sigma/\sigma_R$")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_sigma_nuc(self):
        """ Plots sigma vs a (nuclear component) """
        x = self.a_nuc()
        plt.errorbar(x, self.sigma, yerr=self.sigma_err, label="$E_{CM}$ = " + str(self.E) + " MeV", linestyle="",capsize=3, marker="x", color="tab:red")
        plt.axhline(y=1, color="black", linestyle="-.", linewidth=1)
        plt.grid(alpha=0.4, axis="x")
        plt.xlabel("a (fm)")
        plt.ylabel("$\sigma/\sigma_R$")
        plt.legend()
        plt.tight_layout()
        plt.show()