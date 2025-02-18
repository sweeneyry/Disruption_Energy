import numpy as np
import matplotlib.pyplot as plt
import pdb

class PlasmaProfiles:

    def __init__(self, Radiusq2, MinorRadius, MajorRadius, RadiusVessel, Pressure0, Jz0, BzVac):
        self.dummy=1. 
        self.rq2 = Radiusq2
        self.a = MinorRadius
        self.b = RadiusVessel
        self.R = MajorRadius
        self.L = 2.0 * np.pi * self.R
        self.p0 = Pressure0
        self.Jz0 = Jz0
        self.r = np.linspace(0, self.b, num=50000)
        self.dr = self.r[1] - self.r[0]
        self.indrq2 = np.argmin(np.abs(self.r - self.rq2))
        self.inda = np.argmin(np.abs(self.r - self.a))
        self.indb = -1
        self.mu0 = np.pi * 4e-7
        self.BzVac = BzVac

        # get profiles
        self.pr0 = self.pressure_Profile()
        self.gradpr0 = self.grad_p_Profile()
        self.Jz0r = self.Jz0_Profile()
        self.Btheta0r = self.Btheta0_Profile()
        self.Jtheta0r = self.Jtheta_Profile()
        self.jzr0 = self.jz_Profile()
        self.jthetar0 = self.jtheta_Profile()
        self.deltaBzr0, self.Bzr0 = self.Bz_Profile()
        self.q = self.q_Profile()

        # get scalars
        self.Ip = self.total_Current()
        self.Wmagp = self.total_Poloidal_Magnetic_Energy()
        self.Wmagt = self.total_Toroidal_Magnetic_Energy()
        self.phit = self.toroidal_Flux()
        self.Wth = self.thermal_Energy()

    def pressure_Profile(self):

        pr = np.zeros(np.shape(self.r))
        pr[0:self.indrq2] = self.p0
        pr[self.indrq2:self.inda] = self.p0 * (1.0 - (self.r[self.indrq2:self.inda] - self.rq2) / 
                                               (self.a - self.rq2))
        
        return pr
    
    def grad_p_Profile(self):

        gradpr = np.zeros(np.shape(self.r))
        gradpr[self.indrq2:self.inda] = -self.p0 / (self.a - self.rq2)

        return gradpr
    
    def Jz0_Profile(self):

        Jzr = np.zeros(np.shape(self.r))
        Jzr[0:self.indrq2] = self.Jz0
        Jzr[self.indrq2:self.inda] = self.Jz0 * (1.0 - (self.r[self.indrq2:self.inda] - self.rq2) / 
                                               (self.a - self.rq2))
        
        return Jzr
    
    def Btheta0_Profile(self):


        Btheta0r = np.zeros(np.shape(self.r))
        prefactor = self.mu0 * self.Jz0
        # in the below expression only we have incorporated the 1/r already. 
        Btheta0r[0:self.indrq2] = prefactor * self.r[0:self.indrq2] / 2.0
        Btheta0r[self.indrq2:self.inda] = prefactor * (self.rq2**2 / 2.0 + 
                                           (1.0 + self.rq2 / (self.a - self.rq2)) * (self.r[self.indrq2:self.inda]**2 - self.rq2**2) / 2.0
                                            -(self.r[self.indrq2:self.inda]**3 - self.rq2**3) / (3.0 * (self.a - self.rq2)))
        Btheta0r[self.inda:] = prefactor * (self.rq2**2 / 2.0 + 
                                            (1.0 + self.rq2 / (self.a - self.rq2)) * (self.a**2 - self.rq2**2) / 2.0
                                            -(self.a**3 - self.rq2**3) / (3.0 * (self.a - self.rq2)))
        # now we account for the 1/r for the rest of the domain
        Btheta0r[self.indrq2:] = Btheta0r[self.indrq2:] / self.r[self.indrq2:]

        return Btheta0r
    
    def Jtheta_Profile(self):

        Jthetar = self.Jz0r * self.Btheta0r / self.BzVac

        return Jthetar

    def jz_Profile(self):

        jzr = -self.gradpr0 * (self.Btheta0r / (self.BzVac**2 - self.Btheta0r**2))

        return jzr
    
    def jtheta_Profile(self):

        jthetar = self.gradpr0 * (self.BzVac / (self.BzVac**2 - self.Btheta0r**2))

        return jthetar
    

    def Bz_Profile(self):

        Bzr = np.zeros(np.shape(self.r)) + self.BzVac 
        deltaBzr = np.flip( self.mu0 * np.cumsum( np.flip((self.Jtheta0r[0:self.inda] 
                                                  + self.jthetar0[0:self.inda]) * self.dr )) )
        Bzr[0:self.inda] = Bzr[0:self.inda] + deltaBzr
        
        return deltaBzr, Bzr
    
    def q_Profile(self):

        qr = self.r * self.Bzr0 / (self.R * self.Btheta0r)

        return qr
    
    def total_Current(self):

        Ip = 2 * np.pi * np.sum( (self.Jz0r + self.jzr0) * self.r * self.dr)

        return Ip
    
    def total_Poloidal_Magnetic_Energy(self):

        Wmagp = 2.0 * np.pi * self.L * np.sum( self.Btheta0r**2 / (2.0 * self.mu0) * self.r * self.dr)

        return Wmagp
    
    def total_Toroidal_Magnetic_Energy(self):

        Wmagt = 2.0 * np.pi * self.L * np.sum( self.Bzr0**2 / (2.0 * self.mu0) * self.r * self.dr)

        return Wmagt
    
    def toroidal_Flux(self):

        phit = 2.0 * np.pi * np.sum( self.Bzr0 * self.r * self.dr)

        return phit
    
    def thermal_Energy(self):

        Wth = 2.0 * np.pi * self.L * self.p0 * (self.rq2**2 / 2.0 + 
                                            (1.0 + self.rq2 / (self.a - self.rq2)) * (self.a**2 - self.rq2**2) / 2.0
                                            -(self.a**3 - self.rq2**3) / (3.0 * (self.a - self.rq2)))
        
        return Wth
    
    def plot_Profiles(self):

        fig = plt.figure(figsize=(8, 8))
        ax1 = plt.subplot(5,2,1)
        ax2 = plt.subplot(5,2,2)
        ax3 = plt.subplot(5,2,3)
        ax4 = plt.subplot(5,2,4)
        ax5 = plt.subplot(5,2,5)
        ax6 = plt.subplot(5,2,6)
        ax7 = plt.subplot(5,2,7)
        ax8 = plt.subplot(5,2,8)
        ax9 = plt.subplot(5,2,9)

        ax1.plot(self.r, self.pr0/1e3)
        ax1.set_ylabel('Pressure (kPa)')
        ax1.set_title(str(np.round(self.Wth/1e6, decimals=2)) +  ' MJ')

        ax2.plot(self.r, self.gradpr0)
        ax2.set_ylabel(r'$\nabla p$')

        ax3.plot(self.r, self.Jz0r/1e6)
        ax3.set_ylabel(r'$J_{z}$ (MA/m$^2$)')
        ax3.set_title(str(np.round(self.Ip/1e6, decimals=2)) + ' MA')

        ax4.plot(self.r, self.Jtheta0r/1e6)
        ax4.set_ylabel(r'$J_{\theta}$ (MA/m$^2$)')

        ax5.plot(self.r, self.jzr0/1e6)
        ax5.set_ylabel(r'$j_{z}$ (MA/m$^2$)')        

        ax6.plot(self.r, self.jthetar0/1e6)
        ax6.set_ylabel(r'$j_{\theta}$ (MA/m$^2$)') 

        ax7.plot(self.r, self.Btheta0r)
        ax7.set_ylabel(r'$B_{\theta}$')
        ax7.set_title(str(np.round(self.Wmagp/1e6, decimals=1)) +  ' MJ')

        ax8.plot(self.r, self.Bzr0)
#        ax8.plot(self.r[0:self.inda], self.deltaBzr0)
        ax8.set_ylabel('$B_z$ (T)')
        ax8.set_title(str(np.round(self.Wmagt/1e6, decimals=1)) +  ' MJ, ' + 
                      str(np.round(self.phit, decimals=2)) + ' Wb')
        
        ax9.plot(self.r, self.q)
        ax9.set_ylabel('$q$')
        ax9.set_ylim([0,None])

        plt.tight_layout()
        plt.show()