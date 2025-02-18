import numpy as np
import matplotlib.pyplot as plt
import pdb

b = 0.7 #m
a = 0.4 #m

BT0 = 12.2 # T
Btp0 = 1.0 # T
Btp1 = 0.0 # T
deltaBtp = Btp1 - Btp0 # T
BT1 = BT0 - (a**2 / b**2) * deltaBtp # T
Phi0 = np.pi * b**2 * BT0 + np.pi * a**2 * Btp0
Phi1 = np.pi * b**2 * BT1 + np.pi * a**2 * Btp1

W0 = (BT0**2 * np.pi * (b**2 - a**2) +
      (BT0 + Btp0)**2 * np.pi * a**2)

W1 = (BT1**2 * np.pi * (b**2 - a**2) +
      (BT1 + Btp1)**2 * np.pi * a**2)

print('Phi0', Phi0)
print('Phi1', Phi1)

print('W0', W0)
print('W1', W1)
print('W1/W0', W1/W0*1100)