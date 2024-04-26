# Mason Weiss
# PHYS477 - Physics of Nuclear Energy

# Iodine-135 Constants
lI = 2.9306e-5   # decay constant (1/s)
gI = 0.0629     # fission yield of I-135 & Te-135 from U-235

# Xenon-135 Constants
lX = 2.1066e-5   # decay constant (1/s)
gX = 0.002576    # fission yield of Xe-135 from U-235
sX = 2.65e-18   # microscopic cross-section of absorption of Xe-135

# Fission Constants
sF = 585.1e-24      # microscopic cross-section of fission (cm^2)
nu = 2.43           # average number of neutrons per fission
Ef = 3.2e-11        # average fission output energy (200 MeV) in Joules

# Reactor Constants
N = 4.74e20         # density of U-235 atoms (per cm^3)
NV = 5.13e27        # number of U-235 atoms per core
P = 3000 * 1e6      # power of reactor (thermal) (W)
SF = sF * N         # macroscopic cross-section of fission (1/cm)
phi = P/Ef/sF/NV    # neutron flux

# Physical Constants
avogadro = 6.0221408e+23
mili_avo = avogadro/1000

