from constants import *
import numpy as np
import matplotlib.pyplot as plt

# Mason Weiss
# PHYS477 - Bateman Equations for Xe-135 Concentration
# April 2024

maxhr = 100             # max time in hours
minhr = -10             # min time in hours
t_max = maxhr*60*60     # max time in seconds
timestep = 1            # timestep in seconds

# shorthand for equations
A = gI * SF * phi
B = lI
C = gX * SF * phi
E = lX + sX * phi

# secular equilibrium conditions
i0 = gI * SF * phi / lI
x0 = (gI + gX) * SF * phi / (lX + sX * phi)

# initial concentration values
i0plot = i0 / mili_avo * 1e6
x0plot = x0 / mili_avo * 1e6

# time plotting array in seconds and hours
t_plot = np.arange(0, t_max, timestep)
t_plot_hrs = t_plot / 3600
t_plot_hrs = np.insert(t_plot_hrs, 0, minhr)

################################
#  (1) turn-on to 100% power   #
################################
print("Power from 0% to 100%")
i_plot_1 = i0 * (1-np.exp(-B * t_plot))
x_plot_1 = x0 * (1-np.exp(-B * t_plot)) + C/E * (1-np.exp(-E * t_plot))
x_plot_1 = np.insert(x_plot_1, 0, 0)
i_plot_1 = np.insert(i_plot_1, 0, 0)

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, i_plot_1 / mili_avo * 1e6, 'r', label='Iodine-135')
plt.plot(t_plot_hrs, x_plot_1 / mili_avo * 1e6, 'c', label='Xenon-135')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [i0 / mili_avo * 1e6, i0 / mili_avo * 1e6], 'r--')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [x0 / mili_avo * 1e6, x0 / mili_avo * 1e6], 'c--')
plt.legend()
plt.xlabel(r"time since turn-on $(hours)$")
plt.ylabel(r"concentration $(mmol/m^3)$")
plt.title("Increasing Power from 0% to 100%")
plt.suptitle("Concentration of Fission Products")
plt.xlim([t_plot_hrs[0], t_plot_hrs[-1]])
plt.ylim([-0.1*x0plot, 1.1*i0plot])
plt.savefig("fig2.png")
plt.show()


################################
#  (2) shut down to 0% power   #
################################
print("Power from 100% to 0%")
i_plot_2 = i0 * np.exp(-lI * t_plot)
x_plot_2 = x0 * np.exp(-lX * t_plot) + i0 * lI / (lI - lX) * (np.exp(-lX * t_plot) - np.exp(-lI * t_plot))

tmax_hr= np.argmax(x_plot_2)/3600
xmax = max(x_plot_2)
print("    Maximum Xe-135 concentration of", round(xmax/mili_avo*1e6,5), "mmol/m^3 at", round(tmax_hr,5), "hours")

x_plot_2 = np.insert(x_plot_2, 0, x0)
i_plot_2 = np.insert(i_plot_2, 0, i0)

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, i_plot_2/mili_avo*1e6, 'r', label='Iodine-135')
plt.plot(t_plot_hrs, x_plot_2/mili_avo*1e6, 'c', label='Xenon-135')
plt.plot([tmax_hr, tmax_hr],[-0.1*x0plot, 1.1*i0plot], 'k--', label=r'time of max $^{135}$Xe')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [x0 / mili_avo * 1e6, x0 / mili_avo * 1e6], 'c--')
plt.legend()
plt.xlabel(r"time since shutdown $(hours)$")
plt.ylabel(r"concentration $(mmol/m^3)$")
plt.title("Decreasing Power from 100% to 0%")
plt.suptitle("Concentration of Fission Products")
plt.xlim([t_plot_hrs[0], t_plot_hrs[-1]])
plt.ylim([-0.1*x0plot, 1.1*i0plot])
plt.savefig("fig3.png")
plt.show()

################################
#  (3) reduction to 50% power  #
################################
print("Power from 100% to 50%")
# shorthand for equations
phi = phi * 0.50    # reduction in power is a reduction in neutron flux
A = gI * SF * phi
B = lI
C = gX * SF * phi
E = lX + sX * phi

i_plot_3 = np.zeros_like(t_plot)
x_plot_3 = np.zeros_like(t_plot)
i_plot_3[0] = i0
x_plot_3[0] = x0

# Forward Euler Method
for idx in range(int(t_max/timestep)-1):
    i_plot_3[idx + 1] = i_plot_3[idx] + timestep * (A - B * i_plot_3[idx])
    x_plot_3[idx + 1] = x_plot_3[idx] + timestep * (C + B * i_plot_3[idx] - E * x_plot_3[idx])

tmax_hr = np.argmax(x_plot_3)/3600
xmax = max(x_plot_3)
print("    Maximum Xe-135 concentration of", round(xmax/mili_avo*1e6,5), "mmol/m^3 at", round(tmax_hr,5), "hours")

x_plot_3 = np.insert(x_plot_3, 0, x0)
i_plot_3 = np.insert(i_plot_3, 0, i0)

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, i_plot_3 / mili_avo * 1e6, 'r', label='Iodine-135')
plt.plot(t_plot_hrs, x_plot_3 / mili_avo * 1e6, 'c', label='Xenon-135')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [x0 / mili_avo * 1e6, x0 / mili_avo * 1e6], 'c--')
plt.plot([tmax_hr, tmax_hr],[-0.1*x0plot, 1.1*i0plot], 'k--', label=r'time of max $^{135}$Xe')
plt.legend()
plt.xlabel(r"time since reduction to 50% power $(hours)$")
plt.ylabel(r"concentration $(mmol/m^3)$")
plt.title("Decreasing Power from 100% to 50%")
plt.suptitle("Concentration of Fission Products")
plt.xlim([t_plot_hrs[0], t_plot_hrs[-1]])
plt.ylim([-0.1*x0plot, 1.1*i0plot])
plt.savefig("fig4.png")
plt.show()


################################
#  (4) increase to 100% power  #
################################
print("Power from 50% to 100%")

phi = phi * 2    # increase power back to peak operating level
A = gI * SF * phi
B = lI
C = gX * SF * phi
E = lX + sX * phi

i_plot_4 = np.zeros_like(t_plot)
x_plot_4 = np.zeros_like(t_plot)
i_plot_4[0] = i_plot_3[-1]
x_plot_4[0] = x_plot_3[-1]

# Forward Euler Method
for idx in range(int(t_max/timestep)-1):
    i_plot_4[idx + 1] = i_plot_4[idx] + timestep * (A - B * i_plot_4[idx])
    x_plot_4[idx + 1] = x_plot_4[idx] + timestep * (C + B * i_plot_4[idx] - E * x_plot_4[idx])

tmin_hr = np.argmin(x_plot_4)/3600
xmin = min(x_plot_4)
print("    Minimum Xe-135 concentration of", round(xmin/mili_avo*1e6,5), "mmol/m^3 at", round(tmin_hr,5), "hours")

x_plot_4 = np.insert(x_plot_4, 0, x_plot_3[-1])
i_plot_4 = np.insert(i_plot_4, 0, i_plot_3[-1])

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, i_plot_4/mili_avo*1e6, 'r', label='Iodine-135')
plt.plot(t_plot_hrs, x_plot_4/mili_avo*1e6, 'c', label='Xenon-135')
plt.plot([tmin_hr, tmin_hr],[-0.1*x0plot, 1.1*i0plot], 'k--', label=r'time of min $^{135}$Xe')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [x0 / mili_avo * 1e6, x0 / mili_avo * 1e6], 'c--')
plt.legend()
plt.xlabel(r"time since power increase to 100% $(hours)$")
plt.ylabel(r"concentration $(mmol/m^3)$")
plt.title("Increasing Power from 50% to 100%")
plt.suptitle("Concentration of Fission Products")
plt.xlim([t_plot_hrs[0], t_plot_hrs[-1]])
plt.ylim([-0.1*x0plot, 1.1*i0plot])
plt.savefig("fig5.png")
plt.show()


################################
#    (5) reactivity plot for   #
#   decrease to 0% power (2)   #
################################
prf_factor = -sX / (nu * SF)
print(f"poison reactivity factor: {prf_factor:.5e}")

reactivity_5 = prf_factor * x_plot_2
reactivity_5 = reactivity_5 - reactivity_5[0]

# compute zeros of function to determine normal operation
prod = reactivity_5[1:]*reactivity_5[:-1]
prod = np.where(prod > 0, 0, -1)
prod = np.argwhere(prod)  # returns list of indices
zero_idx = prod[2]+1      # 0, 1, included in list so filter them out

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, reactivity_5, 'g', label='Change in Reactivity')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [0, 0], 'k--', label='Zero Change in Reactivity')
plt.plot([t_plot_hrs[zero_idx], t_plot_hrs[zero_idx]], [-0.024, 0.024], 'c--', label= f'Zero Change in Reactivity at {round(t_plot_hrs[zero_idx][0],2)} hrs')
plt.legend(loc='lower right')
plt.xlabel(r"time since shutdown $(hours)$")
plt.ylabel(r"change in reactivity from $t=0$")
plt.title("Decreasing Power from 100% to 0%")
plt.suptitle("Change in Reactivity")
plt.ylim([-0.024, 0.024])
plt.savefig("fig6.png")
plt.show()


################################
#    (6) reactivity plot for   #
#   decrease to 50% power (3)  #
################################
reactivity_6 = prf_factor * x_plot_3
reactivity_6 = reactivity_6 - reactivity_6[0]

# compute zeros of function to determine normal operation
prod = reactivity_6[1:]*reactivity_6[:-1]
prod = np.where(prod > 0, 0, -1)
prod = np.argwhere(prod)  # returns list of indices
zero_idx = prod[2]+1      # 0, 1, included in list so filter them out

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, reactivity_6, 'g', label='Change in Reactivity')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [0, 0], 'k--', label='Zero Change in Reactivity')
plt.plot([t_plot_hrs[zero_idx], t_plot_hrs[zero_idx]], [-0.024, 0.024], 'c--', label= f'Zero Change in Reactivity at {round(t_plot_hrs[zero_idx][0],2)} hrs')
plt.legend(loc='lower right')
plt.xlabel(r"time since power decrease to 50% $(hours)$")
plt.ylabel(r"change in reactivity from $t=0$")
plt.title("Decreasing Power from 100% to 50%")
plt.suptitle("Change in Reactivity")
plt.ylim([-0.024, 0.024])
plt.savefig("fig7.png")
plt.show()


################################
#   (7) reactivity plot for    #
#  increase to 100% power (4)  #
################################
reactivity_7 = prf_factor * x_plot_4
reactivity_7 = reactivity_7 - reactivity_7[0]

# compute zeros of function to determine normal operation
prod = reactivity_7[1:]*reactivity_7[:-1]
prod = np.where(prod > 0, 0, -1)
prod = np.argwhere(prod)  # returns list of indices
zero_idx = prod[2]+1      # 0, 1, included in list so filter them out

plt.figure()
plt.grid()
plt.plot(t_plot_hrs, reactivity_7, 'g', label='Change in Reactivity')
plt.plot([t_plot_hrs[0], t_plot_hrs[-1]], [0, 0], 'k--', label='Zero Change in Reactivity')
plt.plot([t_plot_hrs[zero_idx], t_plot_hrs[zero_idx]], [-0.024, 0.024], 'c--', label= f'Zero Change in Reactivity at {round(t_plot_hrs[zero_idx][0],2)} hrs')
plt.legend(loc='lower right')
plt.xlabel(r"time since power increase to 100% $(hours)$")
plt.ylabel(r"change in reactivity from $t=0$")
plt.title("Increasing Power from 50% to 100%")
plt.suptitle("Change in Reactivity")
plt.ylim([-0.024, 0.024])
plt.savefig("fig8.png")
plt.show()






