'''
Runs pure diffusion models and produce the spatial patterns.
Corresponds to Figure 3.3 in Section 3.2

Notive Figure 3.3 has three panels: the panels share the same codes except for mu_u value. 
We provide a switch for mu_u (named "mu1_val" below) and you can produce all three panels by changing its value to 1, 25, 50.
'''

from piegy import simulation
from piegy.data_tools import save
import matplotlib as mpl
import matplotlib.pyplot as plt
from random import randrange
from helper_funcs import ave_interval_1D  # provided in the repository

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = plt.rcParams['font.serif']
mpl.rcParams['font.size'] = 13

# set value for mu1 ($\mu_u$ in paper)
# we used 1, 25, 50
mu1_val = 1



# run the simulation
N = 1                   # Number of rows
M = 100                 # Number of cols
maxtime = 200           # how long you want the model to run
record_itv = 0.1        # how often to record data.
sim_time = 1000         # repeat simulation to reduce randomness
boundary = True         # boundary condition.
seed = 36               # seed for random number generation

# initial population for the N x M patches. 
init_popu = [[[444 + randrange(-5, 6), 222 + randrange(-5, 6)] for _ in range(M)] for _ in range(N)]
# flattened payoff matrices, total resource is 0.4, cost of fighting is 0.1
matrices = [[[-1, 4, 0, 2] for _ in range(M)] for _ in range(N)]
# patch parameters
patch_params = [[[mu1_val, 0.2, 0, 0, 0.001, 0.001] for _ in range(M)] for _ in range(N)]

print_pct = 500           # print progress
check_overflow = False

# create and run the simulation object
mod = simulation.model(N, M, maxtime, record_itv, sim_time, boundary, init_popu, matrices, patch_params, 
                        print_pct = print_pct, seed = seed, check_overflow = check_overflow)
simulation.run(mod)
# recommend saving the model as the simulation is expensive
save(mod, 'pure_diffusion/mu1=' + str(patch_params[0][0][0]) + ', mu2=' + str(patch_params[0][0][1]))

# data processing
start = int(0.95 * mod.max_record)
end = int(1.0 * mod.max_record)
xaxis = [i for i in range(mod.M)]
U_end = ave_interval_1D(mod.U, start, end)
V_end = ave_interval_1D(mod.V, start, end)
Upi_end = ave_interval_1D(mod.Hpi, start, end)
Vpi_end = ave_interval_1D(mod.Dpi, start, end)

# plotting
fig, ax = plt.subplots(1, 2, figsize = (15, 5), dpi = 300)
ax[0].scatter(xaxis, U_end, color = 'darkorange', label = 'U', s = 50)
ax[0].scatter(xaxis, V_end, color = 'steelblue', label = 'V', s = 50)
ax[0].set_ylim([-10, 610])
ax[0].tick_params(axis='x', labelsize=16)
ax[0].tick_params(axis='y', labelsize=16)
ax[0].set_xlabel('Patch Location', weight = 400, fontsize = 18)
ax[0].set_ylabel('Population', weight = 400, fontsize = 18)
ax[0].legend(loc = (0.03, 0.03), fontsize = 16)

ax[1].scatter(xaxis, Upi_end, color = 'darkorchid', label = r'$p_H$', s = 50)
ax[1].scatter(xaxis, Vpi_end, marker = 'o', color = 'seagreen', label = r'$p_D$', s = 30)
ax[1].set_ylim([-0.05, 1.45])
ax[1].tick_params(axis='x', labelsize=16)
ax[1].tick_params(axis='y', labelsize=16)
ax[1].set_xlabel('Patch Location', weight = 400, fontsize = 18)
ax[1].set_ylabel('Payoff', weight = 400, fontsize = 18)
ax[1].legend(loc = (0.03, 0.03), fontsize = 16)
fig.suptitle(r'$\mu_u$ = 25, $\mu_v$ = 0.2', fontsize = 20)
fig.savefig("pure_diffusion_" + str(mu1_val) + "_02.png", bbox_inches = "tight")
plt.close(fig)




