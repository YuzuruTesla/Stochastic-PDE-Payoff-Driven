'''
Runs a Demo 1D model.
Corresponds to Figure 3.1 and 3.2 in Section 3.1
'''

from piegy import simulation, figures
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
from random import randrange
from helper_funcs import ave_interval_1D  # provided in the repository

# Configure Matplotlib plot settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = plt.rcParams['font.serif']
mpl.rcParams['font.size'] = 13
mpl.rcParams['savefig.dpi'] = 300


####  Model Initialization  ####

N = 1                   # Number of rows
M = 100                 # Number of cols
maxtime = 100           # how long you want the model to run
record_itv = 0.1        # how often to record data.
sim_time = 1            # repeat simulation to reduce randomness
boundary = True         # boundary condition.
seed = 36               # seed for random number generation
random.seed(seed)       # also seeds the random library

# initial population for the N x M patches. 
init_popu = [[[444 + randrange(-5, 6), 222 + randrange(-5, 6)] for _ in range(M)] for _ in range(N)]

# flattened payoff matrices, total resource is 0.4, cost of fighting is 0.1
matrices = [[[-1, 4, 0, 2] for _ in range(M)] for _ in range(N)]

# patch parameters
patch_params = [[[2, 0.02, 0, 0, 0.001, 0.001] for _ in range(M)] for _ in range(N)]

print_pct = 25           # print progress
check_overflow = False

# create and run the simulation object
mod = simulation.model(N, M, maxtime, record_itv, sim_time, boundary, init_popu, matrices, patch_params, 
                        print_pct = print_pct, seed = seed, check_overflow = check_overflow)
simulation.run(mod)
# You can save the model by uncommenting the line below
# save(mod, "demo_model")



####  Plotting  ####

# averaging over the time interval: 95% ~ 100% of max_time
start = int(0.95 * mod.max_record)
end = int(1.0 * mod.max_record)
xaxis = [i for i in range(mod.M)]
U_end = ave_interval_1D(mod.U, start, end)
V_end = ave_interval_1D(mod.V, start, end)
Hpi_end = ave_interval_1D(mod.Hpi, start, end)
Dpi_end = ave_interval_1D(mod.Dpi, start, end)

fig1, ax1 = plt.subplots(1, 2, figsize = (15, 5))
ax1[0].scatter(xaxis, U_end, color = 'darkorange')
ax1[0].scatter(xaxis, V_end, color = 'steelblue')
ax1[0].set_ylim([-10, 610])
ax1[0].set_xlabel('Patch Location', weight = 400, fontsize = 16)
ax1[0].set_ylabel('Population', weight = 400, fontsize = 16)
ax1[1].scatter(xaxis, Hpi_end, color = 'darkorchid')
ax1[1].scatter(xaxis, Dpi_end, color = 'seagreen')
ax1[1].set_ylim([-0.55, 2.05])
ax1[1].set_xlabel('Patch Location', weight = 400, fontsize = 16)
ax1[1].set_ylabel('Payoff', weight = 400, fontsize = 16)
fig1.savefig("distribution_fig3_2.png", bbox_inches = "tight")
plt.close(fig1)



####  Figure 3.1  ####
# We put Figure 3.1 in the below 3.2 as we need to manually average over the spatial domain

# manually average over the spatial domain (1 x 100)
mod.U /= 100
mod.V /= 100
mod.Hpi /= 100
mod.Dpi /= 100

fig2, ax2 = plt.subplots(figsize = (7.5, 5), dpi = 300)
figures.UV_dyna(mod, ax = ax2)
ax2.set_title("")
ax2.set_ylim(-50, 1050)
ax2.set_xlabel("Time", fontsize = 18)
ax2.set_ylabel("Average Population", fontsize = 18)
ax2.tick_params(axis="both", labelsize = 18) 
for line in ax2.lines:
    line.set_linewidth(4)
ax2.legend(loc = "upper left", fontsize = 16)
ax2.grid()
fig2.subplots_adjust(left=0.18, right=0.99, top=0.95, bottom=0.15)
fig2.savefig("popu_dynamics_fig3_1.png")
plt.close(fig2)

fig3, ax3 = plt.subplots(figsize = (7.5, 5), dpi = 300)
figures.pi_dyna(mod, ax = ax3)
ax3.set_title("")
ax3.set_ylim(0.35, 1.45)
ax3.set_xlabel("Time", fontsize = 18)
ax3.set_ylabel("Average Payoff", fontsize = 18)
ax3.tick_params(axis='both', labelsize = 18) 
ax3.lines[0].set_linewidth(4)
ax3.lines[1].set_linewidth(4)
ax3.lines[1].set_alpha(1)
ax3.lines[1].set_linestyle('--')
ax3.lines[2].set_linewidth(4)
ax3.legend(loc = 'center left', fontsize = 16)
ax3.grid()
fig3.subplots_adjust(left=0.12, right=0.99, top=0.99, bottom=0.15)
fig3.savefig("payoff_dynamics_fig3_1.png")
plt.close(fig3)

