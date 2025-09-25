'''
Runs models with payoff-driven migration and produce the spatial distributions.
Corresponds to Figure 3.5 in Section 3.3

Notice Figure 3.5 has four panels: the panels share the same codes except for w_v value. 
We provide a switch for w_v (a variable named "w2_val" below) 
and you can produce each panel by changing its value to 0, 6, 20, 60.
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


####  Values of w2 ####
# This corresponds to w_v in paper
# We used w_v = 0, 6, 20, 60 when generating the spatial patterns in Figure 3.5
# Change w2_val to the other three values to reproduce each of the panels in Fig 3.5
w2_val = 0


####  Running the Model  ####

N = 1                   # Number of rows
M = 100                 # Number of cols
maxtime = 500           # how long you want the model to run
record_itv = 0.1        # how often to record data.
sim_time = 1000         # repeat simulation to reduce randomness
                        # We repeat 1000 times here to reduce stochasticity to the utmost.
boundary = True         # boundary condition.

# initial population for the N x M patches. 
init_popu = [[[444 + randrange(-5, 6), 222 + randrange(-5, 6)] for _ in range(M)] for _ in range(N)]

# flattened payoff matrices, total resource is 0.4, cost of fighting is 0.1
matrices = [[[-1, 4, 0, 2] for _ in range(M)] for _ in range(N)]

# patch variables
patch_params = [[[1, 1, 0.1, 40, 0.001, 0.001] for _ in range(M)] for _ in range(N)]

print_pct = 100           # print progress
seed = 1000               # seed for random number generation
check_overflow = False

# create a simulation object
mod = simulation.model(N, M, maxtime, record_itv, sim_time, boundary, init_popu, matrices, patch_params, 
                        print_pct = print_pct, seed = seed, check_overflow = check_overflow)
simulation.run(mod)
# recommend saving the model as the simulation is expensive
save(mod, "payoff_driven/w2=" + str(w2_val))


####  Plotting  ####

# Calculate average over the last 5% of time
start = int(0.95 * mod.max_record)
end = int(1.0 * mod.max_record)

xaxis = [i for i in range(mod.M)]
U_end = ave_interval_1D(mod.U, start, end)
V_end = ave_interval_1D(mod.V, start, end)
Upi_end = ave_interval_1D(mod.Hpi, start, end)
Vpi_end = ave_interval_1D(mod.Dpi, start, end)

# Plot
fig, ax = plt.subplots(1, 2, figsize = (15, 5), dpi = 300)
ax[0].scatter(xaxis, U_end, color = 'darkorange', label = 'U', s = 50)
ax[0].scatter(xaxis, V_end, color = 'steelblue', label = 'V', s = 50)
ax[0].set_ylim([-10, 510])
ax[0].tick_params(axis='x', labelsize=16)
ax[0].tick_params(axis='y', labelsize=16)
ax[0].set_xlabel('Patch Location', weight = 400, fontsize = 18)
ax[0].set_ylabel('Population', weight = 400, fontsize = 18)
ax[0].legend(loc = (0.02, 0.02), fontsize = 16)

ax[1].scatter(xaxis, Upi_end, color = 'darkorchid', label = r'$p_H$', s = 50)
# to set smaller dot size, decrease the value of s in the line below. 
# e.g. s = 20
ax[1].scatter(xaxis, Vpi_end, color = 'seagreen', label = r'$p_D$', s = 50)
ax[1].set_ylim([-0.65, 1.05])
ax[1].tick_params(axis='x', labelsize=16)
ax[1].tick_params(axis='y', labelsize=15)
ax[1].set_xlabel('Patch Location', weight = 400, fontsize = 18)
ax[1].set_ylabel('Payoff', weight = 400, fontsize = 18)
ax[1].legend(loc = (0.02, 0.02), fontsize = 16)
fig.suptitle(r"$w_u$ = 0.1, $w_v$ = " + str(w2_val), fontsize = 21)

fig_name = "payoff_driven_pattern_w2_" + str(w2_val) + ".png"
fig.savefig(fig_name)
print("Fig saved: " + fig_name)
