'''
Runs pure diffusion models and produce the correlation between w_v and population/payoff.
Corresponds to Figure 3.6 in Section 3.3

NOT recommend to run directly, see comments on w2_values and mu2_values below.
'''

from piegy import simulation, test_var
import matplotlib as mpl
import matplotlib.pyplot as plt
from random import randrange

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = plt.rcParams['font.serif']
mpl.rcParams['font.size'] = 13


####  Values of w_v  ####
# In practice, we strongly recommend against brutally throwing this huge list of values into the simulation
# as it could not take advantage of multi-processing and will take a very long time.
# We recommend making several copies of this script and split this list into smaller chunks and run seperately
# such as 
# w2_values = 0~5, 5~6, 6~7, ...
# Here w2 is w_v in the paper.
w2_values = [i for i in range(0, 5)] + [round(5 + i * 0.2, 1) for i in range(0, 15)] + [i for i in range(8, 61)]  # finer resolution at 5~8


####  Running the Model  ####

N = 1                   # Number of rows
M = 100                 # Number of cols
maxtime = 2000          # how long you want the model to run
                        # These models take particularly long time to converge.
record_itv = 1          # how often to record data.
sim_time = 100          # repeat simulation to reduce randomness
boundary = True         # boundary condition.

# initial population for the N x M patches. 
init_popu = [[[444 + randrange(-5, 6), 222 + randrange(-5, 6)] for _ in range(M)] for _ in range(N)]
# flattened payoff matrices, total resource is 0.4, cost of fighting is 0.1
matrices = [[[-1, 4, 0, 2] for _ in range(M)] for _ in range(N)]
# patch parameters
patch_params = [[[1, 1, 0.1, None, 0.001, 0.001] for _ in range(M)] for _ in range(N)]

print_pct = 500           # print progress
seed = 36               # seed for random number generation
check_overflow = False

# create a simulation object
mod = simulation.model(N, M, maxtime, record_itv, sim_time, boundary, init_popu, matrices, patch_params, 
                        print_pct = print_pct, seed = seed, check_overflow = check_overflow)


####  Simulation on values of w2  ####
var = 'w2'
values = w2_values
dirs = 'payoff_driven_data'
var_dirs = test_var.test_var1(mod, var, values, dirs, compress_ratio = 5)


####  Plotting  ####
var_dirs = test_var.get_dirs1(var, values, dirs)
fig, ax = plt.subplots(1, 2, figsize = (15, 5), dpi = 300)
test_var.var_UV1(var, values, var_dirs, ax[0], ax[0], color_H = 'dodgerblue', color_D = 'tomato')
test_var.var_pi1(var, values, var_dirs, ax[1], ax[1], color_H = 'dodgerblue', color_D = 'tomato')

ax[0].set_title(None)
for line in ax[0].get_lines():
    line.set_linestyle('None')
    line.set_marker('o')
ax[0].get_lines()[0].set_markersize(7)
handles = ax[0].get_lines()
ax[0].set_ylim([85, 515])
ax[0].legend(handles, ['U', 'V'], loc = 'upper right', borderpad = 0.6)
ax[0].set_xlabel(r'Payoff Sensitivity of Doves $w_v$', weight = 400, fontsize = 16)
ax[0].set_ylabel('Average Population', weight = 400, fontsize = 16)
ax[0].grid()


ax[1].set_title(None)
ax[1].set_title(None)
for line in ax[1].get_lines():
    line.set_linestyle('None')
    line.set_marker('o')
ax[1].get_lines()[0].set_markersize(7)
handles = ax[1].get_lines()
ax[1].set_ylim([-0.23, 0.83])
ax[1].legend(handles, [r'$p_H$', r'$p_D$'], loc = 'upper right', borderpad = 0.6)
ax[1].set_xlabel(r'Payoff Sensitivity of Doves $w_v$', weight = 400, fontsize = 16)
ax[1].set_ylabel('Average Payoff', weight = 400, fontsize = 16)
ax[1].grid()
fig_name = "payoff_driven_w2_01.png"
fig.savefig(fig_name)
print("Fig saved: " + fig_name)

