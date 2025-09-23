'''
Runs pure diffusion models and produce the correlation between mu_u and population/payoff.
Corresponds to Figure 3.4 in Section 3.2

NOT recommend to run directly, see comments on mu1_values and mu2_values below.
'''

from piegy import simulation, test_var
import matplotlib as mpl
import matplotlib.pyplot as plt
from random import randrange

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = plt.rcParams['font.serif']
mpl.rcParams['font.size'] = 13


####  Values of mu_u and mu_v  ####
# In practice, we strongly recommend against brutally feeding this huge list of values into the simulation
# as it could not take advantage of the multi-processing and will take a very long time.
# We recommend making several copies of this script and split this list into smaller chunks and run seperately
# such as 
# mu1_values = [0, 15], [16, 30], ...
# mu2_values = [0.2], [2]

mu2_values = [0.2, 2]
mu1_values = [i for i in range(0, 76)]  # for mu2 = 0.2
mu1_values = [i for i in range(0, 301, 5)]  # for mu2 = 2


N = 1                   # Number of rows
M = 100                 # Number of cols
maxtime = 200           # how long you want the model to run
record_itv = 0.1        # how often to record data.
sim_time = 100            # repeat simulation to reduce randomness
boundary = True         # boundary condition.

# initial population for the N x M patches. 
init_popu = [[[444 + randrange(-5, 6), 222 + randrange(-5, 6)] for _ in range(M)] for _ in range(N)]
# flattened payoff matrices, total resource is 0.4, cost of fighting is 0.1
matrices = [[[-1, 4, 0, 2] for _ in range(M)] for _ in range(N)]
# patch parameters
patch_params = [[[None, None, 0, 0, 0.001, 0.001] for _ in range(M)] for _ in range(N)]

print_pct = 500           # print progress
seed = 36               # seed for random number generation
check_overflow = False

# create a simulation object
mod = simulation.model(N, M, maxtime, record_itv, sim_time, boundary, init_popu, matrices, patch_params, 
                        print_pct = print_pct, seed = seed, check_overflow = check_overflow)


####  w1, w2 values to test  ####
# w1 is w_U, w2 is w_V

var1 = 'mu2'
var2 = 'mu1'
values1 = mu2_values
values2 = mu1_values
dirs = 'pure_diffusion_data'
var_dirs2 = test_var.test_var2(mod, var1, var2, values1, values2, dirs, compress_ratio = 5)



####  Plotting mu2 = 0.2  ####

data_path = dirs

# read values and directories
var1 = 'mu2'
var2 = 'mu1'
val = 0.2
values1 = [val]
values2 = [round(1.0 * i, 3) for i in range(0, 76)]
var_dirs = test_var.get_dirs1(var2, values2, data_path)
for i in range(len(var_dirs)):
    var_dirs[i] = var_dirs[i][:len(data_path)] + '/mu2=' + str(val) + ', ' + var_dirs[i][len(data_path) + 1:]

# plot  
fig, ax = plt.subplots(1, 2, figsize = (15, 5), dpi = 300)
test_var.var_UV1(var2, values2, var_dirs, ax[0], ax[0], color_H = 'dodgerblue', color_D = 'tomato')
test_var.var_pi1(var2, values2, var_dirs, ax[1], ax[1], color_H = 'dodgerblue', color_D = 'tomato')

ax[0].set_title(None)
for line in ax[0].get_lines():
    line.set_linestyle('None')
    line.set_marker('o')
ax[0].get_lines()[0].set_markersize(7)
handles = ax[0].get_lines()
ax[0].set_ylim([190, 510])
ax[0].tick_params(axis='x', labelsize=14)
ax[0].tick_params(axis='y', labelsize=14)
ax[0].legend(handles, ['U', 'V'], loc = (0.02, 0.55), borderpad = 0.5)
ax[0].set_xlabel(r'Diffusivity of Hawks $\mu_u$', weight = 400, fontsize = 16)
ax[0].set_ylabel('Average Population', weight = 400, fontsize = 16)
ax[0].grid()

ax[1].set_title(None)
ax[1].set_title(None)
for line in ax[1].get_lines():
    line.set_linestyle('None')
    line.set_marker('o')
ax[1].get_lines()[0].set_markersize(7)
handles = ax[1].get_lines()
ax[1].set_ylim([0.655, 0.785])
ax[1].tick_params(axis='x', labelsize=13)
ax[1].tick_params(axis='y', labelsize=13)
ax[1].legend(handles, [r'$p_H$', r'$p_D$'], loc = (0.02, 0.55), borderpad = 0.5)
ax[1].set_xlabel(r'Diffusivity of Hawks $\mu_u$', weight = 400, fontsize = 16)
ax[1].set_ylabel('Average Payoff', weight = 400, fontsize = 16)
ax[1].grid()

fig.savefig('/pure_diffusion_,' + var1 + '=' + str(val) + '_Fig3_4.png')
plt.close(fig)



####  Plotting mu2 = 2  ####

data_path = dirs

# read values and directories
var1 = 'mu2'
var2 = 'mu1'
val = 2
values1 = [val]
values2 = [i for i in range(0, 301, 5)]
var_dirs = test_var.get_dirs1(var2, values2, data_path)
for i in range(len(var_dirs)):
    var_dirs[i] = var_dirs[i][:len(data_path)] + '/mu2=' + str(val) + ', ' + var_dirs[i][len(data_path) + 1:]
    
# plot
fig, ax = plt.subplots(1, 2, figsize = (15, 5), dpi = 300)
test_var.var_UV1(var2, values2, var_dirs, ax[0], ax[0], color_H = 'dodgerblue', color_D = 'tomato')
test_var.var_pi1(var2, values2, var_dirs, ax[1], ax[1], color_H = 'dodgerblue', color_D = 'tomato')

ax[0].set_title(None)
for line in ax[0].get_lines():
    line.set_linestyle('None')
    line.set_marker('o')
handles = ax[0].get_lines()
ax[0].set_ylim([190, 510])
ax[0].legend(handles, ['U', 'V'], loc = (0.02, 0.55), borderpad = 0.6)
ax[0].set_xlabel(r'Diffusivity of Hawks $\mu_u$', weight = 400, fontsize = 16)
ax[0].set_ylabel('Average Population', weight = 400, fontsize = 16)
ax[0].grid()

ax[1].set_title(None)
ax[1].set_title(None)
for line in ax[1].get_lines():
    line.set_linestyle('None')
    line.set_marker('o')
handles = ax[1].get_lines()
ax[1].set_ylim([0.655, 0.785])
ax[1].legend(handles, [r'$p_H$', r'$p_D$'], loc = (0.02, 0.55), borderpad = 0.6)
ax[1].set_xlabel(r'Diffusivit of Hawks $\mu_u$', weight = 400, fontsize = 16)
ax[1].set_ylabel('Average Payoff', weight = 400, fontsize = 16)
ax[1].grid()

fig.savefig('/pure_diffusion_,' + var1 + '=' + str(val) + '_Fig3_4.png')
plt.close(fig)

