import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# Read CSV files
data_LF95_df = pd.read_csv('Data/large_tank_data_LF95_3weeks.csv')
data_LF50_df = pd.read_csv('Data/large_tank_data_LF50_3weeks.csv')
data_LF50_wo_wall_df = pd.read_csv('Data/large_tank_data_LF50_3weeks_wo_wall.csv')
data_LF95_wo_wall_df = pd.read_csv('Data/large_tank_data_LF95_3weeks_wo_wall.csv')

# Convert DataFrames to dictionaries to maintain the same syntax
data_LF95 = data_LF95_df.to_dict('series')
data_LF50 = data_LF50_df.to_dict('series')
data_LF50_wo_wall = data_LF50_wo_wall_df.to_dict('series')
data_LF95_wo_wall = data_LF95_wo_wall_df.to_dict('series')

# First subplot - Total, Liquid and Bottom heat transfer rates (normalized)
cmap = plt.get_cmap('inferno', 6)

# Mask to select data starting from a specific time
time_tr = 70
time_to = time_tr + 200
mask1 = data_LF95['Time'] <= time_tr*3600
mask2 = (data_LF95['Time'] >= time_tr*3600) & (data_LF95['Time'] < (time_to)*3600)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=200)

# First subplot
ax1.plot(data_LF95['Time'][mask2]/3600, data_LF95['Q_tot'][mask2]/1000, label=r'$\dot{Q}_{\text{tot}}$', color=cmap(0))
ax1.plot(data_LF95['Time'][mask2]/3600, data_LF95['Q_L'][mask2]/1000, label=r'$\dot{Q}_{\text{L}}$', color=cmap(4))
ax1.plot(data_LF95['Time'][mask2]/3600, data_LF95['Q_b'][mask2]/1000, label=r'$\dot{Q}_{\text{b}}$', color=cmap(3))
ax1.plot(data_LF95['Time'][mask2]/3600, data_LF95['Q_Vw'][mask2]/1000, label=r'$\dot{Q}_{\text{Wi}}$', color=cmap(2))

ax1.set_xlim(time_tr, 240)
ax1.set_ylabel(r'Heat transfer rate / kW', fontsize=fontsize_label)
ax1.set_ylim(0, 50)
ax1.set_xlabel('Time / h', fontsize=fontsize_label)
ax1.tick_params(labelsize=fontsize_ticks)
ax1.legend(fontsize=fontsize_legend, loc='upper right', ncol=1, framealpha=1)

# Second subplot - BOG Temperature
cmap = plt.get_cmap('inferno', 4)

time_tr = 168
init_stationary = 70
end_stationary = init_stationary + 200 

mask_transient_LF50  = data_LF50['Time'] < time_tr *3600
mask_stationary_LF50 = (data_LF50['Time'] >= init_stationary *3600) & (data_LF50['Time'] < (end_stationary) * 3600)
mask_transient_LF95  = data_LF95['Time'] < time_tr *3600
mask_stationary_LF95 = (data_LF95['Time'] >= init_stationary * 3600) & (data_LF95['Time'] < (end_stationary) * 3600)
mask_transient_LF50_wo_wall  = data_LF50_wo_wall['Time'] < time_tr *3600
mask_stationary_LF50_wo_wall = (data_LF50_wo_wall['Time'] >= init_stationary *3600) & (data_LF50_wo_wall['Time'] < (end_stationary) * 3600)
mask_transient_LF95_wo_wall  = data_LF95_wo_wall['Time'] < time_tr *3600
mask_stationary_LF95_wo_wall = (data_LF95_wo_wall['Time'] >= init_stationary *3600) & (data_LF95_wo_wall['Time'] < (end_stationary) * 3600)

# Second subplot - stationary period
ax2.plot(data_LF95['Time'][mask_stationary_LF95]/3600, data_LF95['BOG'][mask_stationary_LF95]*3600, label='Wall model', color=cmap(2))

ax2.plot(data_LF95_wo_wall['Time'][mask_stationary_LF95_wo_wall]/3600, data_LF95_wo_wall['BOG'][mask_stationary_LF95_wo_wall]*3600, 
         label='Simplified', color=cmap(1))

ax2.set_xlabel('Time / h', fontsize=14)
ax2.set_xlim(init_stationary, 240)
ax2.set_ylabel(r'Boil-off gas rate / $kg\ h^{-1}$', fontsize=fontsize_label)
ax2.set_ylim(70, 88)
ax2.tick_params(labelsize=fontsize_ticks)
ax2.legend(loc='upper right', fontsize=fontsize_legend, ncol=2, columnspacing=1.0)
plt.tight_layout()

plt.savefig('Figures/Fig_5.svg', dpi=300, bbox_inches='tight')