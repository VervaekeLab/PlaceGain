clear all
close all

save_dir = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Nora/Paper/Figure1';


fig_1.VIP_GCAMP.load_VIP_GCAMP() % find session data & directory in this file

fig_1.VIP_GCAMP.chance_speed_correlation()

fig_1.VIP_GCAMP.speed_correlation()

fig_1.fig_1_f_g_h_plot_speed_correlation()

fig_1.VIP_GCAMP.start_stop_modulation()

%fig_1.plot_reward_modulation()

