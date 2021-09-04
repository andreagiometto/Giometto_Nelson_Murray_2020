% Written by Andrea Giometto, MIT license
% This function plots the density profiles according to different 
% models parametrized using well-mixed competition experiments between
% killer and sensitive strains. This plot corresponds to Figure 7
% - figure supplement 1 of the manuscript

subplot(7,1,1)
load('data_models/model_output_growth_logistic2D_ideal_inoculum_diffusion_Dfactor_1_Dt_0.003_K_82800000_a1_2.36e-09_a2_9.6e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])
subplot(7,1,2)
load('data_models/model_output_growth_logistic_ideal_inoculum_diffusion_Dfactor_1_Dt_0.003_K_82800000_a1_2.36e-09_a2_9.6e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dz_0.005_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])
subplot(7,1,3)
load('data_models/model_output_growth_logistic_Dgrowth_ideal_inoculum_diffusion_Dfactor_1_Dt_0.003_K_82800000_a1_2.36e-09_a2_9.6e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dz_0.005_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])
subplot(7,1,4)
load('data_models/model_output_growth_logistic_propto_grate_ideal_inoculum_diffusion_Dfactor_1_Dt_0.003_K_82800000_a1_1.57e-09_a2_4.75e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dz_0.005_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])
subplot(7,1,5)
load('data_models/model_output_growth_logistic_propto_grate_Dgrowth_ideal_inoculum_diffusion_Dy_0.003_Dfactor_1_Dt_0.003_K_82800000_a1_1.57e-09_a2_4.75e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dz_0.005_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])
subplot(7,1,6)
load('data_models/model_output_growth_constD_ideal_inoculum_diffusion_cellD_0.001_Dfactor_1_Dg_0.024_Dt_0.003_KS_2e-05_a1_1.76e-09_a2_7.57e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dz_0.005_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])
subplot(7,1,7)
load('data_models/model_output_growth_ideal_inoculum_diffusion_cellD_0.001_Dfactor_1_Dg_0.024_Dt_0.003_KS_2e-05_a1_1.76e-09_a2_7.57e-10_transfers_6_r0_0.5_dilutionF_10000_R_1.5_dr_0.0025_dz_0.005_dt_0.0001_n0IC-n0IM_50000.mat')
plot(r,n1/100,'-m'); hold on; plot(r,n2/100,'-c'); plot(r,n1/100+n2/100,'--k'); hold off
xlabel('R(cm)')
ylabel('Cell density (cells/mm^2)')
xlim([0,1])