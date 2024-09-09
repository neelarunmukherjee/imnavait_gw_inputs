clc; clear;
% %%
% rng('default');
%%
filename = "..\input_data\thickness_riparian_master.xlsx";
D = readtable(filename);
%%  dti/thaw depth
dti_master = table2array(rmmissing(D(:,8)));dti_master(dti_master<0)=[];dti_master(dti_master>98)=[];
dti_mike_nse = rmmissing(D.thaw_nse); 
dti_mike_nse(dti_mike_nse<=0) =[];
dti_mike_nsw =rmmissing(D.thaw_nsw);
dti_mike_nsw(dti_mike_nsw<=0) = [];
dti_mike_br = rmmissing(D.thaw_br);
dti_mike_br(dti_mike_br<=0) = [];

dti = [dti_master;dti_mike_nsw;dti_mike_br;dti_mike_nse]; dti = dti./100; %m
parmhat_dti = lognfit(dti);
figure;
h_thaw = histogram(dti,7,'Normalization','pdf','FaceColor',"#F08080");
hold on;
x_dti = linspace(0,max(dti),100);
dti_param = [-0.67, 0.21];
[dtiM,dtiV] = lognstat(dti_param(1), dti_param(2));
plot(x_dti,lognpdf(x_dti,dti_param(1), dti_param(2)),'r','LineWidth',1);
title('depth to ice (m)'); xlim([0 max(dti)]);
% exportgraphics(gca,'figures_new/dti.pdf','ContentType','vector')
% dti_MC = lognrnd(dti_param(1), dti_param(2),10000000,1);
% figure;
% histogram (dti_MC,100,'FaceColor',"#F08080"); 
% xlim([0 max(dti)]);
% title('depth to ice (m)');

%% dtw
dtw_master = (rmmissing(D.EstimatedDepthToWater));dtw_master(dtw_master<=0)=[];dtw_master(dtw_master>98)=[];
dtw_mike_nse = rmmissing(D.dtw_nse); dtw_mike_nse(dtw_mike_nse<=0) =[];
dtw_mike_nsw = rmmissing(D.dtw_nsw); dtw_mike_nsw(dtw_mike_nsw<=0) =[];
dtw_mike_br = rmmissing(D.dtw_br); dtw_mike_br(dtw_mike_br<=0) = [];
dtw = [dtw_master;dtw_mike_br;dtw_mike_nsw;dtw_mike_nse]; 
dtw = dtw./100;
figure;
h_water = histogram(dtw,13,'Normalization','pdf','FaceColor',"#20b2aa");
hold on; 
x_dtw = linspace(0,max(dtw),100);
dtw_param = [-1.7,1];
[dtwM, dtwV] = lognstat(dtw_param(1), dtw_param(2));
plot(x_dtw, lognpdf(x_dtw, dtw_param(1),dtw_param(2)),'r');
title('depth to water (m)'); xlim([0 max(dtw)]);
% exportgraphics(gca,'figures_new/dtw_master.pdf','ContentType','vector')
% dtw_MC = lognrnd(dtw_param(1), dtw_param(2),10000000,1);
% figure;
% histogram (dtw_MC,1000,'FaceColor',"#20b2aa"); 
% xlim([0 max(dtw)]);
% title('depth to water (m)');
% exportgraphics(gca,'figures_new/dtw_MC.pdf','ContentType','vector')
%% dt Acrotelm 
z_AC = rmmissing(D.dtA); z_AC(z_AC<=0) = []; z_AC = z_AC./100;
parmhat_z_AC = lognfit(z_AC);
figure;
histogram(z_AC,6, 'Normalization','pdf','FaceColor',"#77AC30");
hold on;
z_a = linspace(0,max(z_AC),100);
zAC_param =[-2.15,0.45];% [-2.1,0.4];
[zAC_M, zAC_V]=lognstat(zAC_param(1),zAC_param(2));
plot(z_a,lognpdf(z_a,zAC_param(1), zAC_param(2)),'r','LineWidth',1);
xlim([0 max(z_AC)]);title('Acrotelm thickness (m)');
% exportgraphics(gca,'figures_new/zAC.pdf','ContentType','vector')
% z_AC_MC = lognrnd(zAC_param(1),zAC_param(2),10000000,1); 
% figure;
% histogram (z_AC_MC,100,'FaceColor',"#77AC30"); 
% xlim([0 max(z_AC)]);
% title('Acrotelm thickness (m)');
% exportgraphics(gca,'figures_new/zAC_MC.pdf','ContentType','vector')

%% dt Catotelm
dtC = rmmissing(D.dtC); dtC(dtC<=0) = [];dtC(dtC>98)=[]; dtC = dtC./100;
parmhat_dtC = lognfit(dtC);
figure; 
histogram(dtC,8,'Normalization','pdf','FaceColor',"#ffc20d");
hold on;
z_a = linspace(0,max(dtC),100);
dtC_param = [-1.2,0.26];
[dtC_M, dtC_V] = lognstat(dtC_param(1), dtC_param(2));
plot(z_a,lognpdf(z_a,dtC_param(1), dtC_param(2)),'r','LineWidth',1);
xlim([0 max(dtC)]);title('Depth to Catotelm-MinSoil Bnd  (m)');
% exportgraphics(gca,'figures_new/dtC.pdf','ContentType','vector')
% dtC_MC = lognrnd(dtC_param(1),dtC_param(2),10000000,1); 
% figure;
% histogram (dtC_MC,100,'FaceColor',"#ffc20d"); 
% xlim([0 max(dtC)]);
% title('Depth to Catotelm-MinSoil Bnd (m)');
% exportgraphics(gca,'figures_new/dtC_MC.pdf','ContentType','vector')
%% 
filename = "..\input_data\ksat.xlsx";
K = readtable(filename);
K_AC_raw = rmmissing(table2array(K(:,3)));K_AC_raw = K_AC_raw(K_AC_raw>0);K_AC_raw(K_AC_raw<1e-4) = [];
K_CT_raw = rmmissing(table2array(K(:,6)));K_CT_raw = K_CT_raw(K_CT_raw>0);
K_MN_raw = rmmissing(table2array(K(:,9)));K_MN_raw = K_MN_raw(K_MN_raw>0);
%% acrotelm hydraulic conductivity 
pd = fitdist(K_AC_raw,'Weibull');
figure;
h = histogram(K_AC_raw,5,'Normalization','pdf','FaceColor',"#77AC30");
hold on; 
k_a = linspace(0,max(K_AC_raw),50); 
K_AC_param =[0.00275,3.2];
[K_AC_M, K_AC_V] = wblstat(K_AC_param(1), K_AC_param(2));
plot(k_a,wblpdf(k_a,K_AC_param(1), K_AC_param(2)),'r','LineWidth',1);
title('Acrotelm hydraulic conductivity (m/s): Weibull');
xlim([0 max(K_AC_raw)]);
% exportgraphics(gca,'figures_new/K_AC.pdf','ContentType','vector');
% rng('default');
% K_AC_MC_WBL = wblrnd(K_AC_param(1),K_AC_param(2),10000000,1); 
% figure;
% histogram (K_AC_MC_WBL,100,'FaceColor',"#77AC30");
% xlim([0 max(K_AC_raw)]); 
% title('Acrotelm hydraulic conductivity (m/s)');
% %exportgraphics(gca,'figures_new/K_AC_MC.pdf','ContentType','vector');
% K_AC_MC = K_AC_MC_WBL;
%% catotelm hydraulic conductivity
K_CT_raw(K_CT_raw>3e-4)=[];
parmhat_K_CT = lognfit(K_CT_raw);
figure;
h = histogram (K_CT_raw,8,'Normalization','pdf','FaceColor',"#EDB120");
hold on;
k_a = linspace(0,max(K_CT_raw),500);
K_CT_param = [-10.15,1.25];
[K_CT_M,K_CT_V] = lognstat(K_CT_param(1), K_CT_param(2));
plot(k_a,lognpdf(k_a,K_CT_param(1),K_CT_param(2)),'r','LineWidth',1);
xlim([0 max(K_CT_raw)]);title('Catotelm hydraulic conductivity (m/s)');
%exportgraphics(gca,'figures_new/K_CT.pdf','ContentType','vector');
% rng('default');
% K_CT_MC = lognrnd(K_CT_param(1),K_CT_param(2),10000000,1);
% figure;
% histogram (K_CT_MC,10000,'FaceColor',"#EDB120");
% set(gca,'FontSize',10);   set(gca,'LineWidth',1);xlim([0 max(K_CT_raw)]); 
% title('Catotelm hydraulic conductivity (m/s)');
%exportgraphics(gca,'figures_new/K_CT_MC.pdf','ContentType','vector');
%% Mineral Soil hydraulic conductivity
parmhat_K_MN = lognfit(K_MN_raw);
figure;
histogram(K_MN_raw,63,'Normalization','pdf','FaceColor',"#D95319");
title('Mineral Soil hydraulic conductivity (m/s)'); 
hold on;
k_m = linspace(0,max(K_MN_raw),10000000);
K_MN_param = [-14.55, 0.89];
[K_MN_M,K_MN_V] = lognstat(K_MN_param(1),K_MN_param(2));
plot(k_m,lognpdf(k_m,K_MN_param(1),K_MN_param(2)),'r','LineWidth',1);
xlim([0 0.5e-5]);
% %exportgraphics(gca,'figures_new/K_MN.pdf','ContentType','vector');
% rng('default'); ?????
% K_MN_MC = lognrnd(K_MN_param(1),K_MN_param(2),10000000,1);
% figure;
% histogram (K_MN_MC,1000,'FaceColor',"#D95319");
% set(gca,'FontSize',10);   set(gca,'LineWidth',1);xlim([0 0.5e-5]);
% title('Mineral Soil hydraulic conductivity (m/s)');
% %exportgraphics(gca,'figures_new/K_MN_MC.pdf','ContentType','vector');

%% initialize
KEFF = zeros(10000000,1);
% T = zeros(10000000,1);
% B = zeros(10000000,1);
%%
T1 = zeros(10000000,1);
T2 = zeros(10000000,1);
T3 = zeros(10000000,1);
T4 = zeros(10000000,1);
T5 = zeros(10000000,1);
T6 = zeros(10000000,1);
KEFF1 = zeros(10000000,1);
KEFF2 = zeros(10000000,1);
KEFF3 = zeros(10000000,1);
KEFF4 = zeros(10000000,1);
KEFF5 = zeros(10000000,1);
KEFF6 = zeros(10000000,1);
% DTI = -ones(10000000,1) ;
% DTW = -ones(10000000,1) ;
% DTC = -ones(10000000,1) ;
% zAC = -ones(10000000,1) ;
% K_AC = -ones(10000000,1) ;
% K_CT = -ones(10000000,1) ;
% K_MN = -ones(10000000,1) ;
%%
D = zeros(10000000, 7);
D(:,1)=lognrnd(zAC_param(1),zAC_param(2), 10000000,1);
D(:,2)=lognrnd(dtC_param(1),dtC_param(2), 10000000,1);
D(:,3)=lognrnd(dti_param( 1), dti_param(2), 10000000,1);
D(:,4)=lognrnd(dtw_param(1), dtw_param(2), 10000000,1); 
D(:,5)=wblrnd(K_AC_param(1),K_AC_param(2), 10000000,1);
D(:,6)=lognrnd(K_CT_param(1),K_CT_param(2), 10000000,1);
D(:,7)=lognrnd(K_MN_param(1),K_MN_param(2), 10000000,1);

%% Remove unrelasitic cases (water does not exist)
mask = ((D(:,3)-D(:,4))<0); 
D(mask,:)=[];
T = ones(length(D),1)-1; %initialize
%% Cases when ice table in mineral soil
% caseA1: water table in mineral soil
indxA1 =(D(:,3) >= D(:,2)) & (D(:,4) >= D(:,2));
T(indxA1) = D(indxA1,7) .* ((D(indxA1,3))-D(indxA1,4));

% caseA2: water table in catetelm
indxA2 = (D(:,3) >= D(:,2)) & (D(:,4) >= D(:,1)) & (D(:,4) < D(:,2));
T(indxA2) = D(indxA2, 7).*(D(indxA2, 3)-D(indxA2, 2)) + D(indxA2, 6) .* (D(indxA2, 2)-D(indxA2, 4));

% caseA3: water table in acrotelm
indxA3 = (D(:,3) >= D(:,2)) & (D(:,4) < D(:,1));
T(indxA3) = D(indxA3, 7).*(D(indxA3, 3)-D(indxA3, 2)) + D(indxA3, 6) .* (D(indxA3, 2)-D(indxA3, 1)) +  D(indxA3, 5) .* (D(indxA3, 1)-D(indxA3, 4));

%% Cases when ice table is in catotelm
% caseB1 : water table in catotelm
indxB1 = (D(:,3) < D(:,2)) & (D(:,3) > D(:,1)) & (D(:,4) < D(:,2)) & (D(:,4) > D(:,1));
T(indxB1) = D(indxB1,6).* (D(indxB1,3)-D(indxB1,4));

% case B2 : water table in acrotelm
indxB2 = (D(:,3) < D(:,2)) & (D(:,3) > D(:,1)) & (D(:,4) < D(:,2)) & (D(:,4) <= D(:,1));
T(indxB2) = D(indxB2, 6) .* (D(indxB2, 2) - D(indxB2, 1)) + D(indxB2, 6).* (D(indxB2, 1)-D(indxB2, 4));

%% Cases when ice tabel in acrotelm 
% caseC1 : water tabel in acrotelm
indxC1 = (D(:,3) < D(:,1)) & (D(:,4) < D(:,1));
T(indxC1) = D(indxC1, 5) .* (D(indxC1,3)-D(indxC1,4));

%% 
KEFF = T ./ (D(:,3)-D(:,4));
KEFF1 = KEFF(indxA1); KEFF2 = KEFF(indxA2); KEFF3 = KEFF(indxA3); KEFF4 = KEFF(indxB1); KEFF5 = KEFF(indxB2); KEFF6 = KEFF(indxC1);
KEFF1(KEFF1<=0) = []; figure; histogram(log10(KEFF1)); title('ice minsoil water minsoil');xlim([ -8 , -2]);
%exportgraphics(gca,'keffS\iceMinWaterMin.pdf','ContentType','vector');
KEFF2(KEFF2<=0) = []; figure; histogram(log10(KEFF2)); title('ice minsoil water cato');xlim([ -8 , -2]);
%exportgraphics(gca,'keffS\iceMinWaterCato.pdf','ContentType','vector');
KEFF3(KEFF3<=0) = []; figure; histogram(log10(KEFF3)); title('ice minsoil water acro');xlim([ -8 , -2]);
%exportgraphics(gca,'keffS\iceMinWaterAcro.pdf','ContentType','vector');
KEFF4(KEFF4<=0) = []; figure; histogram(log10(KEFF4)); title('ice cato water cato');xlim([ -8 , -2]);
%exportgraphics(gca,'keffS\iceCatoWaterCato.pdf','ContentType','vector');
KEFF5(KEFF5<=0) = []; figure; histogram(log10(KEFF5)); title('ice cato water acro');xlim([ -8 , -2]);
%exportgraphics(gca,'keffS\iceCatoWaterAcro.pdf','ContentType','vector');
KEFF6(KEFF6<=0) = []; figure; histogram(log10(KEFF6)); title('ice acro water acro');xlim([ -8 , -2]);
%exportgraphics(gca,'keffS\iceAcroWaterAcro.pdf','ContentType','vector');
KEFF(KEFF<0) = [];
figure; histogram(log10(KEFF),100);xlim([-8, -2]);

%% 
B = (D(:,3)-D(:,4));
B(B<=0) = [];
figure; histogram(B); xlim([0 , 1]);
% exportgraphics(gca,'saturated_thickenss_updated_1.pdf','ContentType','vector');
%%
T(T<=0)=[];
% figure; histogram(log10(T),100,'Orientation', 'horizontal');ylim([-10, -3]);
figure; histogram(log10(T),100);xlim([-10, -3])

%%
figure; 
histogram (D(:,1),100,'FaceColor',"#77AC30");xlim([0 max(z_AC)]);title('Acrotelm thickness (m)');
% hold on; plot(z_a,92473.*lognpdf(z_a,zAC_param(1), zAC_param(2)),'r','LineWidth',1); 
%% 
figure; 
histogram (D(:,2),100,'FaceColor',"#ffc20d"); xlim([0 max(dtC)]); title('depth to catotelm mineral soil boundary (m)');
% hold on; 
% z_c = linspace(0,max(dtC),100);
% plot(z_c, 88268*lognpdf(z_c,dtC_param(1), dtC_param(2)),'r','LineWidth',1);
%% 
figure; histogram (D(:,3),100,'FaceColor',"#F08080"); xlim([0 max(dti)]); title('depth to ice  (m)');
%hold on;  plot(x_dti,115748*lognpdf(x_dti,dti_param(1), dti_param(2)),'r','LineWidth',1);
%%
figure; histogram (D(:,4),100,'FaceColor',"#20b2aa");xlim([0 max(dtw)]); title('depth to water (m)');
% hold on; plot(x_dtw,124603 * lognpdf(x_dtw, dtw_param(1),dtw_param(2)),'r');
%% 
figure;histogram (D(:,5),100,'FaceColor',"#77AC30");xlim([0 max(K_AC_raw)]); title('Acrotelm hydraulic conductivity (m/s)');
figure;histogram (D(:,6),10000,'FaceColor',"#EDB120"); xlim([0 max(K_CT_raw)]); title('Catotelm hydraulic conductivity (m/s)');
figure;histogram (D(:,7),1000,'FaceColor',"#D95319");xlim([0 0.5e-5]); title('Mineral Soil hydraulic conductivity (m/s)');

%% Flux calculations
%% Kane weir
cell_size = 0.20524; % in m 
Q_gw_buf_10262 = 261.447; % buffer integrated values from ARC
Q_gw_buf_20524 = 1725.21;  %buffer integrated values from ARC
Buffer_2cell =  Q_gw_buf_10262;  
Q_gw_Kane_model_riparian = 1e3* Buffer_2cell .* T .* (cell_size .*2);
% writematrix(Q_gw_Kane_model_riparian,'Q_gw_Kane_model_revisionB.csv'); % in L/sL/s
%% DOC Kane
filename = "..\input_data\C_C_Hydrozoid.xlsx";
Conc_hydrozoid = readtable(filename);
DOC_hydrozoid_uM = table2array((Conc_hydrozoid(:,1)));
DOC_hydrozoid = DOC_hydrozoid_uM .*1e-6;
rng('default');
DOC_hydrozoid_parmhat = lognfit(DOC_hydrozoid); 
DOC_MC_hydrozoid = lognrnd(DOC_hydrozoid_parmhat(1),DOC_hydrozoid_parmhat(2),length(Q_gw_Kane_model_riparian),1);

figure;
histogram (DOC_hydrozoid,10,'FaceColor',"#B2ADAD");
xlabel('Groundwater DOC concentration (M)');
hold on;
DOC_GW = linspace(0,max(Q_gw_Kane_model_riparian),length(Q_gw_Kane_model_riparian));
plot(DOC_GW,0.26*lognpdf(DOC_GW,DOC_hydrozoid_parmhat(1),DOC_hydrozoid_parmhat(2)),'r','LineWidth',1); 
set(gca,'FontSize',10);   set(gca,'LineWidth',1); xlim([0 8e-3]);
figure; 
histogram (DOC_MC_hydrozoid,100,'FaceColor',"#B2ADAD");
xlabel('Groundwater DOC concentration (M)');
set(gca,'FontSize',10);   set(gca,'LineWidth',1);xlim([0 8e-3]);

DOC_flux = Q_gw_Kane_model_riparian .* DOC_MC_hydrozoid;
%csvwrite('DOC_flux_Kane_Model.csv',DOC_flux);
%% Pool2 Weir
Pool2_buffer_20524 = 573.69; 
Q_gw_Pool2_model_riparian = 1e3* Pool2_buffer_20524 .* T .* (cell_size .*3) ;
% writematrix(Q_gw_Pool2_model_riparian,'Q_gw_pool2_model_revisionB.csv'); % in L/s
%% end