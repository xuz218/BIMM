%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for plotting simulated protein profiles of signature genes SMAD2 
% and TGFβ1, simulated expression level of signature genes FOXO3 and
% TGFβR1, and DEX level over time after one dose of DEX treatment applied  
% to the A549 cells
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
clc; clear; close all;
load('cancer_ga_outputs53.mat', "value", "error"); %  different parameters can be loaded here 
% load('gene_percentile.mat');

Inits = value;

global EGF ERBB2 EGF_ERBB2 MAPK1 FOXO3 TGFB2 TGFBR1 TGFB2_TGFBR1 SMAD23 pSMAD23 MYC DRUG ERBB2_DRUG

% y index for each species
EGF = 1;
ERBB2 = 2;
EGF_ERBB2 = 3;
MAPK1 = 4;
FOXO3 = 5;
TGFB2 = 6;
TGFBR1 = 7;
TGFB2_TGFBR1 = 8;
SMAD23 = 9;
pSMAD23 = 10;
MYC = 11;
ERBB2_DRUG = 12;
DRUG = 13;

t_data = [0 1 3];
FOXO3_data = [1.5403 1.7199 1.7289];
TGFBR1_data = [1.3959 1.6867 1.6087];

FOXO3_scaled = FOXO3_data./FOXO3_data(1);
TGFB2_scaled = [1 1.12];
TGFBR1_scaled = TGFBR1_data./TGFBR1_data(1);
pSMAD23_scaled = [1 1.18];

tstable = [0 1000];
tspan = [0 50];
t_index = t_data + 2 + tstable(end);

y0 = ones(1,13);
y0(ERBB2_DRUG) = 0;
y0(DRUG) = 0;

odefun = @(t, y) model(t, y, Inits);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
[TOUT1, YOUT1] = ode23tb(odefun, tstable, y0, options);

y0 = YOUT1(end,:);
y0(end) = 10;
[TOUT2, YOUT2] = ode23tb(odefun, tspan, y0, options);

TOUT = [TOUT1; TOUT2+1000];
YOUT = [YOUT1; YOUT2];
YOUT2 = [YOUT2 YOUT2(:,TGFBR1)+YOUT2(:,TGFB2_TGFBR1) YOUT2(:,TGFB2)+YOUT2(:,TGFB2_TGFBR1)];
YOUT2_scaled = YOUT2./YOUT2(1,:);

%% general figure 
% figure(1)
% hold on;
% for i=1:13
%     plot(TOUT,YOUT(:,i)) 
% end
% legend('EGF','ERBB2','EGF-ERBB2','MAPK1','FOXO3','TGFB2','TGFBR1','TGFB2-TGFBR1', ...
%     'SMAD23', 'pSMAD23','MYC','ERBB2-DG','DRUG');
extraInputs = {'fontsize',20};
% xlabel('Time',extraInputs{:});
% ylabel('Population',extraInputs{:});
% box on
path ='./result';
% saveas(gcf, fullfile(path, ['ga_out1','.png'] ));

%% figure of protein level of SMAD23 
ytime = 1:1:length(TOUT2);
figure(2)
hold on;
plot(TOUT2(ytime),YOUT2_scaled(ytime, pSMAD23),'k','LineWidth',3)

scatter([0 48], pSMAD23_scaled, 100, 'r', 'filled');
axis([0 50,0.8 2.2])
xlabel('Time (hour)',extraInputs{:});
ylabel('Protein level',extraInputs{:});
ax = gca;ax.FontSize = 20;
r = rectangle('Position',[0 0 50 2.2]');
title('pSMAD2/3');
s_leg =legend('Simulation','Empirical');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['pSMAD23','.png'] ));

%% TGFB1
figure(3)
hold on;
plot(TOUT2(ytime),YOUT2_scaled(ytime, end),'k','LineWidth',3) 
scatter([0 48], TGFB2_scaled, 100, 'r', 'filled');
axis([0 50,0.8 2.2])
xlabel('Time (hour)',extraInputs{:});
ylabel('Protein level',extraInputs{:});
r = rectangle('Position',[0 0 50 2.2]');
ax = gca;ax.FontSize = 20;
title('TGFβ1');
s_leg =legend('Simulation','Empirical');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['TGFB2','.png'] ));

%% FOXO3
figure(4)
hold on;
plot(TOUT2(ytime),YOUT2_scaled(ytime, FOXO3),'k','LineWidth',3) 
scatter(t_data, FOXO3_scaled, 100, 'r', 'filled');
axis([0 50,0.8 1.4])
xlabel('Time (hour)',extraInputs{:});
ylabel('Expression level',extraInputs{:});
r = rectangle('Position',[0 0 50 1.4]');
ax = gca;ax.FontSize = 20;
title('FOXO3');
s_leg =legend('Simulation','Empirical');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['FOXO3','.png'] ));

%% TGFBR1
figure(5)
hold on;
plot(TOUT2(ytime,1),YOUT2_scaled(ytime, TGFBR1),'k','LineWidth',3) 
scatter(t_data, TGFBR1_scaled, 100, 'r', 'filled')
axis([0 50,0.8 1.4])
xlabel('Time (hour)',extraInputs{:});
ylabel('Expression level',extraInputs{:});
r = rectangle('Position',[0 0 50 1.4]');
ax = gca;ax.FontSize = 20;
title('TGFβR1');
s_leg =legend('Simulation','Empirical');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['TGFBR1','.png'] ));

%% DRUG
figure(6)
hold on;
sum = YOUT2(1:length(TOUT2),DRUG) + YOUT2(1:length(TOUT2),ERBB2_DRUG);
plot(TOUT2,sum./10,'k','LineWidth',3); 
axis([0 50,0 1])
xlabel('Time (hour)',extraInputs{:});
ylabel('Drug level',extraInputs{:});
r = rectangle('Position',[0 0 50 1]');
ax = gca;ax.FontSize = 20;
title('DEX');
saveas(gcf, fullfile(path, ['DRUG','.png'] ));
