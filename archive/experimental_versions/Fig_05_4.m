function [] = Fig_05_4
% Convert Arduino readings into OD-vs-distance data + 13x13 map
% Comparison of optical density (OD) for concentrations of 0.25, 1, and 4 mg/mL
% R_L = 200 kΩ, R_LED = 150 Ω, wavelength = 750 nm
% Based on version #3, with a reduced set of figures and using reflectance (R) instead of optical density (OD).

clc
close all
format long g
addpath('/home/kaiser/Dropbox/KSR/University of Houston/Experiment/2.NIRs Project/008 DFP Experiment/code/Characteristics/')

% ---- parameters
R_L    = 200e3;
R_LED  = 150;
lambda = 750;
mu_a = 0.0275;
dlta_d = .17;

[OD_sim_0_25,dist_sim_0_25] = get_OD (0.25,dlta_d);
[OD_sim_1_00,dist_sim_1_00] = get_OD (1.00,dlta_d);

% IMPORTANT: set this to your actual spatial pitch
pitch_cm = 1.0;   % <-- CHANGE if each grid step is not 1 cm

% ---- read
fname = "/home/kaiser/Dropbox/KSR/University of Houston/Experiment/2.NIRs Project/008 DFP Experiment/code/db/OD_grid_2025-12-29-11-40-40-750nm-OD-1gL.txt";     % 1.00 g/L
T = readtable(fname, "FileType","text", "Delimiter","\t"); clearvars fname
db_1_00 = T.Variables; clearvars T fname % assumes fixed column order
fname = "/home/kaiser/Dropbox/KSR/University of Houston/Experiment/2.NIRs Project/008 DFP Experiment/code/db/OD_grid_2025-12-30-12-15-15-750nm-OD-0_25gL.txt";  % 0.25 g/L
T = readtable(fname, "FileType","text", "Delimiter","\t"); clearvars fname
db_0_25 = T.Variables; clearvars T fname % assumes fixed column order

% Columns assumed (based on your code):
% 1=time_ms, 3=voltage_V, 4=i_row, 5=i_col, 6=recording_flag (1=ON)

% ---- extract intensity map
N = 13;
% Preallocate max possible vector length: N*N points
The_Vector = nan(N*N, 1+2);
k = 0;
for i_row = 1:N
    for i_col = 1:N
        sig = db_1_00(db_1_00(:,4)==i_row & db_1_00(:,5)==i_col & db_1_00(:,6)==1, 3);
        sig(sig<eps) = nan;
        if isempty(sig) || all(isnan(sig))
            Io_1 = nan;
        else
            % Vmean = mean(sig, "omitnan");
            Vmean = mean(sig);
            [~, Io_1] = do_FDS100(Vmean, R_L, lambda);
        end
        dist_steps = sqrt((i_row-7).^2 + (i_col-7).^2);
        dist_cm    = dist_steps * pitch_cm;
        clearvars dist_steps sig Vmean

        sig = db_0_25(db_0_25(:,4)==i_row & db_0_25(:,5)==i_col & db_0_25(:,6)==1, 3);
        sig(sig<1e-6) = nan;
        if isempty(sig) || all(isnan(sig))
            Io_25 = nan;
        else
            Vmean = mean(sig, "omitnan");
            [~, Io_25] = do_FDS100(Vmean, R_L, lambda);
        end
        dist_steps = sqrt((i_row-7).^2 + (i_col-7).^2);
        dist_cm    = dist_steps * pitch_cm;
        clearvars dist_steps sig Vmean

        k = k + 1;
        The_Vector(k,:) = [dist_cm, Io_25 Io_1];
        clearvars dist_cm Io_25 Io_1 Io_4
    end
end
The_Vector = The_Vector(1:k,:); % trim
clearvars N k i_row i_col pitch_cm

% ---- source intensity I0
if lambda == 750
    [~, I0, ~] = do_LED750L(R_LED);
else
    [~, I0, ~] = do_LED850LN(R_LED);
end

% Normalize
The_Vector(:,2:3) = The_Vector(:,2:3) ./ I0;
clearvars I0 R_L R_LED lambda



% ---- plots: OD vs distance, and DPF vs distance

OD = The_Vector(:,2:end);
OD = (OD);

r2_0_25 = r_squared(dist_sim_0_25,(OD_sim_0_25),The_Vector(:,1),OD(:,1));
semilogy(dist_sim_0_25, (OD_sim_0_25), '-', 'Color',get_color(3),'LineWidth',2,'LineStyle','-','DisplayName','MC - 0.25 g/L'), hold on
semilogy(The_Vector(:,1), OD(:,1), 'o', ...
    'MarkerFaceColor',get_color(3),'MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2,'DisplayName','data - 0.25 g/L'), hold on
r2_1_00 = r_squared(dist_sim_1_00,(OD_sim_1_00),The_Vector(:,1),OD(:,2));
semilogy(dist_sim_1_00, (OD_sim_1_00), '-', 'Color',get_color(2),'LineWidth',2,'LineStyle','-','DisplayName','MC - 1.00 g/L'), hold on
semilogy(The_Vector(:,1), OD(:,2), 'o', ...
    'MarkerFaceColor',get_color(2),'MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2,'DisplayName','data - 1.00 g/L'), hold on
xlabel('separation (cm)')
ylabel('Reflectance, R (unitless)')
grid on, axis square, hold off, axis([0 8 10^(-6) 10^(-1)])
set(gca,'fontsize',18), legend('show','Location','northeast','NumColumns',1)
title('R vs source-detector separation')
disp(['r^2 for C = 0.25 g/L -> ',num2str(r2_0_25,'%.2f')]), % clearvars r2_0_25
disp(['r^2 for C = 1.00 g/L -> ',num2str(r2_1_00,'%.2f')]), % clearvars r2_1_00
clearvars db_1_00 db_0_25 dist_sim_0_25 dist_sim_1_00 OD_sim_0_25 OD_sim_1_00 OD

end

function [out] = get_color(idx)
if     idx==1
    out = [0.0000 0.4470 0.7410];
elseif idx==2
    out = [0.8500 0.3250 0.0980];
elseif idx==3
    out = [0.9290 0.6940 0.1250];
elseif idx==4
    out = [0.4940 0.1840 0.5560];
elseif idx==5
    out = [0.4660 0.6740 0.1880];
elseif idx==6
    out = [0.3010 0.7450 0.9330];
elseif idx==7
    out = [0.6350 0.0780 0.1840];
else
    out = [0.0000 0.0000 0.0000];
end
end
function [y_bind,x_bind] = get_OD (Cx,dlta_d)
the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_87_concentration_',num2str(Cx,'%.2f'),'.mat'];
Lx = 19.1; Ly = 19.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename); clearvars the_filename

% removing noisy points
u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); u_unique = u_unique(~isnan(u_unique)); freq = nan(size(u_unique));
for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
[~,i_fx] = max(freq);
t_db.x_in = t_db.x_in(u==u_unique(i_fx)); 
t_db.y_in = t_db.y_in(u==u_unique(i_fx)); 
t_db.z_in = t_db.z_in(u==u_unique(i_fx)); 
t_db.x_ot = t_db.x_ot(u==u_unique(i_fx)); 
t_db.y_ot = t_db.y_ot(u==u_unique(i_fx)); 
t_db.z_ot = t_db.z_ot(u==u_unique(i_fx)); 
t_db.s = t_db.s(u==u_unique(i_fx)); 
t_db.w = t_db.w(u==u_unique(i_fx));
u = u(~isnan(u));
disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100 u_unique(i_fx)]))
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
clearvars u u_unique freq i_f i_fx        

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d;
idx = t_db.z_ot<=0;
[~,~,index_in] = histcounts(sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2),d_diff_edges); % d_diffuse bins
clearvars d_diff_edges Lx Ly Lz

% I vs. d
fun_x = @mean; fun_y = @sum;
x_temp = sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2);
y_temp = t_db.w(idx);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = ...
    (            y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d));
end
function [DPF_cnst, DPF_smif, DPF_savg, DPF_slop, DPF_true, DPF_idst, DPF_empr,G] = get_DPF (Cx,dist,dlta_d)
g = 0.93; mua = 0.0275; mus = 55.0.*Cx;

Lx = 19.1; Ly = 19.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_87_concentration_',num2str(Cx,'%.2f'),'.mat'];
t_db = load(the_filename); clearvars the_filename

% removing noisy points
u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); u_unique = u_unique(~isnan(u_unique)); freq = nan(size(u_unique));
for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
[~,i_fx] = max(freq);
t_db.x_in = t_db.x_in(u==u_unique(i_fx)); 
t_db.y_in = t_db.y_in(u==u_unique(i_fx)); 
t_db.z_in = t_db.z_in(u==u_unique(i_fx)); 
t_db.x_ot = t_db.x_ot(u==u_unique(i_fx)); 
t_db.y_ot = t_db.y_ot(u==u_unique(i_fx)); 
t_db.z_ot = t_db.z_ot(u==u_unique(i_fx)); 
t_db.s = t_db.s(u==u_unique(i_fx)); 
t_db.w = t_db.w(u==u_unique(i_fx));
u = u(~isnan(u));
disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100 u_unique(i_fx)]))
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
clearvars u u_unique freq i_f i_fx        

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d;
idx_involved = t_db.z_ot<=0;
[~,~,ind_diff] = histcounts(sqrt( ...
    (t_db.x_in(idx_involved)-t_db.x_ot(idx_involved)).^2 + ...
    (t_db.y_in(idx_involved)-t_db.y_ot(idx_involved)).^2 + ...
    (t_db.z_in(idx_involved)-t_db.z_ot(idx_involved)).^2),d_diff_edges); % d_diffuse bins
clearvars d_diff_edges Lx Ly Lz

% cnst_DPF & smif_DPF
musp = mus*(1-g);
num = sqrt(3*musp).*(dist.*sqrt(3*mua*musp)+0);
den = 2*sqrt(mua) .*(dist.*sqrt(3*mua*musp)+1);
DPF_smif = num./den;
DPF_cnst = (sqrt(3*musp))./(2*sqrt(mua)).*ones(size(dist));
clearvars musp num den

% savg_DPF
fun_x = @mean; fun_y = @mean; index_in = ind_diff; 
x_temp = sqrt(t_db.x_ot(idx_involved).^2+t_db.y_ot(idx_involved).^2+t_db.z_ot(idx_involved).^2); 
y_temp = t_db.s(idx_involved);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
DPF_savg = [x_bind , y_bind./x_bind]; idx = ~isnan(DPF_savg(:,1));
DPF_savg = interp1(DPF_savg(idx,1),DPF_savg(idx,2),dist,'nearest');
clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind idx index_in idx

% slop_DPF
fun_x = @mean; fun_s = @(x)(sum(x.*exp(-mua.*x))./sum(exp(-mua.*x)));
index_in = ind_diff; 
x_temp = sqrt(t_db.x_ot(idx_involved).^2+t_db.y_ot(idx_involved).^2+t_db.z_ot(idx_involved).^2); 
s_temp = t_db.s(idx_involved);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = s_bind./x_bind;                      
DPF_slop = [x_bind , s_bind]; idx = ~isnan(DPF_slop(:,1));
DPF_slop = interp1(DPF_slop(idx,1),DPF_slop(idx,2),dist,'nearest');
clearvars fun_x fun_s index_in TheCode x_temp s_temp x_bind s_bind index_in idx
        
% true_DPF
fun_x = @mean; fun_y = @sum;
index_in = ind_diff; 
x_temp = sqrt(t_db.x_ot(idx_involved).^2+t_db.y_ot(idx_involved).^2+t_db.z_ot(idx_involved).^2);
y_temp = t_db.w(idx_involved);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = -log( ...
    (y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));
G = y_bind(1); % G
y_bind = (y_bind-G)./x_bind./mua;
idx = ~isnan(x_bind);
DPF_true = interp1(x_bind(idx),y_bind(idx),dist,'nearest');
clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind index_in idx

% idst_DPF
fun_x = @mean; fun_y = @sum;
index_in = ind_diff; 
x_temp = sqrt(t_db.x_ot(idx_involved).^2+t_db.y_ot(idx_involved).^2+t_db.z_ot(idx_involved).^2);
y_temp = t_db.w(idx_involved);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = -log( ...
    (y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));

x = x_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
y = y_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);        
Y = (y - G) ./ x;        % response
invX = 1 ./ x;           % predictor 1
logX_over_X = log(x) ./ x; % predictor 2
T = table(Y, invX, logX_over_X, 'VariableNames', {'Y','invX','logX_over_X'});
mdl = fitlm(T, 'Y ~ 1 + invX + logX_over_X');  % A, B, C
clearvars T x y Y invX logX_over_X   
A = mdl.Coefficients.Estimate(1)./mua;
B = mdl.Coefficients.Estimate(2)./mua;
C = mdl.Coefficients.Estimate(3)./mua;
idx = ~isnan(x_bind);
DPF_idst = interp1(x_bind(idx),A+B./x_bind(idx)+C.*log(x_bind(idx))./x_bind(idx),dist,'nearest');
clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind idx index_in A B C mdl

A_coeff = [0.87, -0.51, 0.41];
B_coeff = [0.02, -1.40, 0.72];
C_coeff = [0.03, -1.25, 0.55];
C = A_coeff(1); a = A_coeff(2); b = A_coeff(3); A = C.*(mua.^a).*(mus.^b); clearvars a b c
C = B_coeff(1); a = B_coeff(2); b = B_coeff(3); B = C.*(mua.^a).*(mus.^b); clearvars a b c
C = C_coeff(1); a = C_coeff(2); b = C_coeff(3); C = C.*(mua.^a).*(mus.^b); clearvars a b c
DPF_empr = A + B./dist + C.*log(dist)./dist; 
clearvars A B C a b c mua mus A_coeff B_coeff C_coeff

% figure(7)
% plot(dist,DPF_cnst, 'LineWidth', 2, 'DisplayName', 'const'),      hold on
% plot(dist,DPF_smif, 'LineWidth', 2, 'DisplayName', 'semi-inf'),   hold on
% plot(dist,DPF_savg, 'LineWidth', 2, 'DisplayName', 'mean'),       hold on
% plot(dist,DPF_slop, 'LineWidth', 2, 'DisplayName', 'slope'),      hold on
% plot(dist,DPF_true, 'LineWidth', 2, 'DisplayName', 'true'),       hold on
% plot(dist,DPF_idst, 'LineWidth', 2, 'DisplayName', 'inv. dist.'), hold on
% plot(dist,DPF_empr, 'LineWidth', 2, 'DisplayName', 'empirical'),  hold on
end
function [r2] = r_squared(x,y,x_hat,y_hat)
x = x(~isnan(y));
y = y(~isnan(y));
y_org = interp1(x,y,x_hat,"linear");

idx = ~isnan(y_hat);
ss_res = sum((y_org(idx) - y_hat(idx)).^2);
ss_tot = sum((y_org(idx) - mean(y_org(idx))).^2);
r2 = 1 - ss_res / ss_tot;
end


