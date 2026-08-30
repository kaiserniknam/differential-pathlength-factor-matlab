function [] = Fig_08_1
% FIG_08_1 - Differential MBLL analysis using G(d) = OD(mu_a = 0,d).
%
% Evaluates the recovery of absorption changes using the modified
% Beer-Lambert law (MBLL) in differential mode for several DPF models.
% Unlike Fig_07_1, the geometry/scattering contribution is defined from
% the zero-absorption simulation as
%
%                       G(d) = OD(mu_a = 0,d),
%
% so G is allowed to vary with source-detector distance.
%
% This is the non-assisted inversion: the unknown absorption coefficient
% is estimated by searching over the available mu_a-dependent DPF values.
% The DPF corresponding to the true absorption coefficient is therefore
% not supplied to the inversion.
%
% The script compares constant, semi-infinite, mean-pathlength,
% slope-based, true, inverse-distance, and empirical DPF models and
% evaluates absorption-change recovery versus source-detector distance.

clc
close all
format long g

% ---- parameters
set_of_mua = (0.00:0.01:0.75)+eps;
set_of_mus = 55;
dlta_d = .17;
L = 125;

% ---- variables
set_of________G = nan(1,                 L);
set_of_______OD = nan(length(set_of_mua),L);
set_of_cnst_DPF = nan(length(set_of_mua),L);
set_of_smif_DPF = nan(length(set_of_mua),L);
set_of_savg_DPF = nan(length(set_of_mua),L);
set_of_slop_DPF = nan(length(set_of_mua),L);
set_of_true_DPF = nan(length(set_of_mua),L);
set_of_idst_DPF = nan(length(set_of_mua),L);
set_of_empr_DPF = nan(length(set_of_mua),L);

for i_a = 1:length(set_of_mua)
    mua = set_of_mua(i_a);
    mus = set_of_mus(1);
    [dist, OD, cnst_DPF, smif_DPF, savg_DPF, slop_DPF, true_DPF, idst_DPF, empr_DPF] = get_DPF(mua,mus,dlta_d);
    if (i_a == 1)
        set_of________G(1:length(OD))       = OD;
    end
    set_of_______OD(i_a,1:length(OD))       = OD;
    set_of_cnst_DPF(i_a,1:length(cnst_DPF)) = cnst_DPF;
    set_of_smif_DPF(i_a,1:length(smif_DPF)) = smif_DPF;
    set_of_savg_DPF(i_a,1:length(savg_DPF)) = savg_DPF;
    set_of_slop_DPF(i_a,1:length(slop_DPF)) = slop_DPF;
    set_of_true_DPF(i_a,1:length(true_DPF)) = true_DPF;
    set_of_idst_DPF(i_a,1:length(idst_DPF)) = idst_DPF;
    set_of_empr_DPF(i_a,1:length(empr_DPF)) = empr_DPF;
    clearvars G OD cnst_DPF smif_DPF savg_DPF slop_DPF true_DPF idst_DPF empr_DPF
    clearvars mua mus
end
new_dist = nan(1,L); new_dist(1:length(dist)) = dist; dist = new_dist; clearvars new_dist L i_a

set_of_est_muas = nan(1+7,length(dist),length(set_of_mua));
set_of_r_squard = nan(1+7,length(dist),2);
idx_mua_involved = 2:26; 
for i_d = 1:length(dist)
    disp(['d = ',num2str(dist(i_d)),' cm'])
    for i_od = 1:size(set_of_______OD,1)
        % Calculate the difference in optical density
        ODm = set_of_______OD(i_od,i_d);
        OD0 = set_of_______OD(1   ,i_d);
        dOD = ODm - OD0;
        % true mua
        idx = 1;
        set_of_est_muas(idx,i_d,i_od) = set_of_mua(i_od); clearvars idx
        % constant dpf
        idx = 2;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_cnst_DPF, i_d); clearvars idx
        % semi-infinite dpf
        idx = 3;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_smif_DPF, i_d); clearvars idx
        % mean-pathlength dpf
        idx = 4;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_savg_DPF, i_d); clearvars idx
        % slope-based dpf
        idx = 5;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_slop_DPF, i_d); clearvars idx
        % true dpf
        idx = 6;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_true_DPF, i_d); clearvars idx
        % inverse-distance dpf
        idx = 7;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_idst_DPF, i_d); clearvars idx
        % empirical dpf
        idx = 8;
        set_of_est_muas(idx,i_d,i_od) = solve_equation (dOD, dist, set_of_mua, set_of_empr_DPF, i_d); clearvars idx
        clearvars ODm OD0 dOD
    end
    for idx = 1:8
        set_of_r_squard(idx,i_d,1) = r_squared (set_of_est_muas(1,i_d,idx_mua_involved),set_of_est_muas(idx,i_d,idx_mua_involved));
        set_of_r_squard(idx,i_d,2) = rmse_100  (set_of_est_muas(1,i_d,idx_mua_involved),set_of_est_muas(idx,i_d,idx_mua_involved));
    end
    clearvars i_od idx
end
% save('Fig_05_4.mat','dist','set_of_mua','set_of_est_muas','set_of_mus','set_of_r_squard','idx_mua_involved')
clear set_of________G set_of_______OD set_of_cnst_DPF set_of_smif_DPF set_of_savg_DPF set_of_slop_DPF set_of_true_DPF set_of_idst_DPF set_of_empr_DPF 
clearvars dlta_d i_d

% Plot R-squared
figure(1)
for idx = 1:8
    subplot(1,2,1)
    plot(dist, squeeze(set_of_r_squard(idx,:,1)),     'LineWidth', 2, 'Color', get_color(idx-1), 'LineStyle', '-', 'DisplayName', get_title(idx-1)); hold on
    subplot(1,2,2)
    plot(dist, squeeze(set_of_r_squard(idx,:,2))*100, 'LineWidth', 2, 'Color', get_color(idx-1), 'LineStyle', '-', 'DisplayName', get_title(idx-1)); hold on
end
subplot(1,2,1)
xlabel('Distance (cm)')
ylabel('R^2 of DPF over Distance')
title('Evolution of R^2 over d / not assited')
grid on, axis([min(dist) max(dist) 0 1])
grid on, axis([0 8 0 1])
set(gca, 'FontSize', 16, 'Box', 'on')
legend('show', 'Location', 'northeast', 'NumColumns', 1)
hold off
subplot(1,2,2)
xlabel('Distance (cm)')
ylabel('RMSE of DPF over Distance')
title('Evolution of Abs. RMSE over d / not assited')
grid on, axis([min(dist) max(dist) 0 100])
grid on, axis([0 8 0 250])
set(gca,'ytick',0:50:250)
set(gca, 'FontSize', 16, 'Box', 'on')
legend('show', 'Location', 'northeast', 'NumColumns', 1)
hold off
clearvars idx

% plot 1-to-1 plot
figure(2)
set_of_ids = [12,24,36]; % -> 2, 4, 6
for i_d = 1:length(set_of_ids)
    subplot(1,3,i_d)
    plot([min(set_of_mua(idx_mua_involved)) max(set_of_mua(idx_mua_involved))], [min(set_of_mua(idx_mua_involved)) max(set_of_mua(idx_mua_involved))], 'k:','DisplayName','1-to-1'), hold on
    for idx = 2:7
        plot(squeeze(set_of_est_muas(1,set_of_ids(i_d),:)), squeeze(set_of_est_muas(idx,set_of_ids(i_d),:)), '-o', ...
            'Color',get_color(idx-1),'MarkerFaceColor',get_color(idx-1),'MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2,'DisplayName',[get_title(idx-1),': e = ',num2str(set_of_r_squard(idx,set_of_ids(i_d),2)*100,'%0.2f')]), hold on
    end
    xlabel('True \Delta\mu_a (cm^{-1})')
    ylabel('Estimated \Delta\mu_a (cm^{-1})')
    grid on, axis equal, axis([0 max(set_of_mua(idx_mua_involved)) min(set_of_mua(idx_mua_involved)) 2*max(set_of_mua(idx_mua_involved))]), hold off
    set(gca,'FontSize',12)
    legend('show', 'Location','southeast', 'NumColumns',1)
    title([{'True vs. Estimated \Delta\mu_a / not assited'},{['at d = ',num2str(dist(set_of_ids(i_d))),' cm']}])
    %----- Error table -------------------------------------------------------
    err = squeeze(set_of_r_squard(2:7,set_of_ids(i_d),2))*100;
    tbl = cell(7,2);
    tbl(1,:) = {'Model','e (%)'};
    for k = 1:6
        tbl{k+1,1} = get_title(k);
        tbl{k+1,2} = sprintf('%.2f',err(k));
    end
    annotation_str = sprintf(['Model      e (%%)\n',...
                              'constant   %5.2f\n',...
                              'semi-inf   %5.2f\n',...
                              'mean-path  %5.2f\n',...
                              'slope      %5.2f\n',...
                              'true       %5.2f\n',...
                              'inv-dist   %5.2f'],err);
    text(0.03,0.97,annotation_str,...
        'Units','normalized',...
        'VerticalAlignment','top',...
        'FontName','Courier',...
        'FontSize',10,...
        'BackgroundColor','w',...
        'EdgeColor','k');
    %-----------------------------------------------------------------------
end

fprintf('\n==============================================================\n');
fprintf('Figure 2: Mean Absolute Relative Error (%%)\n');
fprintf('==============================================================\n');
fprintf('%-20s %10s %10s %10s\n', ...
    'Model', ...
    sprintf('%.1f cm',dist(set_of_ids(1))), ...
    sprintf('%.1f cm',dist(set_of_ids(2))), ...
    sprintf('%.1f cm',dist(set_of_ids(3))));
fprintf('--------------------------------------------------------------\n');

for idx = 2:7
    fprintf('%-20s %10.2f %10.2f %10.2f\n', ...
        get_title(idx-1), ...
        100*set_of_r_squard(idx,set_of_ids(1),2), ...
        100*set_of_r_squard(idx,set_of_ids(2),2), ...
        100*set_of_r_squard(idx,set_of_ids(3),2));
end

fprintf('==============================================================\n');
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
function [dist, OD, DPF_cnst, DPF_smif, DPF_savg, DPF_slop, DPF_true, DPF_idst, DPF_empr] = get_DPF (mua, mus, dlta_d)
g = 0.93;

Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_81_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat' ];
t_db = load(the_filename); clearvars the_filename

% removing noisy points
u = round(-log(t_db.w)./t_db.s,4); u_unique = unique(u); freq = nan(size(u_unique));
for i_f = 1:length(u_unique), freq(i_f) = sum(u==u_unique(i_f)); end
[~,i_fx] = max(freq);
t_db.x = t_db.x(u==u_unique(i_fx)); 
t_db.y = t_db.y(u==u_unique(i_fx)); 
t_db.z = t_db.z(u==u_unique(i_fx)); 
t_db.d = t_db.d(u==u_unique(i_fx)); 
t_db.s = t_db.s(u==u_unique(i_fx)); 
t_db.w = t_db.w(u==u_unique(i_fx)); 
t_db.c = t_db.c(u==u_unique(i_fx)); 
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
clearvars u u_unique freq i_f i_fx        

% 1-D sorting
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d; clearvars Lx Ly Lz
[~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_diff_edges); % d_diffuse bins
dist = 1/2.*(d_diff_edges(1:end-1)+d_diff_edges(2:end-0));
clearvars d_diff_edges Lx Ly Lz

% OD
fun_x = @mean; fun_y = @sum;
index_in = ind_diff; TheCode = 0;
x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
y_temp = t_db.w(t_db.c==TheCode);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = -log( ...
    (y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));
idx = ~isnan(x_bind);
OD = interp1(x_bind(idx),y_bind(idx),dist,'nearest','extrap');
clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind index_in idx

% cnst_DPF & smif_DPF
musp = mus*(1-g);
num = sqrt(3*musp).*(dist.*sqrt(3*mua*musp)+0);
den = 2*sqrt(mua) .*(dist.*sqrt(3*mua*musp)+1);
DPF_smif = num./den;
DPF_cnst = (sqrt(3*musp))./(2*sqrt(mua)).*ones(size(dist));
clearvars musp num den

% savg_DPF
fun_x = @mean; fun_y = @mean; index_in = ind_diff; TheCode = 0;
x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
y_temp = t_db.s(t_db.c==TheCode);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
DPF_savg = [x_bind , y_bind./x_bind]; idx = ~isnan(DPF_savg(:,1));
DPF_savg = interp1(DPF_savg(idx,1),DPF_savg(idx,2),dist,'nearest','extrap');
clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind idx index_in idx

% slop_DPF
fun_x = @mean; fun_s = @(x)(sum(x.*exp(-mua.*x))./sum(exp(-mua.*x)));
index_in = ind_diff; TheCode = 0;
x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
s_temp = t_db.s(t_db.c==TheCode);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = s_bind./x_bind;                      
DPF_slop = [x_bind , s_bind]; idx = ~isnan(DPF_slop(:,1));
DPF_slop = interp1(DPF_slop(idx,1),DPF_slop(idx,2),dist,'nearest','extrap');
clearvars fun_x fun_s index_in TheCode x_temp s_temp x_bind s_bind index_in idx
        
% true_DPF
fun_x = @mean; fun_y = @sum;
index_in = ind_diff; TheCode = 0;
x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
y_temp = t_db.w(t_db.c==TheCode);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = -log( ...
    (y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));

persistent cached_mus cached_dlta_d cached_dist_mua_0 cached_OD_mua_0
if mua==0
    G = y_bind;
    cached_mus = mus;
    cached_dlta_d = dlta_d;
    cached_dist_mua_0 = dist;
    cached_OD_mua_0 = OD;
else
    if isempty(cached_OD_mua_0) || cached_mus~=mus || cached_dlta_d~=dlta_d
        [cached_dist_mua_0,cached_OD_mua_0,~,~,~,~,~,~,~] = get_DPF(0,mus,dlta_d);
        cached_mus = mus;
        cached_dlta_d = dlta_d;
    end
    G = interp1(cached_dist_mua_0,cached_OD_mua_0,x_bind,'nearest','extrap');
end

y_bind = (y_bind-G)./x_bind./mua;
idx = ~isnan(x_bind);
DPF_true = interp1(x_bind(idx),y_bind(idx),dist,'nearest','extrap');
clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind index_in idx

% idst_DPF
fun_x = @mean; fun_y = @sum;
index_in = ind_diff; TheCode = 0;
x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
y_temp = t_db.w(t_db.c==TheCode);
x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
y_bind = -log( ...
    (y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));

x = x_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
y = y_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
G = G     (~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
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
DPF_idst = interp1(x_bind(idx),A+B./x_bind(idx)+C.*log(x_bind(idx))./x_bind(idx),dist,'nearest','extrap');
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
function [r2] = r_squared(y,y_hat)
y = y(:); y_hat = y_hat(:);
% align indices where both observed and predicted exist (non-NaN)
idx = ~isnan(y) & ~isnan(y_hat);
% compute r-squared using only indices where both observed and predicted exist (ss_res and ss_tot computed below)
ss_res = sum((y(idx) - y_hat(idx)).^2);
ss_tot = sum((y(idx) - mean(y(idx))).^2);
r2 = 1 - ss_res / ss_tot;
end
function [r2] = rmse_100(y,y_hat)
y = y(:); y_hat = y_hat(:);
% align indices where both observed and predicted exist (non-NaN)
idx = ~isnan(y) & ~isnan(y_hat);
% compute r-squared using only indices where both observed and predicted exist (ss_res and ss_tot computed below)
r2 = mean(abs( (y(idx) - y_hat(idx)) ./ y(idx) ));
end
function [true_mu] = solve_equation (dOD, dist, set_of_mua, set_of_DPF, i_d)
d_error = ((dOD/dist(i_d) + set_of_mua(1)*set_of_DPF(1,i_d))./set_of_DPF(:,i_d)./set_of_mua.'-1);
[~,idx] = min(abs(d_error));
true_mu = set_of_mua(idx);
end
function [out] = get_title(idx)
if     idx==1
    out = 'constant';
elseif idx==2
    out = 'semi-inf';
elseif idx==3
    out = 'mean-pathlength';
elseif idx==4
    out = 'slope-based';
elseif idx==5
    out = 'true';
elseif idx==6
    out = 'inverse-distance';
elseif idx==7
    out = 'empirical';
else
    out = 'undefined';
end
end

