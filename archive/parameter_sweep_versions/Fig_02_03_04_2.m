function [] = Fig_02_03_04_2
% Fig-02/03/04: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - to plot p and q for mu_a's and mu_s's
% based on version #1, with enhanced presentation and readability.

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

n = 1.4; g = 0.95;
set_of_mua = (.00:.05:0.5);
set_of_mus = (000:050:500);
% N_bins = 450;
N_bins = 200;
d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); clearvars Lx Ly Lz
d_cntrs = 1/2*(d_edges(1:end-1)+d_edges(2:end-0));
set_of_cnst_DPF = nan(length(set_of_mua),length(set_of_mus),1);
set_of_idst_DPF = nan(length(set_of_mua),length(set_of_mus),2);
set_of_adst_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_savg_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_slop_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_true_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of________I = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_dpfs_err = nan(length(set_of_mua),length(set_of_mus),8,N_bins);
clearvars N_bins

p_est = nan(size(squeeze(set_of_idst_DPF(:,:,1))));
q_est = nan(size(squeeze(set_of_idst_DPF(:,:,2))));
p_anc = nan(size(squeeze(set_of_idst_DPF(:,:,1))));
q_anc = nan(size(squeeze(set_of_idst_DPF(:,:,2))));
clc

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        % p_est(i_a,i_s) = exp(-0.707 - 0.507*log(set_of_mua(i_a)) + 0.506*log(set_of_mus(i_s)));
        % q_est(i_a,i_s) = exp(+1.135 - 1.370*log(set_of_mua(i_a)) - 0.198*log(set_of_mus(i_s)));
        
        p_coeff = [0.21271, -0.63068, +0.64536];
        q_coeff = [64.4451, -1.07500, -0.76928];
        p_est(i_a,i_s) = p_coeff(1).*((set_of_mua(i_a)).^(p_coeff(2))).*((set_of_mus(i_s)).^(p_coeff(3)));
        q_est(i_a,i_s) = q_coeff(1).*((set_of_mua(i_a)).^(q_coeff(2))).*((set_of_mus(i_s)).^(q_coeff(3)));

        mua  = set_of_mua(i_a);
        musp = set_of_mus(i_s)*(1-g);
        D = 1./(3.*(mua+musp));
        mueff = sqrt(3.*mua.*(mua+musp));
        z0 = 1./(mua+musp);
        Reff = -1.440./n.^2 + 0.710./n + 0.668 + 0.0636.*n;
        A = (1 + Reff)./(1 - Reff);
        zb = 2.*A.*D;
        
        p_anc(i_a,i_s) = mueff/mua;
        q_anc(i_a,i_s) = -log(1/4/pi/D*z0*(z0+2*zb))/mua;
        clearvars mua musp D mueff z0 Reff A zb
    end
end
clearvars i_a i_s

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        % read dbase
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_33_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
        t_db = load(the_filename); clearvars the_filename
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
        
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
        t_db.a = t_db.a(u==u_unique(i_fx));
        disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
        t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
        clearvars u u_unique freq i_f i_fx        
    
        % 1-D sorting
        [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_edges); % d_diffuse bins
        clearvars d_trns_edges d_diff_edges

        % s vs. d
        fun_x = @mean; fun_y = @mean;
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
        y_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
        set_of_savg_DPF(i_a,i_s,1:length(y_bind)) = y_bind./x_bind;
        clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind
        
        % slope DPF
        fun_x = @mean; fun_s = @(x)(sum(x.*exp(-mua.*x))./sum(exp(-mua.*x)));
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
        s_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = s_bind./x_bind;                      
        set_of_slop_DPF(i_a,i_s,1:length(s_bind)) = s_bind;
        clearvars fun_x fun_s index_in TheCode x_temp s_temp x_bind s_bind

        % I vs. d
        fun_x = @mean; fun_y = @sum; fun_s = @(x)(sum(exp(-mua.*x))./t_db.no_of_photons); TheOutFun = @(x)(-log(x));
        index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
        y_temp = t_db.w(t_db.c==TheCode)./t_db.no_of_photons;
        s_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = TheOutFun(y_bind);
        s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = TheOutFun(s_bind)./mua./x_bind;

        mdl = fitlm(1./x_bind(~isnan(x_bind)),y_bind(~isnan(x_bind))./x_bind(~isnan(x_bind)));
        % plot(x_bind,y_bind./x_bind,x_bind,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)./x_bind,'LineWidth',2)
        % title(['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus)])
        % axis ([0 9 0 50]), axis square, pause(0.01)
        p = mdl.Coefficients.pValue(2);
        if p<=0.05
            % set_of_cnst_DPF(i_a,i_s,1) = (x_bind(~isnan(x_bind))\y_bind(~isnan(x_bind)))/mua;
            set_of_idst_DPF(i_a,i_s,1) = mdl.Coefficients.Estimate(1)/mua;
            set_of_idst_DPF(i_a,i_s,2) = mdl.Coefficients.Estimate(2)/mua;
        else
            % set_of_cnst_DPF(i_a,i_s,1) = nan;
            set_of_idst_DPF(i_a,i_s,1) = nan;
            set_of_idst_DPF(i_a,i_s,2) = nan;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        musp = mus*(1-g);
        num = sqrt(3*musp).*(x_bind.*sqrt(3*mua*musp)+0);
        den = 2*sqrt(mua) .*(x_bind.*sqrt(3*mua*musp)+1);
        dpf_tmp = num./den;

        set_of_cnst_DPF(i_a,i_s,1) = (sqrt(3*musp))./(2*sqrt(mua));
        set_of_adst_DPF(i_a,i_s,1:length(dpf_tmp)) = dpf_tmp;
        clearvars musp num den dpf_tmp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set_of________I(i_a,i_s,1:length(y_bind)) = y_bind;
        set_of_true_DPF(i_a,i_s,1:length(s_bind)) = s_bind;
        
        % error
        set_of_dpfs_err(i_a,i_s,1,1:length(x_bind)) = ((y_bind./x_bind./set_of_cnst_DPF(i_a,i_s,1))-mua)./mua; % constant
        set_of_dpfs_err(i_a,i_s,2,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_savg_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % <s>/d
        set_of_dpfs_err(i_a,i_s,3,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_true_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % true
        set_of_dpfs_err(i_a,i_s,4,1:length(x_bind)) = ((y_bind./x_bind./(set_of_idst_DPF(i_a,i_s,1)+set_of_idst_DPF(i_a,i_s,2)./x_bind))-mua)./mua; % iverse-distance
        set_of_dpfs_err(i_a,i_s,5,1:length(x_bind)) = ((y_bind./x_bind./(p_est(i_a,i_s)            +q_est(i_a,i_s)            ./x_bind))-mua)./mua; % empirical
        set_of_dpfs_err(i_a,i_s,6,1:length(x_bind)) = ((y_bind./x_bind./(p_anc(i_a,i_s)            +q_anc(i_a,i_s)            ./x_bind))-mua)./mua; % analytical
        set_of_dpfs_err(i_a,i_s,7,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_slop_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % slope DPF
        set_of_dpfs_err(i_a,i_s,8,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_adst_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % semi-inf DPF

        clearvars mdl p blnFit blnShow fun_x fun_y fun_s i_fig index_in TheCode TheOutFun x_temp y_temp s_temp x_bind y_bind s_bind
        clearvars t_db the_filename mua mus ind_diff ind_trns
    end
end
clearvars i_a i_s

% calc p & q
[params] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_mus(2:end)), log(set_of_idst_DPF(2:end,2:end,1)));
disp(num2str(['C = ',num2str(exp(params(1))),', a = ',num2str(params(2)),', b = ',num2str(params(3))]))
[params] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_mus(2:end)), log(set_of_idst_DPF(2:end,2:end,2)));
disp(num2str(['C = ',num2str(exp(params(1))),', a = ',num2str(params(2)),', b = ',num2str(params(3))]))


figure(1)
i_a = 4; i_s = 4;
subplot(1,2,1)
plot(d_cntrs,squeeze(set_of_cnst_DPF(i_a,i_s,1)).*ones(size(d_cntrs)),      'DisplayName','constant','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_adst_DPF(i_a,i_s,:)),                           'DisplayName','semi-inf','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_savg_DPF(i_a,i_s,:)),                           'DisplayName','mean',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_slop_DPF(i_a,i_s,:)),                           'DisplayName','slope',   'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_true_DPF(i_a,i_s,:)),                           'DisplayName','true',    'LineWidth',2), hold on
plot(d_cntrs,set_of_idst_DPF(i_a,i_s,1)+set_of_idst_DPF(i_a,i_s,2)./d_cntrs,'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 6]),  set(gca,'xtick',0:1:6)
% title(['Diff. DPF definitions for \mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))])
range_1 = [0 +8]; range_2 = [18 40]; the_step = 10;
ylabel('DPF (a.u.)'), ylim([range_1(1) range_2(2)]), set(gca,'ytick',[0 8 20 30 40])
hold off, set(gca,'fontsize',24)
hold off, legend('show','Location','northeast'), set(gca,'fontsize',24), % axis square
breakyaxis([range_1(2) range_2(1)],0.05,0.05); 
clearvars range_1 range_2 the_step
subplot(1,2,2)
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,1,:)).*100,                    'DisplayName','constant','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,8,:)).*100,                    'DisplayName','semi-inf','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,2,:)).*100,                    'DisplayName','mean',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,7,:)).*100,                    'DisplayName','slope',   'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,3,:)).*100,                    'DisplayName','true',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,4,:)).*100,                    'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 6]),   set(gca,'xtick',0:1:6)
range_1 = [-10 +10]; range_2 = [200 350]; the_step = 10;
ylabel('estimation error (%)'), ylim([range_1(1) range_2(2)]), set(gca,'ytick',range_1(1):the_step:range_2(2))
hold off, set(gca,'fontsize',24)
breakyaxis([range_1(2) range_2(1)],0.05,0.05); 

figure(2)
idx = 1<=d_cntrs&d_cntrs<=8;
subplot(2,2,1)
icol = -1:3:50; icol = 0:1:+10; [set_of_lines,my_colormap] = make_colormap(min(icol),max(icol),length(icol)); % my_colormap = 'jet';
contourf(set_of_mua,set_of_mus,abs(squeeze(nanmean(set_of_dpfs_err(:,:,1,idx)*100,4)).'),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
% h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):5:max(icol)]), title('constant DPF'); clim([min(icol) max(icol)])
% set(gca,'fontsize',16), axis square, axis([min(set_of_mua(2)) max(set_of_mua) min(set_of_mus(2)) max(set_of_mus)]); colormap(gca,my_colormap)
subplot(2,2,3)
contourf(set_of_mua,set_of_mus,abs(squeeze(nanmean(set_of_dpfs_err(:,:,4,idx)*100,4)).'),set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):5:max(icol)]), title('inverse DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua(2)) max(set_of_mua) min(set_of_mus(2)) max(set_of_mus)]); colormap(gca,my_colormap)
subplot(2,2,4)
contourf(set_of_mua,set_of_mus,abs(squeeze(nanmean(set_of_dpfs_err(:,:,5,idx)*100,4)).'),set_of_lines,'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',[min(icol):5:max(icol)]), title('empirical DPF'); clim([min(icol) max(icol)])
set(gca,'fontsize',16), axis square, axis([min(set_of_mua(2)) max(set_of_mua) min(set_of_mus(2)) max(set_of_mus)]); colormap(gca,my_colormap)
colormap jet
clearvars h idx icol

figure(3)
subplot(2,1,1), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,1)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,1,1), mesh(set_of_mua,set_of_mus,p_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
% subplot(2,1,1), mesh(set_of_mua,set_of_mus,p_anc.',                          'FaceColor','none','EdgeColor','g','DisplayName','analytical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('p (a.u.)'), set(gca,'ztick',0:25:75), title('p'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
subplot(2,1,2), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,2)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,1,2), mesh(set_of_mua,set_of_mus,q_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
% subplot(2,1,2), mesh(set_of_mua,set_of_mus,q_anc.',                          'FaceColor','none','EdgeColor','g','DisplayName','analytical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('b (a.u.)'), set(gca,'ztick',0:25:75), title('q'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
disp(['r^2 for p = ',num2str(rsquared(set_of_idst_DPF(:,:,1).',p_est.'),'%.2f')])
disp(['r^2 for q = ',num2str(rsquared(set_of_idst_DPF(:,:,2).',q_est.'),'%.2f')])
clearvars p_est q_est



figure(4)
i_a = 7;
i_s = 3;

idx_involved = ~isnan(squeeze(set_of________I(i_a,i_s,:)));
plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1)*d_cntrs(idx_involved)+set_of_idst_DPF(i_a,i_s,2)).*set_of_mua(i_a),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s)),': data'],'LineStyle','-'), hold on
plot(d_cntrs(idx_involved),(set_of_cnst_DPF(i_a,i_s,1)*d_cntrs(idx_involved)                           ).*set_of_mua(i_a),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s)),': data'],'LineStyle','-'), hold on
plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),'Color','k','DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s)),': data']), hold on
xlabel('d (cm)'),     xlim([0 6 ]), set(gca,'xtick',0:1:6)
ylabel('OD (a.u.)'),  ylim([0 25]), set(gca,'ytick',0:5:25)
title(['OD vs. d for \mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))])
hold off, axis square, hold off, legend('show','Location','northeast'), set(gca,'fontsize',16), grid on

% inset_ax = axes('Position',[0.25 0.35 0.3 0.3]); box on
% idx_involved = ~isnan(squeeze(set_of________I(i_a,i_s,:)))&d_cntrs.'<=0.25;
% plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1)*d_cntrs(idx_involved)+set_of_idst_DPF(i_a,i_s,2)).*set_of_mua(i_a),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s)),': data'],'LineStyle','-'), hold on
% plot(d_cntrs(idx_involved),(set_of_cnst_DPF(i_a,i_s,1)*d_cntrs(idx_involved)                           ).*set_of_mua(i_a),'DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s)),': data'],'LineStyle','-'), hold on
% plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),'Color','k','DisplayName',['\mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s)),': data']), hold on
% xlim([0 0.25]), set(gca,'xtick',0:05:0.25)
% ylim([0 5]), set(gca,'ytick',0:1:5)
% hold off
clearvars i_a i_s
end

function [out] = get_color(mua_norm,mus_norm)
col_mua_0_mus_0 = [1 1 0];
col_mua_0_mus_1 = [1 0 0];
col_mua_1_mus_0 = [0 0 1];
col_mua_1_mus_1 = [1 0.0 1];
% out = [mua_norm 0 mus_norm];
out = ...
    (1-mua_norm)*(1-mus_norm)*col_mua_0_mus_0 + ...
    (  mua_norm)*(1-mus_norm)*col_mua_1_mus_0 + ...
    (1-mua_norm)*(  mus_norm)*col_mua_0_mus_1 + ...
    (  mua_norm)*(  mus_norm)*col_mua_1_mus_1;
end
function [set_of_lines,my_colormap] = make_colormap(l_start,l_end,N_lines)
set_of_lines = linspace(l_start,l_end,2*N_lines+1);
my_colormap = [1, 1, 1];
if     l_end<=0
    for idx = 1:2*N_lines
        my_colormap = [...
            [1-idx/2/N_lines, 1-idx/2/N_lines, 1];...
            my_colormap];
    end
elseif l_start>=0
    for idx = 1:2*N_lines
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/2/N_lines, 1-idx/2/N_lines]];
    end
else
    Np = sum(set_of_lines>0);
    Nn = sum(set_of_lines<0);
    for idx = 1:Np
        my_colormap = [...
            my_colormap;...
            [1, 1-idx/Np, 1-idx/Np]];
    end
    for idx = 1:Nn
        my_colormap = [...
            [1-idx/Nn, 1-idx/Nn, 1];...
            my_colormap];
    end
end
end
function [r2] = rsquared(y,yhat)
y    = y(:);
yhat = yhat(:);
idx = ~isnan(y)&~isnan(yhat)&~isinf(y)&~isinf(yhat);
y    = y(idx);
yhat = yhat(idx);
SS_res = sum((y - yhat).^2);    % Residual sum of squares
SS_tot = sum((y - mean(y)).^2); % Total sum of squares
r2 = 1 - (SS_res/SS_tot);
end
function [params] = fitPlaneModel(X, Y, Z)
% Fits Z = C + a*X + b*Y to gridded data

% Create meshgrid (M×N), transpose to match Z
[XX, YY] = meshgrid(X, Y);  % XX, YY: N×M
XX = XX.'; YY = YY.';         % Now M×N

% Flatten everything to column vectors
A = [ones(numel(XX), 1), XX(:), YY(:)];  % Design matrix: [1 X Y]
b = Z(:);                                % Observed values

% Solve for [C; a; b]
params = A \ b;
end



