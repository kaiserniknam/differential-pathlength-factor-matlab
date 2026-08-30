function [] = Fig_02_03_04_3
% Fig-02/03/04: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - to plot p and q for mu_a's and mu_s's
% based on version #2, corrected intensity definition

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

n = 1.4; g = 0.95;
set_of_mua = (.00:.05:0.5);
set_of_mus = (000:050:500);
% N_bins = 450;
N_bins = 120;
d_edges = linspace(0,sqrt((Lx/2).^2+(Ly/2).^2),N_bins+1); clearvars Lx Ly Lz
dlta_d = mean(diff(d_edges));
d_cntrs = 1/2*(d_edges(1:end-1)+d_edges(2:end-0));
set_of_cnst_DPF = nan(length(set_of_mua),length(set_of_mus),1);
set_of_idst_DPF = nan(length(set_of_mua),length(set_of_mus),3);
set_of_adst_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_savg_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_slop_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_true_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of________I = nan(length(set_of_mua),length(set_of_mus),N_bins);
set_of_dpfs_err = nan(length(set_of_mua),length(set_of_mus),8,N_bins);
clearvars N_bins
A_est = nan(size(squeeze(set_of_idst_DPF(:,:,1))));
B_est = nan(size(squeeze(set_of_idst_DPF(:,:,2))));
C_est = nan(size(squeeze(set_of_idst_DPF(:,:,3))));
clc

for i_a = 1:length(set_of_mua)
    for i_s = 2:length(set_of_mus)
        A_coeff = [0.32665, -0.42898, 0.433400];
        B_coeff = [1.66450, -1.14890, 0.018758];
        C_coeff = [0.74579, -1.14810, 0.175400];

        A_est(i_a,i_s) = A_coeff(1).*((set_of_mua(i_a)).^(A_coeff(2))).*((set_of_mus(i_s)).^(A_coeff(3)));
        B_est(i_a,i_s) = B_coeff(1).*((set_of_mua(i_a)).^(B_coeff(2))).*((set_of_mus(i_s)).^(B_coeff(3)));
        C_est(i_a,i_s) = C_coeff(1).*((set_of_mua(i_a)).^(C_coeff(2))).*((set_of_mus(i_s)).^(C_coeff(3)));
        clearvars A_coeff B_coeff C_coeff
    end
end
clearvars i_a i_s A_coeff B_coeff C_coeff

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
        % disp(num2str([length(u) sum(u==u_unique(i_fx)) sum(u==u_unique(i_fx))/length(u)*100]))
        t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));
        clearvars u u_unique freq i_f i_fx        
    
        % 1-D sorting
        [~,~,ind_diff] = histcounts(sqrt(t_db.x(t_db.c==0).^2+t_db.y(t_db.c==0).^2+t_db.z(t_db.c==0).^2),d_edges); % d_diffuse bins
        clearvars d_trns_edges d_diff_edges

        % s vs. d
        fun_x = @mean; fun_y = @mean; index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
        y_temp =      t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan);
        set_of_savg_DPF(i_a,i_s,1:length(y_bind)) = y_bind./x_bind;
        clearvars fun_x fun_y index_in TheCode x_temp y_temp x_bind y_bind
        
        % slope DPF
        fun_x = @mean; fun_s = @(x)(sum(x.*exp(-mua.*x))./sum(exp(-mua.*x))); index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2); 
        s_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        s_bind = accumarray(index_in,s_temp,[],fun_s,nan); s_bind = s_bind./x_bind;                      
        set_of_slop_DPF(i_a,i_s,1:length(s_bind)) = s_bind;
        clearvars fun_x fun_s index_in TheCode x_temp s_temp x_bind s_bind

        % I vs. d
        fun_x = @mean; fun_y = @sum; fun_s = @(x)(sum(exp(-mua.*x))); TheOutFun = @(x)(-log10(x)); index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
        y_temp = t_db.w(t_db.c==TheCode);
        % s_temp = t_db.s(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = TheOutFun( (y_bind./(2.*pi.*x_bind.*dlta_d)) ./ (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)) );
        d_bind = y_bind./x_bind./mua;
        clearvars fun_x fun_y fun_s index_in TheCode x_temp y_temp s_temp
        
        % Least-squares fit
        Zstar = y_bind(~isnan(x_bind)) + log10(x_bind(~isnan(x_bind))); % equals Z - c, but intercept absorbs c anyway
        T = table(Zstar, x_bind(~isnan(x_bind)), log10(x_bind(~isnan(x_bind))), 'VariableNames', {'Zstar','X','logX'});
        try
            mdl = fitlm(T, 'Zstar ~ 1 + X + logX');  % estimates B, A, C
            clearvars Zstar T
            B = mdl.Coefficients.Estimate(1);
            A = mdl.Coefficients.Estimate(2);
            C = mdl.Coefficients.Estimate(3);
            Z_est = A*x_bind + B + C*log10(x_bind);  % (this estimates Zstar; true Z = Zstar + c)
            % plot(x_bind,y_bind,x_bind,Z_est,'LineWidth',2)
            % title(['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus)])
            % set(gca,'fontsize',24)
            % axis ([0 9 0 20]), axis square, pause(0.01)
            p = mdl.Coefficients.pValue(2);
            if p<=0.05
                set_of_idst_DPF(i_a,i_s,1) = A./mua;
                set_of_idst_DPF(i_a,i_s,2) = B./mua;
                set_of_idst_DPF(i_a,i_s,3) = C./mua;
            else
                % set_of_cnst_DPF(i_a,i_s,1) = nan;
                set_of_idst_DPF(i_a,i_s,1) = nan;
                set_of_idst_DPF(i_a,i_s,2) = nan;
                set_of_idst_DPF(i_a,i_s,3) = nan;
            end
            clearvars A B C Z_est
        catch
            set_of_idst_DPF(i_a,i_s,1) = nan;
            set_of_idst_DPF(i_a,i_s,2) = nan;
            set_of_idst_DPF(i_a,i_s,3) = nan;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set_of________I(i_a,i_s,1:length(y_bind)) = y_bind;
        set_of_true_DPF(i_a,i_s,1:length(d_bind)) = d_bind;
        clearvars mdl p blnFit blnShow fun_x fun_y fun_s i_fig index_in TheCode TheOutFun x_temp y_temp s_temp d_bind
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        musp = mus*(1-g);
        num = sqrt(3*musp).*(x_bind.*sqrt(3*mua*musp)+0);
        den = 2*sqrt(mua) .*(x_bind.*sqrt(3*mua*musp)+1);
        dpf_tmp = num./den;

        set_of_cnst_DPF(i_a,i_s,1) = (sqrt(3*musp))./(2*sqrt(mua));
        set_of_adst_DPF(i_a,i_s,1:length(dpf_tmp)) = dpf_tmp;
        clearvars musp num den dpf_tmp
        
        % error
        set_of_dpfs_err(i_a,i_s,1,1:length(x_bind)) = ((y_bind./x_bind./set_of_cnst_DPF(i_a,i_s,1))-mua)./mua; % constant
        set_of_dpfs_err(i_a,i_s,2,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_savg_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % <s>/d
        set_of_dpfs_err(i_a,i_s,3,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_true_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % true
        set_of_dpfs_err(i_a,i_s,4,1:length(x_bind)) = ((y_bind./x_bind./( ...
            set_of_idst_DPF(i_a,i_s,1) + ...
            set_of_idst_DPF(i_a,i_s,2)./x_bind + ...
            set_of_idst_DPF(i_a,i_s,3).*log10(x_bind)./x_bind))-mua)./mua; % inverse-distance
        set_of_dpfs_err(i_a,i_s,5,1:length(x_bind)) = ((y_bind./x_bind./( ...
            A_est(i_a,i_s) + ...
            B_est(i_a,i_s)./x_bind + ...
            C_est(i_a,i_s).*log10(x_bind)./x_bind))-mua)./mua; % empirical
        % set_of_dpfs_err(i_a,i_s,6,1:length(x_bind)) = ((y_bind./x_bind./(p_anc(i_a,i_s)            +q_anc(i_a,i_s)            ./x_bind))-mua)./mua; % analytical
        set_of_dpfs_err(i_a,i_s,7,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_slop_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % slope DPF
        set_of_dpfs_err(i_a,i_s,8,1:length(x_bind)) = ((y_bind./x_bind./squeeze(set_of_adst_DPF(i_a,i_s,1:length(x_bind))))-mua)./mua; % semi-inf DPF
        clearvars x_bind y_bind
        clearvars t_db the_filename mua mus ind_diff ind_trns
    end
end
clearvars i_a i_s

% calc A & B & C
clc
[params,err] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_mus(2:end)), log(set_of_idst_DPF(2:end,2:end,1)));
disp(num2str(['C = ',num2str(exp(params(1)),'%.2f'),', a = ',num2str(params(2),'%.2f'),', b = ',num2str(params(3),'%.2f'),', @ r2 =  ',num2str(err)])), clearvars err params
[params,err] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_mus(2:end)), log(set_of_idst_DPF(2:end,2:end,2)));
disp(num2str(['C = ',num2str(exp(params(1)),'%.2f'),', a = ',num2str(params(2),'%.2f'),', b = ',num2str(params(3),'%.2f'),', @ r2 =  ',num2str(err)])), clearvars err params
[params,err] = fitPlaneModel(log(set_of_mua(2:end)), log(set_of_mus(2:end)), log(set_of_idst_DPF(2:end,2:end,3)));
disp(num2str(['C = ',num2str(exp(params(1)),'%.2f'),', a = ',num2str(params(2),'%.2f'),', b = ',num2str(params(3),'%.2f'),', @ r2 =  ',num2str(err)])), clearvars err params


figure(1)
i_a = 4; i_s = 2;
subplot(1,2,1)
plot(d_cntrs,squeeze(set_of_cnst_DPF(i_a,i_s,1)).*ones(size(d_cntrs)),                                                          'DisplayName','constant','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_adst_DPF(i_a,i_s,:)),                                                                               'DisplayName','semi-inf','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_savg_DPF(i_a,i_s,:)),                                                                               'DisplayName','mean',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_slop_DPF(i_a,i_s,:)),                                                                               'DisplayName','slope',   'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_true_DPF(i_a,i_s,:)),                                                                               'DisplayName','true',    'LineWidth',2), hold on
plot(d_cntrs,set_of_idst_DPF(i_a,i_s,1)+set_of_idst_DPF(i_a,i_s,2)./d_cntrs+set_of_idst_DPF(i_a,i_s,3).*log10(d_cntrs)./d_cntrs,'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 08]), set(gca,'xtick',0:2:8)
ylabel('DPF (a.u.)'), ylim([0 40]), set(gca,'ytick',0:10:40)
title(['Diff. DPF definitions for \mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))])
hold off, set(gca,'fontsize',24)
legend('show','Location','northeast'), set(gca,'fontsize',24), % axis square
subplot(1,2,2)
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,1,:)).*100,                    'DisplayName','constant','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,8,:)).*100,                    'DisplayName','semi-inf','LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,2,:)).*100,                    'DisplayName','mean',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,7,:)).*100,                    'DisplayName','slope',   'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,3,:)).*100,                    'DisplayName','true',    'LineWidth',2), hold on
plot(d_cntrs,squeeze(set_of_dpfs_err(i_a,i_s,4,:)).*100,                    'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 08]), set(gca,'xtick',0:2:8)
ylabel('estimation error (%)'), ylim([-20 500]), set(gca,'ytick',[-10 0 100 300 500])
title(['Diff. DPF estimation error for \mu_a = ',num2str(set_of_mua(i_a)),', \mu_s = ',num2str(set_of_mus(i_s))])
hold off, set(gca,'fontsize',24)

figure(2)
idx = 1<=d_cntrs&d_cntrs<=8;

subplot(2,2,1)
icol = 0:80:+240;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(abs(set_of_dpfs_err(:,:,1,idx)*100),4)).',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('contant DPF');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

subplot(2,2,2)
icol = 0:80:+240;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(abs(set_of_dpfs_err(:,:,1,idx)*100),4)).',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('mean-pathlength DPF');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

subplot(2,2,3)
icol = 0:5:15;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(abs(set_of_dpfs_err(:,:,4,idx)*100),4)).',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('inv. dist DPF');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

subplot(2,2,4)
icol = 0:5:15;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(abs(set_of_dpfs_err(:,:,5,idx)*100),4)).',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('empirical DPF');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h idx icol

figure(3)
subplot(2,2,1), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,1)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,2,1), mesh(set_of_mua,set_of_mus,A_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('A (cm)'), set(gca,'ztick',0:6:18), title('A'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 18])
subplot(2,2,2), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,2)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,2,2), mesh(set_of_mua,set_of_mus,B_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('B (a.u.)'), set(gca,'ztick',0:25:75), title('B'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
subplot(2,2,3), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,3)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,2,3), mesh(set_of_mua,set_of_mus,B_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('C (a.u.)'), set(gca,'ztick',0:25:75), title('C'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 75])
clearvars A_est B_est C_est



figure(4)
set_of_ia = [2,4];
set_of_is = [2,7];
for i_a = set_of_ia
    for i_s = set_of_is
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        TheColor = [1.0,(i_a-min(set_of_ia))./(max(set_of_ia)-min(set_of_ia)),(i_s-min(set_of_is))./(max(set_of_is)-min(set_of_is))];
        if TheColor==[1,1,1], TheColor=[0,0,0]; end
        idx_involved = ~isnan(squeeze(set_of________I(i_a,i_s,:)));
        subplot(2,2,1)
        plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),                                                                                                  'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log10(d_cntrs(idx_involved))).*mua,'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),    xlim([0 08]), set(gca,'xtick',0:2:8)
        ylabel('OD (a.u.)'), ylim([0 10]), set(gca,'ytick',0:2.5:10)
        title('Diff. OD')
        set(gca,'fontsize',16), axis square
        % legend('show','Location','southeast'), % axis square

        subplot(2,2,2)
        plot(d_cntrs(idx_involved),squeeze(set_of_true_DPF(i_a,i_s,idx_involved)),                                                                                                                  'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),set_of_idst_DPF(i_a,i_s,1) + set_of_idst_DPF(i_a,i_s,2)./d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,3).*log10(d_cntrs(idx_involved))./d_cntrs(idx_involved),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),     xlim([0 008]), set(gca,'xtick',0:2:8)
        ylabel('DPF (a.u.)'), ylim([0 120]), set(gca,'ytick',0:30:120)
        title('Diff. DPF definitions')
        set(gca,'fontsize',16), axis square
        legend('show','Location','northeast'), % axis square
        
        subplot(2,2,3)
        plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),                                                                                                  'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log10(d_cntrs(idx_involved))).*mua,'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2)                                                           ).*mua,'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': lin.'],'LineStyle',':' ,'LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),    xlim([0 08]), set(gca,'xtick',0:2:8)
        ylabel('OD (a.u.)'), ylim([0 10]), set(gca,'ytick',0:2.5:10)
        title('Diff. OD')
        set(gca,'fontsize',16), axis square
        legend('show','Location','southeast'), % axis square

        subplot(2,2,4)
        plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),                                                                                                  'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log10(d_cntrs(idx_involved))).*mua,'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(                                                    set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log10(d_cntrs(idx_involved))).*mua,'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': log'], 'LineStyle',':','LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),    xlim([0 08]), set(gca,'xtick',0:2:8)
        ylabel('OD (a.u.)'), ylim([0 10]), set(gca,'ytick',0:2.5:10)
        title('Diff. OD')
        set(gca,'fontsize',16), axis square
        legend('show','Location','southeast'), % axis square

        clearvars mua mus TheColor idx_involved
    end
end
clearvars i_a i_s
end

function [out] = get_color(mua_norm,mus_norm)
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
function [params,err] = fitPlaneModel(X, Y, Z)
% Fits Z = C + a*X + b*Y to gridded data

% Create meshgrid (M×N), transpose to match Z
[XX, YY] = meshgrid(X, Y);  % XX, YY: N×M
XX = XX.'; YY = YY.';         % Now M×N

% Flatten everything to column vectors
A = [ones(numel(XX), 1), XX(:), YY(:)];  % Design matrix: [1 X Y]
b = Z(:);                                % Observed values

% Solve for [C; a; b]
params = A \ b;
err = rsquared(b,A*params);
end



