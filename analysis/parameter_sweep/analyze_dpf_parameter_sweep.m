function [] = analyze_dpf_parameter_sweep
%ANALYZE_DPF_PARAMETER_SWEEP Compare DPF models over optical properties.
%   Loads Photon_33 Monte Carlo histories, forms area-normalized diffuse
%   reflectance and optical density, evaluates six DPF definitions, and fits
%   the three-parameter inverse-distance model used in the paper.

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

cfg = project_config();

% optical & geometry properties
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

n = 1.4; g = 0.95;
set_of_mua = (.00:.05:0.5);
set_of_mus = (000:050:500);
dlta_d = .17;
d_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2); N_bins = length(d_edges)-1; clearvars Lx Ly Lz
d_cntrs = 1/2*(d_edges(1:end-1)+d_edges(2:end-0));
set_of_cnst_DPF = nan(length(set_of_mua),length(set_of_mus),1);       % constant dpf
set_of_idst_DPF = nan(length(set_of_mua),length(set_of_mus),3);       % inverse-distance dpf
set_of_adst_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);  % semi-infinite dpf
set_of_savg_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);  % mean pathlength dpf
set_of_slop_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);  % slope-based dpf
set_of_true_DPF = nan(length(set_of_mua),length(set_of_mus),N_bins);  % true dpf
set_of________I = nan(length(set_of_mua),length(set_of_mus),N_bins);  % OD values

set_of_dpfs_est = nan(length(set_of_mua),length(set_of_mus),8,N_bins);% est. mu_a values
set_of_dpfs_err = nan(length(set_of_mua),length(set_of_mus),8,N_bins);% rel. error values

set_of________G = nan(length(set_of_mua),length(set_of_mus),1);       % G_values
set_of_od_rsqrd = nan(length(set_of_mua), length(set_of_mus));        % R-squared of OD and its approximation
clearvars N_bins
A_est = nan(size(squeeze(set_of_idst_DPF(:,:,1))));
B_est = nan(size(squeeze(set_of_idst_DPF(:,:,2))));
C_est = nan(size(squeeze(set_of_idst_DPF(:,:,3))));
clc

for i_a = 1:length(set_of_mua)
    for i_s = 2:length(set_of_mus)
        % % original
        % A_coeff = [0.87, -0.51, 0.41];
        % B_coeff = [0.02, -1.40, 0.72];
        % C_coeff = [0.03, -1.25, 0.55];

        % based on the output of this code!
        A_coeff = [0.85, -0.52, 0.41];
        B_coeff = [0.03, -1.35, 0.65];
        C_coeff = [0.03, -1.24, 0.54];

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
        the_filename = fullfile(cfg.simulationData, ...
            ['Photon_33_mua_', sprintf('%.2f',mua), ...
             '_mus_', sprintf('%.2f',mus), '.mat']);
        t_db = load(the_filename); clearvars the_filename
        fprintf('mua = %.2f, mus = %.0f, ', mua, mus);
        
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
        fprintf(': %d -> %d meaans %.2f\n', length(u), sum(u == u_unique(i_fx)), sum(u == u_unique(i_fx))/length(u)*100);
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
        fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
        y_temp = t_db.w(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = TheOutFun((y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));
        if isempty(y_bind), G = nan; else, G = y_bind(1); end % G
        set_of________G(i_a,i_s) = G;
        d_bind = (y_bind-G)./x_bind./mua;
        clearvars fun_x fun_y fun_s index_in TheCode x_temp y_temp s_temp
        
        % Least-squares fit
        x = x_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);
        y = y_bind(~isnan(x_bind) & ~isnan(y_bind) & x_bind > 0);        
        Y = (y - G) ./ x;        % response
        invX = 1 ./ x;           % predictor 1
        logX_over_X = log(x) ./ x; % predictor 2
        T = table(Y, invX, logX_over_X, 'VariableNames', {'Y','invX','logX_over_X'});
        try
            mdl = fitlm(T, 'Y ~ 1 + invX + logX_over_X');  % A, B, C
            clearvars T x y Y invX logX_over_X   
            A = mdl.Coefficients.Estimate(1);
            B = mdl.Coefficients.Estimate(2);
            C = mdl.Coefficients.Estimate(3);
            % Z_est = A + B./x_bind + C*log(x_bind)./x_bind;  % (this estimates Zstar; true Z = Zstar + c)
            % plot(x_bind,d_bind.*mua,x_bind,Z_est,'LineWidth',2)
            % title(['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus)])
            % set(gca,'fontsize',24)
            % axis ([0 9 0 20]), axis square, pause(0.01)
            p = mdl.Coefficients.pValue(2);
            set_of_od_rsqrd(i_a,i_s) = mdl.Rsquared.Ordinary;
            if p<=0.05
                set_of_idst_DPF(i_a,i_s,1) = A./mua;
                set_of_idst_DPF(i_a,i_s,2) = B./mua;
                set_of_idst_DPF(i_a,i_s,3) = C./mua;
            else
                set_of_idst_DPF(i_a,i_s,1) = nan;
                set_of_idst_DPF(i_a,i_s,2) = nan;
                set_of_idst_DPF(i_a,i_s,3) = nan;
            end
            clearvars A B C Z_est
        catch e
            disp(e)
            set_of_idst_DPF(i_a,i_s,1) = nan;
            set_of_idst_DPF(i_a,i_s,2) = nan;
            set_of_idst_DPF(i_a,i_s,3) = nan;
            clearvars T x y Y invX logX_over_X               
        end
        clearvars e
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
        set_of_dpfs_est(i_a,i_s,1,1:length(x_bind)) =  ((y_bind-G)./x_bind./set_of_cnst_DPF(i_a,i_s,1));      % constant
        set_of_dpfs_err(i_a,i_s,1,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,1,1:length(x_bind))-mua)./mua; % constant

        set_of_dpfs_est(i_a,i_s,2,1:length(x_bind)) = ((y_bind-G)./x_bind./squeeze(set_of_savg_DPF(i_a,i_s,1:length(x_bind)))); % <s>/d
        set_of_dpfs_err(i_a,i_s,2,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,2,1:length(x_bind))-mua)./mua;                   % <s>/d

        set_of_dpfs_est(i_a,i_s,3,1:length(x_bind)) = ((y_bind-G)./x_bind./squeeze(set_of_true_DPF(i_a,i_s,1:length(x_bind)))); % true
        set_of_dpfs_err(i_a,i_s,3,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,3,1:length(x_bind))-mua)./mua;                   % true

        set_of_dpfs_est(i_a,i_s,4,1:length(x_bind)) = ((y_bind-G)./x_bind./( ...
            set_of_idst_DPF(i_a,i_s,1) + ...
            set_of_idst_DPF(i_a,i_s,2)./x_bind + ...
            set_of_idst_DPF(i_a,i_s,3).*log(x_bind)./x_bind));                                                % inverse-distance
        set_of_dpfs_err(i_a,i_s,4,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,4,1:length(x_bind))-mua)./mua; % inverse-distance


        set_of_dpfs_est(i_a,i_s,5,1:length(x_bind)) = ((y_bind-G)./x_bind./( ...
            A_est(i_a,i_s) + ...
            B_est(i_a,i_s)./x_bind + ...
            C_est(i_a,i_s).*log(x_bind)./x_bind));                                                            % empirical
        set_of_dpfs_err(i_a,i_s,5,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,5,1:length(x_bind))-mua)./mua; % empirical

        % set_of_dpfs_est(i_a,i_s,6,1:length(x_bind)) = (y_bind./x_bind./(p_anc(i_a,i_s)+q_anc(i_a,i_s)./x_bind)); % analytical
        % set_of_dpfs_err(i_a,i_s,6,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,6,1:length(x_bind))-mua)./mua;    % analytical

        set_of_dpfs_est(i_a,i_s,7,1:length(x_bind)) = ((y_bind-G)./x_bind./squeeze(set_of_slop_DPF(i_a,i_s,1:length(x_bind)))); % slope DPF
        set_of_dpfs_err(i_a,i_s,7,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,7,1:length(x_bind))-mua)./mua;                   % slope DPF

        set_of_dpfs_est(i_a,i_s,8,1:length(x_bind)) = ((y_bind-G)./x_bind./squeeze(set_of_adst_DPF(i_a,i_s,1:length(x_bind)))); % semi-inf DPF
        set_of_dpfs_err(i_a,i_s,8,1:length(x_bind)) = (set_of_dpfs_est(i_a,i_s,8,1:length(x_bind))-mua)./mua;                   % semi-inf DPF

        clearvars x_bind y_bind
        clearvars t_db the_filename mua mus ind_diff ind_trns G
    end
end
clearvars i_a i_s
clearvars d_edges

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
plot(d_cntrs,set_of_idst_DPF(i_a,i_s,1)+set_of_idst_DPF(i_a,i_s,2)./d_cntrs+set_of_idst_DPF(i_a,i_s,3).*log(d_cntrs)./d_cntrs,'DisplayName','inverse', 'LineWidth',2), hold on
xlabel('d (cm)'),     xlim([0 08]), set(gca,'xtick',0:2:8)
ylabel('DPF (unitless)'), ylim([0 20]), set(gca,'ytick',0:5:20)
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
clearvars i_a i_s

figure(2)
idx = 1<=d_cntrs&d_cntrs<=8;

subplot(2,2,1)
icol = 200:10:+240;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(abs(set_of_dpfs_err(:,:,1,idx)*100),4)).',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'relative (%)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('constant DPF');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

subplot(2,2,2)
icol = 200:10:+240;
contourf(set_of_mua,set_of_mus,squeeze(nanmean(abs(set_of_dpfs_err(:,:,2,idx)*100),4)).',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
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
zlabel('A (cm)'), set(gca,'ztick',0:20:60), title('A'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 60])
clearvars A_est
subplot(2,2,2), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,2)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,2,2), mesh(set_of_mua,set_of_mus,B_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('B (unitless)'), set(gca,'ztick',0:35:105), title('B'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 105])
clearvars B_est
subplot(2,2,3), mesh(set_of_mua,set_of_mus,squeeze(set_of_idst_DPF(:,:,3)).','FaceColor','none','EdgeColor','b','DisplayName','fitted'), hold on
subplot(2,2,3), mesh(set_of_mua,set_of_mus,C_est.',                          'FaceColor','none','EdgeColor','r','DisplayName','empirical'), hold on
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[min(set_of_mua) mean(set_of_mua) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[min(set_of_mus) mean(set_of_mus) max(set_of_mus)]) 
zlabel('C (unitless)'), set(gca,'ztick',0:12:36), title('C'), hold off; legend
set(gca,'fontsize',16), axis square, axis tight, colormap jet; view([145 22.5])
axis([min(set_of_mua) max(set_of_mua) min(set_of_mus) max(set_of_mus) 0 36])
clearvars C_est



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
        plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),                                                                                                                         'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,3).*log(d_cntrs(idx_involved))).*mua+set_of________G(i_a,i_s),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),    xlim([0 08]), set(gca,'xtick',0:2:8)
        ylabel('OD (unitless)'), ylim([0 24]), set(gca,'ytick',0:6:24)
        title('Diff. OD')
        set(gca,'fontsize',16), axis square
        % legend('show','Location','southeast'), % axis square

        subplot(2,2,2)
        plot(d_cntrs(idx_involved),squeeze(set_of_true_DPF(i_a,i_s,idx_involved)),                                                                                                                'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),set_of_idst_DPF(i_a,i_s,1) + set_of_idst_DPF(i_a,i_s,2)./d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,3).*log(d_cntrs(idx_involved))./d_cntrs(idx_involved),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),     xlim([0 008]), set(gca,'xtick',0:2:8)
        ylabel('DPF (unitless)'), ylim([0 140]), set(gca,'ytick',0:35:140)
        title('Diff. DPF definitions')
        set(gca,'fontsize',16), axis square
        legend('show','Location','northeast'), % axis square
        
        subplot(2,2,3)
        plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),                                                                                                                         'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log(d_cntrs(idx_involved))).*mua+set_of________G(i_a,i_s),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2)                                                         ).*mua+set_of________G(i_a,i_s),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': lin.'],'LineStyle',':' ,'LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),    xlim([0 08]), set(gca,'xtick',0:2:8)
        ylabel('OD (unitless)'), ylim([0 24]), set(gca,'ytick',0:6:24)
        title('Diff. OD')
        set(gca,'fontsize',16), axis square
        legend('show','Location','southeast'), % axis square

        subplot(2,2,4)
        plot(d_cntrs(idx_involved),squeeze(set_of________I(i_a,i_s,idx_involved)),                                                                                                                         'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': true'],'LineStyle','-' ,'LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(set_of_idst_DPF(i_a,i_s,1).*d_cntrs(idx_involved) + set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log(d_cntrs(idx_involved))).*mua+set_of________G(i_a,i_s),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': est.'],'LineStyle','-.','LineWidth',2,'Color',TheColor), hold on
        plot(d_cntrs(idx_involved),(                                                    set_of_idst_DPF(i_a,i_s,2) + set_of_idst_DPF(i_a,i_s,3).*log(d_cntrs(idx_involved))).*mua+set_of________G(i_a,i_s),'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),': log'], 'LineStyle',':','LineWidth',2,'Color',TheColor), hold on
        xlabel('d (cm)'),    xlim([0 08]), set(gca,'xtick',0:2:8)
        ylabel('OD (unitless)'), ylim([0 24]), set(gca,'ytick',0:6:24)
        title('Diff. OD')
        set(gca,'fontsize',16), axis square
        legend('show','Location','southeast'), % axis square

        clearvars mua mus TheColor idx_involved
    end
end
clearvars i_a i_s
clearvars set_of_ia
clearvars set_of_is



figure(5)
icol = 0.85:0.05:1.0;
contourf(set_of_mua,set_of_mus,set_of_od_rsqrd.',linspace(min(icol),max(icol),25),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'r-squared'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('fitting r-squared');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
disp(['R2: mean ± std ~ ',num2str(nanmean(set_of_od_rsqrd(:)),'%.2f'),' ± ',num2str(nanstd(set_of_od_rsqrd(:)),'%.2f')])
clearvars h icol



figure(6)
set_of_mus_to_draw = [2,4,6];
for i_ds = 1:length(set_of_mus_to_draw)
    subplot(1,3,i_ds), hold on
    [MUA,D] = ndgrid(set_of_mua(2:end),d_cntrs);
    x = MUA .* D;
    y = squeeze(set_of________I(2:end,set_of_mus_to_draw(i_ds),:));
    plot(x(:),y(:),'o','Color',get_color(i_ds,length(set_of_mus_to_draw)),'MarkerFaceColor',get_color(i_ds,length(set_of_mus_to_draw)),'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1.5)
    xlabel('\mu_a \times d (unitless)')
    ylabel('OD (unitless)')
    title(['\mu_a \times d at \mu_s = ', ...
        num2str(set_of_mus(set_of_mus_to_draw(i_ds)),'%0.2f'), ...
        ' cm^{-1}'])
    axis square, grid on
    axis([0 1.5 0 20])
    set(gca,'xtick',0:0.5:1.5)
    set(gca,'ytick',0:5:20)
    set(gca,'FontSize',12)
    clearvars MUA D x y
end
clearvars set_of_mus_to_draw i_ds



figure(7)
set_of_mus_to_draw = [2,4,6];
for i_ds = 1:length(set_of_mus_to_draw)
    subplot(1,3,i_ds), hold on
    for i_mua = 2:length(set_of_mua)
        x = set_of_mua(i_mua) .* d_cntrs;
        y = squeeze(set_of________I(i_mua,...
                                    set_of_mus_to_draw(i_ds),:));
        plot(x,y,'o', ...
            'Color',get_color(i_mua-1,length(set_of_mua)-1), ...
            'MarkerFaceColor',get_color(i_mua-1,length(set_of_mua)-1), ...
            'MarkerEdgeColor','k', ...
            'MarkerSize',8, ...
            'LineWidth',1.5, ...
            'DisplayName',['\mu_a = ', ...
            num2str(set_of_mua(i_mua),'%0.3f'),' cm^{-1}'])
        clearvars x y
    end
    xlabel('\mu_a d (unitless)')
    ylabel('OD (unitless)')
    title(['OD vs. \mu_a d at \mu_s = ', ...
        num2str(set_of_mus(set_of_mus_to_draw(i_ds)),'%0.2f'), ...
        ' cm^{-1}'])
    axis square
    grid on
    axis([0 1.5 0 20])
    set(gca,'XTick',0:0.5:1.5)
    set(gca,'YTick',0:5:20)
    set(gca,'FontSize',12)
    legend('show','Location','northwest')
end
clearvars set_of_mus_to_draw i_ds



figure(8)
set_of_mus_to_draw = [2,4,6];
for i_ds = 1:length(set_of_mus_to_draw)
    subplot(1,3,i_ds), hold on
    for i_mua = 1:length(set_of_mua)
        x = d_cntrs;
        y = squeeze(set_of________I(i_mua,...
                                    set_of_mus_to_draw(i_ds),:));
        plot(x,y,'-', ...
            'Color',get_color(i_mua-1,length(set_of_mua)-1), ...
            'MarkerFaceColor',get_color(i_mua-1,length(set_of_mua)-1), ...
            'MarkerEdgeColor','k', ...
            'MarkerSize',8, ...
            'LineWidth',1.5, ...
            'DisplayName',[
            num2str(set_of_mua(i_mua),'%0.3f')])
        clearvars x y
    end
    clearvars i_mua
    xlabel('separation (cm)')
    ylabel('OD (unitless)')
    title(['OD vs. separation distance at \mu_s = ', ...
        num2str(set_of_mus(set_of_mus_to_draw(i_ds)),'%0.2f'),' cm^{-1}'])
    axis square
    grid on
    axis([0 8 0 24])
    set(gca,'XTick',0:4:15)
    set(gca,'YTick',0:8:24)
    set(gca,'FontSize',12)

    lgnd = legend('show','Location','southeast');
    lgnd.Title.String = '\mu_a (cm^{-1})';
    clearvars lgnd
end
clearvars set_of_mus_to_draw i_ds


figure(9)
set_of_mua_to_draw = [2,4,6];
for i_ds = 1:length(set_of_mua_to_draw)
    subplot(1,3,i_ds), hold on
    for i_mus = 1:length(set_of_mus)
        x = d_cntrs;
        y = squeeze(set_of________I(set_of_mua_to_draw(i_ds),i_mus,:));
        plot(x,y,'-','Color',get_color(i_mus,length(set_of_mus)),'MarkerFaceColor',get_color(i_mus,length(set_of_mus)),'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1.5,'DisplayName',['\mu_s = ',num2str(set_of_mus(i_mus),'%0.2f'),' cm^{-1}'])
        clearvars x y
    end
    clearvars i_mus
    xlabel('separation distance (cm)')
    ylabel('OD (unitless)')
    title(['OD vs. separation distance at \mu_a = ',num2str(set_of_mua(set_of_mua_to_draw(i_ds)),'%0.3f'),' cm^{-1}'])
    axis square
    grid on
    axis([0 15 0 24])
    set(gca,'XTick',0:5:15)
    set(gca,'YTick',0:8:24)
    set(gca,'FontSize',12)
    legend('show','Location','southeast')
end
clearvars set_of_mua_to_draw i_ds i_mus


figure(10)
icol = 1:3;
contourf(set_of_mua,set_of_mus,squeeze(set_of________I(:,:,1)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'OD (unitless)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('G (OD_0) vs. \{\mu_a,\mu_s\}');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol


figure(11)
set_of_mua_to_draw = [2, 4, 6];
contour_levels = 0:20;
for i_mua = 1:length(set_of_mua_to_draw)
    subplot(1,3,i_mua), hold on
    contourf(set_of_mus,d_cntrs,squeeze(set_of________I(set_of_mua_to_draw(i_mua),:,:)).',contour_levels,'LineStyle','none')
    xlabel('\mu_s (cm^{-1})')
    ylabel('separation (cm)')
    title(['OD vs. separation & \mu_s at \mu_a = ',num2str(set_of_mua(set_of_mua_to_draw(i_mua)),'%0.2f'),' cm^{-1}'])
    axis square
    set(gca,'FontSize',12)
    h = colorbar; ylabel(h,'OD (unitless)')
    clim([min(contour_levels) max(contour_levels)])
    axis([50 500 min(d_cntrs) 8])
    colormap jet
end
clearvars i_mua set_of_mua_to_draw c_mua contour_levels h


figure(12)
set_of_mus_to_draw = [2, 4, 6];
contour_levels = 0:20;
for i_mus = 1:length(set_of_mus_to_draw)
    subplot(1,3,i_mus), hold on
    contourf(set_of_mua,d_cntrs,squeeze(set_of________I(:,set_of_mus_to_draw(i_mus),:)).', contour_levels, 'LineStyle','none')
    xlabel('\mu_a (cm^{-1})')
    ylabel('separation (cm)')
    title(['OD vs. separation & \mu_a at \mu_s = ',num2str(set_of_mus(set_of_mus_to_draw(i_mus)),'%0.2f'),' cm^{-1}'])
    axis square
    set(gca,'FontSize',12)
    h = colorbar;
    ylabel(h,'OD (unitless)')
    clim([min(contour_levels) max(contour_levels)])
    axis([0.05 0.5 min(d_cntrs) 8])
    colormap jet
end
clearvars i_mus idx_mus set_of_mus_to_draw contour_levels h


figure(13)
set_of___d_to_draw = [12, 24, 36];
contour_levels = 0:20;
for i_d = 1:length(set_of___d_to_draw)
    subplot(1,3,i_d), hold on
    contourf(set_of_mua,set_of_mus,squeeze(set_of________I(:,:,set_of___d_to_draw(i_d))).',contour_levels,'LineStyle','none')
    xlabel('\mu_a (cm^{-1})')
    ylabel('\mu_s (cm^{-1})')
    title(['OD at separation = ',num2str(d_cntrs(set_of___d_to_draw(i_d)),'%0.1f'),' cm'])
    axis square
    set(gca,'FontSize',12)
    h = colorbar;
    ylabel(h,'OD (unitless)')
    clim([min(contour_levels) max(contour_levels)])
    axis([0.05 0.5 50 500])
    colormap jet
end
clearvars i_d idx_d set_of___d_to_draw contour_levels h


figure(14)
contour_levels = 0:20;
[X,Y,Z] = meshgrid(set_of_mus, set_of_mua, d_cntrs);
V = permute(set_of________I,[2 1 3]);
slice( ...
    X, Y, Z, V, ...
    set_of_mus([2 4 6]), ...
    set_of_mua([2 4 6]), ...
    d_cntrs([1 6 13 19 25 37]) )
shading interp
xlabel('\mu_s (cm^{-1})')
ylabel('\mu_a (cm^{-1})')
zlabel('separation (cm)')
h = colorbar;
ylabel(h,'OD (unitless)')
clim([min(contour_levels) max(contour_levels)])
colormap(flipud(jet))
axis([000 500 0.00 0.5 0 3])
view([90+37.50 30]), axis square
set(gca,'FontSize',12)
title('OD as a function of \mu_s, \mu_a, and separation')
clearvars X Y Z V h contour_levels


figure(15)
set_of___d_to_draw = [1, 7, 13, 19, 25, 31];
set_of___d_to_draw = 1:6;
contour_levels = 0:0.5:10;
for i_d = 1:length(set_of___d_to_draw)
    subplot(2,3,i_d), hold on
    contourf(set_of_mua,set_of_mus,squeeze(set_of________I(:,:,set_of___d_to_draw(i_d))).',contour_levels,'LineStyle','none')
    xlabel('\mu_a (cm^{-1})')
    ylabel('\mu_s (cm^{-1})')
    title(['OD at separation = ',num2str(d_cntrs(set_of___d_to_draw(i_d)),'%0.2f'),' cm'])
    axis square
    set(gca,'FontSize',12)
    h = colorbar;
    ylabel(h,'OD (unitless)')
    clim([min(contour_levels) max(contour_levels)])
    axis([0.05 0.5 50 500])
    colormap jet
end
clearvars i_d idx_d set_of___d_to_draw contour_levels h


figure(16)
i_s = 2;
set_of_ids = [12,24,36]; % -> 2, 4, 6
set_of_model_ids = [1,8,2,7,3,4];
for i_d = 1:length(set_of_ids)
    subplot(1,3,i_d)
    plot([min(set_of_mua) max(set_of_mua)], [min(set_of_mua) max(set_of_mua)], 'k:','DisplayName','1-to-1'), hold on
    for i_md = 1:length(set_of_model_ids)
        y = set_of_mua;
        yhat = squeeze(set_of_dpfs_est(:,i_s,set_of_model_ids(i_md),set_of_ids(i_d))).';
        plot(y, yhat, '-o', ...
            'Color',get_color_2(i_md),'MarkerFaceColor',get_color_2(i_md),'MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2,'DisplayName',[get_title(i_md),': e = ',num2str(rmse_100(y,yhat)*100,'%0.2f')]), hold on
        clearvars y yhat
    end
    clearvars i_md
    xlabel('True \mu_a (cm^{-1})')
    ylabel('Estimated \mu_a (cm^{-1})')
    grid on, axis equal, axis([min(set_of_mua) max(set_of_mua) min(set_of_mua) 2*max(set_of_mua)]), hold off
    set(gca,'FontSize',12)
    legend('show', 'Location','southeast', 'NumColumns',1)
    title([{'True vs. Estimated \mu_a Using Absolute DPF Models'},{['at d = ',num2str(d_cntrs(set_of_ids(i_d))),' cm']}])
end
clearvars i_d i_s set_of_ids
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

idx = all(isfinite(A),2) & isfinite(b);
A = A(idx,:);
b = b(idx);

% Solve for [C; a; b]
params = A \ b;
err = rsquared(b,A*params);
end
function out = get_color(idx,N)
if nargin < 2
    N = 12;
end
cmap = turbo(N);     % other options:
% cmap = parula(N);
cmap = jet(N);
% cmap = lines(N);
% cmap = hsv(N);
idx = max(1,min(idx,N));
out = cmap(idx,:);
end

function [r2] = rmse_100(y,y_hat)
y = y(:)+eps; y_hat = y_hat(:)+eps;
% align indices where both observed and predicted exist (non-NaN)
idx = ~isnan(y) & ~isnan(y_hat);
% compute r-squared using only indices where both observed and predicted exist (ss_res and ss_tot computed below)
r2 = mean(abs( (y(idx) - y_hat(idx)) ./ y(idx) ));
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
function [out] = get_color_2(idx)
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

