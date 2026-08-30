function [] = Rev_03
% show how mean pathlengh, and mean-pathlengh, and slope-based DPF's change
% versus mu_a

clc
close all

% Bin width for source-detector separation, in cm
dlta = 0.17;

% Baseline optical properties
set_of_mua = [0.00,0.20];
set_of_mus = [5, 30];

for i_a = 1:length(set_of_mua)
    for i_s = 1:length(set_of_mus)
        mua = set_of_mua(i_a);
        mus = set_of_mus(i_s);
        disp(['mua = ',num2str(mua,'%.2f'),', mus = ',num2str(mus,'%.0f')])
        TheColor = get_color(mua/max(set_of_mua),mus/max(set_of_mus)); % make my own colormap

        % read dbase
        the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_81_mua_',sprintf('%.2f',mua),'_mus_',sprintf('%.2f',mus),'.mat'];
        [y_bind, x_bind, s_bind, e_bind] = get_OD(the_filename, dlta, mua);

        subplot(1,2,1),
        plot(x_bind,s_bind,'LineWidth',2,'LineStyle','-.','Color',TheColor,'DisplayName',['\mu_a = ',num2str(mua),', \mu_s = ',num2str(mus),' cm^{-1}']), hold on
        xlabel('separation (cm)'), ylabel('mean pathlength (cm)')
        title('mean photon pathlength vs. separation')
        set(gca,'FontSize',18)
        axis([0 8 0 45]), legend('show','Orientation','vertical','numcolumns',1,'location','northwest')
        subplot(1,2,2),
        plot(x_bind,s_bind./x_bind,'LineWidth',4,'LineStyle','-.','Color',TheColor,'DisplayName',['DPF_{mean}: \mu_a = ',num2str(mua),', \mu_s = ',num2str(mus)]), hold on
        plot(x_bind,e_bind./x_bind,'LineWidth',1,'LineStyle','-', 'Color',TheColor,'DisplayName',['DPF_{slope}: \mu_a = ',num2str(mua),', \mu_s = ',num2str(mus)]), hold on
        xlabel('separation (cm)'), ylabel('DPF (unitless)')
        title('Various DPF vs. separation')
        set(gca,'FontSize',18)
        axis([0 8 0 9]), legend('show','Orientation','horizontal','numcolumns',2,'location','northeast')
        
        clearvars mua mus TheColor the_filename y_bind x_bind s_bind e_bind
    end
end
end

function [y_bind, x_bind, s_bind, e_bind] = get_OD(the_filename, dlta_d, mua)
% Bin exiting photon weights as a function of diffuse source-detector distance.
Lx = 29.1; Ly = 29.1; % The size of the computational domain (in cm) to prevent reflections from the boundaries.
t_db = load(the_filename);

% 1-D radial binning of photons exiting from the top surface.
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2 + (Ly/2).^2)+dlta_d;
idx = t_db.z <= 0;
x_temp = sqrt(t_db.x(idx).^2 + ...
              t_db.y(idx).^2 + ...
              t_db.z(idx).^2);
y_temp = t_db.w(idx);
s_temp = t_db.s(idx);
[~, ~, index_in] = histcounts(x_temp, d_diff_edges);

% Mean distance and summed detected weight in each radial bin.
x_bind = accumarray(index_in, x_temp,                    [], @mean, nan);
y_bind = accumarray(index_in, y_temp,                    [], @sum,  nan);
s_bind = accumarray(index_in, s_temp,                    [], @mean, nan);
e_bind = accumarray(index_in, s_temp.*exp(-mua.*s_temp), [], @mean, nan);

% Convert summed photon weight to normalized diffuse reflectance.
ring_area = pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d;
source_area = pi.*dlta_d.*dlta_d;
y_bind = (y_bind ./ ring_area) ./ (t_db.no_of_photons ./ source_area);
x_bind = (d_diff_edges(1:end-1)+d_diff_edges(2:end-0))./2; x_bind = x_bind(1:length(y_bind)).';
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



