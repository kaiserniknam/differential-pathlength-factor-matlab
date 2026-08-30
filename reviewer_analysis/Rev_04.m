function Rev_04
% Plot optical density and the spatial distribution of photons exiting the top surface.

clc
close all

% Bin width for source-detector separation, in cm
dlta_d = 0.17;

% Baseline optical properties, in cm^-1
mus = 55;
mua = 0.0275;

% Computational domain size, in cm
Lx = 19.1;
Ly = 19.1;

% TiO2 concentration scaling factors
set_of_musphant = [0.25, 1.00];

figure('Color','w','Position',[100 100 1700 600])

for i_musphant = 1:length(set_of_musphant)

    % Define color changing from blue to red across concentrations
    the_color = [(i_musphant-1)/(length(set_of_musphant)-1), 0, 1-(i_musphant-1)/(length(set_of_musphant)-1)];

    % Load Photon_87 simulation file for the selected TiO2 concentration
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_87_concentration_',num2str(set_of_musphant(i_musphant),'%.2f'),'.mat'];

    % Compute diffuse reflectance and photon exit positions
    [y_bind,x_bind,t_db] = get_OD(the_filename,dlta_d);

    subplot(1,3,1)
    plot(x_bind,-log(y_bind),'LineWidth',2,'LineStyle','-','Color',the_color,'DisplayName',['\mu_a = ',num2str(mua),' cm^{-1}, \mu_s = ',num2str(mus*set_of_musphant(i_musphant)),' cm^{-1}'])
    hold on
    xlabel('Source-detector separation, d (cm)')
    ylabel('Optical density, OD')
    title('Optical Density versus Source-Detector Separation')
    grid on
    axis square
    axis([0 14 2 24])
    set(gca,'FontSize',16)
    legend('show','Location','northwest','NumColumns',1,'Orientation','vertical')

    subplot(1,3,i_musphant+1)
    histogram2(t_db.x_ot,t_db.y_ot,-sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d:5*dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d,-sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d:5*dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d,'DisplayStyle','bar3','FaceColor',the_color)
    set(gca,'ZScale','log')
    xlabel('Photon exit position, x (cm)')
    ylabel('Photon exit position, y (cm)')
    zlabel('Number of detected photons')
    title({['Spatial Distribution of Detected Photons'],['\mu_a = ',num2str(mua),' cm^{-1}, \mu_s = ',num2str(mus*set_of_musphant(i_musphant)),' cm^{-1}']})
    grid on
    axis square
    axis([-10 10 -10 10 1 1e5])
    set(gca,'FontSize',16)

    clearvars x_bind y_bind t_db the_filename the_color
end

sgtitle('Diffuse Reflectance and Photon Exit Distribution for 10^6 Incident Photons','FontSize',18)

end


function [y_bind,x_bind,t_db] = get_OD(the_filename,dlta_d)
% Compute binned diffuse reflectance from Photon_87 simulation data.
% Source-detector separation is calculated from the photon entry and exit positions.

% Computational domain size, in cm
Lx = 19.1;
Ly = 19.1;

% Load simulation database
t_db = load(the_filename);
clearvars the_filename

% Remove noisy or inconsistent photon histories
u = round(-log(t_db.w)./t_db.s,4);
u_unique = unique(u);
u_unique = u_unique(~isnan(u_unique));
freq = nan(size(u_unique));

for i_f = 1:length(u_unique)
    freq(i_f) = sum(u==u_unique(i_f));
end

[~,i_fx] = max(freq);

% Keep photons belonging to the dominant absorption group
t_db.x_in = t_db.x_in(u==u_unique(i_fx));
t_db.y_in = t_db.y_in(u==u_unique(i_fx));
t_db.z_in = t_db.z_in(u==u_unique(i_fx));
t_db.x_ot = t_db.x_ot(u==u_unique(i_fx));
t_db.y_ot = t_db.y_ot(u==u_unique(i_fx));
t_db.z_ot = t_db.z_ot(u==u_unique(i_fx));
t_db.s = t_db.s(u==u_unique(i_fx));
t_db.w = t_db.w(u==u_unique(i_fx));

% Display total photons, retained photons, retained percentage, and dominant mua
u = u(~isnan(u));
disp(num2str([length(u),sum(u==u_unique(i_fx)),sum(u==u_unique(i_fx))/length(u)*100,u_unique(i_fx)]))

% Update photon number after removing noisy photons
t_db.no_of_photons = t_db.no_of_photons-sum(u~=u_unique(i_fx));

clearvars u u_unique freq i_f i_fx

% Define radial bins for source-detector separation
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d;

% Select photons exiting through the top surface
idx = t_db.z_ot<=0;

% Calculate source-detector separation
x_temp = sqrt((t_db.x_in(idx)-t_db.x_ot(idx)).^2+(t_db.y_in(idx)-t_db.y_ot(idx)).^2+(t_db.z_in(idx)-t_db.z_ot(idx)).^2);
y_temp = t_db.w(idx);

% Assign each detected photon to a radial-distance bin
[~,~,index_in] = histcounts(x_temp,d_diff_edges);

% Compute mean separation and total detected weight in each bin
x_bind = accumarray(index_in,x_temp,[],@mean,nan);
y_bind = accumarray(index_in,y_temp,[],@sum,nan);

% Normalize detected weight by annular detection area and incident photon density
y_bind = (y_bind./(pi.*dlta_d.^2+2.*pi.*x_bind.*dlta_d))./(t_db.no_of_photons./(pi.*dlta_d.^2));

end