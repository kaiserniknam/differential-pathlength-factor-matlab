function [] = Rev_02
% Compare optical density (OD) as a function of source-detector separation
% for two TiO2 concentrations. Solid lines show Photon_87 results, and
% dotted lines show reference Photon_58 results with matched optical
% properties.

clc
close all

% Bin width for source-detector separation, in cm
dlta_d = 0.17;

% Baseline optical properties
mus = 55;
mua = 0.0275;

% TiO2 concentration scaling factors
set_of_musphant = [0.25, 1.00];

for i_musphant = 1:length(set_of_musphant)

    % Load Photon_87 simulation file for the selected TiO2 concentration
    the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_87_concentration_', ...
                    num2str(set_of_musphant(i_musphant),'%.2f'),'.mat'];

    % Define color changing from blue to red across concentrations
    the_color = [(i_musphant-1)/(length(set_of_musphant)-1), ...
                 0, ...
                 1-(i_musphant-1)/(length(set_of_musphant)-1)];

    % Compute and plot OD from Photon_87 data
    [y_bind,x_bind] = get_OD(the_filename,dlta_d);
    plot(x_bind,-log(y_bind), ...
        'LineWidth',2, ...
        'LineStyle','-', ...
        'Color',the_color, ...
        'DisplayName',['C_{TiO_2} = ',num2str(set_of_musphant(i_musphant),'%.2f'), ...
        ', # of photons = 10^6']), hold on

    clearvars x_bind y_bind

    % Compute and plot OD from reference Photon_58 data
    % Effective scattering coefficient is scaled by TiO2 concentration
    [y_bind,x_bind] = get_OD_prime(mua, mus*set_of_musphant(i_musphant), dlta_d);
    plot(x_bind,-log(y_bind), ...
        'LineWidth',4, ...
        'LineStyle',':', ...
        'Color',the_color, ...
        'DisplayName',['C_{TiO_2} = ',num2str(set_of_musphant(i_musphant),'%.2f'), ...
        ', # of photons = 10^7']), hold on

    clearvars x_bind y_bind the_filename the_color
end

% Format figure
xlabel('separation (cm)')
ylabel('OD (a.u.)')
grid on
axis square
axis([0 8 2 12])
axis([0 15 2 22])
hold off

set(gca,'fontsize',20)
legend('show','Location','southeast','NumColumns',1,'orientation','horizontal')
title('OD vs source-detector separation')

end


function [y_bind,x_bind] = get_OD(the_filename,dlta_d)

% Compute binned diffuse reflectance and OD from Photon_87 simulation data.
% This version uses input/output photon positions to calculate the lateral
% source-detector separation at the top surface.

% Computational domain size, in cm
Lx = 19.1;
Ly = 19.1;

% Load simulation database
t_db = load(the_filename);
clearvars the_filename

% Remove noisy or inconsistent photon histories.
% Here, absorption coefficient estimates are inferred from -log(w)/s.
% The most frequent rounded value is kept as the reliable photon group.
u = round(-log(t_db.w)./t_db.s,4);
u_unique = unique(u);
u_unique = u_unique(~isnan(u_unique));
freq = nan(size(u_unique));

for i_f = 1:length(u_unique)
    freq(i_f) = sum(u==u_unique(i_f));
end

[~,i_fx] = max(freq);

% Keep only photons belonging to the dominant absorption group
t_db.x_in = t_db.x_in(u==u_unique(i_fx));
t_db.y_in = t_db.y_in(u==u_unique(i_fx));
t_db.z_in = t_db.z_in(u==u_unique(i_fx));
t_db.x_ot = t_db.x_ot(u==u_unique(i_fx));
t_db.y_ot = t_db.y_ot(u==u_unique(i_fx));
t_db.z_ot = t_db.z_ot(u==u_unique(i_fx));
t_db.s    = t_db.s(u==u_unique(i_fx));
t_db.w    = t_db.w(u==u_unique(i_fx));

% Display: total valid photons, kept photons, kept percentage, dominant mua
u = u(~isnan(u));
disp(num2str([length(u), ...
              sum(u==u_unique(i_fx)), ...
              sum(u==u_unique(i_fx))/length(u)*100, ...
              u_unique(i_fx)]))

% Update photon number after removing noisy photons
t_db.no_of_photons = t_db.no_of_photons - sum(u~=u_unique(i_fx));

clearvars u u_unique freq i_f i_fx

% Define radial bins for source-detector separation
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2)+dlta_d;

% Select photons exiting through the top surface
idx = t_db.z_ot <= 0;

% Assign each detected photon to a radial-distance bin
[~,~,index_in] = histcounts(sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2), d_diff_edges);

clearvars d_diff_edges Lx Ly

% Compute mean separation and total detected weight in each bin
fun_x = @mean;
fun_y = @sum;

x_temp = sqrt( ...
    (t_db.x_in(idx)-t_db.x_ot(idx)).^2 + ...
    (t_db.y_in(idx)-t_db.y_ot(idx)).^2 + ...
    (t_db.z_in(idx)-t_db.z_ot(idx)).^2);

y_temp = t_db.w(idx);

x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);

% Normalize detected weight by annular detection area and incident photon density
y_bind = ...
    (y_bind ./ (pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons ./ (pi.*dlta_d.*dlta_d));

end


function [y_bind,x_bind] = get_OD_prime(mua,mus,dlta_d)

% Compute binned diffuse reflectance and OD from Photon_58 reference data.
% The file is selected based on mua and mus.

% Computational domain size, in cm
Lx = 29.1;
Ly = 29.1;

% Load reference simulation file
the_filename = ['/home/kaiser/BackUp/Dropbox/research/2.NIRs Project/code/DB/Photon_58_mua_', ...
                sprintf('%.4f',mua), ...
                '_mus_', ...
                sprintf('%.2f',mus), ...
                '.mat'];

t_db = load(the_filename);
clearvars the_filename

% Remove noisy or inconsistent photon histories.
% The dominant rounded value of -log(w)/s is assumed to represent the
% intended absorption coefficient.
u = round(-log(t_db.w)./t_db.s,4);
u_unique = unique(u);
freq = nan(size(u_unique));

for i_f = 1:length(u_unique)
    freq(i_f) = sum(u==u_unique(i_f));
end

[~,i_fx] = max(freq);

% Keep only photons belonging to the dominant absorption group
t_db.x = t_db.x(u==u_unique(i_fx));
t_db.y = t_db.y(u==u_unique(i_fx));
t_db.z = t_db.z(u==u_unique(i_fx));
t_db.d = t_db.d(u==u_unique(i_fx));
t_db.s = t_db.s(u==u_unique(i_fx));
t_db.w = t_db.w(u==u_unique(i_fx));
t_db.c = t_db.c(u==u_unique(i_fx));

% Update photon number after removing noisy photons
t_db.no_of_photons = t_db.no_of_photons - sum(u~=u_unique(i_fx));

clearvars u u_unique freq i_f i_fx

% Define radial bins for source-detector separation
d_diff_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2);

% Photon exit code used for diffuse reflectance
TheCode = 0;

% Assign photons with the selected exit code to radial-distance bins
[~,~,index_in] = histcounts(sqrt( ...
    t_db.x(t_db.c==TheCode).^2 + ...
    t_db.y(t_db.c==TheCode).^2 + ...
    t_db.z(t_db.c==TheCode).^2), d_diff_edges);

clearvars d_diff_edges Lx Ly

% Compute mean separation and total detected weight in each bin
fun_x = @mean;
fun_y = @sum;

x_temp = sqrt( ...
    t_db.x(t_db.c==TheCode).^2 + ...
    t_db.y(t_db.c==TheCode).^2 + ...
    t_db.z(t_db.c==TheCode).^2);

y_temp = t_db.w(t_db.c==TheCode);

x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
y_bind = accumarray(index_in,y_temp,[],fun_y,nan);

% Normalize detected weight by annular detection area and incident photon density
y_bind = ...
    (y_bind ./ (pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ ...
    (t_db.no_of_photons ./ (pi.*dlta_d.*dlta_d));

end


