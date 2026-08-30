function [] = Fig_02_03_04_6
% Fig-02/03/04: w & s vs. d over mu_s & mu_a: more photons, more combinations
% the analysis: population stat - to plot p and q for mu_a's and mu_s's
% Based on Version 4, extended to investigate the dependence of G
% on the absorption coefficient (mu_a) and reduced scattering coefficient (mu_s').

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% optical & geometry properties
Lx = 29.1; Ly = 29.1; Lz = 05.0; % The size of the computational domain (in cm) to prevent reflections from the boundaries.

n = 1.4; g = 0.95;
set_of_mua = (.00:.05:0.5);
set_of_mus = (000:050:500);
dlta_d = .17;
d_edges = 0:dlta_d:sqrt((Lx/2).^2+(Ly/2).^2); N_bins = length(d_edges)-1; clearvars Lx Ly Lz
set_of________R = nan(length(set_of_mua),length(set_of_mus),N_bins);  % OD values
set_of________G = nan(length(set_of_mua),length(set_of_mus),1);       % G_values
clearvars N_bins
clc

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

        % I vs. d
        fun_x = @mean; fun_y = @sum; TheOutFun = @(x)(-log(x)); index_in = ind_diff; TheCode = 0;
        x_temp = sqrt(t_db.x(t_db.c==TheCode).^2+t_db.y(t_db.c==TheCode).^2+t_db.z(t_db.c==TheCode).^2);
        y_temp = t_db.w(t_db.c==TheCode);
        x_bind = accumarray(index_in,x_temp,[],fun_x,nan);
        y_bind = accumarray(index_in,y_temp,[],fun_y,nan); y_bind = TheOutFun((y_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ (t_db.no_of_photons./(pi.*dlta_d.*dlta_d)));
        r_bind = accumarray(index_in,y_temp,[],fun_y,nan); r_bind =           (r_bind./(pi.*dlta_d.*dlta_d + 2.*pi.*x_bind.*dlta_d)) ./ (t_db.no_of_photons                       );
        if isempty(y_bind), G = nan; else, G = y_bind(1); end % G
        set_of________G(i_a,i_s) = G;
        clearvars fun_x fun_y fun_s index_in TheCode x_temp y_temp s_temp

        set_of________R(i_a,i_s,1:length(y_bind)) = r_bind;
        clearvars mdl p blnFit blnShow fun_x fun_y fun_s i_fig index_in TheCode TheOutFun x_temp y_temp s_temp d_bind
        clearvars musp num den dpf_tmp
        clearvars x_bind y_bind
        clearvars t_db the_filename mua mus ind_diff ind_trns G y r_bind   
    end
end
clearvars i_a i_s n
clearvars d_edges

figure(1), subplot(1,2,1)
icol = 0:4;
contourf(set_of_mua,set_of_mus,squeeze(set_of________R(:,:,1)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'R (cm^{-2})'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('Reflc. vs. \{\mu_a,\mu_s\}');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol
figure(1), subplot(1,2,2)
icol = 1:3;
contourf(set_of_mua,set_of_mus,squeeze(set_of________G(:,:)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'G (unitless)'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('G (OD_0) vs. \{\mu_a,\mu_s\}');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

% Radius of the cylinder and sphere, in cm.
% R0 is not used for the semi-infinite geometry.
R0 = 20;
[set_of_Phi0] = do_Piao (set_of_mua, set_of_mus, g, R0);

figure(2), subplot(2,2,1)
icol = 0:4;
contourf(set_of_mua,set_of_mus,squeeze(set_of________R(:,:,1)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'R(d\rightarrow0) (cm^{-2})'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('Reflectance at d \rightarrow 0: Semi-infinite geometry');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

figure(2), subplot(2,2,2)
icol = 0:20:80;
contourf(set_of_mua,set_of_mus,squeeze(set_of_Phi0(:,:,1)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'\Psi(d\rightarrow0)/S (cm^{-2})'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('Fluence rate at d \rightarrow 0: Semi-infinite geometry');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

figure(2), subplot(2,2,3)
icol = 0:20:80;
contourf(set_of_mua,set_of_mus,squeeze(set_of_Phi0(:,:,2)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'\Psi(d\rightarrow0)/S (cm^{-2})'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('Fluence rate at d \rightarrow 0: Cylindrical geometry');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol

figure(2), subplot(2,2,4)
icol = 0:20:80;
contourf(set_of_mua,set_of_mus,squeeze(set_of_Phi0(:,:,2)).',linspace(min(icol),max(icol),21),'EdgeColor','none'), shading interp
xlabel('\mu_a (cm^{-1})'), set(gca,'xtick',[set_of_mua(2) max(set_of_mua)])
ylabel('\mu_s (cm^{-1})'), set(gca,'ytick',[set_of_mus(2) max(set_of_mus)])
h = colorbar; ylabel(h,'\Psi(d\rightarrow0)/S (cm^{-2})'), set(h,'ylim',[min(icol) max(icol)]), set(h,'ytick',icol), title('Fluence rate at d \rightarrow 0: Spherical geometry');
set(gca,'fontsize',16), axis square, axis([set_of_mua(2) max(set_of_mua) set_of_mus(2) max(set_of_mus)]); colormap jet
clearvars h icol
end

function [set_of_Phi0] = do_Piao (set_of_mua, set_of_mus, g, R0)
% Analytical implementation of the zero-separation fluence rate
% (i.e., detected intensity at d = 0) based on the formulation of
% Piao et al.:
% Piao D, Barbour RL, Graber HL, Lee DC. "On the geometry dependence of
% differential pathlength factor for near-infrared spectroscopy. I.
% Steady-state with homogeneous medium." J Biomed Opt. 2015;20(10):105005.
% doi:10.1117/1.JBO.20.10.105005.

%% Fluence rate at d -> 0
% Based on Table 1 of Piao et al.
%
% Array dimensions:
%   dimension 1 = mu_a
%   dimension 2 = mu_s
%   dimension 3 = geometry
%
% Geometry index:
%   1 = semi-infinite medium
%   2 = infinite cylinder, azimuthal direction
%   3 = sphere
%
% The calculated quantity is normalized fluence rate:
%
%   Phi_0 = Psi(d -> 0)/S
%
% If the optical coefficients are in cm^-1, Phi_0 has units cm^-2.

% Boundary-condition factor.
% A = 1 corresponds to a refractive-index-matched boundary.
A = 1;

% % Radius of the cylinder and sphere, in cm.
% % R0 is not used for the semi-infinite geometry.
% R0 = 5;

%% Create optical-property grid

% Dimensions:
% length(set_of_mua) x length(set_of_mus)

[MUA,MUS] = ndgrid(set_of_mua,set_of_mus);

% Reduced scattering coefficient
MUSP = (1-g).*MUS;

%% Initialize diffusion parameters

D       = nan(size(MUA));
mueff   = nan(size(MUA));
Ra      = nan(size(MUA));
Rb      = nan(size(MUA));

lr0     = nan(size(MUA));
li0     = nan(size(MUA));

eta_azi = nan(size(MUA));
eta_sph = nan(size(MUA));

%% Valid optical-property points

% At mu_s = 0:
%
%   mu_s' = 0
%   D     = infinity
%   Ra    = infinity
%
% Therefore, the diffusion expressions are undefined.

valid_optical = MUSP > 0;

%% Diffusion parameters

% Diffusion coefficient:
%
%   D = 1/(3*mu_s')

D(valid_optical) = ...
    1 ./ (3.*MUSP(valid_optical));

% Effective attenuation coefficient:
%
%   mu_eff = sqrt(3*mu_a*mu_s')

mueff(valid_optical) = ...
    sqrt(3.*MUA(valid_optical).*MUSP(valid_optical));

% Equivalent real-source depth:
%
%   Ra = 1/mu_s'

Ra(valid_optical) = ...
    1 ./ MUSP(valid_optical);

% Extrapolated-boundary distance:
%
%   Rb = 2*A*D

Rb(valid_optical) = ...
    2.*A.*D(valid_optical);

%% Real- and image-source distances at d -> 0

% In the limit d -> 0:
%
%   lr -> Ra
%   li -> Ra + 2*Rb

lr0(valid_optical) = ...
    Ra(valid_optical);

li0(valid_optical) = ...
    Ra(valid_optical) + 2.*Rb(valid_optical);

%% Real-source and image-source contributions

real_term  = nan(size(MUA));
image_term = nan(size(MUA));
prefactor  = nan(size(MUA));

real_term(valid_optical) = ...
    exp(-mueff(valid_optical).*lr0(valid_optical)) ...
    ./ lr0(valid_optical);

image_term(valid_optical) = ...
    exp(-mueff(valid_optical).*li0(valid_optical)) ...
    ./ li0(valid_optical);

% Normalized fluence-rate prefactor:
%
%   Psi/S = 1/(4*pi*D) * (...)

prefactor(valid_optical) = ...
    1 ./ (4.*pi.*D(valid_optical));

%% Initialize normalized zero-separation fluence

set_of_Phi0 = nan( ...
    length(set_of_mua), ...
    length(set_of_mus), ...
    3);

%% Geometry 1: Semi-infinite medium

% For the semi-infinite geometry:
%
%   eta = 1

set_of_Phi0(:,:,1) = ...
    prefactor .* ...
    (real_term - image_term);

%% Validity condition for curved geometries

% The cylindrical and spherical curvature factors contain:
%
%   R0 - Ra
%
% Therefore, R0 must be greater than Ra.

valid_geometry = ...
    valid_optical & ...
    isfinite(Ra) & ...
    (R0 > Ra);

%% Geometry-dependent correction factors

% Infinite cylinder evaluated in the azimuthal direction:
%
%   eta_azi = sqrt[(R0 + Ra + 2*Rb)/(R0 - Ra)]

eta_azi(valid_geometry) = sqrt( ...
    (R0 + Ra(valid_geometry) + 2.*Rb(valid_geometry)) ...
    ./ ...
    (R0 - Ra(valid_geometry)) );

% Sphere:
%
%   eta_sph = (R0 + Ra + 2*Rb)/(R0 - Ra)

eta_sph(valid_geometry) = ...
    (R0 + Ra(valid_geometry) + 2.*Rb(valid_geometry)) ...
    ./ ...
    (R0 - Ra(valid_geometry));

%% Geometry 2: Infinite cylinder, azimuthal direction

Phi_azi = nan(size(MUA));

Phi_azi(valid_geometry) = ...
    prefactor(valid_geometry) .* ...
    ( ...
      real_term(valid_geometry) ...
      - eta_azi(valid_geometry).*image_term(valid_geometry) ...
    );

set_of_Phi0(:,:,2) = Phi_azi;

%% Geometry 3: Sphere

Phi_sph = nan(size(MUA));

Phi_sph(valid_geometry) = ...
    prefactor(valid_geometry) .* ...
    ( ...
      real_term(valid_geometry) ...
      - eta_sph(valid_geometry).*image_term(valid_geometry) ...
    );

set_of_Phi0(:,:,3) = Phi_sph;

%% Remove undefined or nonphysical fluence values

set_of_Phi0(~isfinite(set_of_Phi0)) = nan;
set_of_Phi0(set_of_Phi0 <= 0) = nan;

%% Check the large-radius approximation condition

% Piao et al. suggest that the curved-geometry approximation
% is generally appropriate when:
%
%   mu_eff*R0 >= 10

large_radius_condition = ...
    valid_geometry & ...
    (mueff.*R0 >= 10);

%% Geometry names

geometry_names = { ...
    'Semi-infinite geometry', ...
    'Infinite cylindrical geometry', ...
    'Spherical geometry'};

%% Display numerical ranges

fprintf('\n');
fprintf('Parameters:\n');
fprintf('g  = %.4f\n',g);
fprintf('A  = %.4f\n',A);
fprintf('R0 = %.4f cm\n\n',R0);

for i_geometry = 1:3

    temp_Phi = set_of_Phi0(:,:,i_geometry);
    temp_Phi = temp_Phi(isfinite(temp_Phi));

    fprintf('%s:\n',geometry_names{i_geometry});

    if isempty(temp_Phi)

        fprintf('  No valid fluence values.\n\n');

    else

        fprintf('  Psi(0)/S range = %.6g to %.6g cm^-2\n\n', ...
            min(temp_Phi),max(temp_Phi));

    end

end

clearvars h i_geometry i_s temp_Phi temp_all Z
end

