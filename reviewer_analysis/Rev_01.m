% Visualize the diffusion-theory validity metric:
%
%     sqrt(3 * mua * mus * (1-g))
%
% over a range of absorption and scattering coefficients.
%
% The black contour corresponds to the boundary where:
%
%     sqrt(3 * mua * mus * (1-g)) = 1
%
% which is often used as an approximate diffusion-limit criterion.

function [] = Rev_01

clc
close all

% Anisotropy factor
g = 0.93;

% Range of absorption coefficients (cm^-1)
mua = linspace(0.00,0.5,300);

% Range of scattering coefficients (cm^-1)
mus = linspace(0,500,200);

% Generate 2-D parameter grid
[MUA,MUS] = meshgrid(mua,mus);

% Compute diffusion-validity parameter
val = sqrt(3 .* MUA .* MUS .* (1-g));

figure

% Filled contour plot
contourf(mua,mus,val,25,'LineStyle','none')

% Keep y-axis in standard orientation
set(gca,'YDir','normal')

% Add color scale
colorbar

% Axis labels
xlabel('\mu_a (cm^{-1})')
ylabel('\mu_s (cm^{-1})')

% Figure title
title('\surd(3\mu_a\mu_s(1-g))')

hold on

% Draw the diffusion-limit boundary:
% sqrt(3*mua*mus*(1-g)) = 1
contour(mua,mus,val,[1 1],'k','LineWidth',3)

% Annotate the boundary line
text(0.08,120,' = 1 boundary', ...
    'Color','k', ...
    'FontSize',14)

% Improve figure appearance
set(gca,'FontSize',16)
axis square
colormap hsv

end