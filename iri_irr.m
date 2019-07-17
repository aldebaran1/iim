clc;
clear all;
close all;

plot_vertical_perturbation = 0;
plot_horizontal_perturbation = 1;
vertical = 1;
horizontal = 1;
% Wavelengths
lambda_z = 20;
lambda_x = 3;
lambda_y = 0;
% Dimensions and resolution
X = 10; % km
Y = 10; % km
dx = 0.1; % km resolution
dy = 0.1; % km resolution
dz = 1; % km
% Amplittude of irregularity in Z- direction
A_z = 0.1;
% Horizontal baselines
x = 0 : dx : X;
y = 0 : dy : Y;

% IRI import
irifn = 'iri_baseline.txt';
IRI = readtable(irifn);
altkm = table2array(IRI(:,1));
Ne = table2array(IRI(:,2));
% Vertical baseline after the IRI profile
z = altkm(1) : dz : altkm(end);
Nei = interp1(altkm, Ne, z, 'spline');

% 3D empty boxes
Ne_box = ones(length(z), length(x), length(y));
dNe_box = ones(length(z), length(x), length(y));

if vertical
    %% Vertical Irregularity Parameters
    % Initial parameters
    z0 = 150; % km
    z1 = 350; % km
    idh = (z >= z0) & (z <= z1);
    zi_n = z(idh)';

    k_z = 2 * pi / lambda_z; % km
     
    % Taper?
    taper = hanning(length(zi_n))';
    Zh = zeros(1, length(z));
    Zh(idh) = Zh(idh) + taper;
    % Sinusioidal peturbation
    P_z = ((A_z * Nei) .* sin(k_z * z)) .* Zh;
    Ne_irr_z = Nei + P_z;
    % Plot sin
    if plot_vertical_perturbation
        % Plot
        figure()
        hold on
        semilogx(Ne, altkm, 'b')
        semilogx(Ne_irr_z, z, 'r')
        axis([0, 5*1e11, 1e2, 600])
        xlabel('Electron denisty Ne [el/m3]');
        ylabel('Altitude [km]');
    end
    
else
    Ne_irr_z = Nei;
end
%% Horizontal perturbations
% Initial profile is accounted for vertical perturbations if any
Ne_x = ones(1, length(x)) .* Ne_irr_z';

if horizontal
    k_x = 2 * pi / lambda_x;
    A_x = 0.1 * Ne_irr_z;
    % Sinusoidal Perturbation
    P_x = A_x' .* sin(k_x * x);
    Ne_irr_x = Ne_x + P_x;
    
    % Fill the Ne box
    for i=1 : length(y)
        Ne_box(:, :, i) = Ne_irr_x;
        dNe_box(:, :, i) = P_x;
    end
    % Compte TEC and dTEC
    % Plot 3D TEC field
    tec = zeros(length(x), length(y));
    dtec = zeros(length(x), length(y));
    for i = 1 : length(x)
        for j = 1 : length(y)
            tec(i, j) = trapz(z, Ne_box(:, i, j) * 1e3 / 1e16);
            dtec(i, j) = trapz(z, dNe_box(:, i, j) * 1e3 / 1e16);
        end
    end
    % Make plots
    if plot_horizontal_perturbation
    % Plot Perturbations in one horizontal direction
        figure()
        h = pcolor(x, z, P_x/1e10);
        colormap jet
        cbar = colorbar; ylabel(cbar, 'Ne perturbation [1e10 el/m3]');
        set(h, 'EdgeColor', 'none')
        xlabel('Horizontal distance x [km]'); ylabel('Altitude z [km]')
        % Plot Ionosphere profile with perturbations
        figure()
        ax1 = subplot(2, 1, 1);
        h = pcolor(x, z, Ne_irr_x);
        colormap jet
        cbar = colorbar; ylabel(cbar, 'Total Ne [el/m3]');
        xlabel('Horizontal distance x [km]'); ylabel('Altitude z [km]')
        set(h, 'EdgeColor', 'none')
        ax2 = subplot(2, 1, 2);
        plot(x, trapz(z, Ne_irr_x,1)*1e3/1e16, 'b')
        xlabel('Horizontal distance x [km]'); ylabel('TEC [TECu]')
        linkaxes([ax1,ax2],'x')
        
        figure()
        subplot(1, 2, 1)
        h=pcolor(x, y, tec');
        colormap jet
        cbar = colorbar; ylabel(cbar, 'TEC [TECu]')
        set(h, 'EdgeColor', 'none')
        xlabel('Horizontal distance x [km]'); ylabel('Horizontal distance y [km]')
        subplot(1, 2, 2)
        h=pcolor(x, y, dtec');
        colormap jet
        cbar = colorbar; ylabel(cbar, 'dTEC [TECu]')
        set(h, 'EdgeColor', 'none')
        xlabel('Horizontal distance x [km]'); ylabel('Horizontal distance y [km]')
    end
else 
    Ne_irr_x = 0; % ??
end

%% Spectral analysis
Y = fftn(dNe_box);
figure()
imagesc(log(squeeze(sum(abs(Y),1))))