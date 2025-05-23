% Given Parameters
f_min = 1200e6; % Minimum frequency
f_max = 2520e6; % Maximum frequency
directivity_dB = 8.7; % Directivity in dB
tau = 0.925; % Geometric ratio
sigma = 0.1738; % Spacing factor
Z0 = 50; % Characteristic impedance
c = 3e8; % Speed of light
lambda_min = c / f_max; % Min wavelength
lambda_max = c / f_min; % Max wavelength

% Frequency range
f = linspace(f_min, f_max, 1000); % Frequency array
lambda = c ./ f; % Wavelength array

% Apex half angle
alpha = atan((1 - tau) / (4 * sigma));

% Calculated Parameters for Log-Periodic Design
B = f_max / f_min; % Bandwidth
B_ar = 1.1 + 7.7 * (1 - tau)^2 * cot(alpha); % Bandwidth of Active Region
B_s = B * B_ar; % Designed Bandwidth
L = (lambda_max / 4) * (1 - 1 / B_s) * cot(alpha); % Antenna length

% Number of Elements in the Array
N = ceil(1 + log(B_s) / log(1 / tau));

% Average Characteristic Impedance
d_n = 0.01; % Dipole element diameter
Z_a = 120 * (log(lambda_max / d_n) - 2.25);

% Spacing Between Feeder Line Conductors
s = d_n * cosh(Z0 / 120);

% Calculated Parameters
fprintf('Antenna Parameters:\n');
fprintf('Bandwidth of Active Region, B_ar: %.2f\n', B_ar);
fprintf('Designed Bandwidth, B_s: %.2f\n', B_s);
fprintf('Antenna Length, L: %.2f m\n', L);
fprintf('Number of Elements, N: %d\n', N);
fprintf('Average Characteristic Impedance, Z_a: %.2f Ohms\n', Z_a);
fprintf('Spacing Between Conductors, s: %.2f mm\n', s * 1000);

% Logarithmic Impedance Variation and Return Loss
k = 0.2; % Rate of change of impedance with frequency
Rin = Z0 .* (1 + k * log(f ./ f_min)); % Input impedance across frequency

% Reflection Coefficient and Return Loss
Gamma = (Rin - Z0) ./ (Rin + Z0); % Reflection coefficient
RL = 20 * log10(abs(Gamma)); % Return loss in dB

% Return Loss across Frequency Range Plot
figure;
plot(f / 1e6, RL, 'LineWidth', 1.5);
xlabel('Frequency (MHz)');
ylabel('Return Loss (dB)');
title('Return Loss across Frequency Range');
grid on;

% Antenna Structure Plot
l_max = 0.5 * lambda_max; % Max element length
l_min = 0.5 * lambda_min; % Min element length
lengths = l_max * tau.^((0:N-1)); % Element lengths
positions = cumsum([0 lengths(1:end-1)] * sigma); % Element positions
Dipole_Number = (1:N)'; % Dipole number
Length_mm = lengths' * 1000; % Length of each dipole
Spacing_mm = [0; diff(positions)'] * 1000; % Spacing between elements
Dipole_Table = table(Dipole_Number, Length_mm, Spacing_mm);
disp('Dipole Lengths and Spacings:');
disp(Dipole_Table);
figure;
hold on;
for i = 1:length(lengths)
plot([-lengths(i)/2, lengths(i)/2], [positions(i), positions(i)], 'k', 'LineWidth', 2);
end
xlabel('Length (m)');
ylabel('Position along the boom (m)');
title('Log-Periodic Dipole Antenna Structure');
grid on;
axis equal;
hold off;

% Angles for the Radiation Pattern
theta = linspace(0, 2*pi, 360);
phi = linspace(0, pi, 180);

% 2D Radiation Pattern
D_lin = 10^(directivity_dB / 10); % Convert directivity to linear scale
U = D_lin * cos(theta).^2; % Radiation intensity

% Radiation Pattern in Polar Coordinates Plot
figure;
polarplot(theta, U, 'LineWidth', 1.5);
title('Radiation Pattern in Polar Coordinates');

% 3D Radiation Pattern
[Theta, Phi] = meshgrid(theta, phi);
R = D_lin * (cos(Theta).^2 .* sin(Phi)); % Radiation intensity

% Polar to Cartesian Coordinates
X = R .* sin(Phi) .* cos(Theta);
Y = R .* sin(Phi) .* sin(Theta);
Z = R .* cos(Phi);
figure;
surf(X, Y, Z, R, 'EdgeColor', 'none');
colorbar;
title('3D Radiation Pattern');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
axis equal;
view([45 45]);
grid on;

% Gain vs Frequency Plot
gain = 10 * log10(D_lin) * ones(size(f)); % Gain across frequency range
figure;
plot(f / 1e6, gain, 'LineWidth', 1.5);
xlabel('Frequency (MHz)');
ylabel('Gain (dB)');
title('Gain vs Frequency');
grid on;

% VSWR vs Frequency Plot
VSWR = (1 + abs(Gamma)) ./ (1 - abs(Gamma));
figure;
plot(f / 1e6, VSWR, 'LineWidth', 1.5);
xlabel('Frequency (MHz)');
ylabel('VSWR');
title('VSWR across Frequency Range');
grid on;

% Impedance vs Frequency Plot
figure;
plot(f / 1e6, real(Rin), 'LineWidth', 1.5); hold on;
plot(f / 1e6, imag(Rin), 'LineWidth', 1.5);
xlabel('Frequency (MHz)');
ylabel('Impedance (Ohms)');
title('Impedance vs Frequency');
legend('Real(Z)', 'Imag(Z)');
grid on;

% Directivity vs Frequency Plot
figure;
plot(f / 1e6, directivity_dB * ones(size(f)), 'LineWidth', 1.5);
xlabel('Frequency (MHz)');
ylabel('Directivity (dB)');
title('Directivity vs Frequency');
grid on;