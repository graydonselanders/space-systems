%% Question 5
c = 299792458;          % Speed of light (m/s)
h = 6.62607015e-34;     % Planck's constant (J s)

R_km = 1200;
R = R_km * 1000; % Range in meters
lambda_nm = 1550;
lambda = lambda_nm * 1e-9; % Wavelength in meters
freq = c / lambda; % Frequency in Hz (approx 193.4 THz)

P_tx_watts = 5;
P_tx_dBm = 10 * log10(P_tx_watts * 1000); % Convert 5W to dBm (~37 dBm)

D_sat = 0.1; % Satellite Telescope Diameter (m)
e_A_sat = 0.81; % Satellite Aperture Efficiency
L_Tx = 3; % Satellite Transmit Path Loss (dB)
Delta_theta_Tx = 3e-6; % Mispointing Error (radians)

L_atmo = 30; % Atmospheric Loss (dB)

D_gnd = 0.5; % Ground Station Telescope Diameter (m)
e_A_gnd = 0.81; % Ground Aperture Efficiency
L_Rx = 8; % Ground Receive Path Loss (dB)
N_f = 3; % Ground Receiver Noise Figure (dB)

BER_req = 1e-10; % Target Bit Error Rate
Margin = 3; % Required System Margin (dB)

% Gain G = (pi * D / lambda)^2 * efficiency
% Optical antennas are diffraction limited
G_Tx_lin = ((pi * D_sat) / lambda)^2 * e_A_sat;
G_Tx_dB = 10 * log10(G_Tx_lin);

G_Rx_lin = ((pi * D_gnd) / lambda)^2 * e_A_gnd;
G_Rx_dB = 10 * log10(G_Rx_lin);

% Free Space Path Loss: L_FS = (4 * pi * R / lambda)^2
L_FS_lin = ((4 * pi * R) / lambda)^2;
L_FS_dB = 10 * log10(L_FS_lin);

% Pointing Loss (Gaussian Beam Approximation)
% Beam waist w0 is approx D/2 for optimal truncation
w0 = D_sat / 2;
% Divergence half-angle theta_div
theta_div = lambda / (pi * w0);
% Loss formula: L(dB) = 8.686 * (Error / Divergence)^2
L_P_dB = 8.686 * (Delta_theta_Tx / theta_div)^2;

% P_Rx = P_Tx + G_Tx - L_Tx - L_FS - L_P - L_atmo + G_Rx - L_Rx - Margin
P_Rx_dBm = P_tx_dBm + G_Tx_dB - L_Tx - L_FS_dB - L_P_dB - L_atmo +...
           G_Rx_dB - L_Rx - Margin;

% BER = Q(sqrt(2 * SNR_lin))
% invert the Q-function. In MATLAB, Q(x) = 0.5 * erfc(x/sqrt(2))
% Therefore: x = sqrt(2) * erfcinv(2 * BER)
% And since x = sqrt(2 * SNR), then SNR = x^2 / 2
x_req = sqrt(2) * erfcinv(2 * BER_req);
SNR_lin_req = (x_req^2) / 2;
SNR_dB_req = 10 * log10(SNR_lin_req);

% Using the optical sensitivity equation:
% SNR_dB = P_Rx_dBm + 159 - 10*log10(Rs) - N_f
% Solving for Rs (Symbol Rate):
Log_Rs = P_Rx_dBm + 159 - N_f - SNR_dB_req;
Rs = 10^(Log_Rs / 10);

% For BPSK, Data Rate = Symbol Rate
Data_Rate_Gbps = Rs / 1e9;

fprintf('\n---------- Question 5 ----------\n');
fprintf('Wavelength: %.2f nm | Frequency: %.2f THz\n', lambda_nm, freq/1e12);
fprintf('Range: %.0f km\n', R_km);
fprintf('\nGAIN & LOSS ANALYSIS:\n');
fprintf('  Tx Gain:          %8.2f dB\n', G_Tx_dB);
fprintf('  Rx Gain:          %8.2f dB\n', G_Rx_dB);
fprintf('  Free Space Loss:  %8.2f dB\n', L_FS_dB);
fprintf('  Pointing Loss:    %8.2f dB (Error: %.1f urad / Div: %.1f urad)\n',...
        L_P_dB, Delta_theta_Tx*1e6, theta_div*1e6);
fprintf('  Atmospheric Loss: %8.2f dB\n', L_atmo);
fprintf('\nPOWER BUDGET:\n');
fprintf('  Tx Power:         %8.2f dBm\n', P_tx_dBm);
fprintf('  Rx Power (Net):   %8.2f dBm\n', P_Rx_dBm);
fprintf('\nPERFORMANCE METRICS:\n');
fprintf('  Required SNR:     %8.2f dB (for BER 1e-10)\n', SNR_dB_req);
fprintf('  Max Symbol Rate:  %8.2f Gsym/s\n', Rs/1e9);
fprintf('------------------------------------------------\n');
fprintf('  MAX DATA RATE:    %8.4f Gbps\n', Data_Rate_Gbps);
fprintf('------------------------------------------------\n');
