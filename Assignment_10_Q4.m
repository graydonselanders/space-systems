%% Question 4
% Given values
f = 30e9; % Frequency in Hz (30 GHz)
lambda = 0.01; % Wavelength (m)

P_tx_watts = 5;  
D_sat = 0.4; % Antenna Diameter 
e_sat = 0.65; % Aperture efficiency
L_tx = 3; % Tx System loss 
theta_mis_tx= 0.2; % Tx mis pointing (degrees)

R_km = 1200; % Range (km)
R = R_km * 1000;  % Range (m)
L_atmo = 3;  % Atmospheric Loss (dB)

% Receiver (Ground Station)
D_gnd = 3.0; % Antenna Diameter (m)
e_gnd = 0.65; % Aperture Efficiency
L_rx = 3; % Rx System Loss (dB)
theta_mis_rx = 0.05; % Rx Mispointing (degrees)
Nf = 4; % Receiver Noise Figure (dB)

% Requirements
BER_req = 1e-4; % Target Bit Error Rate (10^-4)
Margin  = 3; % Required System Margin (dB)

% Power & Gains
P_tx_dBm = 10 * log10(P_tx_watts * 1000); % Convert Watts to dBm

% Antenna Gains
G_tx_lin = e_sat * (pi * D_sat / lambda)^2;
G_tx_dBi = 10 * log10(G_tx_lin);

G_rx_lin = e_gnd * (pi * D_gnd / lambda)^2;
G_rx_dBi = 10 * log10(G_rx_lin);

% Losses 
% Free Space Path Loss: L_fs = 20*log10(4*pi*R/lambda)
L_fs_dB = 20 * log10(4 * pi * R / lambda);

% Pointing Losses
theta_3dB_tx = 70 * lambda / D_sat;
theta_3dB_rx = 70 * lambda / D_gnd;

% Loss
L_point_tx = 12 * (theta_mis_tx / theta_3dB_tx)^2;
L_point_rx = 12 * (theta_mis_rx / theta_3dB_rx)^2;

% Total System Losses
L_total = L_tx + L_rx + L_fs_dB + L_atmo + L_point_tx + L_point_rx + Margin;

% Received Power
P_rx_dBm = P_tx_dBm + G_tx_dBi + G_rx_dBi - L_total;

% Noise Characterization
N0_dBm_Hz = -174; % Standard -174 dBm/Hz

% Effective System Noise Density (with Noise Figure)
N_sys_dBm_Hz = N0_dBm_Hz + Nf;

% Available Carrier-to-Noise Density (C/N0)
CN0_dBHz = P_rx_dBm - N_sys_dBm_Hz;

% Required SNR Calculation
sqrt_EbN0   = 3.719; %erfcinv(BER_req); %weirdly doesnt align with tables...
EbN0_req_lin = 1/2 * sqrt_EbN0^2;
EbN0_req_dB  = 10 * log10(EbN0_req_lin);

% Data Rate Solution
Rb_dBHz = CN0_dBHz - EbN0_req_dB;

% Converting to linear bits per second
Rb_bps = 10^(Rb_dBHz / 10);
Rb_Gbps = Rb_bps / 1e9;

% outputs
fprintf('\n---------- Question 4 ----------\n');
fprintf('Frequency:              %.2f GHz\n', f/1e9);
fprintf('Slant Range:            %.0f km\n', R_km);
fprintf('--------------------------------------------------\n');
fprintf('Tx Antenna Gain:        %.2f dBi (Beamwidth: %.2f deg)\n', G_tx_dBi, theta_3dB_tx);
fprintf('Rx Antenna Gain:        %.2f dBi (Beamwidth: %.2f deg)\n', G_rx_dBi, theta_3dB_rx);
fprintf('--------------------------------------------------\n');
fprintf('Free Space Path Loss:   %.2f dB\n', L_fs_dB);
fprintf('Pointing Loss (Tx+Rx):  %.2f dB\n', L_point_tx + L_point_rx);
fprintf('Atmospheric Loss:       %.2f dB\n', L_atmo);
fprintf('--------------------------------------------------\n');
fprintf('Received Power (Prx):   %.2f dBm\n', P_rx_dBm);
fprintf('System Noise Floor:     %.2f dBm/Hz\n', N_sys_dBm_Hz);
fprintf('Available C/N0:         %.2f dB-Hz\n', CN0_dBHz);
fprintf('--------------------------------------------------\n');
fprintf('Required Eb/N0 (BPSK):  %.2f dB (for BER %.0e)\n', EbN0_req_dB, BER_req);
fprintf('--------------------------------------------------\n');
fprintf('MAX ACHIEVABLE DATA RATE: %.4f Gbps\n\n', Rb_Gbps);
