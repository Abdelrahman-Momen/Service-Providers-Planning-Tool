function [N, cells, R_cell, A_cell, A_sect, Ptx, Pr, d] = ...
    planning_tool(gos, area, user_density, sir_min_dB, sectors)

%% constant parameters:
total_users = user_density * area;
total_channels = 340;
Au = 0.025;             % traffic per user
At = total_users * Au;  % total traffic
f = 900;                % freq. in MHz
ht = 20;                % Tx hight
hr = 1.5;               % Rx hight
Pr_min = -95;           % Rx_sensitivity in dBm
n = 4;                  % path loss exponent

% import json standard cluster sizes (N)
N_std = jsondecode(fileread('cluster_sizes.json'));
% disp(N_std)

% import modified csv erlangB table (contains numbers only)
erlangb_table = readmatrix('ErlangB_table.csv');
traffics = erlangb_table(2:end, 2:end); % erlang values
gos_table = erlangb_table(1, 2:end);    % gos values
gos_idx = find(gos_table == gos, 1);    % input gos index
if isempty(gos_idx)
    error('GOS value %.2f is not found in the ErlangB table.', gos);
end

%% cluster size based on SIR:
sir_linear = 10^(sir_min_dB/10);

% define the no. of interferers w.r.t no. of sectors
ni_map = containers.Map([1, 3, 6], [6, 2, 1]);
ni = ni_map(sectors);
N_req = (1/3) * ((ni * sir_linear)^(1/n)+1)^2; % fraction -> refused
N = N_std(find(N_std >= N_req, 1));            % std. cluster size


%% cell / sector parameters:
cell_channels = floor(total_channels/N);        % channels per cell
sect_channels = floor(cell_channels/sectors);   % channels per sector

A_sect = traffics(sect_channels, gos_idx);      % sector traffic intensity
A_cell = A_sect * sectors;                      % cell traffic intensity

cells = ceil(At / A_cell);                      % total no. of cells
cell_area = area / cells;                       % cell area in km^2
R_cell = sqrt(cell_area / (1.5 * sqrt(3)));     % cell radius

%% Hata model:
d = linspace(0.1, R_cell, 100);

a_hr_dB = (1.1 * log10(f) - 0.7) * hr - (1.56 * log10(f) - 0.8);

L_dB = 69.55 + 26.16 * log10(f) - 13.82 * log10(ht) - a_hr_dB + ...
    (44.9 - 6.55 * log10(ht)) * log10(d);

Ptx = Pr_min + L_dB(end);

Pr = Ptx - L_dB;
end