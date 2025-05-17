clc; clear; close all;

%% user input parameters:
gos = input('Enter GOS in %: ');
area = input('Enter city area in km^2: ');
user_density = input('Enter user density (users/km^2): ');
sir_min_dB = input('Enter SIR min in dB: ');
sectors = input('Enter sectorization method (1, 3, or 6): ');

% insure valide sectorisation method:
while sectors ~= 1 && sectors ~= 3 && sectors ~= 6
    sectors = input('Invalid sectorization. Enter 1, 3, or 6: ');
end

%% planning tool outputs:
[N, cells, R_cell, A_cell, A_sector, Ptx, Pr, d] = ...
    planning_tool(gos, area, user_density, sir_min_dB, sectors);

fprintf('\n----- Planning Tool Outputs -----\n')
fprintf('1. Cluster Size N = %d\n', N);
fprintf('2. Total Cells = %d\n', cells);
fprintf('3. Cell Radius = %.2f km\n', R_cell);
fprintf('4. Traffic per Cell = %.2f Erlang\n', A_cell);
fprintf('   Traffic per Sector = %.2f Erlang\n', A_sector);
fprintf('5. BS Transmit Power = %.2f dBm\n', Ptx);

figure;
plot(d, Pr);
xlabel('Distance from BS (km)');
ylabel('Received Power (dBm)');
title('MS Received Power vs Distance');
grid on;
