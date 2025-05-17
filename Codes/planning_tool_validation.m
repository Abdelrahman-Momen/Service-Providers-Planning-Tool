%% initial parameters:
clc; clear; close all;

area = 100;             % city area
user_density = 1400;    % user density/km^2
sectorisation = [1, 3, 6];  % no. of sectors w.r.t sectorisation method
sectorisation_methods = {'omni', '120°', '60°'};

%% point 1: cluster size vs SIR min
sir_min_range = 1:1:30; % sir min in db
gos = 0.2;              % assumed gos (not given)

% 2d array for N values w.r.t each sectorisation method
N_cells = zeros(length(sectorisation), length(sir_min_range)); 

for i = 1:length(sectorisation)
    sectors = sectorisation(i);
    for j = 1:length(sir_min_range)
        sir_min_dB = sir_min_range(j);
        try
            [N, ~, ~, ~, ~, ~, ~, ~] = planning_tool(gos, area, ...
                user_density, sir_min_dB, sectors);
            N_cells(i, j) = N;
        catch
            N_cells(i, j) = NaN;    % put nan if calc. fails
            fprintf('error calculating N');
        end
    end
end

figure;
hold on;
for i = 1:length(sectorisation)
    plot(sir_min_range, N_cells(i, :));
end
hold off;
xlabel('SIR_{min} (dB)');
ylabel('Cluster Size (N)');
title('Cluster Size vs SIR_{min}');
legend(sectorisation_methods);
grid on;

%% point 2: no. of cells and cell traffic vs gos @ sir = 19
sir_min_dB = 19;
gos_range = 1:0.1:30;
N_cells = zeros(length(sectorisation), length(gos_range));
A_cell_vals = zeros(length(sectorisation), length(gos_range));

for i = 1:length(sectorisation)
    sectors = sectorisation(i);
    for j = 1:length(gos_range)
        gos = gos_range(j);
        try
            [~, cells, ~, A_cell, ~, ~, ~, ~] = planning_tool(gos, area, ...
                user_density, sir_min_dB, sectors);
            N_cells(i, j) = cells;
            A_cell_vals(i, j) = A_cell;
        catch
            N_cells(i, j) = NaN;    % put nan if calc. fails
            A_cell_vals(i, j) = NaN;
            %fprintf('error finding gos');
        end
    end
end

valid_idx = find(~isnan(A_cell_vals(1,:))); % valid indices to ignore nans

figure;
hold on;
for i = 1:length(sectorisation)
    plot(gos_range(valid_idx), N_cells(i, valid_idx));
end
hold off;
xlabel('Grade of Service (GOS) (%)');
ylabel('Total Number of Cells (N_{cell})');
title('N_{cell} vs. GOS (SIR_{min}=19 dB, Density=1400)');
legend(sectorisation_methods);
grid on;

figure;
hold on;
for i = 1:length(sectorisation)
    plot(gos_range(valid_idx), A_cell_vals(i, valid_idx));
end
hold off;
xlabel('Grade of Service (GOS) (%)');
ylabel('Traffic Intensity per Cell (A_{cell}) (Erlangs)');
title('A_{cell} vs. GOS (SIR_{min}=19 dB, Density=1400)');
legend(sectorisation_methods);
grid on;

%% point 3: no. of cells and cell traffic vs gos @ sir = 14
sir_min_dB = 14;
gos_range = 1:0.1:30;
N_cells = zeros(length(sectorisation), length(gos_range));
A_cell_vals = zeros(length(sectorisation), length(gos_range));

for i = 1:length(sectorisation)
    sectors = sectorisation(i);
    for j = 1:length(gos_range)
        gos = gos_range(j);
        try
            [~, cells, ~, A_cell, ~, ~, ~, ~] = planning_tool(gos, area, ...
                user_density, sir_min_dB, sectors);
            N_cells(i, j) = cells;
            A_cell_vals(i, j) = A_cell;
        catch
            N_cells(i, j) = NaN;    % put nan if calc. fails
            A_cell_vals(i, j) = NaN;
            %fprintf('error finding gos');
        end
    end
end

figure;
hold on;
for i = 1:length(sectorisation)
    plot(gos_range(valid_idx), N_cells(i, valid_idx));
end
hold off;
xlabel('Grade of Service (GOS) (%)');
ylabel('Total Number of Cells (N_{cell})');
title('N_{cell} vs. GOS (SIR_{min}=14 dB, Density=1400)');
legend(sectorisation_methods);
grid on;

figure;
hold on;
for i = 1:length(sectorisation)
    plot(gos_range(valid_idx), A_cell_vals(i, valid_idx));
end
hold off;
xlabel('Grade of Service (GOS) (%)');
ylabel('Traffic Intensity per Cell (A_{cell}) (Erlangs)');
title('A_{cell} vs. GOS (SIR_{min}=14 dB, Density=1400)');
legend(sectorisation_methods);
grid on;

%% point 4: no. of cells and cell radius vs density @ sir = 14, gos = 2
sir_min_dB = 14;
gos = 2;
density_range = 100:100:2000; % users/km^2
N_cells = zeros(length(sectorisation), length(density_range));
R_cell_vals = zeros(length(sectorisation), length(density_range));

for i = 1:length(sectorisation)
    sectors = sectorisation(i);
    for j = 1:length(density_range)
        user_density = density_range(j);
        try
            [~, cells, R_cell, ~, ~, ~, ~, ~] = planning_tool(gos, area, ...
                user_density, sir_min_dB, sectors);
            N_cells(i, j) = cells;
            R_cell_vals(i, j) = R_cell;
        catch
            N_cells(i, j) = NaN;    % put nan if calc. fails
            R_cell_vals(i, j) = NaN;
            %fprintf('error -> density');
        end
    end
end

figure;
hold on;
for i = 1:length(sectorisation)
    plot(density_range, N_cells(i, :));
end
hold off;
xlabel('User Density (users/km^2)');
ylabel('Total Number of Cells (N_{cell})');
title('N_{cell} vs. User Density (SIR_{min}=14 dB, GOS=2%)');
legend(sectorisation_methods);
grid on;

figure;
hold on;
for i = 1:length(sectorisation)
    plot(density_range, R_cell_vals(i, :));
end
hold off;
xlabel('User Density (users/km^2)');
ylabel('Cell Radius (R) (km)');
title('Cell Radius vs. User Density (SIR_{min}=14 dB, GOS=2%)');
legend(sectorisation_methods);
grid on;

%% point 5: no. of cells and cell radius vs density @ sir = 19, gos = 2
sir_min_dB = 19;
gos = 2;
density_range = 100:100:2000; % users/km^2
N_cells = zeros(length(sectorisation), length(density_range));
R_cell_vals = zeros(length(sectorisation), length(density_range));

for i = 1:length(sectorisation)
    sectors = sectorisation(i);
    for j = 1:length(density_range)
        user_density = density_range(j);
        try
            [~, cells, R_cell, ~, ~, ~, ~, ~] = planning_tool(gos, area, ...
                user_density, sir_min_dB, sectors);
            N_cells(i, j) = cells;
            R_cell_vals(i, j) = R_cell;
        catch
            N_cells(i, j) = NaN;    % put nan if calc. fails
            R_cell_vals(i, j) = NaN;
            %fprintf('error -> density');
        end
    end
end

figure;
hold on;
for i = 1:length(sectorisation)
    plot(density_range, N_cells(i, :));
end
hold off;
xlabel('User Density (users/km^2)');
ylabel('Total Number of Cells (N_{cell})');
title('N_{cell} vs. User Density (SIR_{min}=19 dB, GOS=2%)');
legend(sectorisation_methods);
grid on;

figure;
hold on;
for i = 1:length(sectorisation)
    plot(density_range, R_cell_vals(i, :));
end
hold off;
xlabel('User Density (users/km^2)');
ylabel('Cell Radius (R) (km)');
title('Cell Radius vs. User Density (SIR_{min}=19 dB, GOS=2%)');
legend(sectorisation_methods);
grid on;

