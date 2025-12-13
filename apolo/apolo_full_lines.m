clear; close all; clc; rng(0);

load('gardenMap.mat');
map = binaryOccupancyMap(imrotate(garden(:,:,1), 180), "Resolution", 10);
trajectories = load('trayectorias.mat');

traj = trajectories.trajDirect;
start = traj(1,:);

load('map/garden_lines.mat', 'map_lines');
num_walls = size(map_lines, 1);

%% Initialize Robot Interface
robot = Apolo("MazeRunner.xml", "Dafne", "LMS100");

%% Vehicle and Simulation Parameters
dt = 0.1;                       % Simulation time step [s]
control_dt = 0.1;                % Control time step [s]
control_rate = control_dt / dt;  % Control update every N steps

controller = NMPCCBFController(...
    'HorizonLength', 20, ...
    'TimeStep', control_dt, ...
    'StateWeights', [1000, 1000, 100], ...
    'ControlWeights', [1, 1], ...
    'VelocityLimits', [0, 1], ...
    'AngularLimits', [-2, 2], ...
    'SafetyRadius', 0.5, ...
    'AlphaCBF', 0.1, ...
    'ScanDownsample', 20, ...
    'ConstraintRange', 3.0, ...
    'MaxIterations', 100, ...
    'UseSlack', true, ...
    'SlackPenalty', 100);

%% Interpolate Trajectory
% Calculate trajectory length by calculating the euclidean distance
% between original points
traj_sparse = traj; traj_length = 0;
for i = 1:size(traj_sparse, 1)-1
    traj_length = traj_length + norm(traj_sparse(i+1, 1:2) - traj_sparse(i, 1:2));
end

% Estimate traversal time based on desired average speed
pause_time = 0;                  % For "realism" use dt [s]
avg_speed = 1.0;  % Desired average speed [m/s]
traj_time = traj_length / avg_speed;
simulation_time = traj_time;
num_steps = round(simulation_time / dt);

% Create dense time vector at EKF rate
num_traj_points = round(traj_time / dt);
t_sparse = linspace(0, 1, size(traj_sparse, 1));  % Normalized time for sparse traj
t_dense = linspace(0, 1, num_traj_points);        % Normalized time for dense traj

% Interpolate trajectory using pchip
% https://www.mathworks.com/help/matlab/ref/pchip.html
traj = zeros(num_traj_points, 3);
traj(:, 1) = interp1(t_sparse, traj_sparse(:, 1), t_dense, 'pchip');  % x
traj(:, 2) = interp1(t_sparse, traj_sparse(:, 2), t_dense, 'pchip');  % y
traj(:, 3) = interp1(t_sparse, traj_sparse(:, 3), t_dense, 'pchip');  % θ

fprintf('Trajectory interpolation:\n');
fprintf('  Sparse waypoints:  %d points\n', size(traj_sparse, 1));
fprintf('  Dense trajectory:  %d points\n', size(traj, 1));
fprintf('  Trajectory length: %.2f m\n', traj_length);
fprintf('  Est. traversal:    %.1f s @ %.1f m/s avg\n\n', traj_time, avg_speed);

% Start at first trajectory point
pathPoint = 1;

%% Beacon Positions
beacons = robot.getBeaconPositions();
if isempty(beacons)
    disp("Beacons are not present")
    return;
end
num_beacons = size(beacons, 1);
num_measurements = 2 * num_beacons;  % Range + bearing per beacon

%% Initial Conditions
robot.reset(start);
x = robot.getState();

% For EKF (in odometry space: Δd, Δβ)
process_noise_d = 9.5003e-05;      % Distance increment noise [m]
process_noise_beta = 3.9080e-05;   % Heading increment noise [rad]
Q = diag([process_noise_d^2, process_noise_beta^2]);

% Precompute standard deviations
Q_std = sqrt(diag(Q));

% Range and bearing to beacons noise (extracted from Apolo sensor calibration)
measurement_noise_range = 0.018085189925279;   % Range measurement noise [m]
measurement_noise_bearing = 0.023174091647608;  % Bearing measurement noise [rad]

% Measurement covariance matrix is 2x2
% [σ^2_r 0; 0 σ^2_α]
R = diag([measurement_noise_range^2, measurement_noise_bearing^2]);

% Precompute standard deviations 
R_std = sqrt(diag(R));

% Initialize estimated state with certain deviations
x0 = x + [0.05; 0.05; 0.01];

% Initialize EKF
P = diag([0.1^2, 0.1^2, 0.01^2]);
% Create EKF instance
chi2_threshold = 9.21;  % 99% confidence for 2 DOF
ekf = EKFLines(x0, P, map_lines, Q, chi2_threshold);

%% Variable storage
true_trajectory = zeros(3, num_steps);
estimated_trajectory = zeros(3, num_steps);
deviation_history = zeros(3, num_steps);
control_history = zeros(2, num_steps);

%% Main Simulation Loop
fprintf('Starting simulation with multi-rate control:\n');
fprintf('  EKF rate:     %.0f Hz (dt = %.3f s)\n', 1/dt, dt);
fprintf('  Control rate: %.0f Hz (dt = %.3f s)\n', 1/control_dt, control_dt);
fprintf('  Total time:   %.1f s (%d steps)\n\n', simulation_time, num_steps);

% Initialize control to zero
u_current = [0; 0];

for k = 1:num_steps
    %% Get Lidar scan
    scan = robot.getLaserScan();

    %% Control Update
    if mod(k-1, control_rate) == 0
        % Compute control using remaining trajectory
        remaining_traj = traj(k:end, :);
        u_current = controller.compute(ekf.x, remaining_traj, scan);
    end

    %% Robot Movement
    % Keep applying the current control signal if no
    % new value has been calculated (zero-order hold)
    robot.move(u_current, dt, pause_time);
    x = robot.getState();
    true_trajectory(:, k) = x;

    %% EKF Prediction
    [delta_d_hat, delta_beta_hat] = robot.getOdometry();
    ekf.predict(delta_d_hat, delta_beta_hat);

    %% EKF Update
    % Extract lines from scan using RANSAC
    lines_observed = ransac_lines(scan, 0.015, 3);
    ekf.update(lines_observed);

    %% Store Results
    estimated_trajectory(:, k) = ekf.x;
    deviation_history(:, k) = sqrt(diag(ekf.P));
    control_history(:, k) = u_current;
end

% Adapt storage in case we have fewer samples
if k < num_steps
    true_trajectory = true_trajectory(:, 1:k);
    estimated_trajectory = estimated_trajectory(:, 1:k);
    deviation_history = deviation_history(:, 1:k);
    control_history = control_history(:, 1:k);
end

fprintf('\nSimulation complete! Final step: %d/%d\n', k, num_steps);

%% Create Time Vector
t = (0:k-1) * dt;  % Time vector at simulation rate

%% Compute Estimation Errors
position_error = sqrt(sum((true_trajectory(1:2,:) - estimated_trajectory(1:2,:)).^2, 1));
angle_error = abs(true_trajectory(3,:) - estimated_trajectory(3,:));

%% TODO: Compute Tracking Error

%% Visualization
figure('Name', 'Maze Solver Full Stack', 'Position', [50 50 1400 900]);

% Plot 1: 2D Trajectory with Beacons
subplot(2,3,1);
hold on; grid on; axis equal;
plot(true_trajectory(1,:), true_trajectory(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'True');
plot(estimated_trajectory(1,:), estimated_trajectory(2,:), 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated');
plot(true_trajectory(1,1), true_trajectory(2,1), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
xlabel('X [m]'); ylabel('Y [m]');
title('Trajectory');
legend('Ground truth', 'Estimated', 'Initial Position', 'Location', 'best');

% Plot 2: X Position
subplot(2,3,2);
plot(t, true_trajectory(1,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, estimated_trajectory(1,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('X [m]');
title('X Position'); legend('Ground truth', 'Estimated');
grid on;

% Plot 3: Y Position
subplot(2,3,3);
plot(t, true_trajectory(2,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, estimated_trajectory(2,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Y [m]');
title('Y Position'); legend('Ground truth', 'Estimated');
grid on;

% Plot 4: Heading Angle
subplot(2,3,4);
plot(t, rad2deg(true_trajectory(3,:)), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, rad2deg(estimated_trajectory(3,:)), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Heading [deg]');
title('Heading Angle'); legend('Ground truth', 'Estimated');
grid on;

% Plot 5: Estimation and Tracking Errors
subplot(2,3,5);
plot(t, position_error, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Estimation Error (EKF vs True)');
hold on;
xlabel('Time [s]'); ylabel('Error [m]');
title('Position Errors');
legend('Location', 'best');
grid on;

% Plot 6: Covariance
subplot(2,3,6);
semilogy(t, deviation_history(1,:), 'r-', 'LineWidth', 1.5); hold on;
semilogy(t, deviation_history(2,:), 'g-', 'LineWidth', 1.5);
semilogy(t, deviation_history(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Standard Deviation');
title('State Uncertainty');
legend('\sigma_x', '\sigma_y', '\sigma_\theta');
grid on;

% Plot 7: Control Inputs
figure('Name', 'Control Signals', 'Position', [50 50 800 600]);
subplot(2,1,1);
plot(t, control_history(1,:), 'b-', 'LineWidth', 1.5);
ylabel('v [m/s]'); xlabel('Time [s]');
title('Linear Velocity');
grid on;

subplot(2,1,2);
plot(t, control_history(2,:), 'r-', 'LineWidth', 1.5);
ylabel('\omega [rad/s]'); xlabel('Time [s]');
title('Angular Velocity');
grid on;

figure("Name","Trajectories");
hold on; grid on; axis equal;
show(map);
plot(traj(:,1), traj(:,2), ':pentagramy', 'LineWidth', 2, 'DisplayName', 'Planned');
plot(true_trajectory(1,:), true_trajectory(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'True');
plot(estimated_trajectory(1,:), estimated_trajectory(2,:), 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated');
plot(true_trajectory(1,1), true_trajectory(2,1), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot(waypoints(:,1), waypoints(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'cyan');
xlabel('X [m]'); ylabel('Y [m]');
title('Trajectory');
legend('Planned', 'Ground truth', 'Estimated', 'Initial Position', 'Waypoints', 'Location', 'best');

