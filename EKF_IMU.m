close all

%% Load data
truth_data = load('Datasets\euler(xaxismove).log');
raw_data = load('Datasets\raw(xaxismove).log');
acc = raw_data(:, 2:4);
gyro = raw_data(:, 5:7);
mag = raw_data(:, 8:10);
mag  = mag / norm(mag);
timestamp = raw_data(:, 1);
Ts = timestamp(2) - timestamp(1);
num_samples = length(timestamp);

%% Define process and measurement noises
Qd = Ts^2*diag([zeros(1, 4), 1e2*ones(1, 3), 4e-8*ones(1, 3), 1e-2*ones(1, 3), 9.95e-10*ones(1, 3)]);
Rd = diag([[0.02, 0.005, 0.1], [0.002, 0.005, 5.0], [0.00005, 0.0002, 0.00001]]);

%% Create containers for data storage
state = zeros(num_samples, 6);
state(:, 4:6) = truth_data(:, 2:4)*pi/180;
state(1, 1:3) = state(1, 4:6);

%% Apply Extended Kalman Filter
initial = eul2quat([state(1, 6), state(1, 5), state(1, 4)]);
if initial(1) < 0
    initial = -initial;
end
mu_corr = [[initial(2); initial(3); initial(4); initial(1)]; zeros(12, 1)];
sigma_corr = 1e-5*eye(16);
q_corr = mu_corr(1:4);
for i = 2:num_samples
    % Perform prediction
    mu_pred = process_model(mu_corr, Ts);
    F = jacobian_process(mu_corr, Ts);
    sigma_pred = F*sigma_corr*F' + Qd;

    % Perform correction
    z = [acc(i, :), mag(i, :), gyro(i, :)]';
    H = jacobian_measurement(mu_pred);
    K = sigma_pred*H'/(H*sigma_pred*H' + Rd);
    mu_corr = mu_pred + K*(z - measurement_model(mu_pred));
    mu_corr(1:4) = mu_corr(1:4)/norm(mu_corr(1:4));
    if mu_corr(4) < 0
        mu_corr(1:4) = -mu_corr(1:4);
    end
    sigma_corr = sigma_pred - K*H*sigma_pred;
    
    % Trick to boost accuracy
    phi = atan2(acc(i, 2), acc(i, 3));
    theta = atan2(-acc(i, 1), norm(acc(i, 2:3)));
    psi = atan2(-mag(i, 2)*cos(phi) + mag(i, 3)*sin(phi), ...
                 mag(i, 1)*cos(theta) + mag(i, 2)*sin(theta)*sin(phi)+ mag(i, 3)*sin(theta)*cos(phi));
    q_tilde = eul2quat([psi, theta, phi]);
    q_tilde = [q_tilde(2:4) q_tilde(1)]';
    if q_tilde(4) < 0
        q_tilde = -q_tilde;
    end
    mu_corr(1:4) = mu_pred(1:4) + 0.009*(q_tilde - q_corr);
    mu_corr(1:4) = mu_corr(1:4)/norm(mu_corr(1:4));
    if mu_corr(1:4)<0
        mu_corr(1:4) = -mu_corr(1:4);
    end
    q_corr = mu_corr(1:4);
    
    % Save data for evaluation
    euler_angles = quat2eul([mu_corr(4) mu_corr(1:3)']);
    state(i, 1) = euler_angles(3);
    state(i, 2) = euler_angles(2);
    state(i, 3) = euler_angles(1);
end

%% Visualize results
figure(1)
figure(2)
for i = 1:3
    % Get angle's name
    switch i
        case 1
            angle = 'roll';
        case 2
            angle = 'pitch';
        case 3
            angle = 'yaw';
    end
    
    % Plot state of estimated value and ground truth
    figure(1)
    subplot(1, 3, i);
    plot(timestamp, state(:, i)*180/pi, 'b')
    hold on
    plot(timestamp, state(:, i + 3)*180/pi, 'r')
    title(strcat(angle, ' (deg)'))
    legend('Estimated', 'True')
    
    % Plot error
    figure(2)
    subplot(1, 3, i)
    error = normalize_angle(state(:, i) - state(:, i + 3))*180/pi;
    plot(timestamp, error, 'b');
    title(strcat(angle, ' (deg)'))
    rmse = rms(error);
    title(strcat(angle, ', rmse = ', num2str(rmse), ' deg'))
    disp(['rmse_', angle, ' = ', num2str(rmse)])
end