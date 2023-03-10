imu_data_bag = rosbag('LocationD.bag');
imu_data_topic = select(imu_data_bag,'Topic','/vectornav');
imu_data = readMessages(imu_data_topic,'DataFormat','struct');

%% Take gyro data
for i = 1:length(imu_data)
    gyro_dat = split(imu_data{i}.Data, ',');
    %disp(imu_data{i});
    %disp(i);
    if length(gyro_dat) == 13
        gyro_x = str2double(gyro_dat{11});
        gyro_y = str2double(gyro_dat{12});
        temp_gyro_z = split(gyro_dat{13}, '*');
        gyro_z = str2double(temp_gyro_z{1});
        omega(i,:) = [gyro_x, gyro_y, gyro_z];
    end
     if  length(gyro_dat)<=13
         gyro_dat = zeros(13,1);
     end
end
%%%Gyro Time-Series Plot


figure('Name',"Time Series Plot of Gyro", 'NumberTitle', 'off')
subplot(2,2,1)
plot(linspace(1,length(imu_data)/40,length(imu_data)), omega(:,1),'b+')
legend ('Gyro_X')	
title('Gyro Plot X')	
xlabel('Time (seconds)')
ylabel('Gyro (rad/s)')
hold on
subplot(2,2,2)	
plot(linspace(1,length(imu_data)/40,length(imu_data)), omega(:,2),'g+')
legend ('Gyro_Y')	
title('Gyro Plot Y')	
xlabel('Time (seconds)')
ylabel('Gyro (rad/s)')
hold on
subplot(2,2,3)
plot(linspace(1,length(imu_data)/40,length(imu_data)), omega(:,3),'r+')
legend ('Gyro_Z')
title('Gyro Plot Z')	
xlabel('Time (seconds)')
ylabel('Gyro (rad/s)')
hold on
subplot(2,2,4)
plot(linspace(1,length(imu_data)/40,length(imu_data)), omega(:,1),'b+')
hold on
plot(linspace(1,length(imu_data)/40,length(imu_data)), omega(:,2),'g+')
hold on
plot(linspace(1,length(imu_data)/40,length(imu_data)), omega(:,3),'r+')
grid on
legend ('Gyro_X, Gyro_Y, Gyro_Z')
title('Time Series Gyro Plot X Y X')	
xlabel('Time (seconds)')
ylabel('Gyro (rad/s)')
hold off


% Compute the Allan variance for gyro x
Fs = 40; % Sample rate in Hz
t0 = 1/Fs;
theta = cumsum(omega(:,1), 1)*t0;
maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.
tau = m*t0;
[avarFromFunc, tauFromFunc] = allanvar(omega(:,1), m, Fs);
adevFromFunc = sqrt(avarFromFunc);
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N_GyroX = 10^logN
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K_GyroX = 10^logK
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tauFromFunc); %tau = tauFromFunc
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B_GyroX = 10^logB
% Plot the results.
tauN = 1;
lineN = N_GyroX ./ sqrt(tau);
% Plot the results.
tauK = 3;
lineK = K_GyroX .* sqrt(tau/3);
% Plot the results.
tauB = tau(i);
lineB = B_GyroX * scfB * ones(size(tau));
% Plot results
figure('Name',"Allan Deviation of Gyro X", 'NumberTitle', 'off')
loglog(tauFromFunc, adevFromFunc, tauFromFunc, lineN, '--', tauN, N_GyroX, 'o', tauFromFunc, lineK, '--', tauK, K_GyroX, 'o', tau, lineB, '--', tauB, scfB*B_GyroX, 'o');
legend('\sigma (rad/s)', '\sigma_N ((rad/s)sqrt{Hz})', '','\sigma_K ((rad/s)sqrt{Hz})','','\sigma_B (rad/s)')
title('Allan Deviation of Gyro X')
xlabel('\tau')
ylabel('\sigma(\tau)')
text(tauN, N_GyroX, 'N')
text(tauK, K_GyroX, 'K')
text(tauB, scfB*B_GyroX, '0.664B')
grid on
axis equal



%% 
% Compute the Allan variance for gyro y
Fs = 40; % Sample rate in Hz
t0 = 1/Fs;
theta = cumsum(omega(:,2), 1)*t0;
maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.
tau = m*t0;
[avarFromFunc, tauFromFunc] = allanvar(omega(:,2), m, Fs);
adevFromFunc = sqrt(avarFromFunc);
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N_GyroY = 10^logN
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K_GyroY = 10^logK
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tauFromFunc); %tau = tauFromFunc
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B_GyroY = 10^logB
% Plot the results.
tauN = 1;
lineN = N_GyroY ./ sqrt(tau);
% Plot the results.
tauK = 3;
lineK = K_GyroY .* sqrt(tau/3);
% Plot the results.
tauB = tau(i);
lineB = B_GyroY * scfB * ones(size(tau));
% Plot results
figure('Name',"Allan Deviation of Gyro Y", 'NumberTitle', 'off')
loglog(tauFromFunc, adevFromFunc, tauFromFunc, lineN, '--', tauN, N_GyroY, 'o', tauFromFunc, lineK, '--', tauK, K_GyroY, 'o', tau, lineB, '--', tauB, scfB*B_GyroY, 'o');
legend('\sigma (rad/s)', '\sigma_N ((rad/s)sqrt{Hz})', '','\sigma_K ((rad/s)sqrt{Hz})','','\sigma_B (rad/s)')
title('Allan Deviation of Gyro Y')
xlabel('\tau')
ylabel('\sigma(\tau)')
text(tauN, N_GyroY, 'N')
text(tauK, K_GyroY, 'K')
text(tauB, scfB*B_GyroY, '0.664B')
grid on
axis equal



%% 
% Compute the Allan variance for gyro z
Fs = 40; % Sample rate in Hz
t0 = 1/Fs;
theta = cumsum(omega(:,3), 1)*t0;
%{
%for loop
for i = 1:length(omega(:,3))
    disp(omega(i,3))
    if class(omega(i,3)) == "imag"
        error('class not double')
    end
end
%}
maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.
tau = m*t0;
%class(omega(:,3))
omega(402372,3) = double(0.0);
omega(402371,3) = double(0.0);
%{
for i = 1: length(omega(:,3))
    if isnan(omega(i,3))
        disp(i);
        error('error in value')
    end
    disp(i)
end
%}
[avarFromFunc, tauFromFunc] = allanvar(omega(:,3), m, Fs);
adevFromFunc = sqrt(avarFromFunc);
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N_GyroZ = 10^logN
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tauFromFunc);
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K_GyroZ = 10^logK
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tauFromFunc); %tau = tauFromFunc
logadev = log10(adevFromFunc);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B_GyroZ = 10^logB
% Plot the results.
tauN = 1;
lineN = N_GyroZ ./ sqrt(tau);
% Plot the results.
tauK = 3;
lineK = K_GyroZ .* sqrt(tau/3);
% Plot the results.
tauB = tau(i);
lineB = B_GyroZ * scfB * ones(size(tau));
% Plot results
figure('Name',"Allan Deviation of Gyro Z", 'NumberTitle', 'off')
loglog(tauFromFunc, adevFromFunc, tauFromFunc, lineN, '--', tauN, N_GyroZ, 'o', tauFromFunc, lineK, '--', tauK, K_GyroZ, 'o', tau, lineB, '--', tauB, scfB*B_GyroZ, 'o');
legend('\sigma (rad/s)', '\sigma_N ((rad/s)sqrt{Hz})', '','\sigma_K ((rad/s)sqrt{Hz})','','\sigma_B (rad/s)')
title('Allan Deviation of Gyro Z')
xlabel('\tau')
ylabel('\sigma(\tau)')
text(tauN, N_GyroZ, 'N')
text(tauK, K_GyroZ, 'K')
text(tauB, scfB*B_GyroZ, '0.664B')
grid on
axis equal