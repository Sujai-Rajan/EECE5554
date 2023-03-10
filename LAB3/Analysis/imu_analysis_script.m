% IMU stationary data analysis
	
	
imu_stationary_data = rosbag('imu_stn_data.bag');
	
imu_data_Topic = select(imu_stationary_data, "Topic", '/imu');
	
imu_data_read = readMessages(imu_data_Topic, 'DataFormat','struct');
	
Angular_velocity = zeros(length(imu_data_read),3);
	
Linear_Acceleration = zeros(length(imu_data_read),3);
	
Quaternion = zeros(length(imu_data_read),4);
	
Magnetic_Field = zeros(length(imu_data_read),3);
	
for s = 1:length(imu_data_read)
	
    %Angular_velocity(s,1) = imu_data_read{s}.geometry_msgs/Vector3;
	
    Angular_velocity(s,1) = imu_data_read{s}.Imu.AngularVelocity.X;
	
    Angular_velocity(s,2) = imu_data_read{s}.Imu.AngularVelocity.Y;
	
    Angular_velocity(s,3) = imu_data_read{s}.Imu.AngularVelocity.Z;
	
    Linear_Acceleration(s,1) = imu_data_read{s}.Imu.LinearAcceleration.X;
	
    Linear_Acceleration(s,2) = imu_data_read{s}.Imu.LinearAcceleration.Y;
	
    Linear_Acceleration(s,3) = imu_data_read{s}.Imu.LinearAcceleration.Z;
	
    Quaternion(s,1) = imu_data_read{s}.Imu.Orientation.W;
	
    Quaternion(s,2) = imu_data_read{s}.Imu.Orientation.X;
	
    Quaternion(s,3) = imu_data_read{s}.Imu.Orientation.Y;
	
    Quaternion(s,4) = imu_data_read{s}.Imu.Orientation.Z;
	
    Magnetic_Field(s,1) = imu_data_read{s}.MagField.MagneticField_.X;
	
    Magnetic_Field(s,2) = imu_data_read{s}.MagField.MagneticField_.Y;
	
    Magnetic_Field(s,3) = imu_data_read{s}.MagField.MagneticField_.Z;
	
end 
	
	
%Gyro Plot
	
figure('Name',"Stationary Data, Gyro Plot", 'NumberTitle', 'off')

subplot(2,2,1)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Angular_velocity(:,1),'b+')

legend ('Gyro_X')
	
hold on

subplot(2,2,2)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Angular_velocity(:,2),'g+')

legend ('Gyro_Y')
	
hold on

subplot(2,2,3)

plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Angular_velocity(:,3),'r+')

legend ('Gyro_Z')
	
hold on

subplot(2,2,4)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Angular_velocity(:,1),'b+')
	
hold on
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Angular_velocity(:,2),'g+')
	
hold on
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Angular_velocity(:,3),'r+')
	
legend ('Gyro_X, Gyro_Y, Gyro_Z')

grid on
	
title('Stationary Data, Gyro Plot')
	
xlabel('Time (seconds)')
	
ylabel('Gyro (rad/s)')
	
hold off
	

%% 
	
	
%Acceleration Plot
	
figure('Name',"Stationary Data, Acceleration Plot", 'NumberTitle', 'off')
	
subplot(2,2,1)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Linear_Acceleration(:,1),'b+')
	
legend ('Accel_X')
	
hold on
	
subplot(2,2,2)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Linear_Acceleration(:,2),'g+')
	
legend ('Accel_Y')
	
hold on
	
subplot(2,2,3)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Linear_Acceleration(:,3),'r+')
	
legend ('Accel_Z')
	
hold on
	
subplot(2,2,4)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Linear_Acceleration(:,1),'b+')
	
hold on
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Linear_Acceleration(:,2),'g+')
	
hold on
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Linear_Acceleration(:,3),'r+')
	
legend ('Accel_X, Accel_Y, Accel_Z')
	
grid on
	
title('Stationary Data, Acceleration Plot')
	
%legend ('Accel_X','Accel_Y','Accel_Z')
	
xlabel('Time (seconds)')
	
ylabel('Acceleration (m/s^2)')
	
hold off
	
%% 
	
	
Euler_Angles = quat2eul(Quaternion);
	
%Euler Angle Plot
	
figure('Name',"Stationary Data, Euler Angle Plot", 'NumberTitle', 'off')
	
subplot(2,2,1)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Euler_Angles(:,1),'b+')
	
legend ('Ori_X')
	
hold on
	
subplot(2,2,2)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Euler_Angles(:,2),'g+')
	
legend ('Ori_Y')
	
hold on
	
subplot(2,2,3)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Euler_Angles(:,3),'r+')
	
legend ('Ori_Z')
	
hold on
	
subplot(2,2,4)
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Euler_Angles(:,1),'b+')
	
hold on
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Euler_Angles(:,2),'r+')
	
hold on
	
plot(linspace(1,length(imu_data_read)/40,length(imu_data_read)), Euler_Angles(:,3),'g+')
	
grid on
	
title('Stationary Data, Euler Angle Plot')
	
legend ('Orientation_Z','Orientation_Y','Orientation_X')
	
xlabel('Time (seconds)')
	
ylabel('Orientation')
	
hold off
	
mean_accel_x = mean(Linear_Acceleration(:,1));
	
median_accel_x = median(Linear_Acceleration(:,1));
	
standard_deviatiion_x = std(Linear_Acceleration(:,1))
	
fprintf('Mean accel x: %.2f m/s^2\n', mean_accel_x);
	
fprintf('Median accel x: %.2f m/s^2\n', median_accel_x);


%Histogram Plot
	
figure ('Name', "Histogram of Accel X", 'NumberTitle', 'off')
	
histogram(Linear_Acceleration(:,1), 10)

hold on
	
% Overlay the fitted Gaussian distribution on the histogram
	
pd = fitdist(Linear_Acceleration(:,1),'Normal');
	
x_values = linspace(min(Linear_Acceleration(:,1)),max(Linear_Acceleration(:,1)),100);
	
y_values = pdf(pd,x_values);
	
area_under_curve = trapz(x_values, y_values); % compute area under the curve
	
bin_edges = linspace(min(Linear_Acceleration(:,1)), max(Linear_Acceleration(:,1)), 11);
	
bin_width = bin_edges(2) - bin_edges(1);
	
scaling_factor = bin_width * numel(Linear_Acceleration(:,1));
	
y_values_scaled = y_values * scaling_factor / area_under_curve; % scale Gaussian function
	
plot(x_values,y_values_scaled,'LineWidth',2)
	
legend('Histogram Accel_x', 'gaussian distribution')
	
hold on
	
title('Sensor Output Distribution for Accel X')
	
xlabel('Acceleration (m/s^2)')
	
ylabel('Frequency')

hold off
	
