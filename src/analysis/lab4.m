CircleData = rosbag('data_going_in_circles.bag');
MovingData = rosbag('data_driving.bag');

%%%%%%% GPS Data %%%%%%%%%
%% Circle %%
CircleDataGPS_TopicData = select(CircleData,'Topic','/gps');
msgStructs_CircleDataGPS = readMessages(CircleDataGPS_TopicData,'DataFormat','struct');

EastingOffset_Circle = cellfun(@(m) double(m.UTMEasting), msgStructs_CircleDataGPS(1));
NorthingOffset_Circle = cellfun(@(m) double(m.UTMNorthing), msgStructs_CircleDataGPS(1));
Circle_X = cellfun(@(m) double(m.UTMEasting),msgStructs_CircleDataGPS)-EastingOffset_Circle;
Circle_Y = cellfun(@(m) double(m.UTMNorthing),msgStructs_CircleDataGPS)-NorthingOffset_Circle;

GPS_Time1_c = cellfun(@(m) double(m.Header.Stamp.Sec),msgStructs_CircleDataGPS);
GPS_Time_c = GPS_Time1_c-GPS_Time1_c(1,1);

%%% Moving data %%%
MovingDataGPS_TopicData = select(MovingData,'Topic','/gps');
msgStructs_MovingDataGPS = readMessages(MovingDataGPS_TopicData,'DataFormat','struct');

EastingOffset_Moving =cellfun(@(m) double(m.UTMEasting), msgStructs_MovingDataGPS(1));
NorthingOffset_Moving =cellfun(@(m) double(m.UTMNorthing), msgStructs_MovingDataGPS(1));
Moving_Easting = cellfun(@(m) double(m.UTMEasting),msgStructs_MovingDataGPS)-EastingOffset_Moving;
Moving_Northing = cellfun(@(m) double(m.UTMNorthing),msgStructs_MovingDataGPS)-NorthingOffset_Moving;

GPS_Time1 = cellfun(@(m) double(m.Header.Stamp.Sec),msgStructs_MovingDataGPS);
GPS_Time = GPS_Time1-GPS_Time1(1,1);

%%%%%%% Imu Data %%%%%%%%%
CircleDataImu_TopicData = select(CircleData,'Topic','/imu');
msgStructs_CircleDataImu = readMessages(CircleDataImu_TopicData,'DataFormat','struct');

Circle_Orientation_X = cellfun(@(m) double(m.Imu.Orientation.X),msgStructs_CircleDataImu);
Circle_Orientation_Y = cellfun(@(m) double(m.Imu.Orientation.Y),msgStructs_CircleDataImu);
Circle_Orientation_Z = cellfun(@(m) double(m.Imu.Orientation.Z),msgStructs_CircleDataImu);
Circle_Orientation_W = cellfun(@(m) double(m.Imu.Orientation.W),msgStructs_CircleDataImu);

Circle_AngularVelocity_X = cellfun(@(m) double(m.Imu.AngularVelocity.X),msgStructs_CircleDataImu);
Circle_AngularVelocity_Y = cellfun(@(m) double(m.Imu.AngularVelocity.Y),msgStructs_CircleDataImu);
Circle_AngularVelocity_Z = cellfun(@(m) double(m.Imu.AngularVelocity.Z),msgStructs_CircleDataImu);

Circle_LinearAcceleration_X = cellfun(@(m) double(m.Imu.LinearAcceleration.X),msgStructs_CircleDataImu);
Circle_LinearAcceleration_Y = cellfun(@(m) double(m.Imu.LinearAcceleration.Y),msgStructs_CircleDataImu);
Circle_LinearAcceleration_Z = cellfun(@(m) double(m.Imu.LinearAcceleration.Z),msgStructs_CircleDataImu);

% Circle_MagneticField_X_initial =cellfun(@(m) double(m.MagField.MagneticField_.X),msgStructs_CircleDataImu);
% Circle_MagneticField_Y_initial =cellfun(@(m) double(m.MagField.MagneticField_.Y),msgStructs_CircleDataImu);
% Circle_MagneticField_Z_initial =cellfun(@(m) double(m.MagField.MagneticField_.Z),msgStructs_CircleDataImu);
% Circle_MagneticField_offset_X = (min(Circle_MagneticField_X_initial)+ max(Circle_MagneticField_X_initial))/2;
% Circle_MagneticField_offset_Y = (min(Circle_MagneticField_Y_initial)+ max(Circle_MagneticField_Y_initial))/2;
% Circle_MagneticField_offset_Z = (min(Circle_MagneticField_Z_initial)+ max(Circle_MagneticField_Z_initial))/2;

Circle_MagneticField_X = cellfun(@(m) double(m.MagField.MagneticField_.X),msgStructs_CircleDataImu);
Circle_MagneticField_Y = cellfun(@(m) double(m.MagField.MagneticField_.Y),msgStructs_CircleDataImu);
Circle_MagneticField_Z = cellfun(@(m) double(m.MagField.MagneticField_.Z),msgStructs_CircleDataImu);



%%% Moving Data
MovingDataImu_TopicData = select(MovingData,'Topic','/imu');
msgStructs_MovingDataImu = readMessages(MovingDataImu_TopicData,'DataFormat','struct');

Moving_Orientation_X = cellfun(@(m) double(m.Imu.Orientation.X),msgStructs_MovingDataImu);
Moving_Orientation_Y = cellfun(@(m) double(m.Imu.Orientation.Y),msgStructs_MovingDataImu);
Moving_Orientation_Z = cellfun(@(m) double(m.Imu.Orientation.Z),msgStructs_MovingDataImu);
Moving_Orientation_W = cellfun(@(m) double(m.Imu.Orientation.W),msgStructs_MovingDataImu);

Moving_AngularVelocity_X = cellfun(@(m) double(m.Imu.AngularVelocity.X),msgStructs_MovingDataImu);
Moving_AngularVelocity_Y = cellfun(@(m) double(m.Imu.AngularVelocity.Y),msgStructs_MovingDataImu);
Moving_AngularVelocity_Z = cellfun(@(m) double(m.Imu.AngularVelocity.Z),msgStructs_MovingDataImu);

Moving_LinearAcceleration_X = cellfun(@(m) double(m.Imu.LinearAcceleration.X),msgStructs_MovingDataImu);
Moving_LinearAcceleration_Y = cellfun(@(m) double(m.Imu.LinearAcceleration.Y),msgStructs_MovingDataImu);
Moving_LinearAcceleration_Z = cellfun(@(m) double(m.Imu.LinearAcceleration.Z),msgStructs_MovingDataImu);

Moving_MagneticField_X = cellfun(@(m) double(m.MagField.MagneticField_.X),msgStructs_MovingDataImu);
Moving_MagneticField_Y = cellfun(@(m) double(m.MagField.MagneticField_.Y),msgStructs_MovingDataImu);
Moving_MagneticField_Z = cellfun(@(m) double(m.MagField.MagneticField_.Z),msgStructs_MovingDataImu);

%% 

TimeOffset_Circle = double(msgStructs_CircleDataImu{1}.Header.Stamp.Sec);
Time_seconds_Circle = cellfun(@(m) double(m.Header.Stamp.Sec),msgStructs_CircleDataImu)-TimeOffset_Circle;
Time_nseconds_Circle = cellfun(@(m) double(m.Header.Stamp.Nsec),msgStructs_CircleDataImu);

radEul_Circle = quat2eul([Circle_Orientation_W Circle_Orientation_X Circle_Orientation_Y Circle_Orientation_Z]);
eul_Circle = rad2deg(radEul_Circle);

Time_Circle = Time_seconds_Circle;
for i= 1:size(Time_seconds_Circle,1)
    B = fix(abs(log10(abs(Time_nseconds_Circle(1,1)))))+1;
    Time_Circle(i,1) = Time_seconds_Circle(i,1) + (Time_nseconds_Circle(i,1)/10.^B);
end
%%

%%% Moving Data
TimeOffset_Moving = double(msgStructs_MovingDataImu{1}.Header.Stamp.Sec);
Time_seconds_Moving = cellfun(@(m) double(m.Header.Stamp.Sec),msgStructs_MovingDataImu)-TimeOffset_Moving;
Time_nseconds_Moving = cellfun(@(m) double(m.Header.Stamp.Nsec),msgStructs_MovingDataImu);

radEul_Moving = quat2eul([Moving_Orientation_W Moving_Orientation_X Moving_Orientation_Y Moving_Orientation_Z]);
eul_Moving = rad2deg(radEul_Moving);

Time_Moving = Time_seconds_Moving;
for i= 1:size(Time_seconds_Moving,1)
    B = fix(abs(log10(abs(Time_nseconds_Moving(1,1)))))+1;
    Time_Moving(i,1) = Time_seconds_Moving(i,1) + (Time_nseconds_Moving(i,1)/10.^B);
end

%% %%%%  PLOTS   %%%%%%%%%
figure
plot(Circle_X,Circle_Y);
title("UTM Easting Northing Data - Circle")
xlabel({'UTM Easting','(in m)'})
ylabel({'UTM Northing','(in m)'})
%
figure
plot(Moving_Easting,Moving_Northing);
title("UTM Easting Northing Data - Moving")
xlabel({'UTM Easting','(in m)'})
ylabel({'UTM Northing','(in m)'})
%%
figure
plot(Circle_MagneticField_X,Circle_MagneticField_Y,'o')
hold on
ellipse = fit_ellipse(Circle_MagneticField_X,Circle_MagneticField_Y);
plot(ellipse.X0_in,ellipse.Y0_in,'rx');
hold on
theta_ellipse = linspace(0, 2*pi);
x_ellipse = ellipse.X0_in + ellipse.a*cos(theta_ellipse)*cos(ellipse.phi) - ellipse.b*sin(theta_ellipse)*sin(ellipse.phi);
y_ellipse = ellipse.Y0_in + ellipse.a*cos(theta_ellipse)*sin(ellipse.phi) + ellipse.b*sin(theta_ellipse)*cos(ellipse.phi);
plot(x_ellipse, y_ellipse, 'm-', 'LineWidth', 2);
grid on
title("Magnetometer Data Before Soft and Hard Iron Correction")
xlabel({'Magnetometer X','(in Gauss)'})
ylabel({'Magnetometer Y','(in Gauss)'})
legend('Data','Center','Fit');
hold off


%% Calibration
Rotation_matrix = [cos(-ellipse.phi) -sin(-ellipse.phi);sin(-ellipse.phi) cos(-ellipse.phi)];
Input_matrix = [Circle_MagneticField_X Circle_MagneticField_Y];
SoftIronCorrection = Input_matrix*Rotation_matrix;
ellipse2 = fit_ellipse(SoftIronCorrection(:,1),SoftIronCorrection(:,2));
sigma = ellipse2.short_axis/ellipse2.long_axis;
%%
ellipse2 = fit_ellipse(SoftIronCorrection(:,1)*sigma,SoftIronCorrection(:,2));
%
figure
plot(Circle_MagneticField_X,Circle_MagneticField_Y); 
hold on
plot(ellipse.X0_in,ellipse.Y0_in,'bx');
hold on
theta_ellipse = linspace(0, 2*pi);
x_ellipse = ellipse.X0_in + ellipse.a*cos(theta_ellipse)*cos(ellipse.phi) - ellipse.b*sin(theta_ellipse)*sin(ellipse.phi);
y_ellipse = ellipse.Y0_in + ellipse.a*cos(theta_ellipse)*sin(ellipse.phi) + ellipse.b*sin(theta_ellipse)*cos(ellipse.phi);
plot(x_ellipse, y_ellipse, 'b-', 'LineWidth', 2);
hold on
plot((SoftIronCorrection(:,1)*sigma)-ellipse2.X0_in,SoftIronCorrection(:,2)-ellipse2.Y0_in);
hold on
ellipse2 = fit_ellipse((SoftIronCorrection(:,1)*sigma)-ellipse2.X0_in,(SoftIronCorrection(:,2))-ellipse2.Y0_in);
plot(ellipse2.X0_in,ellipse2.Y0_in,'rx');
hold on
theta_ellipse2 = linspace(0, 2*pi);
x_ellipse2 = ellipse2.X0_in + ellipse2.a*cos(theta_ellipse2)*cos(ellipse2.phi) - ellipse2.b*sin(theta_ellipse2)*sin(ellipse2.phi);
y_ellipse2 = ellipse2.Y0_in + ellipse2.a*cos(theta_ellipse2)*sin(ellipse2.phi) + ellipse2.b*sin(theta_ellipse2)*cos(ellipse2.phi);
plot(x_ellipse2, y_ellipse2, 'r-', 'LineWidth', 2);
hold off
axis equal
grid on
title("Magnetometer Data Before and After Soft and Hard Iron Calibration")
xlabel({'Magnetometer X','(in Gauss)'})
ylabel({'Magnetometer Y','(in Gauss)'})
legend('Data Before','Center Before','Fit Before', 'Data After','Center After','Fit After');


%%
figure
plot(Time_Circle,Circle_MagneticField_X)
grid on
title("Magnetometer Data vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Magnetometer X','(in Gauss)'})
legend('Data');
%
figure
plot(Time_Circle,Circle_MagneticField_Y)
grid on
title("Magnetometer Data vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Magnetometer Y','(in Gauss)'})
legend('Data');

xmag = Circle_MagneticField_X;
ymag = Circle_MagneticField_Y;
zmag = Circle_MagneticField_Z;

Circle_pitch = atan2(Circle_LinearAcceleration_X,sqrt(Circle_LinearAcceleration_Y.^2+Circle_LinearAcceleration_Z.^2));
Circle_roll = atan2(Circle_LinearAcceleration_Y,sqrt(Circle_LinearAcceleration_X.^2+Circle_LinearAcceleration_Z.^2));

yaw = atan2((-ymag.*cos(Circle_roll)+zmag.*sin(Circle_roll)),(xmag.*cos(Circle_pitch)+ymag.*sin(Circle_pitch).*sin(Circle_roll)+zmag.*sin(Circle_pitch).*cos(Circle_roll)));

%% MOVING DATA CORRECTION

 Moving_MagneticField_X_HardIron = Moving_MagneticField_X-ellipse2.X0_in;
 Moving_MagneticField_Y_HardIron = Moving_MagneticField_Y-ellipse2.Y0_in;
 
 Moving_Input_matrix = [Moving_MagneticField_X_HardIron Moving_MagneticField_Y_HardIron];

 Moving_SoftIronCorrection1 = Moving_Input_matrix*Rotation_matrix;

 Moving_SoftIronCorrection(:,1) = Moving_SoftIronCorrection1(:,1)*sigma;
 Moving_SoftIronCorrection(:,2) = Moving_SoftIronCorrection1(:,2)*sigma;

 Moving_xmag = Moving_MagneticField_X;
 Moving_ymag = Moving_MagneticField_Y;
 Moving_zmag = Moving_MagneticField_Z;

Corrected_xmag = Moving_SoftIronCorrection(:,1);
Corrected_ymag = Moving_SoftIronCorrection(:,2);
Corrected_zmag = Moving_MagneticField_Z;
movingYaw2 = (atan2(-Corrected_ymag,Corrected_xmag)-2);


%%
Moving_MagneticField (:,1) = Moving_MagneticField_X;
Moving_MagneticField (:,2) = Moving_MagneticField_Y;

%%
figure
plot(Time_Moving,Moving_MagneticField_X)
hold on
plot(Time_Moving,Moving_SoftIronCorrection(:,1))
hold off
grid on
title("Magentometer X Before and After correction vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Yaw','(in Gauss)'})
legend('Magnetometer Reading Before Correction','Magnetometer Reading After Correction','Location','best');
%%
figure
plot(Time_Moving,Moving_MagneticField_Y)
hold on
plot(Time_Moving,Moving_SoftIronCorrection(:,2))
hold off
grid on
title("Magentometer Y Before and After correction vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Yaw','(in Gauss)'})
legend('Magnetometer Reading Before Correction','Magnetometer Reading After Correction','Location','best');
%%
figure
Moving_yaw = drawYaw(Moving_xmag,Moving_ymag,Moving_zmag,Moving_LinearAcceleration_X,Moving_LinearAcceleration_Y,Moving_LinearAcceleration_Z,Time_Moving);
hold on
Corrected_Moving_yaw = drawYaw(Corrected_xmag,Corrected_ymag,Corrected_zmag,Moving_LinearAcceleration_X,Moving_LinearAcceleration_Y,Moving_LinearAcceleration_Z,Time_Moving);
hold off
grid on
title("Magentometer Yaw and Corrected Magnetometer Yaw vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Yaw','(in Gauss)'})
legend('Magnetometer Yaw','Corrected Magnetometer Yaw','Location','best');

%%

yaw_angle_vel = cumtrapz(Time_Moving,Moving_AngularVelocity_Z);
figure
plot(Time_Moving,unwrap(yaw_angle_vel),'r')
hold on
plot(Time_Moving,unwrap(movingYaw2),'k')
hold off
grid on
title('Yaw Angle from two methods (Integrated from Gyro & Corrected Magnetometer)')
legend('Yaw integrated from gyro','yaw from magnetometer')


%% Complementory Filter%%
alpha1 = 0.79;
alpha_low = 3;
alpha_high = 0.0003;



mag_lowpass = lowpass(unwrap(movingYaw2),alpha_low,40);
gyro_highpass = highpass(unwrap(yaw_angle_vel),alpha_high,40);
complimentary_filter = unwrap(alpha1*mag_lowpass+(1-alpha1)*gyro_highpass);

%% PLOTS %%
figure
plot(Time_Moving,mag_lowpass,'r','LineWidth',0.5)
hold on
plot(Time_Moving,gyro_highpass,'g')
hold on
plot(Time_Moving,complimentary_filter,'blue','LineWidth',0.5);
hold off
grid on
title("LPF, HPF & CF plots vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Yaw','(in Rad)'})
legend('Lowpass filter','Highpass filter','Complementary filter');
%%
figure
plot(Time_Moving,complimentary_filter,'blue')
hold on 
plot(Time_Moving,unwrap(Moving_yaw))
hold off
grid on
title("Complementary filter Yaw and Yaw angle from Imu sensor vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Yaw','(in Rad)'})
legend('Complementary filter','Yaw from Imu','location','best');
%%
figure
plot(Time_Moving,Corrected_Moving_yaw)
hold on
plot(Time_Moving,unwrap(Corrected_Moving_yaw))
grid on
title("Plot of Magnetic yaw before and after using unwrap function vs Time")
xlabel({'Time','(in Seconds)'})
ylabel({'Yaw','(in Rad)'})
legend('Data','Unwrapped Data','location','best');
%%
figure
plot(Time_Moving,unwrap(Moving_yaw),'-y')
hold on
plot(Time_Moving,unwrap(yaw_angle_vel),'-r')
hold on
plot(Time_Moving,complimentary_filter,'-b')
hold on
plot(Time_Moving,unwrap(movingYaw2))
hold off
grid on
title('Yaw Comparison between four methods')
legend('Magnetometer','Yaw integrated from gyro','Complementary Filter','Yaw angle computed from IMU')
%%
figure
plot(Time_Moving, Moving_LinearAcceleration_X,'b')
hold on
plot(Time_Moving, Moving_LinearAcceleration_Y,'g')
hold on
plot(Time_Moving, Moving_LinearAcceleration_Z,'r')
hold off
grid on 
title('Comarison on Linear Acceleration on the IMU')
legend('Linear Acc in X','Linear Acc in Y','Linear Acc in Z')

%% 







%% Estimation of Forward velocity !!
% integrate forward acc to get forward velocity
Imu_forward_velocity = cumtrapz(Time_Moving,Moving_LinearAcceleration_X);

GPS_data(:,1) = Moving_Easting;
GPS_data(:,2) = Moving_Northing;
GPS_velocity = zeros(size(Moving_Easting));

for i = 1:size(GPS_data,1)-1
     distance = GPS_data(i+1,:)-GPS_data(i,:);
     GPS_velocity(i) = ((distance(1)^2 +distance(2)^2)^0.5)/((Time_Moving(i+1)-Time_Moving(i)));
end

%%
figure
plot(GPS_Time,GPS_velocity/10)
hold on
plot(Time_Moving,Imu_forward_velocity)
hold off
title("Plot of GPS velocity and Imu velocity vs Time before adjustment")
xlabel({'Time','(in Seconds)'})
ylabel({'Velocity','(in m/s)'})
legend('GPS Velocity','Imu Velocity','location','best');
%%
Jerk = zeros(size(Moving_LinearAcceleration_X));

for i = 1:size(Moving_LinearAcceleration_X,1)-1
     distance = Moving_LinearAcceleration_X(i+1,:)-Moving_LinearAcceleration_X(i,:);
     Jerk(i) = distance/((Time_Moving(i+1)-Time_Moving(i)));
end
figure
plot(Time_Moving,Jerk)
max(Jerk)
min(Jerk)

%%

final_velocity = zeros(size(Jerk));
window = 100;
for n = 1:size(Jerk)
    %disp(Jerk(n,1))
    if Jerk(n,1) < 7.8 && Jerk(n,1) > 4.4
        final_velocity(n,1) = 0;
        
    else
        final_velocity(n,1) = Imu_forward_velocity(n,1);
        disp(final_velocity(n,1))
    end 
        
end

%%
count = zeros(size(Jerk));
for n = 1:size(Jerk)-window
    
    for m = 1:window

        if final_velocity(n+m,1) == 0 
            count(n,1) = 81.45;
        end
    end
end
velocity2 = zeros(size(Jerk));
for n = 1:size(Jerk)-window
    if count(n,1) == 0
        for m = n:size(Jerk)-window
            velocity2(m,1) = Imu_forward_velocity(m,1)-Imu_forward_velocity(n,1);
            
        end
    else
        velocity2(n,1) = Imu_forward_velocity(n,1);
    end
end

figure
plot(Time_Moving,count)
hold on
plot(Time_Moving,Jerk)
plot(Time_Moving,velocity2)
hold off


%%
velocity3 = zeros(size(Jerk));
for n = 2:size(Jerk)
    if velocity2(n-1,1) == 0
        for m = n:size(Jerk)
            velocity3(m,1) = velocity2(m,1) - velocity2(n,1);
        end
    end
end
for n = 2:size(Jerk)
    if velocity3(n,1) < 0
        velocity3(n,1) = 0;
    end
end
figure
plot(Time_Moving,velocity3)
hold on
plot(GPS_Time,GPS_velocity./6.11)
hold off
legend('Imu velocity','gps velocity')

%%
p = fittype('Poly1');

[fitresult, gof] = fit(Time_Moving,velocity3, p);

figure
plot(fitresult,Time_Moving,velocity3,'predoba')

%%
i = 0;
m = 1;
limit = size(Jerk);
while m < limit(1,1)
        if count(m,1) == 81.45
            i = i+1;
            sections(i,1) = m;
            while count(m,1) == 81.45
                m = m +1;
            end
            sections(i,2) = m-1;
        else
            m = m + 1;
        end

end

%%
limit2 = size(sections);
p = fittype('Poly1');
figure
velocity7 = size(Jerk);
for n = 1:limit2(1,1)
    start = sections(n,1);
    ending = sections(n,2);
    [fitresult, gof] = fit(Time_Moving(start:ending,1),velocity3(start:ending,1), p);

    %plot(fitresult,Time_Moving(start:ending,1),velocity3(start:ending,1),'predoba')
    hold on
    %plot(Time_Moving(start:ending,1),velocity3(start:ending,1))
    for x = start:ending
        %distance1(x,1) = GetPointDistance(fitresult.p1,fitresult.p2,start,ending,Time_Moving);
        distance1(x,1) = fitresult.p1.*(Time_Moving(x,1)) + fitresult.p2;
        %velocity7(x,1) = velocity3(x,1) - distance1;
    end
    %plot(Time_Moving,distance)

end
%%
distance2 = distance1 - distance1*10/100;
plot(Time_Moving(1:size(distance1)),distance2)
hold off
distance2(size(distance1):size(Jerk),1) = 0;
velocity7 = zeros(size(velocity3));
velocity7 = velocity3(1:size(distance2),1) - distance2;

for n = 1:size(Jerk)
    if velocity7(n,1) < 0
        velocity7(n,1) = 0;
    end
end

figure 
plot(Time_Moving(1:size(Jerk)-101),velocity7(1:size(Jerk)-101))
hold on
plot(GPS_Time,GPS_velocity./60.11)
hold off
title("Plot of GPS velocity and Imu velocity vs Time after adjustment")
xlabel({'Time','(in Seconds)'})
ylabel({'Velocity','(in m/s)'})
legend('GPS Velocity','Imu Velocity','location','best');


%% Dead Reckoning !!

Imu_displacement = cumtrapz(Time_Moving,velocity7);
figure 
plot(Time_Moving,Imu_displacement)  
hold on
plot(GPS_Time,Moving_Northing)



%%

t = linspace(1, 1500, 58791);

Imu_acc_mv2 = Moving_LinearAcceleration_X; 

vel_car_X = velocity7;

w = Moving_AngularVelocity_Z;
value = w .* vel_car_X./3;

acc_car_Y = Moving_AngularVelocity_Y;

figure;
plot(t, value, 'color', 'b');
hold on 
xlabel('time series (second)'); 
title('Compare w * X dot AND acc car y');
grid on
plot(t, acc_car_Y);
legend('w * Xdot', 'acc car y');

%% 
yaw_frommag_shift = Corrected_Moving_yaw;

fwd_vel_from_GPS = GPS_velocity;
new_GPS_vel = GPS_velocity;%Get_New_GPS_Vel(fwd_vel_from_GPS);

vn = [];
ve = [];

yaw_from_gyro = Moving_AngularVelocity_Z;
yaw_from_gyro_ture = Moving_AngularVelocity_Z;
yaw_from_filter = Moving_AngularVelocity_Z;
yaw_from_Imu_shift = Moving_AngularVelocity_Z;

for ii = 1:length(yaw_frommag_shift)
    angle = yaw_frommag_shift(ii);
    
    vn = [vn; vel_car_X(ii) * cos(angle)]; 
    ve = [ve; vel_car_X(ii) * sin(angle)];
end
v_head = [ve, vn];

xe = cumtrapz(t, ve);
xn = cumtrapz(t, vn);
figure(2);
plot(xe, xn);
xlabel('East'); 
ylabel('North');
title('The estimated trajectory before adjustment');
grid on

%% plotting the GPS track and estimated trajectory on the same plot
gps_data_utm(:,1) = Moving_Easting; 
gps_data_utm(:,2) = Moving_Northing; 


%%

xe = xe + gps_data_utm(1, 1);
xn = xn + gps_data_utm(1, 2);

GPS_point_A = [gps_data_utm(27, 1), gps_data_utm(27, 2)];
GPS_point_B = [gps_data_utm(840, 1), gps_data_utm(840, 2)];
slope_A = (GPS_point_B(2) - GPS_point_A(2)) / (GPS_point_B(1) - GPS_point_A(1));
degree_A = rad2deg(atan(slope_A) + pi);

estimated_pointA = [152.469, -241];
estimated_pointB = [xe(5748), xn(5748)];
slope_B = (estimated_pointB(2) - estimated_pointA(2)) / (estimated_pointB(1) - estimated_pointA(1));
degree_B = rad2deg(atan(slope_B) + pi);
theta = -deg2rad(degree_A - degree_B);
theta = 0.05;
R_matrix = [cos(theta), sin(theta);
            -sin(theta), cos(theta)];
new_x = R_matrix * [xe, xn]';
new_x = new_x';
        
figure
plot(gps_data_utm(:, 1), gps_data_utm(:, 2), 'linewidth', 0.5);
xlabel('UTM easting'); 
ylabel('UTM northing');
title('GPS UTM trajectory and estimated trajectory w/o adjust');
grid off
hold on

plot(xe, -xn);
legend('GPS trajectory', 'estimated trajectory');
hold off

%%
% after adjustment
figure(4);
plot(gps_data_utm(:, 1), gps_data_utm(:, 2), 'linewidth', 0.5);
xlabel('UTM easting'); 
ylabel('UTM northing');
title('GPS UTM trajectory and adjusted estimated trajectory');
grid on
hold on
plot(new_x(:, 1), new_x(:, 2));
legend('GPS trajectory', 'adjusted trajectory');
hold off

%% 3. Estimate x_c
gps_acc = [(new_GPS_vel(1)) / t(1)];
for ii = 2:length(new_GPS_vel)-168
    a_temp = (new_GPS_vel(ii)-new_GPS_vel(ii-1)) / (t(ii)-t(ii-1));
    gps_acc = [gps_acc; a_temp];
end

w_dot = [(w(1)) / t(1)];
for ii = 2:length(w)
    w_dot_temp = (w(ii)-w(ii-1)) / (t(ii)-t(ii-1));
    w_dot = [w_dot; w_dot_temp];
end

Imu_acc = Moving_LinearAcceleration_X;
zzz = w_dot + w.*w;
Imu_acc_drift = mean(Imu_acc(1:650, 1));
x_c = (Imu_acc(:, 1)-Imu_acc_drift - Imu_acc_mv2 - value) ./ zzz
mean(x_c)
%%
% Compute forward velocity by integrating forward acceleration
forward_acceleration = Moving_LinearAcceleration_X; % Assuming forward acceleration is in X-axis
forward_velocity = cumtrapz(Time_Moving, forward_acceleration);

% Compute heading from magnetometer data
heading_rad = atan2(Moving_MagneticField_Y_HardIron, Moving_MagneticField_X_HardIron);
heading_deg = rad2deg(heading_rad);

% Rotate forward velocity to obtain Easting and Northing components
ve = velocity7 .* sin(heading_rad);
vn = velocity7 .* cos(heading_rad);

% Integrate ve and vn to estimate the trajectory of the vehicle (xe, xn)
xe = cumtrapz(Time_Moving, -ve);
xn = cumtrapz(Time_Moving, -vn);

% Plot IMU-based trajectory and GPS-based trajectory
figure;
plot(xe, xn, 'r', 'LineWidth', 1.5);
hold on;
plot(gps_data_utm(:,1), gps_data_utm(:,2), 'b', 'LineWidth', 1.5);
xlabel('Easting (m)');
ylabel('Northing (m)');
legend('Dead Reckoning (IMU)', 'GPS Trajectory');
title('Dead Reckoning vs GPS Trajectory');
grid on;
hold off;

%%
% Define the desired starting point
xe_start = 150; % Replace with the desired starting value for Easting
xn_start = -239; % Replace with the desired starting value for Northing

% Compute the offsets to adjust the starting point
xe_offset = xe_start - xe(1);
xn_offset = xn_start - xn(1);

% Adjust the trajectory with the computed offsets
xe_adjusted = xe + xe_offset;
xn_adjusted = xn + xn_offset;

% Plot the adjusted IMU-based trajectory and GPS-based trajectory
figure;
plot(xe_adjusted, xn_adjusted, 'r', 'LineWidth', 1.5);
hold on;
plot(Moving_Easting(1:1450), Moving_Northing(1:1450), 'b', 'LineWidth', 1.5);
xlabel('Easting (m)');
ylabel('Northing (m)');
legend('Dead Reckoning (IMU)', 'GPS Trajectory');
title('Dead Reckoning vs GPS Trajectory');
grid on;
hold off;

%%





















%% Function for Yaw calculation
function Moving_yaw = drawYaw(Moving_xmag1,Moving_ymag1,Moving_zmag1,Moving_LinearAcceleration_X,Moving_LinearAcceleration_Y,Moving_LinearAcceleration_Z,Time_Moving)

 mag_norm = sqrt((Moving_xmag1.*Moving_xmag1)+(Moving_ymag1.*Moving_ymag1)+(Moving_zmag1.*Moving_zmag1));
 
 Moving_xmag = Moving_xmag1./mag_norm;
 Moving_ymag = Moving_ymag1./mag_norm;
 Moving_zmag = Moving_zmag1./mag_norm;

 Moving_pitch = atan2(Moving_LinearAcceleration_X,sqrt(Moving_LinearAcceleration_Y.^2+Moving_LinearAcceleration_Z.^2));
 Moving_roll = atan2(Moving_LinearAcceleration_Y,sqrt(Moving_LinearAcceleration_X.^2+Moving_LinearAcceleration_Z.^2));

 Moving_yaw = (atan2((-Moving_ymag.*cos(Moving_roll)+Moving_zmag.*sin(Moving_roll)),(Moving_xmag.*cos(Moving_pitch)+Moving_ymag.*sin(Moving_pitch).*sin(Moving_roll)+Moving_zmag.*sin(Moving_pitch).*cos(Moving_roll))))-2;
 
 plot(Time_Moving,unwrap(Moving_yaw))

end
%%
