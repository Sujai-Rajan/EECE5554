
stn_data = rosbag("stn.bag");
lne_data = rosbag("lne.bag");

stn_data_topic = select(stn_data,'Topic', '/gps');
lne_data_topic = select(lne_data,'Topic', '/gps');

stn_struct = readMessages(stn_data_topic, 'DataFormat', 'struct');
lne_struct = readMessages(lne_data_topic, 'DataFormat', 'struct');


stn_known_position = [x, y]; % Replace x and y with your known position

stn_struct{1};
lne_struct{1};


%% Stationary Data

scale_value_stn_x = stn_struct{1}.UTMEasting;
scale_value_stn_y = stn_struct{1}.UTMNorthing;

stn_x = cellfun(@(m) double(m.UTMEasting),stn_struct) - scale_value_stn_x;
stn_y = cellfun(@(m) double(m.UTMNorthing),stn_struct) - scale_value_stn_y;

stn_mean_x = mean(cellfun(@(m) double(m.UTMEasting),stn_struct) - scale_value_stn_x);
stn_mean_y = mean(cellfun(@(m) double(m.UTMNorthing),stn_struct) - scale_value_stn_y);

time = cellfun(@(m) double(m.Header.Stamp.Sec),stn_struct);
err_x = (cellfun(@(m) double(m.UTMEasting),stn_struct))-stn_mean_x;
err_y = (cellfun(@(m) double(m.UTMNorthing),stn_struct))-stn_mean_y;

figure
scatter(stn_x,stn_y)
hold on
plot(stn_mean_x,stn_mean_y,'--ro')
hold off
title('Stationary - GPS Data')
xlabel({'UTM_Easting'})
ylabel({'UTM_Northing'})
legend('Data','Mean')


%% Stationary Altitude

stn_alt = (cellfun(@(m) double(m.Altitude),stn_struct));
stn_time = (cellfun(@(m) double(m.Header.Stamp.Sec),stn_struct)-double(stn_struct{1}.Header.Stamp.Sec));
syn_atl_mean = mean((cellfun(@(m) double(m.Altitude),stn_struct)));

figure
scatter(stn_time,stn_alt)
hold on
title('Stationary - Altitude Data')
plot(stn_time,syn_atl_mean,'--ro','LineWidth',0.5)
xlabel({'Time','(in seconds)'})
ylabel({'Altitude','(in meters)'})
legend('Altitude Data','Mean Line')
hold off


for i = 1:length(stn_x)
        % Calculate the error from known position to measured position
        measured_position = [stn_x(i), stn_y(i)];
        error(i) = norm(measured_position - stn_known_position);
end


% Plot the histogram of the error values
figure;
histogram(error);
xlabel('Error (in meters)');
ylabel('Frequency');
title('Error plot - Histogram');



scale_value_lne_x = lne_struct{1}.UTMEasting;
scale_value_lne_y = lne_struct{1}.UTMNorthing;



%% Straight Line - GPS Data

lne_x = cellfun(@(m) double(m.UTMEasting),lne_struct) - scale_value_lne_x;
lne_y = cellfun(@(m) double(m.UTMNorthing),lne_struct) - scale_value_lne_y;

lne_err = point_dist(lne_x, lne_y);
p = fittype('Poly1');

[lne_fit, gof ] = fit( lne_x,lne_y, p);

rmse = sqrt(mean(lne_err)^2);
disp(rmse)
disp(lne_fit)

figure
plot(lne_fit, lne_x, lne_y, 'predoba')
title('Straight Line - GPS Data')
xlabel({'UTM_Easting'})
ylabel({'UTM_Northing'})


%% Straight Line Altitude

lne_alt = (cellfun(@(m) double(m.Altitude),lne_struct));
lne_time = (cellfun(@(m) double(m.Header.Stamp.Sec),lne_struct)-double(lne_struct{1}.Header.Stamp.Sec));
lne_atl_mean = mean((cellfun(@(m) double(m.Altitude),lne_struct)));

figure
scatter(lne_time,lne_alt)
hold on
plot(lne_time,lne_atl_mean,'--ro','LineWidth',0.5)
hold off
title('Straight Line Altitude Data Plot')
xlabel({'Time','(in seconds)'})
ylabel({'Altitude','(in meters)'})
legend('Altitude Data','Mean Line')

figure
plot(lne_time,lne_err)
title('Straight Line Error Plot')
xlabel({'Time','(in seconds)'})
ylabel({'Error','(in meters)'})
legend('Spread of Residual in Straigth Line data')


%% defined fuctions 

function dist = point_dist(x3, y3)
    x1 = 18.3492;
    y1 = 125.416;
    x2 = 113.852;
    y2 = 105.53;
    num =abs((x2-x1)*(y1-y3)-(x1-x3)*(y2-y1));
    den = sqrt((x2-x1)^2 + (y2-y1)^2);
    dist = num ./ den;
    return
end
