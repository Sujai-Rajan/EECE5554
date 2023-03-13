stn_opn_data = rosbag('stationary_open.bag');
mov_opn_data = rosbag('walking_open.bag');
stn_opn_data_occ = rosbag('stationary_occ.bag');
mov_opn_data_occ = rosbag('walking_occ.bag');

stn_opn_data_TopicData = select(stn_opn_data,'Topic','/gps');
mov_opn_data_TopicData = select(mov_opn_data,'Topic','/gps');
stn_opn_data_occ_TopicData = select(stn_opn_data_occ,'Topic','/gps');
mov_opn_data_occ_TopicData = select(mov_opn_data_occ,'Topic','/gps');

stn_opn = readMessages(stn_opn_data_TopicData,'DataFormat','struct');
mov_opn = readMessages(mov_opn_data_TopicData,'DataFormat','struct');
stn_opn_occ = readMessages(stn_opn_data_occ_TopicData,'DataFormat','struct');
mov_opn_occ = readMessages(mov_opn_data_occ_TopicData,'DataFormat','struct');

%Stationary Data Open
stn_opn_x_off = stn_opn{1}.UTMEasting;
stn_opn_y_off = stn_opn{1}.UTMNorthing;
stn_opn_x = cellfun(@(m) double(m.UTMEasting),stn_opn)-stn_opn_x_off;
stn_opn_y = cellfun(@(m) double(m.UTMNorthing),stn_opn)-stn_opn_y_off;
stn_opn_quality = cellfun(@(m) double(m.FixQuality),stn_opn);
stn_opn_x_Gearth = stn_opn_x_off-stn_opn_x_off;
stn_opn_y_Gearth = stn_opn_y_off-stn_opn_y_off;
stn_opn_xdeg = cellfun(@(m) double(m.Latitude),stn_opn)-0;
stn_opn_ydeg = cellfun(@(m) double(m.Longitude),stn_opn)-0;
Sx = std(stn_opn_x);
Sy = std(stn_opn_y);
CEP = 0.59*(Sx+Sy);
DRMS = (sqrt(Sx^2+Sy^2))*2;
stn_opn_x_mean = (mean(cellfun(@(m) double(m.UTMEasting),stn_opn)))-stn_opn_x_off;
stn_opn_y_mean = (mean(cellfun(@(m) double(m.UTMNorthing),stn_opn)))-stn_opn_y_off;
stn_opn_alt = cellfun(@(m) double(m.Altitude),stn_opn);

figure
histfit(stn_opn_alt,50);
title('Histogram for Stationary Altitude Data in Open Space')
xlabel({'Altitude in meters'})
ylabel({'Frequency','n'})
legend('data','curve');

figure
gscatter(stn_opn_x,stn_opn_y,stn_opn_quality,'br','xo')
hold on
scatter(stn_opn_x_mean,stn_opn_y_mean,'black+')
gscatter(stn_opn_x_Gearth,stn_opn_y_Gearth)
viscircles([stn_opn_x_mean stn_opn_y_mean],CEP,"Color",'green');
viscircles([stn_opn_x_mean stn_opn_y_mean],DRMS,"Color",'cyan');
hold off
grid on
title('Stationary GPS Data in Open Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'})
legend('Quality=4 ','Quality=5 ','mean','google earth true value','Location','northwest');

stn_opn_dataError = GetPointDistance(stn_opn_x_Gearth,stn_opn_y_Gearth,stn_opn_x,stn_opn_y);
stn_opn_data_RMSE = sqrt(mean(stn_opn_dataError)^2);
disp('RMSE for Stationary Open Data');
disp(stn_opn_data_RMSE);

figure
hh = histogram2(stn_opn_x+stn_opn_x_off,stn_opn_y+stn_opn_y_off);
title('Histogram for Stationary UTMEasting and UTMNorthing in Open Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'});
zlabel({'Frequency','n'})
legend('data');
counts = hh.BinCounts;
values = hh.Values;

figure
histfit(stn_opn_x+stn_opn_x_off,50);
title('Histogram for Stationary UTMEasting in Open Space')
xlabel({'UTMEasting in meters'})
ylabel({'Frequency','n'})
legend('data','curve');

figure
histfit(stn_opn_y+stn_opn_y_off,50);
title('Histogram for Stationary UTMNorthing in Open Space')
xlabel({'UTMNorthing in meters'})
ylabel({'Frequency','n'})
legend('data','curve');
wm = webmap('World Imagery');
s = geoshape(stn_opn_xdeg,stn_opn_ydeg);
hold on
wmline(s,'Color','red','Width',10);
hold off

%Moving Data Open
movopn_x_off = mov_opn{1}.UTMEasting;
movopn_y_off = mov_opn{1}.UTMNorthing;
movopn_x = cellfun(@(m) double(m.UTMEasting),mov_opn)-movopn_x_off;
movopn_y = cellfun(@(m) double(m.UTMNorthing),mov_opn)-movopn_y_off;
movopn_q = cellfun(@(m) double(m.FixQuality),mov_opn);
movopn_xdeg = cellfun(@(m) double(m.Latitude),mov_opn)-0;
movopn_ydeg = cellfun(@(m) double(m.Longitude),mov_opn)-0;
movopn_alt = cellfun(@(m) double(m.Altitude),mov_opn);


figure
histfit(movopn_alt,50);
title('Histogram for Moving Altitude Data in Open Space')
xlabel({'Altitude in meters'})
ylabel({'Frequency','n'})
legend('data','curve');

figure
gscatter(movopn_x,mov_opn_data_Y,movopn_q,'br','oo')
grid on
title('Moving GPS Data in Open Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'})
legend('Quality=4 ','Quality=5 ','Location','northwest');

p = fittype('Poly1');
[fitresult1, gof1] = fit(movopn_x(1:40),mov_opn_data_Y(1:40), p);
[fitresult2, gof2] = fit(movopn_x(41:48),mov_opn_data_Y(41:48), p);
[fitresult3, gof3] = fit(movopn_x(49:86),mov_opn_data_Y(49:86), p);
[fitresult4, gof4] = fit(movopn_x(87:95),mov_opn_data_Y(87:95), p);
[fitresult5, gof5] = fit(movopn_x(96:103),mov_opn_data_Y(96:103), p);

figure
plot(fitresult1,movopn_x(1:40),mov_opn_data_Y(1:40),'predoba')
hold on
plot(fitresult2,movopn_x(41:48),mov_opn_data_Y(41:48),'predoba')
plot(fitresult3,movopn_x(49:86),mov_opn_data_Y(49:86),'predoba')
plot(fitresult4,movopn_x(87:95),mov_opn_data_Y(87:95),'predoba')
plot(fitresult5,movopn_x(96:103),mov_opn_data_Y(96:103),'predoba')
hold off
title('Best fit plot for Moving Data in Open Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'})
legend('Data','Best Fit','Prediction Bounds','Location','northwest');
mov_opn_data_RMSE = (gof1.rmse + gof2.rmse + gof3.rmse + gof4.rmse + gof5.rmse)/5;
disp('RMSE for Moving Open Data');
disp(mov_opn_data_RMSE);
wm2 = webmap('World Imagery');
s = geoshape(movopn_xdeg,movopn_ydeg);
hold on
wmline(s,'Color','red','Width',3);
hold off

%Stationary Data Occluded
stn_opn_data_occ_X_Offset = stn_opn_occ{1}.UTMEasting;
stn_opn_data_occ_Y_Offset = stn_opn_occ{1}.UTMNorthing;
stn_opn_data_occ_X = cellfun(@(m) double(m.UTMEasting),stn_opn_occ)-stn_opn_data_occ_X_Offset;
stn_opn_data_occ_Y = cellfun(@(m) double(m.UTMNorthing),stn_opn_occ)-stn_opn_data_occ_Y_Offset;
stn_opn_data_occ_Quality = cellfun(@(m) double(m.FixQuality),stn_opn_occ);

stn_opn_data_occ_X_Gearth = stn_opn_data_occ_X_Offset - stn_opn_data_occ_X_Offset;
stn_opn_data_occ_Y_Gearth = stn_opn_data_occ_Y_Offset - stn_opn_data_occ_Y_Offset;
stn_opn_data_occ_Xdeg = cellfun(@(m) double(m.Latitude),stn_opn_occ)-0;
stn_opn_data_occ_Ydeg = cellfun(@(m) double(m.Longitude),stn_opn_occ)-0;
stn_opn_data_occ_Sx = std(stn_opn_data_occ_X);
stn_opn_data_occ_Sy = std(stn_opn_data_occ_Y);
disp('std :');
disp(stn_opn_data_occ_Sx);
disp(stn_opn_data_occ_Sy);
stn_opn_data_occ_CEP = 0.59*(stn_opn_data_occ_Sx+stn_opn_data_occ_Sy);
disp(stn_opn_data_occ_CEP);
stn_opn_data_occ_DRMS = (sqrt(stn_opn_data_occ_Sx^2+stn_opn_data_occ_Sy^2))*2;
stn_opn_data_occ_MeanX = (mean(cellfun(@(m) double(m.UTMEasting),stn_opn_occ)))-stn_opn_data_occ_X_Offset;
stn_opn_data_occ_MeanY = (mean(cellfun(@(m) double(m.UTMNorthing),stn_opn_occ)))-stn_opn_data_occ_Y_Offset;
stn_opn_data_occ_Altitude = cellfun(@(m) double(m.Altitude),stn_opn_occ);

figure
histfit(stn_opn_data_occ_Altitude,50);
title('Histogram for Stationary Altitude Data in Occluded Space')
xlabel({'Altitude in meters'})
ylabel({'Frequency','n'})
legend('data','curve');

figure
gscatter(stn_opn_data_occ_X,stn_opn_data_occ_Y,stn_opn_data_occ_Quality,'grk','*ox')
hold on
scatter(stn_opn_data_occ_MeanX,stn_opn_data_occ_MeanY,'cyan+')
gscatter(stn_opn_data_occ_X_Gearth,stn_opn_data_occ_Y_Gearth)
viscircles([stn_opn_data_occ_MeanX stn_opn_data_occ_MeanY],stn_opn_data_occ_CEP,"Color",'green');
viscircles([stn_opn_data_occ_MeanX stn_opn_data_occ_MeanY],stn_opn_data_occ_DRMS,"Color",'cyan');
hold off
grid on
title('Stationary GPS Data in Occluded Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'})
legend('Quality =4 ','Quality = 5','mean','google earth true value','Location','northeast');
stn_opn_data_occError = GetPointDistance(stn_opn_data_occ_X_Gearth,stn_opn_data_occ_Y_Gearth,stn_opn_data_occ_X,stn_opn_data_occ_Y);
stn_opn_data_occ_RMSE = sqrt(mean(stn_opn_data_occError)^2);
disp('RMSE for Stationary Occluded Data');
disp(stn_opn_data_occ_RMSE);

figure
hh = histogram2(stn_opn_data_occ_X+stn_opn_data_occ_X_Offset,stn_opn_data_occ_Y+stn_opn_data_occ_Y_Offset);
title('Histogram for Stationary UTMEasting and UTMNorthing in Occluded Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'});
zlabel({'Frequency','n'})
legend('data');

figure
histfit(stn_opn_data_occ_X+stn_opn_data_occ_X_Offset,50);
title('Histogram for Stationary UTMEasting in Occluded Space')
xlabel({'UTMEasting in meters'})
ylabel({'Frequency','n'})
legend('data','curve');

figure
histfit(stn_opn_data_occ_Y+stn_opn_data_occ_Y_Offset,50);
title('Histogram for Stationary UTMNorthing in Occluded Space')
xlabel({'UTMNorthing in meters'})
ylabel({'Frequency','n'})
legend('data','curve');
wm3 = webmap('World Imagery');
s = geoshape(stn_opn_data_occ_Xdeg,stn_opn_data_occ_Ydeg);
hold on
wmline(s,'Color','red','Width',10);
hold off

%Moving Data Occluded
movopn_occ_x_off = mov_opn_occ{1}.UTMEasting;
movopn_occ_y_off = mov_opn_occ{1}.UTMNorthing;
movopn_occ_x = cellfun(@(m) double(m.UTMEasting),mov_opn_occ)-movopn_occ_x_off;
movopn_occ_y = cellfun(@(m) double(m.UTMNorthing),mov_opn_occ)-movopn_occ_y_off;
movopn_occ_quality = cellfun(@(m) double(m.FixQuality),mov_opn_occ);
movopn_occ_x_deg = cellfun(@(m) double(m.Latitude),mov_opn_occ)-0;
movopn_occ_ydeg = cellfun(@(m) double(m.Longitude),mov_opn_occ)-0;
movopn_occ_alt = cellfun(@(m) double(m.Altitude),mov_opn_occ);


figure
histfit(movopn_occ_alt,50);
title('Histogram for Moving Altitude Data in Occluded Space')
xlabel({'Altitude in meters'})
ylabel({'Frequency','n'})
legend('data','curve');

figure
gscatter(mov_opn_data_occ_X,movopn_occ_y,movopn_occ_quality,'br','oo')
grid on
title('Moving GPS Data in Occluded Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'})
legend('Quality=4','Location','northeast');
p = fittype('Poly1');
[Occluded_fitresult1, Occluded_gof1] = fit(mov_opn_data_occ_X(1:20),movopn_occ_y(1:20), p);
[Occluded_fitresult2, Occluded_gof2] = fit(mov_opn_data_occ_X(21:30),movopn_occ_y(21:30), p);
[Occluded_fitresult3, Occluded_gof3] = fit(mov_opn_data_occ_X(31:51),movopn_occ_y(31:51), p);
[Occluded_fitresult4, Occluded_gof4] = fit(mov_opn_data_occ_X(52:62),movopn_occ_y(52:62), p);
[Occluded_fitresult5, Occluded_gof5] = fit(mov_opn_data_occ_X(63:65),movopn_occ_y(63:65), p);

figure
plot(Occluded_fitresult1,mov_opn_data_occ_X(1:20),movopn_occ_y(1:20),'predoba')
hold on
plot(Occluded_fitresult2,mov_opn_data_occ_X(21:30),movopn_occ_y(21:30),'predoba')
plot(Occluded_fitresult3,mov_opn_data_occ_X(31:51),movopn_occ_y(31:51),'predoba')
plot(Occluded_fitresult4,mov_opn_data_occ_X(52:62),movopn_occ_y(52:62),'predoba')
plot(Occluded_fitresult5,mov_opn_data_occ_X(63:65),movopn_occ_y(63:65),'predoba')
hold off
title('Best fit plot for Moving Data in Occluded Space')
xlabel({'UTMEasting in meters'})
ylabel({'UTMNorthing in meters'})
grid on
legend('Data','Best Fit','Prediction Bounds','Location','northeast');
mov_opn_data_occ_RMSE = (Occluded_gof1.rmse + Occluded_gof2.rmse + Occluded_gof3.rmse + Occluded_gof4.rmse + Occluded_gof5.rmse)/5;
disp('RMSE for Moving Occluded Data');
disp(mov_opn_data_occ_RMSE);
wm4 = webmap('World Imagery');
s = geoshape(movopn_occ_x_deg,movopn_occ_ydeg);
hold on
wmline(s,'Color','red','Width',3);
hold off

function distance = GetPointDistance(x1,y1,x2,y2)
    distance = sqrt((x2-x1).^2 + (y2-y1).^2);
    return
end
