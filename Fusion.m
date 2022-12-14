mag_fl = readtable('fusion-mag_message.csv');
imu_fl = readtable('fusion-imu_message.csv');
gps_fl = readtable('fusion-gps_message.csv');

cl_mag_x = table2array(mag_fl(:,6));
cl_mag_y = table2array(mag_fl(:,7));
cl_time = table2array(mag_fl(:,1));
cl_vel_z = table2array(imu_fl(:,13));
cl_vel_x = table2array(imu_fl(:,11));
cl_orientation_x = table2array(imu_fl(:,6));
cl_orientation_y = table2array(imu_fl(:,7));
cl_orientation_z = table2array(imu_fl(:,8));
cl_orientation_w = table2array(imu_fl(:,9));
cl_acc_x = table2array(imu_fl(:,15));
cl_acc_y = table2array(imu_fl(:,16));
cl_utm_easting = table2array(gps_fl(:,6));
cl_utm_northing = table2array(gps_fl(:,7));
gps_time = table2array(gps_fl(:,1));

mag_x = cl_mag_x(7040:11280,:);
mag_y = cl_mag_y(7040:11280,:);
time_val = cl_time(7040:11280,:);
vel_z = cl_vel_z(7040:11280,:);
vel_x = cl_vel_x(7040:11280,:);
orientation_x = cl_orientation_x(7040:11280,:);
orientation_y = cl_orientation_y(7040:11280,:);
orientation_z = cl_orientation_z(7040:11280,:);
orientation_w = cl_orientation_w(7040:11280,:);
%% 

%Hard iron factor elimination
alpha = (max(mag_x)+min(mag_x))/2;
beta = (max(mag_y)+min(mag_y))/2;

hard_corrected_mag_value = [(mag_x-alpha), (mag_y-beta)];

%Soft iron factor elimination
temp_mag_value = cat(2,hard_corrected_mag_value(:,1),hard_corrected_mag_value(:,2));
%% 
[z,a,b,al] = fitellipse(temp_mag_value);
%% 
z2 = z;
a2 = a;
b2 = b;
al2 = al;
R = [cos(al2+pi), sin(al2+pi); -sin(al2+pi), cos(al2+pi)];
% R = [cos(al), -sin(al); sin(al), cos(al)];
temp2 = cat(2,hard_corrected_mag_value(:,1),hard_corrected_mag_value(:,2));
tilt_correct_temp_mag = R*temp2';
tilt_correct_mag=tilt_correct_temp_mag';
tilt_correct_mag_x = tilt_correct_mag(:,1);
tilt_correct_mag_y = tilt_correct_mag(:,2);
%% 

[z1,a1,b1,al1] = fitellipse(tilt_correct_mag);  %to check for the major axis and
%minior axis values are correct or not.
%% 
z3 = z1;
a3 = a1;
b3 = b1;
al3 = al1;
gamma = a2/b2;
soft_corrected_mag_x = (tilt_correct_mag_x/gamma);
soft_corrected_mag_y = (tilt_correct_mag_y);
%% 

figure;
scatter(mag_x,mag_y,"red")
hold on;
plotellipse(z,a,b,al)
xlabel("Original Mag-x data")
ylabel("Original Mag-y data")
title('Original Magnatometer Data')
grid on
daspect([1 1 1]);
hold off;

figure;
scatter(hard_corrected_mag_value(:,1),hard_corrected_mag_value(:,2))
hold on;
plotellipse(z,a,b,al)
xlabel("Hard Iron corrected Mag-X")
ylabel("Hard Iron corrected Mag-y")
title('Hard Iron Corrected Data')
grid on
daspect([1 1 1])
hold off;

figure;
scatter(soft_corrected_mag_x,soft_corrected_mag_y,"green")
hold on;
plotellipse(z1,a1,b1,al1)
xlabel("soft-iron corrected mag-x")
ylabel("soft iron corrected mag-y")
title("soft iron corrected")
grid on
daspect([1 1 1])
hold off;

%% 
hard_corrected = [(cl_mag_x-alpha), (cl_mag_y-beta)];
temp3 = cat(2,hard_corrected(:,1),hard_corrected(:,2));
tilt_corrected_temp = R*temp3';
tilt_corrected = tilt_corrected_temp';
tilt_corrected_x = tilt_corrected(:,1);
tilt_corrected_y = tilt_corrected(:,2);

soft_corrected_x = (tilt_corrected_x/gamma);
soft_corrected_y = (tilt_corrected_y);


%% 
%Calculate magnetometer yaw angle
yaw_angle = atan2(soft_corrected_y,soft_corrected_x);
yaw_angle_unwraped = unwrap(yaw_angle);
shifted_yaw_angle = yaw_angle_unwraped - min(yaw_angle_unwraped);
deg_yaw_angle = rad2deg(shifted_yaw_angle);

figure;
plot(cl_time,deg_yaw_angle)
xlabel("Time")
ylabel("Magnitude Yaw angle")
title("Yaw angle plot of magnetometer with time")
grid on

%Calculate yaw angle for gyro

vel_yaw_angle = -(cumtrapz(cl_time,cl_vel_z));
 
figure;
plot(cl_time,vel_yaw_angle)
% hold on
% plot(time_val,deg_mag_yaw_angle)
xlabel("Time")
ylabel("Yaw angle")
title("Yaw angle plot of gyroscope with time")
grid on
% hold off

%Using Lowpass Filter for magnetometer readings
filtered_yaw_angle1 = lowpass(deg_yaw_angle,0.01,40);
%filtered_mag_yaw_angle2 = lowpass(deg_mag_yaw_angle,0.1,40);

%subplot(1,2,1);
% figure;
% plot(cl_time,filtered_yaw_angle1)
% hold on
% plot(cl_time,deg_yaw_angle)
% xlabel('Time')
% ylabel("Lowpass filtered yaw angle")
% title('Yaw angle plot of magnetometer with time')
% grid on
% daspect([1 1 1]);
% hold off

% subplot(1,2,2);
% plot(time_val,filtered_mag_yaw_angle2)
% xlabel('Time')
% ylabel("Lowpass filtered yaw angle")
% title('Yaw angle plot of magnetometer with time')
% grid on
% daspect([1 1 1]);

%Using Highpass filter for gyro readings
filtered_vel_yaw_angle = highpass(vel_yaw_angle,0.01,40);

% figure;
% plot(cl_time,filtered_vel_yaw_angle)
% hold on
% plot(cl_time,vel_yaw_angle)
% xlabel('Time')
% ylabel("Lowpass filtered yaw angle")
% title('Yaw angle plot of gyroscope with time')
% grid on
% daspect([1 1 1]);
% hold off

added_signal = filtered_yaw_angle1 + filtered_vel_yaw_angle;

figure;
plot(cl_time,added_signal)
hold on
plot(cl_time,deg_yaw_angle)
plot(cl_time,vel_yaw_angle)
xlabel('Time')
ylabel("Complementry filtered signal")
title('Corrected Signal')
legend('Complimentory filtered signal','Magnetometer reading','Velocity reading')
grid on
% daspect([1 1 1]);
hold off

eul = quat2eul([cl_orientation_w,cl_orientation_x,cl_orientation_y,cl_orientation_z]);
yaw_eul = -unwrap(eul(:,1));

figure;
plot(cl_time,added_signal)
hold on
plot(cl_time,yaw_eul)
xlabel('Time')
ylabel("Complementry filtered signal")
title('Filtered Signal vs IMU yaw angle data plot')
legend('Complimentory filtered signal','Euler Yaw angle data')
grid on
% daspect([1 1 1]);
hold off


%% 
%Lab4 part-2
%Calculating velocity from linear acceleration
init_lnr_vl = cl_acc_x - mean(cl_acc_x);
lnr_vel = cumtrapz(cl_time,init_lnr_vl);

figure;
plot(cl_time,lnr_vel)
xlabel("Time")
ylabel("Linear Velocity")
title("Integrated Velocity Plot")
grid on
%daspect([1 1 1]);

%Calculating Velocity from gps plot
len = length(cl_utm_easting);
for i = 1: len-1
    sq_diff1 = (cl_utm_easting(i+1,:)-cl_utm_easting(i,:))^2;
    sq_diff2 = (cl_utm_northing(i+1,:)-cl_utm_northing(i,:))^2;
    dist = sqrt(sq_diff1 + sq_diff2);
    gps_velocity(i) = dist;
end

real_gps_vel = gps_velocity';
temp_gps_vel = real_gps_vel;
real_gps_vel(end+1) = 1;

figure;
plot(gps_time,real_gps_vel)
xlabel('time')
ylabel('Velocity from GPS data')
title("Velocity Calculated from GPS data")
grid on

scale1 = cl_time(end)/gps_time(end);
temp_gps_time = gps_time*scale1;

scale2 = (max(lnr_vel)/max(real_gps_vel));
temp_gps_vel = real_gps_vel*scale2;

figure;
plot(cl_time,lnr_vel);
hold on
plot(temp_gps_time,temp_gps_vel)
legend('Integrated Velocity','Velocity from GPS data')
xlabel('Time')
ylabel('Linear Velocity')
title("Velocity calculated from accelaration and gps data")
grid on
hold off
%% 

for i = 1: 6475
    cl_acc_x(i) = 0;  
end

m1 = mean(cl_acc_x(6476:11952));

for i = 6476: 11952
    cl_acc_x(i) = (cl_acc_x(i) - m1);
end

for i = 11952: 13924
    cl_acc_x(i) = 0;
end

m2 = mean(cl_acc_x(13925:16309));

for i = 13925: 16369
    cl_acc_x(i) = (cl_acc_x(i) - m2);
end

for i = 16370: 17109
    cl_acc_x(i) = 0;
end

m3 = mean(cl_acc_x(17110:20107));

for i = 17110:20107
    cl_acc_x(i) = (cl_acc_x(i) - m3);
end

for i = 20108:20679
    cl_acc_x(i) = 0;
end

m4 = mean(cl_acc_x(20108:22083));

for i = 20108:22083
    cl_acc_x(i) = (cl_acc_x(i) - m4);
end

for i = 22084:25135
    cl_acc_x(i) = 0;
end

m5 = mean(cl_acc_x(25136:27249));

for i = 25136:27249
    cl_acc_x(i) = (cl_acc_x(i) - m5);
end

for i = 27249:28538
    cl_acc_x(i) = 0;
end

m6 = mean(cl_acc_x(28539:33003));

for i = 28539:33003
    cl_acc_x(i) = (cl_acc_x(i) - m6);
end

for i = 33003:34856
    cl_acc_x(i) = 0;
end

m7 = mean(cl_acc_x(34857:35546));

for i = 34857:35546
    cl_acc_x(i) = (cl_acc_x(i) - m7);
end

m8 = mean(cl_acc_x(35547:36090));

for i = 35547:36090
    cl_acc_x(i) = (cl_acc_x(i) - m8);
end

m9 = mean(cl_acc_x(36091:36881));

for i = 36091:36881
    cl_acc_x(i) = (cl_acc_x(i) - m9);
end

for i = 36882:37322
    cl_acc_x(i) = 0;
end

m10 = mean(cl_acc_x(37323:43731));

for i = 37323:43731
    cl_acc_x(i) = (cl_acc_x(i) - m10);
end

for i = 43732:44891
    cl_acc_x(i) = 0;
end

m11 = mean(cl_acc_x(44892:50114));

for i = 44584:50114
    cl_acc_x(i) = (cl_acc_x(i) - m11);
end

for i = 50114: 51247
    cl_acc_x(i) = 0;
end

figure;
plot(cl_time,cl_acc_x)
xlabel('time')
ylabel('Accelaration')
grid on

corrected_lnr_vel = cumtrapz(cl_time,cl_acc_x);

figure;
plot(cl_time,corrected_lnr_vel);
hold on
plot(temp_gps_time,temp_gps_vel)
legend('Integrated Velocity','Velocity from GPS data')
xlabel('Time')
ylabel('Linear Velocity')
title("Velocity calculated from accelaration and gps data")
grid on
hold off
%% 
%Lab4 - Part-3
x_dot = lnr_vel;
omega = cl_vel_z;
outpt = omega.*x_dot;

% new_acc_y = (cl_acc_y - mean(cl_acc_y));
y_double_dot = cl_acc_y;

figure;
plot(cl_time,outpt);
hold on
plot(cl_time,50*y_double_dot)
legend('Multiplied result','Linear Acceleration')
xlabel('Time')
ylabel('Function')
title('Multiplied result vs linear acceleration')
grid on
hold off

%% 

xcordi = corrected_lnr_vel.*(cos(deg2rad(yaw_eul)));
ycordi = corrected_lnr_vel.*(sin(deg2rad(yaw_eul)));
x_disp = cumtrapz(cl_time,xcordi);
y_disp = cumtrapz(cl_time,ycordi);

figure
plot(x_disp,y_disp)
xlabel('Displacement in X')
ylabel('Displacement in Y')
title('Calculated Displacement')
grid on

thet = 66;
rot = [cos(thet), sin(thet); -sin(thet), cos(thet)];

new_cordi = cat(2,x_disp,y_disp);
temp_new_disp = rot*new_cordi';
new_disp = temp_new_disp';

new_x_disp = new_disp(:,1);
new_y_disp = new_disp(:,2);

figure
plot(new_x_disp,new_y_disp)
xlabel('Displacement in X')
ylabel('Displacement in Y')
title('Calculated Displacement')
grid on

figure
plot(cl_utm_easting,cl_utm_northing)
xlabel('utm_easting')
ylabel('utm_northing')
title('GPS plot')
grid on
%%
temp_lnr_vel = lnr_vel;
m_y_double_dot = y_double_dot - mean(y_double_dot);
n = 1;
for i = 1:51247
    if corrected_lnr_vel(i) < 0.4
        continue
    else
        xc(n) = ((m_y_double_dot(i) - outpt(i))./corrected_lnr_vel(i));
        n = n + 1;
    end
end
disp(mean(xc))
