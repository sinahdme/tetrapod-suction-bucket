clc;
clear all;
close all;
%% %%%%%%%%%% Geometric parameters %%%%%%%%%%%%
do_b = 6.5; %bucket outer diameter
thickness = 0.025; %bucket thickness
l_b = 6.5; %length of the bucket skirt
Z_lid =0.0001; % thickness of the lid of bucket
n = 32; %number of strips in a single circumference
n_l =32; % number of rings in the bottom
n_q =12;
eccentricity = 33;
th1= 45;
th2 = 135;
th3 = 225;
th4 = 315;
Deadload=5000; %kN
l = 19.5;
%% %%%%%%%%%% Soil parameters %%%%%%%%%%%%
gamma = 9.21; %kN/m3
S_gamma = 0.5;
e_max =0;
e_min =0;
D_r = 0.75;
m_c = 161.42*D_r^2+199.8*D_r+36.877;
A_c = 0.3474*D_r^2+0.4222*D_r+0.328;
fi_crit = 30;
m =0;
delta = 22.6;
%%
ro_b = do_b/2; %outer radious of the bucket
di_b =do_b-2*thickness; %bucket inner diameter
ri_b = di_b/2; %outer radious of the bucket
r_pile = di_b/2 + thickness/2; %from the center of buckt to the center of the thickness
r_b = abs(di_b-do_b)/4; %radius of the pile
penetration =l_b; %fully penetrated
circle_split = linspace(0,360,n+1);
R_split = linspace(0,r_pile,n_q);
penetration1 =l_b; %fully penetrated
penetration2 =l_b; %fully penetrated
penetration3 =l_b; %fully penetrated
penetration4 =l_b; %fully penetrated
%%
ro_b = do_b/2; %outer radious of the bucket
di_b =do_b-2*thickness; %bucket inner diameter
ri_b = di_b/2; %outer radious of the bucket
r_pile = ri_b/2 + thickness/2; %from the center of buckt to the center of the thickness
r_b = abs(di_b-do_b)/4; %radius of the pile

circle_split = linspace(0,360,n+1);
R_split = linspace(0,r_pile,n_q);
R_sec = R_split(2)-R_split(1);

%% %%%%%%%%%% Depth sweep parameters %%%%%%%%%%%%
depth_values = [6.0, 6.1, 6.2, 6.35, 6.5];  % penetration depths to sweep
num_depths = length(depth_values);

%% %%%%%%%%%% Tilt direction sweep parameters %%%%%%%%%%%%
tilt_directions = 0:30:330;  % 12 tilt directions in degrees
dd_steps = 0.25*pi/180:0.25*pi/180:1.5*pi/180;  % 0.25 degree steps up to 1.5 degrees
num_dirs = length(tilt_directions);
num_steps = length(dd_steps);
M_polar = zeros(num_dirs, num_steps, num_depths);  % Store moment for each (direction, step, depth)
F_H_polar = zeros(num_dirs, num_steps, num_depths); % Store horizontal force

%% Code initiation
kkk =0;
abc=0;
cc1=0;
% full penetration action
    abc=abc+1;
 cc1=cc1+1;

%% $$$$$$$$$$$$$$$$$$$$$$$$$ DEPTH LOOP $$$$$$$$$$$$$$$$$$$$$$$$$ %%
for i_depth = 1:num_depths
    l_b = depth_values(i_depth);
    penetration = l_b;
    penetration1 = l_b;
    penetration2 = l_b;
    penetration3 = l_b;
    penetration4 = l_b;
    fprintf('\n\n########## DEPTH: %.1f m (%d/%d) ##########\n', l_b, i_depth, num_depths);

%%  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ dealload effect with no tilt angle $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  %%
dz=0;
SUMz=0;
QQ = zeros(1,10000);
while true
    if SUMz>Deadload
    dz=dz-0.00001;
    else
     dz=dz+0.00001;
    end
for dd=0
   q=[0 -dd];
phi = q(1);     % Rotation around x-axis
theta = q(2);   % Rotation around y-axis
psi = 0;     % Rotation around z-axis
j=0;
for k=0:l_b/n_l:l_b-(l_b/n_l)
for i=2:n+1
    j=j+1;
    first_s_end(:,j) = [l*cosd(th1)+ro_b*cosd(circle_split(i)) l*sind(th1)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame
    second_s_end(:,j) = [l*cosd(th2)+ro_b*cosd(circle_split(i)) l*sind(th2)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame
    third_s_end(:,j) = [l*cosd(th3)+ro_b*cosd(circle_split(i)) l*sind(th3)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame
    fourth_s_end(:,j) = [l*cosd(th4)+ro_b*cosd(circle_split(i)) l*sind(th4)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame

    % notice: the stars are in the middle of layers.
end
end
% w. r. t.  fixed frame
A1 = zeros(size(first_s_end));
A2 = zeros(size(second_s_end));
A3 = zeros(size(third_s_end));
A4 = zeros(size(fourth_s_end));
A1(3,:) = -dz;
A2(3,:) = -dz;
A3(3,:) = -dz;
A4(3,:) = -dz;
first_s_end_fixed       = A1 + first_s_end;
second_s_end_fixed = A2+ second_s_end;
third_s_end_fixed     = A3 + third_s_end;
fourth_s_end_fixed   = A4 + fourth_s_end;

j=0;
for i=2:n+1
    j=j+1;
    Pile_first_s_end(:,j) = [l*cosd(th1)+r_pile*cosd(circle_split(i)) l*sind(th1)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
        Pile_second_s_end(:,j) = [l*cosd(th2)+r_pile*cosd(circle_split(i)) l*sind(th2)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
            Pile_third_s_end(:,j) = [l*cosd(th3)+r_pile*cosd(circle_split(i)) l*sind(th3)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
                Pile_fourth_s_end(:,j) = [l*cosd(th4)+r_pile*cosd(circle_split(i)) l*sind(th4)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
end
B1 = zeros(size(Pile_first_s_end));
B2 = zeros(size(Pile_second_s_end));
B3 = zeros(size(Pile_third_s_end));
B4 = zeros(size(Pile_fourth_s_end));

B1(3,:) = -dz;
B2(3,:) = -dz;
B3(3,:) = -dz;
B4(3,:) = -dz;

Pile_first_s_end_fixed = Pile_first_s_end + B1;
Pile_second_s_end_fixed = Pile_second_s_end + B2;
Pile_third_s_end_fixed = Pile_third_s_end + B3;
Pile_fourth_s_end_fixed = Pile_fourth_s_end + B4;

j=0;
jj = 0;
co =0;
R_sec = R_split(2)-R_split(1);
for jj =2:size(R_split,2)
j=0;
% w.r.t body frame
for i=2:n+1
    j=j+1;
r_Pile_first_s_end(:,j+co*n) = [l*cosd(th1)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th1)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';
r_Pile_second_s_end(:,j+co*n) = [l*cosd(th2)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th2)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';
r_Pile_third_s_end(:,j+co*n) = [l*cosd(th3)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th3)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';
r_Pile_fourth_s_end(:,j+co*n) = [l*cosd(th4)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th4)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';

end
co = co+1;
end
j=0;
C1 = zeros(size(r_Pile_first_s_end));
C2 = zeros(size(r_Pile_second_s_end));
C3 = zeros(size(r_Pile_third_s_end));
C4 = zeros(size(r_Pile_fourth_s_end));

C1(3,:) = -dz;
C2(3,:) = -dz;
C3(3,:) = -dz;
C4(3,:) = -dz;

r_Pile_first_s_end_fixed = r_Pile_first_s_end+C1;
r_Pile_second_s_end_fixed = r_Pile_second_s_end+C2;
r_Pile_third_s_end_fixed = r_Pile_third_s_end+C3;
r_Pile_fourth_s_end_fixed = r_Pile_fourth_s_end+C4;

active_z_values1 = Zvalues(l_b,n,n_l,penetration1);
active_z_values2 = Zvalues(l_b,n,n_l,penetration2);
active_z_values3 = Zvalues(l_b,n,n_l,penetration3);
active_z_values4 = Zvalues(l_b,n,n_l,penetration4);

FI_peak1 = fi_finder(gamma,active_z_values1,D_r,fi_crit,m);
FI_peak2 = fi_finder(gamma,active_z_values2,D_r,fi_crit,m);
FI_peak3 = fi_finder(gamma,active_z_values3,D_r,fi_crit,m);
FI_peak4 = fi_finder(gamma,active_z_values4,D_r,fi_crit,m);

fi1 = FI_peak1;
fi2= FI_peak2;
fi3 = FI_peak3;
fi4 = FI_peak4;

fi_end1 = fi1(end);
fi_end2 = fi2(end);
fi_end3 = fi3(end);
fi_end4 = fi4(end);

displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3;
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3;
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3;
displacement4 = (fourth_s_end(3, :) - fourth_s_end_fixed(3, :)) * 1e3;

%% py and tz forces self-weight
[~,~,P_i_out1] = coeffinder_deadload(first_s_end,first_s_end_fixed,active_z_values1,l_b,n_l,fi1,n,do_b,gamma); %KN
[~,~,P_i_out2] = coeffinder_deadload(second_s_end,second_s_end_fixed,active_z_values2,l_b,n_l,fi2,n,do_b,gamma); %KN
[~,~,P_i_out3] = coeffinder_deadload(third_s_end,third_s_end_fixed,active_z_values3,l_b,n_l,fi3,n,do_b,gamma); %KN
[~,~,P_i_out4] = coeffinder_deadload(fourth_s_end,fourth_s_end_fixed,active_z_values4,l_b,n_l,fi4,n,do_b,gamma); %KN

[tz_in1,tz_out1]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,tz_out2]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,tz_out3]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in4,tz_out4]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

KKz1 = tz_in1+tz_out1;
KKz2 = tz_in2+tz_out2;
KKz3 = tz_in3+tz_out3;
KKz4 = tz_in4+tz_out4;

KKz_in = tz_in1+tz_in2+tz_in3+tz_in4;
KKz_out = tz_out1+tz_out2+tz_out3+tz_out4;
KKz =KKz_in+KKz_out;
%% qz forces self-weight
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_4] = qzcurve2(fi_end4,Pile_fourth_s_end,Pile_fourth_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi1(1),Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi2(1),Z_lid,r_Pile_second_s_end, r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi3(1),Z_lid,r_Pile_third_s_end, r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_4] = qzcurve4(fi4(1),Z_lid,r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz1= sum(KKz1)+sum(qb1_1)+sum(qb2_1);
SUMz2= sum(KKz2)+sum(qb1_2)+sum(qb2_2);
SUMz3= sum(KKz3)+sum(qb1_3)+sum(qb2_3);
SUMz4= sum(KKz4)+sum(qb1_4)+sum(qb2_4);
SUMz = SUMz1+SUMz2+SUMz3+SUMz4;
end
% (Deadload - SUMz)
 if abs(Deadload - SUMz)<Deadload
    break
 end

end
%% defining the local frams and relative
 %dz=0
mask1 = (active_z_values1 == 0);
mask2 = (active_z_values2 == 0);
mask3 = (active_z_values3 == 0);
mask4 = (active_z_values4 == 0);

fprintf('\n===== Deadload equilibrium found. dz = %f =====\n', dz);
fprintf('===== Starting polar sweep: %d directions x %d rotation steps =====\n', num_dirs, num_steps);

%%  $$$$$$$$$$$$$$$$$$$$$$$$$   LOOP OVER TILT DIRECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$ %%
for i_dir = 1:num_dirs
    tilt_dir = tilt_directions(i_dir);  % current tilt direction in degrees
    fprintf('\n\n========== TILT DIRECTION: %d degrees (%d/%d) ==========\n', tilt_dir, i_dir, num_dirs);

    Q = [0,0,0]';     %initial center of rotation
    count =0;
    QQ = zeros(1,10000);
    cc1=0;

%%  $$$$$$$$$$$$$$$$$$$$$$$$$   bucket rotation angle change step by step $$$$$$$$$$$$$$$$$$$$$$$$$ %%
for dd= dd_steps
% Build rotation matrix via Rodrigues formula (replaces Euler decomposition)
R_tilt = rodrigues_rotation(dd, tilt_dir + 180);  % original convention: tilt_dir=0 -> toward +x

first_s_lo = [l*cosd(th1), l*sind(th1), 0];
second_s_lo = [l*cosd(th2), l*sind(th2), 0];
third_s_lo = [l*cosd(th3), l*sind(th3), 0];
fourth_s_lo = [l*cosd(th4), l*sind(th4), 0];
translation_vector = [l,0,0];

%%%%%%%%%%%Suction arrangement angle and geometrica vectors%%%%%%%%%%%
tower_tip = [0,0,eccentricity]';
tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
SUMz =0;
count = count + 1;
ddd(count) = dd*180/pi;
%%%%%%%%%%%

active_z_values1 = Zvalues(l_b,n,n_l,penetration1);
active_z_values2 = Zvalues(l_b,n,n_l,penetration2);
active_z_values3 = Zvalues(l_b,n,n_l,penetration3);
active_z_values4 = Zvalues(l_b,n,n_l,penetration4);
ee=0;
ee_tar = eccentricity;
nnn=0;
mmm=0;
QQQ = zeros(1,10000);
ccc =0;
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  Outer While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%
while  abs((ee)-ee_tar)>0.002*ee_tar %while condition to find the z-axis of the center of rotation
        tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
ccc = ccc+1;

 if ((ee)-ee_tar)>0.0
 mmm=mmm-0.015;
 else
 mmm=mmm+0.015;
 end
if ccc==1
 Q = [nnn*cosd(angle_indicator), nnn*sind(angle_indicator), -(1.0*l_b)+mmm]';
else
 Q = [(nnn)*cosd(angle_indicator), nnn*sind(angle_indicator), -(1.0*l_b)+mmm]';
end
if rem(ccc,5)==0
fprintf('\rdir: %d | degree: %0.3f  | ((ee)-ee_tar): %0.2f |(SUMz-Deadload): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', tilt_dir, dd*180/pi, CON2,CON1,Q(1),Q(2),Q(3));
end
QQQ(ccc+2)=Q(3);
if abs(abs(QQQ(ccc+2))-abs(QQQ(ccc)))==0
     break
end

count2 =0;
j=0;
fi1=FI_peak1;
fi2=FI_peak2;
fi3=FI_peak3;
fi4=FI_peak4;
cc=0;
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%
 while abs(SUMz-Deadload)>(Deadload/4000)
CON1 = abs(SUMz-Deadload); %track the condition
tower_tip = [0,0,eccentricity]';
tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
cc=cc+1;
if SUMz-Deadload>0
nnn=nnn+0.02;
else
nnn=nnn-0.02;
end
Q = [(0+nnn)*cosd(angle_indicator), (0+nnn)*sind(angle_indicator), Q(3)]';% iterative increasing or decreasing the x-axis of the rotation center
fprintf('\rdir: %d | degree: %0.3f  | ((ee)-ee_tar): %0.3f |(Deadload - SUMz): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', tilt_dir, dd*180/pi,((ee)-ee_tar) ,CON1,Q(1),Q(2),Q(3));

QQ(cc+2)=Q(1);
if abs(abs(QQ(cc+2))-abs(QQ(cc)))==0
     break
end

[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Pile_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Pile_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end, r_Pile_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[fourth_s_end,fourth_s_end_in,fourth_s_end_fixed,fourth_s_end_fixed_in,Pile_fourth_s_end,Pile_fourth_s_end_fixed,r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed]=points(R_split,th4,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
%% relative points calculation inner while
[relative_first_s_end,relative_first_s_end_fixed,relative_first_s_end_in,relative_first_s_end_fixed_in]...
    =relative_points(first_s_end, first_s_end_in, first_s_end_fixed, first_s_end_fixed_in, Pile_first_s_end, Pile_first_s_end_fixed, r_Pile_first_s_end, r_Pile_first_s_end_fixed, th1,l,dz,Q,R_tilt);

[relative_second_s_end,relative_second_s_end_fixed,relative_second_s_end_in,relative_second_s_end_fixed_in]...
    =relative_points(second_s_end, second_s_end_in, second_s_end_fixed, second_s_end_fixed_in, Pile_second_s_end, Pile_second_s_end_fixed, r_Pile_second_s_end, r_Pile_second_s_end_fixed, th2,l,dz,Q,R_tilt);

[relative_third_s_end,relative_third_s_end_fixed,relative_third_s_end_in,relative_third_s_end_fixed_in]...
    =relative_points(third_s_end, third_s_end_in, third_s_end_fixed, third_s_end_fixed_in, Pile_third_s_end, Pile_third_s_end_fixed, r_Pile_third_s_end, r_Pile_third_s_end_fixed, th3,l,dz,Q,R_tilt);

[relative_fourth_s_end,relative_fourth_s_end_fixed,relative_fourth_s_end_in,relative_fourth_s_end_fixed_in]...
    =relative_points(fourth_s_end, fourth_s_end_in, fourth_s_end_fixed, fourth_s_end_fixed_in, Pile_fourth_s_end, Pile_fourth_s_end_fixed, r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed, th4,l,dz,Q,R_tilt);


%% py and tz forces inner while
displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement4 = (fourth_s_end(3, :) - fourth_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in4 = (fourth_s_end_in(3, :) - fourth_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[~,~,P_i_out2] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[~,~,P_i_out3] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
[~,~,P_i_out4] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4); %KN

P_i_out1(mask1) = 0;
P_i_out2(mask2) = 0;
P_i_out3(mask3) = 0;
P_i_out4(mask4) = 0;

[~,tz_out1]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out2]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out3]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out4]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

[tz_in1,~]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement_in1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,~]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement_in2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,~]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement_in3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in4,~]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement_in4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

tz_out1(mask1) = 0;
tz_out2(mask2) = 0;
tz_out3(mask3) = 0;
tz_out4(mask4) = 0;

tz_in1(mask1) = 0;
tz_in2(mask2) = 0;
tz_in3(mask3) = 0;
tz_in4(mask4) = 0;

KKz_in1 = tz_in1;
KKz_in2 = tz_in2;
KKz_in3 = tz_in3;
KKz_in4 = tz_in4;

KKz_out1 = tz_out1;
KKz_out2 = tz_out2;
KKz_out3 = tz_out3;
KKz_out4 = tz_out4;

KKz1 =KKz_in1+KKz_out1;
KKz2 =KKz_in2+KKz_out2;
KKz3 =KKz_in3+KKz_out3;
KKz4 =KKz_in4+KKz_out4;

%% qz forces inner while
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_4] = qzcurve2(fi_end4,Pile_fourth_s_end,Pile_fourth_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi_end1,Z_lid,r_Pile_first_s_end_fixed, r_Pile_first_s_end_fixed,gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi_end2,Z_lid,r_Pile_second_s_end,r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi_end3,Z_lid,r_Pile_third_s_end,r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_4] = qzcurve4(fi_end4,Z_lid,r_Pile_fourth_s_end,r_Pile_fourth_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1)+...
             sum(KKz2)+sum(qb1_2)+sum(qb2_2)+...
             sum(KKz3)+sum(qb1_3)+sum(qb2_3)+...
             sum(KKz4)+sum(qb1_4)+sum(qb2_4);


 end
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  End of inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%
%% points  calculation outer while
[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Pile_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Pile_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end,r_Pile_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[fourth_s_end,fourth_s_end_in,fourth_s_end_fixed,fourth_s_end_fixed_in,Pile_fourth_s_end,Pile_fourth_s_end_fixed,r_Pile_fourth_s_end,r_Pile_fourth_s_end_fixed]=points(R_split,th4,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);

%% py and tz forces outer while
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[KKx2,KKy2,~] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[KKx3,KKy3,~] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
[KKx4,KKy4,~] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4); %KN

KKx1(mask1) = 0;
KKy1(mask1) = 0;
KKx2(mask2) = 0;
KKy2(mask2) = 0;
KKx3(mask3) = 0;
KKy3(mask3) = 0;
KKx4(mask4) = 0;
KKy4(mask4) = 0;

SUMx = sum(KKx1)+sum(KKx2)+sum(KKx3)+sum(KKx4);
SUMy = sum(KKy1)+sum(KKy2)+sum(KKy3)+sum(KKy4);

displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement4 = (fourth_s_end(3, :) - fourth_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in4 = (fourth_s_end_in(3, :) - fourth_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[~,~,P_i_out2] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[~,~,P_i_out3] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
[~,~,P_i_out4] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4); %KN

P_i_out1(mask1) = 0;
P_i_out2(mask2) = 0;
P_i_out3(mask3) = 0;
P_i_out4(mask4) = 0;

[~,tz_out1]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out2]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out3]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out4]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

[tz_in1,~]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,~]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,~]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in4,~]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);


tz_out1(mask1) = 0;
tz_out2(mask2) = 0;
tz_out3(mask3) = 0;
tz_out4(mask4) = 0;

tz_in1(mask1) = 0;
tz_in2(mask2) = 0;
tz_in3(mask3) = 0;
tz_in4(mask4) = 0;

KKz_in1 = tz_in1;
KKz_in2 = tz_in2;
KKz_in3 = tz_in3;
KKz_in4 = tz_in4;

KKz_out1 = tz_out1;
KKz_out2 = tz_out2;
KKz_out3 = tz_out3;
KKz_out4 = tz_out4;

KKz1 =KKz_in1+KKz_out1;
KKz2 =KKz_in2+KKz_out2;
KKz3 =KKz_in3+KKz_out3;
KKz4 =KKz_in4+KKz_out4;

%% qz forces outer while
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_4] = qzcurve2(fi_end4,Pile_fourth_s_end,Pile_fourth_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi_end1,Z_lid,r_Pile_first_s_end,r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi_end2,Z_lid,r_Pile_second_s_end,r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi_end3,Z_lid,r_Pile_third_s_end,r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_4] = qzcurve4(fi_end4,Z_lid,r_Pile_fourth_s_end,r_Pile_fourth_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1)+...
             sum(KKz2)+sum(qb1_2)+sum(qb2_2)+...
             sum(KKz3)+sum(qb1_3)+sum(qb2_3)+...
             sum(KKz4)+sum(qb1_4)+sum(qb2_4);
%% py moments outer while
forces_py1 = [KKx1;KKy1;0*KKz1]';
forces_py2 = [KKx2;KKy2;0*KKz2]';
forces_py3 = [KKx3;KKy3;0*KKz3]';
forces_py4 = [KKx4;KKy4;0*KKz4]';

positions1 = first_s_end_fixed';
positions2 = second_s_end_fixed';
positions3 = third_s_end_fixed';
positions4 = fourth_s_end_fixed';

moments_py1 = cross(positions1, forces_py1);
moments_py2 = cross(positions2, forces_py2);
moments_py3 = cross(positions3, forces_py3);
moments_py4 = cross(positions4, forces_py4);

total_moment_py1 = sum(moments_py1);
total_moment_py2 = sum(moments_py2);
total_moment_py3 = sum(moments_py3);
total_moment_py4 = sum(moments_py4);

total_moment_py_vec = total_moment_py1+total_moment_py2+total_moment_py3+total_moment_py4;
total_moment_py = sqrt(total_moment_py_vec(1)^2+total_moment_py_vec(2)^2+total_moment_py_vec(3)^2);

for i=1:n_l
    pc(i)= n*i;
end

%% tz moments outer while
forces_tz_out1 = [0*KKx1;0*KKy1;tz_out1]';
forces_tz_out2 = [0*KKx2;0*KKy2;tz_out2]';
forces_tz_out3 = [0*KKx3;0*KKy3;tz_out3]';
forces_tz_out4 = [0*KKx4;0*KKy4;tz_out4]';

forces_tz_in1 = [0*KKx1;0*KKy1;tz_in1]';
forces_tz_in2 = [0*KKx2;0*KKy2;tz_in2]';
forces_tz_in3 = [0*KKx3;0*KKy3;tz_in3]';
forces_tz_in4 = [0*KKx4;0*KKy4;tz_in4]';

positions_out1 = first_s_end_fixed';
positions_out2 = second_s_end_fixed';
positions_out3 = third_s_end_fixed';
positions_out4 = fourth_s_end_fixed';

positions_in1 = first_s_end_fixed_in';
positions_in2 = second_s_end_fixed_in';
positions_in3 = third_s_end_fixed_in';
positions_in4 = fourth_s_end_fixed_in';

moments_tz_out1 = cross(positions_out1, forces_tz_out1);
moments_tz_out2 = cross(positions_out2, forces_tz_out2);
moments_tz_out3 = cross(positions_out3, forces_tz_out3);
moments_tz_out4 = cross(positions_out4, forces_tz_out4);

moments_tz_in1 = cross(positions_in1, forces_tz_in1);
moments_tz_in2= cross(positions_in2, forces_tz_in2);
moments_tz_in3 = cross(positions_in3, forces_tz_in3);
moments_tz_in4 = cross(positions_in4, forces_tz_in4);

total_moment_tz1 = sum(moments_tz_in1)+sum(moments_tz_out1);
total_moment_tz2 = sum(moments_tz_in2)+sum(moments_tz_out2);
total_moment_tz3 = sum(moments_tz_in3)+sum(moments_tz_out3);
total_moment_tz4 = sum(moments_tz_in4)+sum(moments_tz_out4);

total_moment_tz_vec = total_moment_tz1+total_moment_tz2+total_moment_tz3+total_moment_tz4;
total_moment_tz = sqrt(total_moment_tz_vec(1)^2+total_moment_tz_vec(2)^2+total_moment_tz_vec(3)^2);
%% qz lid moments outer while

positions_pile1 = Pile_first_s_end_fixed';
positions_pile2 = Pile_second_s_end_fixed';
positions_pile3 = Pile_third_s_end_fixed';
positions_pile4 = Pile_fourth_s_end_fixed';

XX = zeros(size(qb1_1));
YY = zeros(size(qb1_1));

force_vectors1 = [XX;YY;qb1_1]';
force_vectors2 = [XX;YY;qb1_2]';
force_vectors3 = [XX;YY;qb1_3]';
force_vectors4 = [XX;YY;qb1_4]';

moments_qz_1 = cross(positions_pile1, force_vectors1);
moments_qz_2 = cross(positions_pile2, force_vectors2);
moments_qz_3 = cross(positions_pile3, force_vectors3);
moments_qz_4 = cross(positions_pile4, force_vectors4);

total_moment_qz_1 = sum(moments_qz_1);
total_moment_qz_2 = sum(moments_qz_2);
total_moment_qz_3 = sum(moments_qz_3);
total_moment_qz_4 = sum(moments_qz_4);

total_moment_qz_tip_vec = total_moment_qz_1+total_moment_qz_2+total_moment_qz_3+total_moment_qz_4;
total_moment_tip = sqrt(total_moment_qz_tip_vec(1)^2+total_moment_qz_tip_vec(2)^2+total_moment_qz_tip_vec(3)^2);

%% qz tip moments outer while
r_positions_pile1 = r_Pile_first_s_end_fixed';
r_positions_pile2 = r_Pile_second_s_end_fixed';
r_positions_pile3 = r_Pile_third_s_end_fixed';
r_positions_pile4 = r_Pile_fourth_s_end_fixed';

XXX = zeros(size(qb2_1));
YYY = zeros(size(qb2_1));

force_vectors_2_1 = [XXX;YYY;qb2_1]';
force_vectors_2_2 = [XXX;YYY;qb2_2]';
force_vectors_2_3 = [XXX;YYY;qb2_3]';
force_vectors_2_4 = [XXX;YYY;qb2_4]';

moments_qz_2_1 = cross(r_positions_pile1, force_vectors_2_1);
moments_qz_2_2 = cross(r_positions_pile2, force_vectors_2_2);
moments_qz_2_3 = cross(r_positions_pile3, force_vectors_2_3);
moments_qz_2_4 = cross(r_positions_pile4, force_vectors_2_4);

total_moment_qz_2_1 = sum(moments_qz_2_1);
total_moment_qz_2_2 = sum(moments_qz_2_2);
total_moment_qz_2_3 = sum(moments_qz_2_3);
total_moment_qz_2_4 = sum(moments_qz_2_4);

total_moment_qz_lid_vec = total_moment_qz_2_1+total_moment_qz_2_2+total_moment_qz_2_3+total_moment_qz_2_4;
total_moment_lid = sqrt(total_moment_qz_lid_vec(1)^2+total_moment_qz_lid_vec(2)^2+total_moment_qz_lid_vec(3)^2);
%% outer While condition calculation
MOMENT =total_moment_lid + total_moment_tip + total_moment_tz + total_moment_py;

F_Hx = SUMx;
F_Hy = SUMy;
F_HH  = sqrt (F_Hx^2+F_Hy^2);
ee = abs(MOMENT/F_HH);

CON2=abs((ee)-ee_tar);

end
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  End of outer While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%

j=0;
fi1=FI_peak1;
fi2=FI_peak2;
fi3=FI_peak3;
fi4=FI_peak4;
%% points  calculation final
[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Piles_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Piles_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end,r_Piles_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[fourth_s_end,fourth_s_end_in,fourth_s_end_fixed,fourth_s_end_fixed_in,Pile_fourth_s_end,Pile_fourth_s_end_fixed,r_Pile_fourth_s_end,r_Piles_fourth_s_end_fixed]=points(R_split,th4,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);

%% py and tz forces final
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[KKx2,KKy2,~] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[KKx3,KKy3,~] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
[KKx4,KKy4,~] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4); %KN

KKx1(mask1) = 0;
KKy1(mask1) = 0;
KKx2(mask2) = 0;
KKy2(mask2) = 0;
KKx3(mask3) = 0;
KKy3(mask3) = 0;
KKx4(mask4) = 0;
KKy4(mask4) = 0;

SUMx1= sum(KKx1);
SUMx2 = sum(KKx2);
SUMx3 = sum(KKx3);
SUMx4 = sum(KKx4);
SUMx = SUMx1+SUMx2+SUMx3+SUMx4;

SUMy1= sum(KKy1);
SUMy2 = sum(KKy2);
SUMy3 = sum(KKy3);
SUMy4 = sum(KKy4);
SUMy = SUMy1+SUMy2+SUMy3+SUMy4;

displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement4 = (fourth_s_end(3, :) - fourth_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in4 = (fourth_s_end_in(3, :) - fourth_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[~,~,P_i_out2] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[~,~,P_i_out3] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
[~,~,P_i_out4] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4); %KN

P_i_out1(mask1) = 0;
P_i_out2(mask2) = 0;
P_i_out3(mask3) = 0;
P_i_out4(mask4) = 0;

[~,tz_out1]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out2]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out3]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out4]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

[tz_in1,~]=tzcurve(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement_in1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,~]=tzcurve(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement_in2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,~]=tzcurve(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement_in3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in4,~]=tzcurve(P_i_out4, fi4, gamma, do_b, di_b, (fourth_s_end_fixed(3, :)), displacement_in4, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

tz_out1(mask1) = 0;
tz_out2(mask2) = 0;
tz_out3(mask3) = 0;
tz_out4(mask4) = 0;

tz_in1(mask1) = 0;
tz_in2(mask2) = 0;
tz_in3(mask3) = 0;
tz_in4(mask4) = 0;

KKz_in1 = tz_in1;
KKz_in2 = tz_in2;
KKz_in3 = tz_in3;
KKz_in4 = tz_in4;

KKz_out1 = tz_out1;
KKz_out2 = tz_out2;
KKz_out3 = tz_out3;
KKz_out4 = tz_out4;

KKz1 =KKz_in1+KKz_out1;
KKz2 =KKz_in2+KKz_out2;
KKz3 =KKz_in3+KKz_out3;
KKz4 =KKz_in4+KKz_out4;

%% qz forces final
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_4] = qzcurve2(fi_end4,Pile_fourth_s_end,Pile_fourth_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi_end1,Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi_end2,Z_lid,r_Pile_second_s_end, r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi_end3,Z_lid,r_Pile_third_s_end, r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_4] = qzcurve4(fi_end4,Z_lid,r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1)+...
             sum(KKz2)+sum(qb1_2)+sum(qb2_2)+...
             sum(KKz3)+sum(qb1_3)+sum(qb2_3)+...
             sum(KKz4)+sum(qb1_4)+sum(qb2_4);
%% py moments final
forces_py1 = [KKx1;KKy1;0*KKz1]';
forces_py2 = [KKx2;KKy2;0*KKz2]';
forces_py3 = [KKx3;KKy3;0*KKz3]';
forces_py4 = [KKx4;KKy4;0*KKz4]';

positions1 = first_s_end_fixed';
positions2 = second_s_end_fixed';
positions3 = third_s_end_fixed';
positions4 = fourth_s_end_fixed';

moments_py1 = cross(positions1, forces_py1);
moments_py2 = cross(positions2, forces_py2);
moments_py3 = cross(positions3, forces_py3);
moments_py4 = cross(positions4, forces_py4);

total_moment_py1 = sum(moments_py1);
total_moment_py2 = sum(moments_py2);
total_moment_py3 = sum(moments_py3);
total_moment_py4 = sum(moments_py4);

total_moment_py_vec = total_moment_py1+total_moment_py2+total_moment_py3+total_moment_py4;
total_moment_py = sqrt(total_moment_py_vec(1)^2+total_moment_py_vec(2)^2+total_moment_py_vec(3)^2);

for i=1:n_l
    pc(i)= n*i;
end

%% tz moments final
forces_tz_out1 = [0*KKx1;0*KKy1;tz_out1]';
forces_tz_out2 = [0*KKx2;0*KKy2;tz_out2]';
forces_tz_out3 = [0*KKx3;0*KKy3;tz_out3]';
forces_tz_out4 = [0*KKx4;0*KKy4;tz_out4]';

forces_tz_in1 = [0*KKx1;0*KKy1;tz_in1]';
forces_tz_in2 = [0*KKx2;0*KKy2;tz_in2]';
forces_tz_in3 = [0*KKx3;0*KKy3;tz_in3]';
forces_tz_in4 = [0*KKx4;0*KKy4;tz_in4]';

positions_out1 = first_s_end_fixed';
positions_out2 = second_s_end_fixed';
positions_out3 = third_s_end_fixed';
positions_out4 = fourth_s_end_fixed';

positions_in1 = first_s_end_fixed_in';
positions_in2 = second_s_end_fixed_in';
positions_in3 = third_s_end_fixed_in';
positions_in4 = fourth_s_end_fixed_in';

moments_tz_out1 = cross(positions_out1, forces_tz_out1);
moments_tz_out2 = cross(positions_out2, forces_tz_out2);
moments_tz_out3 = cross(positions_out3, forces_tz_out3);
moments_tz_out4 = cross(positions_out4, forces_tz_out4);

moments_tz_in1 = cross(positions_in1, forces_tz_in1);
moments_tz_in2= cross(positions_in2, forces_tz_in2);
moments_tz_in3 = cross(positions_in3, forces_tz_in3);
moments_tz_in4 = cross(positions_in4, forces_tz_in4);

total_moment_tz1 = sum(moments_tz_in1)+sum(moments_tz_out1);
total_moment_tz2 = sum(moments_tz_in2)+sum(moments_tz_out2);
total_moment_tz3 = sum(moments_tz_in3)+sum(moments_tz_out3);
total_moment_tz4 = sum(moments_tz_in4)+sum(moments_tz_out4);

total_moment_tz_vec = total_moment_tz1+total_moment_tz2+total_moment_tz3+total_moment_tz4;
total_moment_tz = sqrt(total_moment_tz_vec(1)^2+total_moment_tz_vec(2)^2+total_moment_tz_vec(3)^2);
%% qz lid moments final

positions_pile1 = Pile_first_s_end_fixed';
positions_pile2 = Pile_second_s_end_fixed';
positions_pile3 = Pile_third_s_end_fixed';
positions_pile4 = Pile_fourth_s_end_fixed';

XX = zeros(size(qb1_1));
YY = zeros(size(qb1_1));

force_vectors1 = [XX;YY;qb1_1]';
force_vectors2 = [XX;YY;qb1_2]';
force_vectors3 = [XX;YY;qb1_3]';
force_vectors4 = [XX;YY;qb1_4]';

moments_qz_1 = cross(positions_pile1, force_vectors1);
moments_qz_2 = cross(positions_pile2, force_vectors2);
moments_qz_3 = cross(positions_pile3, force_vectors3);
moments_qz_4 = cross(positions_pile4, force_vectors4);

total_moment_qz_1 = sum(moments_qz_1);
total_moment_qz_2 = sum(moments_qz_2);
total_moment_qz_3 = sum(moments_qz_3);
total_moment_qz_4 = sum(moments_qz_4);

total_moment_qz_tip_vec = total_moment_qz_1+total_moment_qz_2+total_moment_qz_3+total_moment_qz_4;
total_moment_qz_tip = sqrt(total_moment_qz_tip_vec(1)^2+total_moment_qz_tip_vec(2)^2+total_moment_qz_tip_vec(3)^2);
%% qz tip moments final
r_positions_pile1 = r_Pile_first_s_end_fixed';
r_positions_pile2 = r_Pile_second_s_end_fixed';
r_positions_pile3 = r_Pile_third_s_end_fixed';
r_positions_pile4 = r_Pile_fourth_s_end_fixed';

XXX = zeros(size(qb2_1));
YYY = zeros(size(qb2_1));

force_vectors_2_1 = [XXX;YYY;qb2_1]';
force_vectors_2_2 = [XXX;YYY;qb2_2]';
force_vectors_2_3 = [XXX;YYY;qb2_3]';
force_vectors_2_4 = [XXX;YYY;qb2_4]';

moments_qz_2_1 = cross(r_positions_pile1, force_vectors_2_1);
moments_qz_2_2 = cross(r_positions_pile2, force_vectors_2_2);
moments_qz_2_3 = cross(r_positions_pile3, force_vectors_2_3);
moments_qz_2_4 = cross(r_positions_pile4, force_vectors_2_4);

total_moment_qz_2_1 = sum(moments_qz_2_1);
total_moment_qz_2_2 = sum(moments_qz_2_2);
total_moment_qz_2_3 = sum(moments_qz_2_3);
total_moment_qz_2_4 = sum(moments_qz_2_4);

total_moment_qz_lid_vec = total_moment_qz_2_1+total_moment_qz_2_2+total_moment_qz_2_3+total_moment_qz_2_4;
total_moment_qz_lid = sqrt(total_moment_qz_lid_vec(1)^2+total_moment_qz_lid_vec(2)^2+total_moment_qz_lid_vec(3)^2);
%% Moment Calculation final
MOMENT =total_moment_qz_tip+ total_moment_qz_lid + total_moment_tz + total_moment_py;

F_Hx = sum(KKx1)+sum(KKx2)+sum(KKx3)+sum(KKx4);
F_Hy = sum(KKy1)+sum(KKy2)+sum(KKy3)+sum(KKy4);
F_HH  = sqrt (F_Hx^2+F_Hy^2);

%% Store results for polar plot
M_polar(i_dir, count, i_depth) = MOMENT;
F_H_polar(i_dir, count, i_depth) = F_HH;

fprintf('\n  >> dir=%d deg, dd=%.1f deg => MOMENT=%.2f, F_HH=%.2f\n', tilt_dir, dd*180/pi, MOMENT, F_HH);

end  % end of dd loop (rotation steps)

end  % end of i_dir loop (tilt directions)

end  % end of i_depth loop (depth values)

%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  POLAR PLOTS  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%
fprintf('\n\n===== Generating bi1/bi2 style polar plots =====\n');

% Close tilt_directions to make full circle for plotting
tilt_dirs_closed = [tilt_directions, 360];
tilt_dirs_rad_closed = tilt_dirs_closed * pi / 180;

% Colors for different depths (deepest=red to shallowest=purple)
colors_depth = [
    0.85, 0.13, 0.10;  % Red (z=6.5m)
    0.93, 0.59, 0.13;  % Orange (z=6m)
    0.17, 0.67, 0.19;  % Green (z=5m)
    0.00, 0.25, 0.74;  % Blue (z=4m)
    0.56, 0.18, 0.56;  % Purple (z=3.5m)
];

%% Figure 1: bi1-style - Polar subplots of M vs Tilt Direction (one per tilt angle, curves per depth)
figure('Position', [50 50 1400 900]);
sgtitle({'Moment Resistance M [MN-m] vs Tilt Direction', '(Polar Plot for Each Tilt Angle)'}, ...
    'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

for i_step = 1:num_steps
    subplot(2, 3, i_step);
    h_lines = [];
    for i_depth = num_depths:-1:1
        M_closed = [M_polar(:, i_step, i_depth)', M_polar(1, i_step, i_depth)] / 1000;
        h = polarplot(tilt_dirs_rad_closed, M_closed, '-o', 'LineWidth', 1.5, ...
            'Color', colors_depth(i_depth,:), 'MarkerSize', 4, 'MarkerFaceColor', colors_depth(i_depth,:));
        hold on;
        h_lines(end+1) = h;
    end
    title(sprintf('\\theta_{tilt} = %.2f\\circ', dd_steps(i_step)*180/pi), ...
        'FontSize', 13, 'FontWeight', 'bold');
    ax = gca;
    ax.ThetaDir = 'counterclockwise';
    ax.ThetaZeroLocation = 'right';
    set(ax, 'FontSize', 9, 'FontWeight', 'bold');
end

% Add shared legend in last subplot area
legend_labels_depth = cell(1, num_depths);
for i_depth = num_depths:-1:1
    legend_labels_depth{num_depths - i_depth + 1} = sprintf('z = %.1f m', depth_values(i_depth));
end
legend(h_lines, legend_labels_depth, 'FontSize', 11, 'Location', 'bestoutside', 'Interpreter', 'tex');

%% Figure 2: bi2-style - Cylindrical Surface M(theta, z)
figure('Position', [50 50 1500 500]);
sgtitle({'Cylindrical Surface: M(\theta, z)', 'Radius = Moment Resistance, Height = Penetration Depth'}, ...
    'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

bi2_steps = [2, 4, 6];  % indices for 0.5, 1.0, 1.5 degrees

% Fine mesh for smooth surface
theta_fine = linspace(0, 2*pi, 73);
depth_fine = linspace(min(depth_values), max(depth_values), 30);
[THETA_c, DEPTH_c] = meshgrid(tilt_dirs_rad_closed, depth_values);
[THETA_f, DEPTH_f] = meshgrid(theta_fine, depth_fine);

for idx = 1:3
    i_step = bi2_steps(idx);
    subplot(1, 3, idx);

    % Build coarse surface data (num_depths x n_theta_closed)
    M_coarse = zeros(num_depths, length(tilt_dirs_rad_closed));
    for i_depth = 1:num_depths
        M_coarse(i_depth, 1:num_dirs) = M_polar(:, i_step, i_depth)' / 1000;
        M_coarse(i_depth, end) = M_polar(1, i_step, i_depth) / 1000;
    end

    % Interpolate to fine mesh
    M_fine = interp2(THETA_c, DEPTH_c, M_coarse, THETA_f, DEPTH_f, 'spline');

    X_f = M_fine .* cos(THETA_f);
    Y_f = M_fine .* sin(THETA_f);
    Z_f = DEPTH_f;

    surf(X_f, Y_f, Z_f, M_fine, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
    hold on;

    % Plot data points as red dots
    for i_depth = 1:num_depths
        for i_dir = 1:num_dirs
            M_val = M_polar(i_dir, i_step, i_depth) / 1000;
            th_val = tilt_directions(i_dir) * pi / 180;
            plot3(M_val*cos(th_val), M_val*sin(th_val), depth_values(i_depth), ...
                'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
        end
    end

    colormap(gca, parula);
    cb = colorbar;
    ylabel(cb, 'M [MN-m]', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('M\cdotcos(\theta) [MN-m]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('M\cdotsin(\theta) [MN-m]', 'FontSize', 11, 'FontWeight', 'bold');
    zlabel('Depth [m]', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Tilt Angle = %.1f\\circ', dd_steps(i_step)*180/pi), ...
        'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'ZDir', 'reverse');
    grid on;
    view([-37.5, 30]);
    hold off;
end

%% Figure 3: Polar Moment-Rotation Diagram (all steps, deepest depth only)
figure('Position', [50 50 800 700]);
legend_labels = {};
colors_polar = lines(num_steps);
for i_step = 1:num_steps
    M_closed = [M_polar(:, i_step, end)', M_polar(1, i_step, end)] / 1000;
    polarplot(tilt_dirs_rad_closed, M_closed, '-o', 'LineWidth', 2, ...
        'Color', colors_polar(i_step,:), 'MarkerSize', 5);
    hold on;
    legend_labels{i_step} = sprintf('\\theta = %.2f\\circ', dd_steps(i_step)*180/pi);
end
title('Polar Moment-Rotation Diagram (Four Pod, z = 6.5 m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(legend_labels, 'fontsize', 12, 'interpreter', 'tex', 'location', 'bestoutside');
set(gca, 'fontweight', 'bold', 'fontsize', 14);
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'right';

%% Figure 4: Moment vs Rotation for all directions (Cartesian, deepest depth)
figure('Position', [50 50 900 600]);
colors_dir = lines(num_dirs);
for i_dir = 1:num_dirs
    plot(dd_steps*180/pi, squeeze(M_polar(i_dir, 1:num_steps, end))/1000, '-o', 'LineWidth', 2, ...
        'Color', colors_dir(i_dir,:), 'MarkerSize', 5);
    hold on;
end
grid on;
set(gca,'fontweight','bold','fontsize',16);
xlabel('Rotation (degree)', 'FontSize', 20, 'FontWeight', 'bold','FontName','Times New Roman');
ylabel('Moment (MN-m)', 'FontSize', 20, 'FontWeight', 'bold','FontName','Times New Roman');
title('Moment vs Rotation for All Tilt Directions (z = 6.5 m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend_labels4 = {};
for i_dir = 1:num_dirs
    legend_labels4{i_dir} = sprintf('Dir = %d\\circ', tilt_directions(i_dir));
end
legend(legend_labels4, 'fontsize', 10, 'interpreter', 'tex', 'location', 'bestoutside');

fprintf('\n===== DONE =====\n');
