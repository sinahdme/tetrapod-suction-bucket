clc;
clear all;
close all;
%% %%%%%%%%%% Geometric parameters %%%%%%%%%%%%
%% VALIDATION: Liu et al. (2024) tetrapod - L/D=1.0, S_cc/D=1.0
%% Paper: "Tilt adjustment strategy and applicability of four-bucket
%%         jacket foundation for offshore wind turbines"
%%         Ocean Engineering 302, 117642
%% Scale: 1:100 (1g model test) - ALL VALUES SCALED TO PROTOTYPE (x100)
%% Soil: Fine sand, Dr_compact=0.69, phi=34.46 deg
%% WARNING: 1g test - stress levels NOT corrected for prototype
do_b = 15.0; %bucket outer diameter [m] - Liu2024: 150mm x100
thickness = 0.100; %bucket wall thickness [m] - Liu2024: 1mm x100
l_b = 15.0; %length of the bucket skirt [m] - Liu2024: 150mm x100, L/D=1.0
Z_lid =0.0001; % thickness of the lid of bucket
n = 32; %number of strips in a single circumference
n_l =32; % number of rings in the bottom
n_q =10;
eccentricity = 76; % [m] - estimated (~5D, proportional to Barari)
th1= 45;
th2 = 135;
th3 = 225;
th4 = 315;
Deadload=8000; %kN - estimated (non-dim Deadload/(gamma'*D^3)~2.0, same as Barari)
l = 10.607; % [m] - Liu2024: R_jacket = S_cc*sqrt(2)/2 = 15.0*0.7071
            % S_cc=15.0m (adjacent bucket spacing, 150mm x100), S_cc/D=1.0
%% Group interaction p-multiplier (parametric calibration)
SD_ratio = l / do_b;  % spacing-to-diameter ratio (center-to-bucket / D = 0.707)
delta_ref = 0.05;  % reference displacement (m) for full f_g reduction
bucket_angles = [th1, th2, th3, th4];  % [45, 135, 225, 315]
%% %%%%%%%%%% Soil parameters %%%%%%%%%%%%
%% Fine sand (Table 2 in Liu et al. 2024)
gamma = 11.77; % buoyant unit weight [kN/m3] - Liu2024: gamma_sat=21.58, gamma_w=9.81
S_gamma = 0.5;
e_max =0;
e_min =0;
D_r = 0.69; % relative compaction - Liu2024: 0.69
m_c = 161.42*D_r^2+199.8*D_r+36.877; % 80-150 soft sand % 150-250 medium sand % 250-400 stiff sand  (rewrite based on dr)
A_c = 0.3474*D_r^2+0.4222*D_r+0.328;
fi_crit = 30; % critical state friction angle [deg] - estimated from Bolton: phi_peak(34.46)-psi~30
m =0; % for finding the FI_preak %m is either 0 or 3
delta = 23.0; % interface friction [deg] - estimated: ~2/3 * phi_peak(34.46)
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
%% Code initiation
QQ = zeros(1,10000);
kkk =0;
abc=0;
cc1=0;
% full penetration action
    abc=abc+1;
 cc1=cc1+1;
%    if abs(l_b - penetration)<=0.01
%        FP =1;
%    else
%         FP=0;d
%    end
count =0;

%%  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ dealload effect with no tilt angle $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  %%
dz=0;
SUMz=0;
while true
    if SUMz>Deadload
    dz=dz-0.00001;
    else
     dz=dz+0.00001;
    end
for dd=0
   q=[-dd -dd];
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
(Deadload - SUMz)
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
%% Pre-compute bucket adjacency (passive wedge overlap)
passive_extent = do_b/2 + l_b*tand(60);
overlap_threshold = 2 * passive_extent;
n_buckets = length(bucket_angles);
bucket_centers = l * [cosd(bucket_angles); sind(bucket_angles)];
adj_matrix = false(n_buckets);
for ii = 1:n_buckets
    for jj = ii+1:n_buckets
        d_ij = norm(bucket_centers(:,ii) - bucket_centers(:,jj));
        adj_matrix(ii,jj) = d_ij < overlap_threshold;
        adj_matrix(jj,ii) = adj_matrix(ii,jj);
    end
end
%%  $$$$$$$$$$$$$$$$$$$$$$$$$   bucket rotation angle change step by step $$$$$$$$$$$$$$$$$$$$$$$$$ %%
tilt_angles = [0,45];  % tilt directions to simulate (degrees)
tilt_colors = {'b', 'r'};  % colors for plotting: 0deg=blue, 45deg=red
for tilt_idx = 1:length(tilt_angles)
    tilt_angle = tilt_angles(tilt_idx);
    fprintf('\n--- Running tilt direction: %d deg ---\n', tilt_angle);
    count = 0;
Q = [0,0,0]';     %initial acenter of rotation
for dd=0.25*pi/180:0.25*pi/180:3.25*pi/180
R_tilt = rodrigues_rotation(dd, 180 - tilt_angle);  % Liu2024 convention: tilt_angle=0 -> between pods (+x)
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
cc=0;
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
 mmm=mmm-0.015;%*abs(((ee)-ee_tar));
 else
 mmm=mmm+0.015;%*abs(((ee)-ee_tar));
 end
if cc==1
 Q = [nnn*cosd(angle_indicator), nnn*sind(angle_indicator), -(1.0*l_b)+mmm]';
else
 Q = [(nnn)*cosd(angle_indicator), nnn*sind(angle_indicator), -(1.0*l_b)+mmm]';
end
if rem(ccc,5)==0
fprintf('\rdegree: %0.3f  | ((ee)-ee_tar): %0.2f |(SUMz-Deadload): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', dd*180/pi, CON2,CON1,Q(1),Q(2),Q(3));
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
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %%
 while abs(SUMz-Deadload)>(Deadload/800)
CON1 = abs(SUMz-Deadload); %track the condition
tower_tip = [0,0,eccentricity]';
tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
cc=cc+1;
if SUMz-Deadload>0
nnn=nnn+0.015;
else
nnn=nnn-0.015;
end
Q = [(0+nnn)*cosd(angle_indicator), (0+nnn)*sind(angle_indicator), Q(3)]';% iterative increasing or decreasing the x-axis of the rotation center
fprintf('\rdegree: %0.3f  | ((ee)-ee_tar): %0.3f |(Deadload - SUMz): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', dd*180/pi,((ee)-ee_tar) ,CON1,Q(1),Q(2),Q(3));

QQ(cc+2)=Q(1);
if abs(abs(QQ(cc+2))-abs(QQ(cc)))==0
     break
 end

[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Pile_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Pile_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end, r_Pile_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[fourth_s_end,fourth_s_end_in,fourth_s_end_fixed,fourth_s_end_fixed_in,Pile_fourth_s_end,Pile_fourth_s_end_fixed,r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed]=points(R_split,th4,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
%% relative points calculation outer while
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

KKZ1(count,:) = tz_out1+tz_in1;
KKZ2(count,:) = tz_out2+tz_in2;
KKZ3(count,:) = tz_out3+tz_in3;
KKZ4(count,:) = tz_out4+tz_in4;
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
%% relative points calculation outer while
% [relative_first_s_end,relative_first_s_end_fixed,relative_first_s_end_in,relative_first_s_end_fixed_in]...
%     =relative_points(first_s_end, first_s_end_in, first_s_end_fixed, first_s_end_fixed_in, Pile_first_s_end, Pile_first_s_end_fixed, r_Pile_first_s_end, r_Pile_first_s_end_fixed, th1,l,dz,Q,R_tilt);
%
% [relative_second_s_end,relative_second_s_end_fixed,relative_second_s_end_in,relative_second_s_end_fixed_in]...
%     =relative_points(second_s_end, second_s_end_in, second_s_end_fixed, second_s_end_fixed_in, Pile_second_s_end, Pile_second_s_end_fixed, r_Pile_second_s_end, r_Pile_second_s_end_fixed, th2,l,dz,Q,R_tilt);
%
% [relative_third_s_end,relative_third_s_end_fixed,relative_third_s_end_in,relative_third_s_end_fixed_in]...
%     =relative_points(third_s_end, third_s_end_in, third_s_end_fixed, third_s_end_fixed_in, Pile_third_s_end, Pile_third_s_end_fixed, r_Pile_third_s_end, r_Pile_third_s_end_fixed, th3,l,dz,Q,R_tilt);
%
% [relative_fourth_s_end,relative_fourth_s_end_fixed,relative_fourth_s_end_in,relative_fourth_s_end_fixed_in]...
%     =relative_points(fourth_s_end, fourth_s_end_in, fourth_s_end_fixed, fourth_s_end_fixed_in, Pile_fourth_s_end, Pile_fourth_s_end_fixed, r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed, th4,l,dz,Q,R_tilt);

%% py and tz forces outer while
% Detect per-bucket f_g (compression + adjacency)
fg_vec = detect_fg_per_bucket( ...
    {first_s_end,second_s_end,third_s_end,fourth_s_end}, ...
    {first_s_end_fixed,second_s_end_fixed,third_s_end_fixed,fourth_s_end_fixed}, ...
    adj_matrix, SD_ratio, delta_ref, l_b/do_b, 'tetra');
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1,fg_vec(1)); %KN
[KKx2,KKy2,~] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2,fg_vec(2)); %KN
[KKx3,KKy3,~] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3,fg_vec(3)); %KN
[KKx4,KKy4,~] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4,fg_vec(4)); %KN

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

[~,~,P_i_out1] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1,fg_vec(1)); %KN
[~,~,P_i_out2] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2,fg_vec(2)); %KN
[~,~,P_i_out3] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3,fg_vec(3)); %KN
[~,~,P_i_out4] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4,fg_vec(4)); %KN

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

KKZ1(count,:) = tz_out1+tz_in1;
KKZ2(count,:) = tz_out2+tz_in2;
KKZ3(count,:) = tz_out3+tz_in3;
KKZ4(count,:) = tz_out4+tz_in4;
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
%   Q=[0,0,0]';
fi1=FI_peak1;
fi2=FI_peak2;
fi3=FI_peak3;
fi4=FI_peak4;
%% points  calculation outer while
[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Piles_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Piles_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end,r_Piles_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[fourth_s_end,fourth_s_end_in,fourth_s_end_fixed,fourth_s_end_fixed_in,Pile_fourth_s_end,Pile_fourth_s_end_fixed,r_Pile_fourth_s_end,r_Piles_fourth_s_end_fixed]=points(R_split,th4,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
%% relative points calculation final
% [relative_first_s_end,relative_first_s_end_fixed,relative_first_s_end_in,relative_first_s_end_fixed_in]...
%     =relative_points(first_s_end, first_s_end_in, first_s_end_fixed, first_s_end_fixed_in, Pile_first_s_end, Pile_first_s_end_fixed, r_Pile_first_s_end, r_Pile_first_s_end_fixed, th1,l,dz,Q,R_tilt);
%
% [relative_second_s_end,relative_second_s_end_fixed,relative_second_s_end_in,relative_second_s_end_fixed_in]...
%     =relative_points(second_s_end, second_s_end_in, second_s_end_fixed, second_s_end_fixed_in, Pile_second_s_end, Pile_second_s_end_fixed, r_Pile_second_s_end, r_Pile_second_s_end_fixed, th2,l,dz,Q,R_tilt);
%
% [relative_third_s_end,relative_third_s_end_fixed,relative_third_s_end_in,relative_third_s_end_fixed_in]...
%     =relative_points(third_s_end, third_s_end_in, third_s_end_fixed, third_s_end_fixed_in, Pile_third_s_end, Pile_third_s_end_fixed, r_Pile_third_s_end, r_Pile_third_s_end_fixed, th3,l,dz,Q,R_tilt);
%
% [relative_fourth_s_end,relative_fourth_s_end_fixed,relative_fourth_s_end_in,relative_fourth_s_end_fixed_in]...
%     =relative_points(fourth_s_end, fourth_s_end_in, fourth_s_end_fixed, fourth_s_end_fixed_in, Pile_fourth_s_end, Pile_fourth_s_end_fixed, r_Pile_fourth_s_end, r_Pile_fourth_s_end_fixed, th4,l,dz,Q,R_tilt);

%% py and tz forces final
% Detect per-bucket f_g (compression + adjacency)
fg_vec = detect_fg_per_bucket( ...
    {first_s_end,second_s_end,third_s_end,fourth_s_end}, ...
    {first_s_end_fixed,second_s_end_fixed,third_s_end_fixed,fourth_s_end_fixed}, ...
    adj_matrix, SD_ratio, delta_ref, l_b/do_b, 'tetra');
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1,fg_vec(1)); %KN
[KKx2,KKy2,~] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2,fg_vec(2)); %KN
[KKx3,KKy3,~] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3,fg_vec(3)); %KN
[KKx4,KKy4,~] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4,fg_vec(4)); %KN

KKx1(mask1) = 0;
KKy1(mask1) = 0;
KKx2(mask2) = 0;
KKy2(mask2) = 0;
KKx3(mask3) = 0;
KKy3(mask3) = 0;
KKx4(mask4) = 0;
KKy4(mask4) = 0;

% KKX1(count,:) = KKx1;
% KKY1(count,:) = KKy1;
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

Reshaped_first_s_end = reshape(first_s_end(3, :),[n,n_l]);
Reshaped_second_s_end = reshape(second_s_end(3, :),[n,n_l]);
Reshaped_third_s_end = reshape(third_s_end(3, :),[n,n_l]);
Reshaped_fourth_s_end = reshape(fourth_s_end(3, :),[n,n_l]);

Reshaped_first_s_end_fixed = reshape(first_s_end_fixed(3, :) ,[n,n_l]);
Reshaped_second_s_end_fixed  = reshape(second_s_end_fixed(3, :) ,[n,n_l]);
Reshaped_third_s_end_fixed  = reshape(third_s_end_fixed(3, :) ,[n,n_l]);
Reshaped_fourth_s_end_fixed  = reshape(fourth_s_end_fixed(3, :) ,[n,n_l]);

Reshaped_displacement1 = reshape(displacement1,[n,n_l]);
Reshaped_displacement2 = reshape(displacement2,[n,n_l]);
Reshaped_displacement3 = reshape(displacement3,[n,n_l]);
Reshaped_displacement4 = reshape(displacement4,[n,n_l]);

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in4 = (fourth_s_end_in(3, :) - fourth_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1,fg_vec(1)); %KN
[~,~,P_i_out2] = coeffinder(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2,fg_vec(2)); %KN
[~,~,P_i_out3] = coeffinder(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3,fg_vec(3)); %KN
[~,~,P_i_out4] = coeffinder(fourth_s_end,fourth_s_end_fixed,abs(fourth_s_end_fixed(3, :)),l_b,n_l,fi4,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th4,fg_vec(4)); %KN


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

Reshaped_tz_out1 = reshape(tz_out1,[n,n_l]);
Reshaped_tz_out2 = reshape(tz_out2,[n,n_l]);
Reshaped_tz_out3 = reshape(tz_out3,[n,n_l]);
Reshaped_tz_out4 = reshape(tz_out4,[n,n_l]);
Reshaped_tz_in1 = reshape(tz_in1,[n,n_l]);
Reshaped_tz_in2 = reshape(tz_in2,[n,n_l]);
Reshaped_tz_in3 = reshape(tz_in3,[n,n_l]);
Reshaped_tz_in4 = reshape(tz_in4,[n,n_l]);

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

KKZ1(count,:) = tz_out1+tz_in1;
KKZ2(count,:) = tz_out2+tz_in2;
KKZ3(count,:) = tz_out3+tz_in3;
KKZ4(count,:) = tz_out4+tz_in4;
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

MOMENTT(count) =total_moment_qz_tip+ total_moment_qz_lid + total_moment_tz + total_moment_py;

F_Hx = sum(KKx1)+sum(KKx2)+sum(KKx3)+sum(KKx4);
F_Hy = sum(KKy1)+sum(KKy2)+sum(KKy3)+sum(KKy4);
F_HH  = sqrt (F_Hx^2+F_Hy^2);
M_1(count)  = F_HH  *eccentricity;
M_2(count)  = MOMENT;

colors2 = [
    0.111, 0.111, 0.111;    % Dark Grey
    0.8500, 0.3250, 0.0980; % Orange
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330; % Light Blue
    0.9290, 0.6940, 0.1250; % Yellow
    0.6350, 0.0780, 0.1840  % Dark Red
];

figure(4)
plot(dd*180/pi,M_2(count)/1000,['*' tilt_colors{tilt_idx}],'LineWidth',5)

grid on
hold on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',22)
xlabel('Rotation (degree)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Moment Load (MNm)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
figure(107)
scatter3(Q(1),Q(2),Q(3),'o','LineWidth',8)
xlabel('X', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Y', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Z', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')

hold on
hold on
grid on;
xlabel('Angle (degrees)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
ylabel('Skirt wall friction (F_t_z/\gamma^\prime/D^3) ','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
legend({'Inner right side friction','Inner left side friction','Outer left side friction','Outer right side friction'},'fontsize',16,'interpreter','latex','location','northeast')

 CENTE_OF_ROTATION(:,count) = Q;

end  % dd loop

% Store results for this tilt direction
ddd_all{tilt_idx} = ddd(1:count);
M_2_all{tilt_idx} = M_2(1:count);
M_1_all{tilt_idx} = M_1(1:count);

%% Print results for this direction
fprintf('\n\n====== LIU2024 VALIDATION RESULTS (%d deg tilt) ======\n', tilt_angle);
fprintf('Rotation(deg) | M_tower(MNm) | M_soil(MNm)\n');
for ii=1:count
    fprintf('%8.3f      | %12.3f | %12.3f\n', ddd(ii), M_1(ii)/1000, M_2(ii)/1000);
end
fprintf('=======================================\n\n');

end  % tilt_idx loop

%% Comparison plot: 0deg vs 45deg
figure(200)
hold on; grid on;
for tilt_idx = 1:length(tilt_angles)
    plot(ddd_all{tilt_idx}, M_2_all{tilt_idx}/1000, ['-o' tilt_colors{tilt_idx}], ...
        'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('%d deg tilt', tilt_angles(tilt_idx)));
end
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',22)
xlabel('Rotation (degree)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Moment (MNm)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
legend('show', 'FontSize', 16, 'Location', 'northwest')
title('Liu2024 Tetrapod (prototype scale x100): 0 vs 45 deg tilt', 'FontSize', 20)

%% Print directional gap
fprintf('\n====== DIRECTIONAL GAP ANALYSIS ======\n');
fprintf('Rotation(deg) | M_0deg(MNm) | M_45deg(MNm) | Gap(%%)\n');
n_common = min(length(ddd_all{1}), length(ddd_all{2}));
for ii = 1:n_common
    gap_pct = (M_2_all{2}(ii) - M_2_all{1}(ii)) / M_2_all{2}(ii) * 100;
    fprintf('%8.3f      | %11.3f | %12.3f | %6.1f\n', ...
        ddd_all{1}(ii), M_2_all{1}(ii)/1000, M_2_all{2}(ii)/1000, gap_pct);
end
fprintf('==========================================\n\n');

M_2_f = M_2_all{1}'/1000;  % use 0deg for backward compatibility

figure(15)
 hold on
  scatter3(Q(1),Q(2),Q(3),'*green', 'LineWidth', 10)
  hold on

hold on;  % Hold the current plot

% Number of points in each segment
segment_size = n;

% Loop through the points in segments
for i = 1:segment_size:size(first_s_end, 2)
    % Get the indices for the current segment
    indices = i:min(i+segment_size-1, size(first_s_end, 2));

    % Ensure that we have a full segment of points
    if numel(indices) == segment_size
        % Get the points for the current segment
        segment_points1 = first_s_end(:, indices);
        segment_points2 = second_s_end(:, indices);
        segment_points3 = third_s_end(:, indices);
        segment_points4 = fourth_s_end(:, indices);

        % Add the first point to the end to close the circle
        segment_points1 = [segment_points1, segment_points1(:, 1)];
        segment_points2 = [segment_points2, segment_points2(:, 1)];
        segment_points3 = [segment_points3, segment_points3(:, 1)];
        segment_points4 = [segment_points4, segment_points4(:, 1)];
        % Plot the closed circle for the current segment
        plot3(segment_points1(1,:), segment_points1(2,:), segment_points1(3,:),'-black', 'LineWidth', 1);
        plot3(segment_points2(1,:), segment_points2(2,:), segment_points2(3,:),'-black', 'LineWidth', 1);
        plot3(segment_points3(1,:), segment_points3(2,:), segment_points3(3,:),'-black', 'LineWidth', 1);
        plot3(segment_points4(1,:), segment_points4(2,:), segment_points4(3,:),'-black', 'LineWidth', 1);

    end
end
% Adding labels to each point
num_points = size(first_s_end, 2); % Number of points in the matrix
offset = [0.2; 0.2; 0.2];  % Define an offset vector [x_offset; y_offset; z_offset]
FontSize = 12;  % Define the desired font size
%% put the number outside the bucket

% Define the origin of the cylinder
cylinder_origin = mean(first_s_end, 2);  % Assuming the cylinder is centered around the mean of the points

% Define the offset scale (how far to move the numbers outside the cylinder)
offset_scale = 0.38;  % Adjust this value as needed

% Define the font size
FontSize = 12;
FontColor = 'blue';
FontWeight = 'bold';
FontName = 'Times New Roman';

%%
 X =sqrt(first_s_end_fixed(1,:).^2 + first_s_end_fixed(2,:).^2)
indices_p = find(X(:) >= do_b/2)'; % nodes that pressurize the soil from outside
indices_Not_p = find(X(:) < do_b/2)'; % nodes that NOT pressurize the soil from outside
hold on
 scatter3(first_s_end_fixed(1,find(KKx1<0)),first_s_end_fixed(2,find(KKx1<0)),first_s_end_fixed(3,find(KKx1<0)),'*g', 'LineWidth', 5)
 hold on
 scatter3(first_s_end_fixed(1,find(KKx1==0)),first_s_end_fixed(2,find(KKx1==0)),first_s_end_fixed(3,find(KKx1==0)),'*k', 'LineWidth', 5)
 hold on
scatter3(first_s_end_fixed(1,find(KKx1>0)),first_s_end_fixed(2,find(KKx1>0)),first_s_end_fixed(3,find(KKx1>0)),'*r', 'LineWidth', 5)
hold on


%%
 scatter3(second_s_end_fixed(1,find(KKx2<0)),second_s_end_fixed(2,find(KKx2<0)),second_s_end_fixed(3,find(KKx2<0)),'*g', 'LineWidth', 5)
 hold on
 scatter3(second_s_end_fixed(1,find(KKx2==0)),second_s_end_fixed(2,find(KKx2==0)),second_s_end_fixed(3,find(KKx2==0)),'*k', 'LineWidth', 5)
 hold on
scatter3(second_s_end_fixed(1,find(KKx2>0)),second_s_end_fixed(2,find(KKx2>0)),second_s_end_fixed(3,find(KKx2>0)),'*r', 'LineWidth', 5)
hold on
%%
 scatter3(third_s_end_fixed(1,find(KKx3<0)),third_s_end_fixed(2,find(KKx3<0)),third_s_end_fixed(3,find(KKx3<0)),'*g', 'LineWidth', 5)
 hold on
 scatter3(third_s_end_fixed(1,find(KKx3==0)),third_s_end_fixed(2,find(KKx3==0)),third_s_end_fixed(3,find(KKx3==0)),'*k', 'LineWidth', 5)
 hold on
scatter3(third_s_end_fixed(1,find(KKx3>0)),third_s_end_fixed(2,find(KKx3>0)),third_s_end_fixed(3,find(KKx3>0)),'*r', 'LineWidth', 5)
hold on
%%
 scatter3(fourth_s_end_fixed(1,find(KKx4<0)),fourth_s_end_fixed(2,find(KKx4<0)),fourth_s_end_fixed(3,find(KKx4<0)),'*g', 'LineWidth', 5)
 hold on
 scatter3(fourth_s_end_fixed(1,find(KKx4==0)),fourth_s_end_fixed(2,find(KKx4==0)),fourth_s_end_fixed(3,find(KKx4==0)),'*k', 'LineWidth', 5)
 hold on
scatter3(fourth_s_end_fixed(1,find(KKx4>0)),fourth_s_end_fixed(2,find(KKx4>0)),fourth_s_end_fixed(3,find(KKx4>0)),'*r', 'LineWidth', 5)
hold on
%% make a local frame
% Define the frame size
frame_size = 15; % Scaled for prototype

% Define the origin
origin = [0, 0, 0];

% Define the endpoints of the axes
x_axis = [frame_size, 0, 0];
y_axis = [0, frame_size, 0];
z_axis = [0, 0, frame_size];


hold on;
grid on;
axis equal;

% Plot the x-axis
quiver3(origin(1), origin(2), origin(3), x_axis(1), x_axis(2), x_axis(3), ...
    'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Plot the y-axis
quiver3(origin(1), origin(2), origin(3), y_axis(1), y_axis(2), y_axis(3), ...
    'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Plot the z-axis
quiver3(origin(1), origin(2), origin(3), z_axis(1), z_axis(2), z_axis(3), ...
    'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% Add labels to the axes
text(x_axis(1) + offset(1), x_axis(2) + offset(2), x_axis(3) + offset(3), 'X', 'FontSize', 16, 'Color', 'r');
text(y_axis(1) + offset(1), y_axis(2) + offset(2), y_axis(3) + offset(3), 'Y', 'FontSize', 16, 'Color', 'r');
text(z_axis(1) + offset(1), z_axis(2) + offset(2), z_axis(3) + offset(3), 'Z', 'FontSize', 16, 'Color', 'r');

% Set the view and labels
set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
hold off;
%%


figure(16);
hold on;
grid on;
axis equal;
X = r_Pile_first_s_end';
V =qb2_1 ;
hold on
plot3(X(:,1),X(:,2),X(:,3), '*r')
hold on
axis([-20, 20, -20, 20, 0, 20]); % scaled for prototype
grid
X = r_Pile_first_s_end';
t = delaunay(X(:,1),X(:,2));
defaultFaceColor  = [0.6875 0.8750 0.8984];
trisurf(t,X(:,1),X(:,2),V, 'FaceColor', ...
    defaultFaceColor, 'FaceAlpha',0.9)
hold off
view(122.5, 30);

figure(17);
scatter3(r_Pile_first_s_end(1,:), r_Pile_first_s_end(2,:), r_Pile_first_s_end(3,:), 'ko', 'LineWidth', 1);
hold on
scatter3(r_Pile_first_s_end_fixed(1,:), r_Pile_first_s_end_fixed(2,:), r_Pile_first_s_end_fixed(3,:), 'bo', 'LineWidth', 1);
hold on
