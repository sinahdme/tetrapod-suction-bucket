function [qb2] = qzcurve3(fi_end, r_tip_end, r_tip_end_fixed, gamma, n, nq, R_split,Do)
    % qzcurve2 calculates the resisting moment components qb1 and qb2 for a suction bucket foundation.
properties_Nq = [15, 8; 20, 12; 25, 20; 30, 40; 35, 50];

input_Nq = properties_Nq(:, 1);
output_Nq = properties_Nq(:, 2);
% Generate a set of input values (xi) for interpolation
% This spans from the minimum to the maximum input values in your dataset

% Perform linear interpolation to find the output values (yi) at the input points xi
Nq = interp1(input_Nq, output_Nq, atand(2/3*tand(fi_end)), 'linear');

    % Calculate w; where w is the displacement of the end-bearin
 w =abs(r_tip_end_fixed(3,:))-abs(r_tip_end(3,:));
    
    for i=1:size(w,2)
        if w(i)<=0
            w(i)=0;
        end
    end
    
    WoverD = w / Do;
    Reshaped_w = reshape(w,[n,nq-1]);

    %the nodes which are not under pressure , their q-z force = ZERO
properties = [0, 0; 0.002, 0.25; 0.013, 0.5; 0.042, 0.75; 0.073, 0.9; 0.1, 1];
% Extract the input and output values from the properties matrix
input = properties(:, 1);
output = properties(:, 2);
% Generate a set of input values (xi) for interpolation
% This spans from the minimum to the maximum input values in your dataset

% Perform linear interpolation to find the output values (yi) at the input points xi
qend_over_qb = interp1(input, output, WoverD, 'linear');


    
    % Calculate soil stiffness
%     sigma_p = gamma * tip_end(3,:);
    sigma_p2 = -gamma * r_tip_end(3,:);
        Reshaped_sigma_p2 = reshape(sigma_p2,[n,nq-1]);

%     Initialize the area matrix
    area = zeros(n, nq - 1);

    k=0;
    for i = 1:n
        for k = 2:nq
            area(i, k - 1) = pi * (R_split(k)^2 - R_split(k-1)^2) / n;
        end
    end
    indices_wN = find(w <= 0)';
   area =  reshape(area,[1,n*(nq - 1)]);
 qb_lid = area .* sigma_p2*Nq;
    qb2 = qb_lid .* qend_over_qb;
    qb2(indices_wN)=0;
end
