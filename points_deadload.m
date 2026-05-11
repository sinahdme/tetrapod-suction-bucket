
function [skirt_pts, tip_pts, lid_ring_pts] = points_deadload(R_split, th_deg, n, n_l, Z_lid, l_b, l, ro_b, ri_b, r_pile)
	circle_split = linspace(0,360,n+1);

	% Skirt nodes: n per ring, n_l rings => n*n_l
	j = 0;
	for k = 0:l_b/n_l:l_b-(l_b/n_l)
		for i = 2:n+1
			j = j + 1;
			skirt_pts(:,j) = [l*cosd(th_deg)+ro_b*cosd(circle_split(i)), ...
			                  l*sind(th_deg)+ro_b*sind(circle_split(i)), ...
			                 -k-(l_b/n_l)/2]';
		end
	end

	% Tip nodes: n points
	j = 0;
	for i = 2:n+1
		j = j + 1;
		tip_pts(:,j) = [l*cosd(th_deg)+r_pile*cosd(circle_split(i)), ...
		                l*sind(th_deg)+r_pile*sind(circle_split(i)), ...
		               -l_b]';
	end

	% Lid-ring nodes: (n_q-1) rings, n per ring => n*(n_q-1)
	R_sec = R_split(2) - R_split(1);
	lid_ring_pts = zeros(3, n*(numel(R_split)-1));
	for jj = 2:numel(R_split)
		j = 0;
		for i = 2:n+1
			j = j + 1;
			idx = (jj-2)*n + j;
			lid_ring_pts(:,idx) = [l*cosd(th_deg)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)), ...
			                       l*sind(th_deg)+(R_split(jj)-R_sec/2)*sind(circle_split(i)), ...
			                      -Z_lid]';
		end
	end
end

