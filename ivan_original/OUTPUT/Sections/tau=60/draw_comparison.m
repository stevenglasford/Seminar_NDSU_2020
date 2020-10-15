graphics_toolkit gnuplot;


clc;
clear;


tau = 60.0;

T = 60.0;

k = 1.0 / 6.0;
alpha = 0.2;
beta = 0.1;
gamma = 0.06;
mu = 0.03;

n_1 = 9.0;
n_2 = 1.0;

N = 501;
N_full = 750;

M_1_min = 0.0;
M_1_max = 1.002;

M_2_min = 0.0;
M_2_max = 1.002;

M_1_left = 0.0;
M_1_right = 1.0;

M_2_left = 0.0;
M_2_right = 1.0;

hat_M_aster = 0.35596;

VF = load( "VF240.dat" );


exp_minus_mu_t = exp(-mu * (T - tau));

dM_1 = (M_1_max - M_1_min) / N;
dM_2 = (M_2_max - M_2_min) / N;

two_dM_1 = 2.0 * dM_1;
two_dM_2 = 2.0 * dM_2;

M_1_arr = M_1_min : dM_1 : M_1_max;
M_2_arr = M_2_min : dM_2 : M_2_max;


value_func_full = zeros(N_full + 1, N_full + 1);

k = 1;

for j = 1 : (N_full + 1)

	for i = 1 : (N_full + 1)
	
		value_func_full(i, j) = VF(k, 3);
		
		k = k + 1;
	
	end
	
end
	
clear VF;
	
value_func = value_func_full( 1 : (N + 1), 1 : (N + 1) );
	
clear value_func_full;


u_1 = zeros(N + 1, N + 1);
u_2 = zeros(N + 1, N + 1);


for i = 1 : (N + 1)

	for j = 1 : (N + 1)
	
		if ((i > 1) && (i < N + 1))
		
			if ( exp_minus_mu_t + (value_func(i + 1, j) - value_func(i - 1, j)) / two_dM_1 < 0.0 )
				u_1(i, j) = 0.0;
			else
				u_1(i, j) = 1.0;
			end
		
		end
		
		if ((j > 1) && (j < N + 1))
		
			if ( exp_minus_mu_t - (value_func(i, j + 1) - value_func(i, j - 1)) / two_dM_2 < 0.0 )
				u_2(i, j) = 0.0;
			else
				u_2(i, j) = 1.0;
			end
		
		end
	
	end

end


for j = 1 : (N + 1)

	u_1(1, j) = u_1(2, j);
	u_1(N + 1, j) = u_1(N, j);
	
	u_2(1, j) = u_2(2, j);
	u_2(N + 1, j) = u_2(N, j);
 
end


for i = 1 : (N + 1)

	u_1(i, 1) = u_1(i, 2);
	u_1(i, N + 1) = u_1(i, N);
	
	u_2(i, 1) = u_2(i, 2);
	u_2(i, N + 1) = u_2(i, N);
	
end


u_1_sw_surf = zeros(1, N + 1);
u_2_sw_surf = zeros(1, N + 1);


for i = 1 : (N + 1)


	ind_arr_1 = find( u_1(:, i) > 0.5 );
	ind_sw_1 = ind_arr_1(1);
	
	if (ind_sw_1 > 1)
		coord_sw_1 = 0.5 * (M_1_arr(ind_sw_1 - 1) + M_1_arr(ind_sw_1));
	else
		coord_sw_1 = M_1_arr(ind_sw_1);
	end
	
	u_1_sw_surf(i) = coord_sw_1;
	
	
	ind_arr_2 = find( u_2(i, :) > 0.5 );
	ind_sw_2 = ind_arr_2(1);
	
	if (ind_sw_2 > 1)
		coord_sw_2 = 0.5 * (M_2_arr(ind_sw_2 - 1) + M_2_arr(ind_sw_2));
	else
		coord_sw_2 = M_2_arr(ind_sw_2);
	end
	
	u_2_sw_surf(i) = coord_sw_2;


end


clear u_1;
clear u_2;


traj_1 = load( "trajectory_regime_1.txt" );
traj_2 = load( "trajectory_regime_2.txt" );
traj_3 = load( "trajectory_regime_3.txt" );
traj_4 = load( "trajectory_regime_4.txt" );


plot( traj_1(:, 2), traj_1(:, 3), '-', 'LineWidth', 25, 'Color', [0.5, 0.5, 0.5] );

hold on;

plot( [-6.0, -6.0], [-5.0, -5.0], '--k', 'LineWidth', 14 );

hold on;

plot( traj_2(:, 2), traj_2(:, 3), '-', 'LineWidth', 25, 'Color', [0.5, 0.5, 0.5] );

hold on;

plot( traj_3(:, 2), traj_3(:, 3), '-', 'LineWidth', 25, 'Color', [0.5, 0.5, 0.5] );

hold on;

plot( traj_4(:, 2), traj_4(:, 3), '-', 'LineWidth', 25, 'Color', [0.5, 0.5, 0.5] );

hold on;

plot( M_1_arr, u_2_sw_surf, '--k', 'LineWidth', 14 );

hold on;

plot( u_1_sw_surf, M_2_arr, '--k', 'LineWidth', 14 );

hold on;

plot( hat_M_aster, hat_M_aster, 'o', 'MarkerSize', 13, 'MarkerFaceColor', 'Black', 'MarkerEdgeColor', 'Black' );


xlabel( 'm_1', 'fontsize', 21, 'fontname', 'Franklin Gothic Book' );
ylabel( 'm_2', 'fontsize', 21, 'fontname', 'Franklin Gothic Book' );

axis( [0.0, 1.0, 0.0, 1.0] );

set( gca, 'fontsize', 18, 'fontname', 'Franklin Gothic Book' );

grid on;

set( gca, 'XTick', [0, 0.2, hat_M_aster, 0.6, 0.8, 1] );
set( gca, 'XTickLabel', {'0', '0.2', 'm^{**}', '0.6', '0.8', '1'} );

set( gca, 'YTick', [0, 0.2, hat_M_aster, 0.6, 0.8, 1] );
set( gca, 'YTickLabel', {'0', '0.2', 'm^{**}', '0.6', '0.8', '1'} );

leg = legend( '\Sigma_1-\Sigma_4 (Fig. 3)', 'Fig. 5d' );
set( leg, 'fontsize', 18, 'fontname', 'Franklin Gothic Book' );


print -color sw_surf_comparison.eps
