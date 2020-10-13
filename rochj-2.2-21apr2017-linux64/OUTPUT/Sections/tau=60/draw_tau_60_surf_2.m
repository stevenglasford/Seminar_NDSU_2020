graphics_toolkit gnuplot;


clc;
clear;


draw_surf = 1;
draw_controls = 0;


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


%surf_axis_limits = [0.0,  1.0, 0.0, 1.0, -0.8, 1.0];

%surf_ind_M_1_left = 1;
%surf_ind_M_1_right = 501;

%surf_ind_M_2_left = 1;
%surf_ind_M_2_right = 501;

surf_axis_limits = [0.1, 0.7, 0.1, 0.7, -0.7, 0.8];

surf_ind_M_1_left = 51;
surf_ind_M_1_right = 351;

surf_ind_M_2_left = 51;
surf_ind_M_2_right = 351;


surf_ind_step = 2;


rhs_1  =  @(m_1, m_2, contr_1)  m_1 * ( (1.0 - contr_1) * (1.0 / (1.0 + beta * (n_1 * m_1 + n_2 * m_2))) * alpha / (m_1 + k) - gamma );
rhs_2  =  @(m_1, m_2, contr_2)  m_2 * ( (1.0 - contr_2) * (1.0 / (1.0 + beta * (n_1 * m_1 + n_2 * m_2))) * alpha / (m_2 + k) - gamma );


vel_field_ind_arr = 10 : 37 : 491;
N_vel_field_ind_arr = size( vel_field_ind_arr, 2 );


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


vel_field_1 = zeros( 1, N_vel_field_ind_arr );
vel_field_2 = zeros( 1, N_vel_field_ind_arr );


for i = 1 : N_vel_field_ind_arr

	m_1 = M_1_arr( vel_field_ind_arr(i) );

	for j = 1 : N_vel_field_ind_arr
	
		m_2 = M_2_arr( vel_field_ind_arr(j) );
		
		contr_1 = u_1( vel_field_ind_arr(i), vel_field_ind_arr(j) );
		contr_2 = u_2( vel_field_ind_arr(i), vel_field_ind_arr(j) );
		
		vel_field_1(i, j) = rhs_1(m_1, m_2, contr_1);
		vel_field_2(i, j) = rhs_2(m_1, m_2, contr_2);
	
	end

end


clear u_1;
clear u_2;


if ( draw_surf )
	
	M_1_arr_part = M_1_arr( surf_ind_M_1_left : surf_ind_step : surf_ind_M_1_right );
	M_2_arr_part = M_2_arr( surf_ind_M_2_left : surf_ind_step : surf_ind_M_2_right );
	value_func_part = value_func( surf_ind_M_1_left : surf_ind_step : surf_ind_M_1_right,...
								  surf_ind_M_2_left : surf_ind_step : surf_ind_M_2_right );
	
	[Y_part, X_part] = meshgrid( M_2_arr_part, M_1_arr_part );
	
	surf( X_part, Y_part, value_func_part, 'edgecolor', 'none' );
	
	title( 'V \cdot 10^{-4} ,  (d) \tau = T = 60 (t = 0)', 'fontsize', 25, 'fontname', 'Franklin Gothic Book' );
	
	xlabel( 'm_1', 'fontsize', 20, 'fontname', 'Franklin Gothic Book' );
    ylabel( 'm_2', 'fontsize', 20, 'fontname', 'Franklin Gothic Book' );
    
    axis( surf_axis_limits );
    
    grid on;
    
    set( gca, 'fontsize', 18, 'fontname', 'Franklin Gothic Book' );
    
    colormap('gray');
    map = colormap;
    colormap( map(1 : 55, :) );
    
    view(120, 40);
    
    print -color "-S600,400" 'Value_function_tau=60.eps'
	
	clf;
	
	clear X_part;
	clear Y_part;
	clear value_func_part;

end


clear value_func;


if (draw_controls)


	M_1_arr_part = M_1_arr( vel_field_ind_arr );
	M_2_arr_part = M_2_arr( vel_field_ind_arr );
	
	[Y_part, X_part] = meshgrid( M_2_arr_part, M_1_arr_part );


	ind_right = ceil( (M_2_right - M_2_left) / dM_2 ) + 1;
	
	area( [u_1_sw_surf(1 : ind_right), M_1_right], [M_2_arr(1 : ind_right), M_2_arr(ind_right)],...
		  'BaseValue', M_1_left, 'EdgeColor', 'Black', 'LineWidth', 1, 'FaceColor', [0.5, 0.5, 0.5] );
	
	hold on;
	
	plot( M_1_arr, M_1_arr, '--k', 'LineWidth', 1 );
	
	hold on;
	
	plot( hat_M_aster, hat_M_aster, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'Black', 'MarkerEdgeColor', 'Black' );
	
	hold on;
	
	plot( [hat_M_aster, hat_M_aster, 0], [0, hat_M_aster, hat_M_aster], '--k', 'LineWidth', 1 );
	
	hold on;
	
	plot( M_1_arr, u_2_sw_surf, '--k', 'LineWidth', 1 );
	
	hold on;
	
	quiver( X_part, Y_part, vel_field_1, vel_field_2, 0.6, '-k', 'LineWidth', 2 );
	
	title( 'u_1 ,  (d) \tau = T = 60 (t = 0)', 'fontsize', 25, 'fontname', 'Franklin Gothic Book' );
	
	xlabel( 'm_1', 'fontsize', 20, 'fontname', 'Franklin Gothic Book' );
	ylabel( 'm_2', 'fontsize', 20, 'fontname', 'Franklin Gothic Book' );
	
	axis( [M_1_left, M_1_right, M_2_left, M_2_right] );
	
	set( gca, 'fontsize', 18, 'fontname', 'Franklin Gothic Book' );
	
	set( gca, 'XTick', [0, 0.2, hat_M_aster, 0.6, 0.8, 1] );
	set( gca, 'XTickLabel', {'0', '0.2', 'm^{**}', '0.6', '0.8', '1'} );
	
	set( gca, 'YTick', [0, 0.2, hat_M_aster, 0.6, 0.8, 1] );
	set( gca, 'YTickLabel', {'0', '0.2', 'm^{**}', '0.6', '0.8', '1'} );
	
	print -color "-S600,400" 'u_1_tau=60.eps'
	
	clf;
	
	
	ind_right = ceil( (M_1_right - M_1_left) / dM_1 ) + 1;
	
	area( M_1_arr(1 : ind_right), u_2_sw_surf(1 : ind_right),...
		  'BaseValue', M_1_right, 'EdgeColor', 'Black', 'LineWidth', 1, 'FaceColor', [0.5, 0.5, 0.5] );
	
	hold on;
	
	plot( M_1_arr, M_1_arr, '--k', 'LineWidth', 1 );
	
	hold on;
	
	plot( hat_M_aster, hat_M_aster, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'Black', 'MarkerEdgeColor', 'Black' );
	
	hold on;
	
	plot( [hat_M_aster, hat_M_aster, 0], [0, hat_M_aster, hat_M_aster], '--k', 'LineWidth', 1 );
	
	hold on;
	
	plot( u_1_sw_surf, M_2_arr, '--k', 'LineWidth', 1 );
	
	hold on;
	
	quiver( X_part, Y_part, vel_field_1, vel_field_2, 0.6, '-k', 'LineWidth', 2 );
	
	title( 'u_2 ,  (d) \tau = T = 60 (t = 0)', 'fontsize', 25, 'fontname', 'Franklin Gothic Book' );
	
	xlabel( 'm_1', 'fontsize', 20, 'fontname', 'Franklin Gothic Book' );
	ylabel( 'm_2', 'fontsize', 20, 'fontname', 'Franklin Gothic Book' );
	
	axis( [M_1_left, M_1_right, M_2_left, M_2_right] );
	
	set( gca, 'fontsize', 18, 'fontname', 'Franklin Gothic Book' );
	
	set( gca, 'XTick', [0, 0.2, hat_M_aster, 0.6, 0.8, 1] );
	set( gca, 'XTickLabel', {'0', '0.2', 'm^{**}', '0.6', '0.8', '1'} );
	
	set( gca, 'YTick', [0, 0.2, hat_M_aster, 0.6, 0.8, 1] );
	set( gca, 'YTickLabel', {'0', '0.2', 'm^{**}', '0.6', '0.8', '1'} );
	
	print -color "-S600,400" 'u_2_tau=60.eps'
	
	clf;


end
