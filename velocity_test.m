leg_params = struct();
leg_params.num_vertices = 7;
leg_params.num_linkages = 10;
leg_params.link_to_vertex_list = [ 1, 3;3, 4;2, 3;2, 4;4, 5;2, 6;1, 6;5, 6;5, 7;6, 7];
leg_params.link_lengths = [ 50.0,55.8,41.5,40.1,39.4,39.3,61.9,36.7,65.7,49.0];
leg_params.crank_length = 15.0;
leg_params.vertex_pos0 = [0;0];
leg_params.vertex_pos2 = [-38.0;-7.8];
vertex_coords_guess = [[ 0; 50];[ -50; 0];[ -50; 50];[-100; 0];[-100; -50];[ -50; -50];[ -50; -100]];



vertex_coords = compute_coords(vertex_coords_guess, leg_params, 0);

%create 4x4 identity matrix and 4x10 of zeros
I = eye(4);
ZM = zeros(4,10);
IZM = [I,ZM];

%create 1x10 of zeros
ZB = zeros(10,1);

%calculate jacobian of link length error
J = approximate_jacobian(@(X) link_length_error_func(X, leg_params),vertex_coords_guess);

%calculate partial derivative of X1 Y1 X2 Y2
partial_coord_errors = partial_fixed_coord_error_func(vertex_coords, leg_params, theta);

leg_drawing = initialize_leg_drawing(leg_params);
axis square; axis([-120, 20, -100, 40]);
for theta = 0:0.1:100
    vertex_coords = compute_coords(vertex_coords, leg_params, theta);
    J = approximate_jacobian(@(X) link_length_error_func(X, leg_params),vertex_coords);
    M = [IZM;J];
    partial_coord_errors = partial_fixed_coord_error_func(vertex_coords, leg_params, theta);
    B = [partial_coord_errors;ZB];
    %use matrix division to calculate the DV/Dtheta
    dv_dtheta = M\B

    update_leg_drawing(vertex_coords, leg_drawing, leg_params);
    plot()
 
    pause(0.01);
end



function coords_out = column_to_matrix(coords_in)
    num_coords = length(coords_in);
    coords_out = [coords_in(1:2:(num_coords-1)),coords_in(2:2:num_coords)];
end
function coords_out = matrix_to_column(coords_in)
    num_coords = 2*size(coords_in,1);
    coords_out = zeros(num_coords,1);
    coords_out(1:2:(num_coords-1)) = coords_in(:,1);
    coords_out(2:2:num_coords) = coords_in(:,2);
end
function length_errors = link_length_error_func(vertex_coords, leg_params)
    length_errors = zeros(leg_params.num_linkages, 1);
    vertex_matrix = column_to_matrix(vertex_coords);
    for i = 1:leg_params.num_linkages
        xa = vertex_matrix(leg_params.link_to_vertex_list(i, 1), 1);
        ya = vertex_matrix(leg_params.link_to_vertex_list(i, 1), 2);
        xb = vertex_matrix(leg_params.link_to_vertex_list(i, 2), 1);
        yb = vertex_matrix(leg_params.link_to_vertex_list(i, 2), 2);
        d_i = leg_params.link_lengths(i);
        length_errors(i) = (xb - xa)^2 + (yb - ya)^2 - d_i^2;
    end
end
function coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta)
    x1_error = vertex_coords(1) - leg_params.vertex_pos0(1) + leg_params.crank_length * cos(theta);
    y1_error = vertex_coords(2) - leg_params.vertex_pos0(2) + leg_params.crank_length * sin(theta);
    x2_error = vertex_coords(3) - leg_params.vertex_pos2(1);
    y2_error = vertex_coords(4) - leg_params.vertex_pos2(2);
    coord_errors = [x1_error; y1_error; x2_error; y2_error];
end
%function for the partial derivitaves of the defined coordinates
function partial_coord_errors = partial_fixed_coord_error_func(vertex_coords, leg_params, theta)
    x1_error = vertex_coords(1) - leg_params.vertex_pos0(1) + leg_params.crank_length * -sin(theta);
    y1_error = vertex_coords(2) - leg_params.vertex_pos0(2) + leg_params.crank_length * cos(theta);
    x2_error = 0;
    y2_error = 0;
    partial_coord_errors = [x1_error; y1_error; x2_error; y2_error];
end
function error_vec = linkage_error_func(vertex_coords, leg_params, theta)
    distance_errors = link_length_error_func(vertex_coords, leg_params);
    coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta);
    error_vec = [distance_errors;coord_errors];
end
function vertex_coords_root = compute_coords(vertex_coords_guess, leg_params, theta)
    fun = @(vertex_coords) linkage_error_func(vertex_coords, leg_params, theta);
    vertex_coords_root = multi_newton_solver(fun, vertex_coords_guess, struct());
end
function leg_drawing = initialize_leg_drawing(leg_params)
    leg_drawing = struct();
    leg_drawing.linkages = cell(leg_params.num_linkages,1);
    for linkage_index = 1:leg_params.num_linkages
        leg_drawing.linkages{linkage_index} = line([0,0],[0,0],'color','k','linewidth',2);
    end
    leg_drawing.crank = line([0,0],[0,0],'color','k','linewidth',1.5);
    leg_drawing.vertices = cell(leg_params.num_vertices,1);
    for vertex_index = 1:leg_params.num_vertices
        leg_drawing.vertices{vertex_index} = line([0],[0],'marker',...
            'o','markerfacecolor','r','markeredgecolor','r','markersize',8);
    end
end
function update_leg_drawing(complete_vertex_coords, leg_drawing, leg_params)
    vertex_matrix = column_to_matrix(complete_vertex_coords);
    for linkage_index = 1:leg_params.num_linkages
        xa = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 1), 1);
        ya = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 1), 2);
        xb = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 2), 1);
        yb = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 2), 2);
        line_x = [xa, xb];
        line_y = [ya, yb];
        set(leg_drawing.linkages{linkage_index},'xdata',line_x,'ydata',line_y);
    end
    for vertex_index = 1:leg_params.num_vertices
        dot_x = vertex_matrix(vertex_index, 1);
        dot_y = vertex_matrix(vertex_index, 2);
        set(leg_drawing.vertices{vertex_index},'xdata',dot_x,'ydata',dot_y);
    end
    crank_x = [leg_params.vertex_pos0(1), vertex_matrix(1, 1)];
    crank_y = [leg_params.vertex_pos0(2), vertex_matrix(1, 2)];
    set(leg_drawing.crank,'xdata',crank_x,'ydata',crank_y);
end



% Computes the theta derivatives of each vertex coordinate for the Jansen linkage
% INPUTS:
% vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% these are assumed to be legal values that are roots of the error funcs!
% leg_params: a struct containing the parameters that describe the linkage
% theta: the current angle of the crank
% OUTPUTS:
% dVdtheta: a column vector containing the theta derivates of each vertex coord
% function dVdtheta = compute_velocities(vertex_coords, leg_params, theta)
% 
% end