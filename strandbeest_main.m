%initialize leg_params structure
leg_params = struct();
%number of vertices in linkage
leg_params.num_vertices = 7;
%number of links in linkage
leg_params.num_linkages = 10;
%matrix relating links to vertices
leg_params.link_to_vertex_list = ...
[ 1, 3;... %link 1 adjacency
3, 4;... %link 2 adjacency
2, 3;... %link 3 adjacency
2, 4;... %link 4 adjacency
4, 5;... %link 5 adjacency
2, 6;... %link 6 adjacency
1, 6;... %link 7 adjacency
5, 6;... %link 8 adjacency
5, 7;... %link 9 adjacency
6, 7 ... %link 10 adjacency
];

%list of lengths for each link
%in the leg mechanism
leg_params.link_lengths = ...
[ 50.0,... %link 1 length
55.8,... %link 2 length
41.5,... %link 3 length
40.1,... %link 4 length
39.4,... %link 5 length
39.3,... %link 6 length
61.9,... %link 7 length
36.7,... %link 8 length
65.7,... %link 9 length
49.0 ... %link 10 length
];

%length of crank shaft
leg_params.crank_length = 15.0;
%fixed position coords of vertex 0
leg_params.vertex_pos0 = [0;0];
%fixed position coords of vertex 2
leg_params.vertex_pos2 = [-38.0;-7.8];

%column vector of initial guesses
%for each vertex location.
%in form: [x1;y1;x2;y2;...;xn;yn]
vertex_coords_guess = [...
[ 0; 50];... %vertex 1 guess
[ -50; 0];... %vertex 2 guess
[ -50; 50];... %vertex 3 guess
[-100; 0];... %vertex 4 guess
[-100; -50];... %vertex 5 guess
[ -50; -50];... %vertex 6 guess
[ -50; -100]... %vertex 7 guess
];

vertex_coords = compute_coords(vertex_coords_guess, leg_params, 0);
leg_drawing = initialize_leg_drawing(leg_params);
axis square; axis([-120, 20, -100, 40]);
for theta = 0:0.1:100
    vertex_coords = compute_coords(vertex_coords, leg_params, theta);
    update_leg_drawing(vertex_coords, leg_drawing, leg_params);
    pause(0.01);
end


%Converts from the column vector form of the coordinates to a
%friendlier matrix form
%INPUTS:
%coords_in = [x1;y1;x2;y2;...;xn;yn] (2n x 1 column vector)
%OUTPUTS:
%coords_out = [x1,y1;x2,y2;...;xn,yn] (n x 2 matrix)
function coords_out = column_to_matrix(coords_in)
    num_coords = length(coords_in);
    coords_out = [coords_in(1:2:(num_coords-1)),coords_in(2:2:num_coords)];
end

%Converts from the matrix form of the coordinates back to the
%original column vector form
%INPUTS:
%coords_in = [x1,y1;x2,y2;...;xn,yn] (n x 2 matrix)
%OUTPUTS:
%coords_out = [x1;y1;x2;y2;...;xn;yn] (2n x 1 column vector)
function coords_out = matrix_to_column(coords_in)
    num_coords = 2*size(coords_in,1);
    coords_out = zeros(num_coords,1);
    coords_out(1:2:(num_coords-1)) = coords_in(:,1);
    coords_out(2:2:num_coords) = coords_in(:,2);
end

%Error function that encodes the link length constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
%             in the linkage. There are two ways that I would recommend stacking
%             the coordinates. You could alternate between x and y coordinates:
%             i.e. vertex_coords = [x1;y1;x2;y2;...;xn;y_n], or alternatively
%             you could do all the x's first followed by all of the y's
%             i.e. vertex_coords = [x1;x2;...xn;y1;y2;...;yn]. You could also do
%             something else entirely, the choice is up to you.
%leg_params: a struct containing the parameters that describe the linkage
%          importantly, leg_params.link_lengths is a list of linakge lengths
%          and leg_params.link_to_vertex_list is a two column matrix where
%          leg_params.link_to_vertex_list(i,1) and
%          leg_params.link_to_vertex_list(i,2) are the pair of vertices connected
%          by the ith link in the mechanism
%OUTPUTS:
%length_errors: a column vector describing the current distance error of the ith
%               link specifically, length_errors(i) = (xb-xa)ˆ2 + (yb-ya)ˆ2 - d_iˆ2
%               where (xa,ya) and (xb,yb) are the coordinates of the vertices that
%               are connected by the ith link, and d_i is the length of the ith link
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

%Error function that encodes the fixed vertex constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
%               same input as link_length_error_func
%leg_params: a struct containing the parameters that describe the linkage
%            importantly, leg_params.crank_length is the length of the crank
%            and leg_params.vertex_pos0 and leg_params.vertex_pos2 are the
%            fixed positions of the crank rotation center and vertex 2.
%theta: the current angle of the crank
%OUTPUTS:
%coord_errors: a column vector of height four corresponding to the differences
%              between the current values of (x1,y1),(x2,y2) and
%              the fixed values that they should be
function coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta)
    x1_error = vertex_coords(1) - leg_params.vertex_pos0(1) + leg_params.crank_length * cos(theta);
    y1_error = vertex_coords(2) - leg_params.vertex_pos0(2) + leg_params.crank_length * sin(theta);
    x2_error = vertex_coords(3) - leg_params.vertex_pos2(1);
    y2_error = vertex_coords(4) - leg_params.vertex_pos2(2);
    coord_errors = [x1_error; y1_error; x2_error; y2_error];
end

%Error function that encodes all necessary linkage constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%error_vec: a vector describing each constraint on the linkage
%           when error_vec is all zeros, the constraints are satisfied
function error_vec = linkage_error_func(vertex_coords, leg_params, theta)
    distance_errors = link_length_error_func(vertex_coords, leg_params);
    coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta);
    error_vec = [distance_errors;coord_errors];
end

%Computes the vertex coordinates that describe a legal linkage configuration
%INPUTS:
%vertex_coords_guess: a column vector containing the (x,y) coordinates of every vertex
%                      these coords are just a GUESS! It's used to seed Newton's method
%leg_params: a struct containing the parameters that describe the linkage
%theta: the desired angle of the crank
%OUTPUTS:
%vertex_coords_root: a column vector containing the (x,y) coordinates of every vertex
%                    these coords satisfy all the kinematic constraints!
function vertex_coords_root = compute_coords(vertex_coords_guess, leg_params, theta)
    fun = @(vertex_coords) linkage_error_func(vertex_coords, leg_params, theta);
    vertex_coords_root = multi_newton_solver(fun, vertex_coords_guess, struct());
end

%Creates a set of plotting objects to keep track of each link drawing
%each vertex drawing, and the crank drawing
%INPUTS:
%leg_params: a struct containing the parameters that describe the linkage
%OUTPUTS:
%leg_drawing: a struct containing all the plotting objects for the linkage
% leg_drawing.linkages is a cell array, where each element corresponds
% to a plot of a single link (excluding the crank)
% leg_drawing.crank is a plot of the crank link
% leg_drawing.vertices is a cell array, where each element corresponds
% to a plot of one of the vertices in the linkage
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

%Updates the plot objects that visualize the leg linkage
%for the current leg configuration
%INPUTS:
%complete_vertex_coords: a column vector containing the (x,y) coordinates of every vertex
%leg_drawing: a struct containing all the plotting objects for the linkage
% leg_drawing.linkages is a cell array, where each element corresponds
% to a plot of a single link (excluding the crank)
% leg_drawing.crank is a plot of the crank link
% leg_drawing.vertices is a cell array, where each element corresponds
% to a plot of one of the vertices in the linkage
function update_leg_drawing(complete_vertex_coords, leg_drawing, leg_params)
    vertex_matrix = column_to_matrix(complete_vertex_coords);
    %iterate through each link, and update corresponding link plot
    for linkage_index = 1:leg_params.num_linkages
        %linkage_index is the label of the current link
        xa = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 1), 1);
        ya = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 1), 2);
        xb = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 2), 1);
        yb = vertex_matrix(leg_params.link_to_vertex_list(linkage_index, 2), 2);
        %line_x and line_y should both be two element arrays containing
        %the x and y coordinates of the line segment describing the current link
        line_x = [xa, xb];
        line_y = [ya, yb];
        set(leg_drawing.linkages{linkage_index},'xdata',line_x,'ydata',line_y);
    end
    %iterate through each vertex, and update corresponding vertex plot
    for vertex_index = 1:leg_params.num_vertices
        %vertex_index is the label of the current vertex
        %dot_x and dot_y should both be scalars
        %specifically the x and y coordinates of the corresponding vertex
        dot_x = vertex_matrix(vertex_index, 1);
        dot_y = vertex_matrix(vertex_index, 2);
        set(leg_drawing.vertices{vertex_index},'xdata',dot_x,'ydata',dot_y);
    end
    %crank_x and crank_y should both be two element arrays
    %containing the x and y coordinates of the line segment describing the crank
    crank_x = [leg_params.vertex_pos0(1), vertex_matrix(1, 1)];
    crank_y = [leg_params.vertex_pos0(2), vertex_matrix(1, 2)];
    set(leg_drawing.crank,'xdata',crank_x,'ydata',crank_y);
end