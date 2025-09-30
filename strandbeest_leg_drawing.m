function strandbeest_leg_drawing(vertex_coords_guess, leg_params)
    vertex_coords = strandbeest_compute_coords(vertex_coords_guess, leg_params, 0);
    leg_drawing = initialize_leg_drawing(leg_params);
    axis square; axis([-120, 20, -100, 40]);
    for theta = 0:0.1*pi:100*pi
        vertex_coords = strandbeest_compute_coords(vertex_coords, leg_params, theta);
        update_leg_drawing(vertex_coords, leg_drawing, leg_params);
        pause(0.01);
    end
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