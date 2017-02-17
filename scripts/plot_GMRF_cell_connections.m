clear;
clc;
close all;

cell_size = 1;

% Load cell connections (from ROS pkg)
relations = dlmread('/home/jgmonroy/GMRF.txt');
%relations = relations + 1;  %To correct Matlab index starting 1

N = max (relations(:,1)) +1;
n_x_cells = relations(2,2);
n_y_cells = N/n_x_cells;
relations = sortrows(relations);

figure();
hold on;
for i=1:size(relations,1)
    a = relations(i,1);
    b = relations(i,2);
    if a < b
        plot([mod(a,n_x_cells) mod(b,n_x_cells)],[floor(a/n_x_cells) floor(b/n_x_cells)],'Color',[0 0 1],'LineWidth',2);
    end;
end;