clear variables

input_path = "../build/mat"

% load([input_path '/p001.mat'])
load([input_path '/h001.mat'])

x = vertices(1, :)';
y = vertices(2, :)';

x_m = max(x);
y_m = max(y);

tri = double(faces');
z = data;

% graphics_toolkit gnuplot
figure
trisurf(tri, x, y, z, 'FaceColor', 'interp')

% figure
% hold off
% trimesh(tri, x, y, z * 0, z, 'EdgeAlpha', 0, 'FaceColor', 'interp')
% hold on
% trimesh(tri, -x, y, z * 0, z, 'EdgeAlpha', 0, 'FaceColor', 'interp')
% view(2)

colormap jet
colorbar

xlabel('r [nm]')
ylabel('z [nm]')

xlim([0 x_m])
ylim([0 y_m])
