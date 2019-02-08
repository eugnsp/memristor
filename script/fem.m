clear variables

input_path = "../build/mat"

%graphics_toolkit fltk

k = 5;
file = [input_path "/p" num2str(k, "%.3d") ".mat"];
load(file)

x = vertices(1, :)';
y = vertices(2, :)';

x_m = max(x);
y_m = max(y);

tri = double(faces');
z = data;

figure
hold off
trisurf(tri, x, y, z, 'FaceColor', 'interp')
%hold on
%trisurf(tri, -x, y, z, 'FaceColor', 'interp')
%view(2)

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

xlim([-x_m x_m])
%ylim([0 40])
%caxis([300 650])
  
title('Temperature distribution')
