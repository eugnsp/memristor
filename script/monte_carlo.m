clear variables
graphics_toolkit gnuplot

input_path = "../build/mat"
output_path = "../build/png"

load([input_path "/iv.mat"])

k = 0;
while (true)
  file = [input_path "/m" num2str(k, "%.3d") ".mat"];
  if ~exist(file)
    break
  endif
  
  load(file)

  d = .25;  
  x = double(x) * d;
  y = double(y) * d;
  z = double(z) * d;
  
  x_mid = max(x) / 2;
  y_mid = max(y) / 2;
  mask = (y == y_mid);
  
  hf = figure('visible', 'off');
  scatter(x(mask) - x_mid, z(mask), 'k')
  
  hold on
  plot([-1 -1 1 1], [0 5 5 0], 'r', 'linewidth', 2)
  plot(1.5 * cos(0:.1:6.4), 20 + 3 * sin(0:.1:6.4), 'r', 'linewidth', 2)
  
  axis equal
  xlim([-x_mid x_mid])
  ylim([0 max(z)])
  
  title(['Vacancy distribution, bias = ' num2str(v(k + 1), "%.2f") ' V'])
  xlabel('x [nm]')
  ylabel('y [nm]')
  
  print(hf, [output_path "/m" num2str(k, "%.3d") ".png"]);
  close(hf)
  
  k
  k = k + 1;
end

disp("Done.")
