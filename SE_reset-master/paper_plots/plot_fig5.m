function plot_fig5
% Plotting for Figure 5
% Show relative RMS error and runtime as a function of P

  load('data/fig5.mat');

  markers = {'^', 'd', 's', 'o'};
  width = 1;

  sfigure(1); clf
  col = get(gca, 'ColorOrder');
  colors = [0 0 0; col(1,:); col(2,:); col(3,:)];
  for i=1:numel(per)
    for w=1:numel(windows)
      P_err = squeeze(err(i,w,:));
      semilogy(Pvec, P_err, sprintf('%s-', markers{i}), 'Color', colors(w,:), ...
               'LineWidth', width, 'DisplayName', ...
               sprintf('%s %dP', replace(windows{w}, '_', ' '), per(i)));
      hold on
    end
  end
  xlabel('P')
  ylabel('RMS error (rel.)')
  legend('show', 'Location', 'NorthEast')
  grid on
  title('Error vs P')
  ylim([1e-15 1e0])
end

% tim = time.total - exclude_pre*(time.pre + time.prefft);
% tim2 = time.grid + time.int;
% exclude_pre = true; % exclude precomputation from timing
