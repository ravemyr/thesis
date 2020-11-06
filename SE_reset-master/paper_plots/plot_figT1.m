function plot_figT1
% Plotting for Figure T1
% Show relative RMS error and runtime as a function of P

  load('data/figT1.mat');

  markers = {'^', 'd', 's', 'o'};
  width = 1;

  sfigure(2); clf
  col = get(gca, 'ColorOrder');
  colors = [0 0 0; col(1,:); col(2,:); col(3,:)];
  for i=1:numel(per)
    sfigure(i+1); clf
    for w=1:numel(windows)
      for k=1:numel(Pvec)
        M_err = squeeze(err(i,w,k,:));
        semilogy(pi*Mvec/xi, M_err, sprintf('%s-', markers{i}), 'Color', colors(w,:), ...
                 'LineWidth', width, 'DisplayName', ...
                 sprintf('%s %dP (P=%d)', replace(windows{w}, '_', ' '), per(i), Pvec(k)));
        hold on
      end
    end
    xlabel('kinf/xi')
    ylabel('RMS error (rel.)')
    legend('show', 'Location', 'best')
    grid on
    title(sprintf('Error vs M (%dP)', per(i)))
    ylim([1e-15 1e0])
  end
end

% tim = time.total - exclude_pre*(time.pre + time.prefft);
% tim2 = time.grid + time.int;
% exclude_pre = true; % exclude precomputation from timing
