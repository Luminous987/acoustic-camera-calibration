function closed_form_error = closed_form_calib_func(input)
% This file performs least square SLAM

%% parameters

% gonna estimate clock drift?
est_drift_on = 0;
% gonna estimate starting time delay?
est_delay_on = 0;
% display starting time delay estimation error?
display_delay_error_on = 0;
% display norm(dx) for each iteration?
display_norm_dx_on = 0;

% the maximum number of iterations
numIterations = 50;

% maximum allowed dx
EPSILON = input.eps;%1.5/1e-2;

% Error
err = 0;

% load the graph into the variable "g"
load(input.graph_file);
mic_position = linearize_and_solve_v2(g);
%% start slSLAM
g.x(1:24) = mic_position;
% carry out the iterations
for i = 1:numIterations
%   disp(['Performing iteration ', num2str(i)]);
  Mic_pos_err = compute_RMS_error(g);
  disp(Mic_pos_err);
  % solve the dx
  [dx,H] = linearize_and_solve_with_H_new(g);

  % TODO: apply the solution to the state vector g.x
  g.x(1:3*g.M,1) = g.x(1:3*g.M,1) + dx;


  % compute current error
  % err = compute_global_error(g);

  % TODO: implement termination criterion as suggested on the sheet
  if display_norm_dx_on>0
    disp(['norm(dx) = ' num2str(norm(dx))]);
  end
  
  if (norm(dx)<EPSILON)
    break;
  end


end
closed_form_error = Mic_pos_err;
% 
% plot the current state of the graph
figure;
closed_form_plot_graph_with_cov(g, i, H,mic_position);
title(input.fig.title);
view(input.fig.view_a, input.fig.view_e);

end

