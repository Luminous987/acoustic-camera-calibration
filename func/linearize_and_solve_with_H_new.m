% performs one iteration of the Gauss-Newton algorithm
% each constraint is linearized and added to the Hessian

function [dx,H] = linearize_and_solve_with_H_new(g)

% nnz = nnz_of_graph(g);

% allocate the sparse H and the vector b
% H = spalloc(length(g.x), length(g.x), nnz);
% b = zeros(length(g.x), 1);
H = zeros(3*g.M, 3*g.M);
b = zeros(3*g.M, 1);
needToAddPrior = true;

% compute the addend term to H and b for each of our constraints
for eid = 1:length(g.edges)
    edge = g.edges(eid);

    x1 = g.x(edge.toIdx:edge.toIdx+2);   % the robot pose
    x2 = g.x(edge.fromIdx:edge.fromIdx+(3*g.M-1));     % the landmark

    % Computing the error and the Jacobians
    % e the error vector
    % A Jacobian wrt x1
    % B Jacobian wrt x2
    [e, A] = linearize_pose_landmark_constraint_new(x1, x2, edge.measurement,g);

    b(edge.fromIdx:edge.fromIdx+(3*g.M-1)) = (b(edge.fromIdx:edge.fromIdx+(3*g.M-1))' + (e')*edge.information*A)';
    H(edge.fromIdx:edge.fromIdx+(3*g.M-1),edge.fromIdx:edge.fromIdx+(3*g.M-1)) = H(edge.fromIdx:edge.fromIdx+(3*g.M-1),edge.fromIdx:edge.fromIdx+(3*g.M-1)) + A'*edge.information*A;
%     lambda_l = 1 * 1e-9;
%     max_error = 0.001;
%     Cov_e_l = zeros(7,7);
%     for i =1 : 7
%         % if abs(e(i,1)) > max_error
%         %     e(i,1) = sign(e(i,1)) * max_error;
%         % end
%     Cov_e_l(i,i) = e(i,1) * e(i,1);
%     end
% %   Cov_e_l = (e * e');
%     Cov_e_inv_l = inv(Cov_e_l + lambda_l * eye(size(Cov_e_l))) ;
%     edge.information = 2 * Cov_e_inv_l / 10^-0;

end

if (needToAddPrior)
  % add the prior for one pose of this edge
  % This fixes one node to remain at its current location

  % 1st mic
  H(1:3,1:3) = eye(3);
 
% solve the linear system, whereas the solution should be stored in dx
% Remember to use the backslash operator instead of inverting H

dx = H\(-b);

end
