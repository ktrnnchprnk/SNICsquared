function [U,eval,v1,v2] = compute_fp(f,IC1,IC2,delta1,delta2,n)

% Create a grid of initial guesses for fsolve
X1 = IC1(1):delta1:IC1(2);
X2 = IC2(1):delta2:IC2(2);
[E, I] = meshgrid(X1, X2);
gridpts = [E(:), I(:)];

% Initialize storage for fixed points and solver exit flags
FixedPTS = NaN(length(gridpts), 2);
FixedPTS2 = NaN(length(gridpts), 2);
Exit = NaN(length(gridpts), 1);

% Set fsolve options for high precision and no display
options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);

% Try to find fixed points from each initial guess
for i = 1:length(gridpts)
    [FixedPTS(i,:), fval, Exit(i)] = fsolve(@(x) f(x), gridpts(i,:), options);
    if Exit(i) == 1  % Convergence successful
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

% Remove unsuccessful or NaN results
FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2), 2), :);

% Round fixed points to eliminate numerical duplicates, then find unique ones
UniqueEig = unique(round(FixedPTS3, n), 'rows');

% For each unique fixed point, compute Jacobian and its eigenvalues/vectors
for i = 1:size(UniqueEig, 1)
    U(i,1) = UniqueEig(i,1);
    U(i,2) = UniqueEig(i,2);
    
    A = jacobianest(@(x) f(x), UniqueEig(i,:));  % Estimate Jacobian matrix
    [v, evalue] = eig(A);                        % Eigen-decomposition

    eval(i,:) = [evalue(1,1); evalue(2,2)];      % Store eigenvalues
    v1(i,:) = v(:,1);                            % Store first eigenvector
    v2(i,:) = v(:,2);                            % Store second eigenvector
end

end
