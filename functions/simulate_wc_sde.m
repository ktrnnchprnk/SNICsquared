function [t, X] = simulate_wc_sde(par, x0, dt, T)

% Time vector
t = 0:dt:T;
N = length(t);

% Initialize state matrix
X = zeros(2, N);
X(:,1) = x0;

% Extract parameters
cee = par.cee; 
cie = par.cie; 
cei = par.cei; 
ai = par.ai; 
ae = par.ae; 
beta1 = par.beta1; 
thetai = par.thetai; 
thetae = par.thetae; 
Kp = par.Kp;  
alpha = par.alpha; 
sigma = par.noise_std;  % noise intensity

% Sigmoid function
phi = @(a,F,theta) 1./(1 + exp(-a.*(F - theta))) - 1./(1 + exp(a.*theta));

% Euler–Maruyama integration
for k = 1:N-1
    x = X(:,k);
    
    %stochastic noise
    dW = sqrt(dt) * randn;  % Wiener increment
    
    Fe = cee*x(1) - (1-beta1)*cie*x(2) + alpha*Kp ;
    Fi = cei*x(1) + (1-alpha)*Kp;
    
    dx1 = -x(1) + (1 - x(1)) * phi(ae, Fe, thetae);
    dx2 = -x(2) + (1 - x(2)) * phi(ai, Fi, thetai);

    X(1,k+1) = x(1) + dx1 * dt+sigma * dW;
    X(2,k+1) = x(2) + dx2 * dt+sigma * dW;
end

end
