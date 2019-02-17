function P = solve_num(delta,x1,E,n)
% Solves the TISE d^2(P)/dx^2=f(x)P using the numerov solution
% Example run: solve_num(0.05,6.5,3,1)
%Input:
% delta: spacing between grid points for numerical solution
% x1: upper limit for integration
% E = energy eigenvalue (E = 2n+1)
% n = E = 2n+1 eigenstte

% Output: 
% P: wavefunction

%Calculating grid points P values for numerical solution

x = 0:delta:x1; %spaces grid points lineraly in x for numerov calculation
analyticP = zeros(size(x)); %Creates zero vector to put on grid points for analytical solution
P = zeros(x1/delta,1); %Does the same for numerical solution
f = @(x) x.^2-E; %This defines the function f from the schrodinger equation

%Calculating value of P at second grid point for odd and even n
if mod(n,2) == 0 % divisible by 2
    P0 = 1;
    dP0=0; % Inputs boundary conditions
    P(1) = P0;
    P(2) = P0 + P0*(delta^2/2)*f(0) + delta^4/24*(2*P0+f(0)^2);
    % Taylor expansion
end
if mod(n,2) == 1 % not divisible by 2
    P0 = 0;
    dP0=1; 
    P(1) = P0;
    P(2) = dP0*delta + dP0*delta^3/6*f(0)+1/120*delta^5*(6+f(0)^2);
    %Also Taylor expansion
end

% iteration formula to calculate value of P at next grid point
% is just the Taylor expansion of the Schrodinger equation
for j = 2:x1/delta
   P(j+1) = ((2+(5/6)*delta^2*f(x(j)))*P(j)-(1-(delta^2*f(x(j-1)))/12)*P(j-1))/(1-delta^2*f(x(j+1))/12);
end

analyticP = hermiteH(n,x).*exp(-x.^2/2); %built in hermite polynomials 
P = P.*analyticP(3)/P(3);

%plots
plot(x,analyticP); % plots the analytic solution
hold on
plot(x,P); % Plots the numerical solution
title('Quantum harmonic oscillator wavefunction value')

legend('Analytical Solution','Numerical Solution');
xlabel('x')
ylabel('\psi(x)')

end