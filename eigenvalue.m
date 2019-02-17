function E = oscillator(E,n)
%Looks for energy eigenvalue close to input energy
%Input:
% E0: Input energy guess

%Output:
% E: Actual energy eigenvalue close to input

% E0 has to be above zero and real
 
x1=5; %upper limit
delta = 0.05;
P = solve_num(delta,x1,E,n); %starts the numerical function
a = abs(P(101)); %absolute value of P with E0
i=0.025; %incriment which the eigenvalue is changed by
e = 0.005; %small absolute value of P 

E=E-i; %Looks at how P differs 
P = solve_num(delta,x1,E,n); %recalculates the solution with new value of E
if abs(P(101)) < a %investigates if subtracting i moves the function closer to the minimum
    i=-i;
else %if taking away i takes P further away from the root
    i=i; %this means i is added instead of subtracted
end
   while abs(P(101)) > e %find E when P is very small(~minimum)
        E=E+i %Repeatedly removes or adds e from E 
        P = solve_num(delta,x1,E,n); %Calculates P again
    end
end