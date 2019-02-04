clear

tau = 1;
zeta = 0.5;
a = 1;
f = -1;
V = -(6*zeta*f)/(a*tau);

steps = 100;
x = linspace(-10,10,steps);
t = linspace(0,1,steps);

for n=1:steps
    for i=1:steps
        phi(i,n) = 1/2*(1 - tanh(a*(x(i)-V*t(n))/(2*zeta)));
    end
end

for n=1:steps
    plot(x,phi(:,n))
    pause(0.2)
end