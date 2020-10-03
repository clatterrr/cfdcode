% LBM 1D diffusion equation D1Q3
clear
m = 101;
dx = 1.0;
rho = zeros(m);%这玩意就是温度T
w0 = 4./6.;
w1 = 1./6.;
w2 = w1;
c2 = 1./3.;

alpha = 0.25;
omega = 1/(3 * alpha + 0.5);
f0 = zeros(m);
f1 = zeros(m);
flux = zeros(m);
x = zeros(m);
x(1) = 0.0;

for i = 1:m-1
    x(i+1)=x(i)+dx;
end

twall = 1.0;%左边界温度为1.0
nstep = 200;

for i = 1:m
    f0(i) = 0.5*rho(i);
    f1(i) = 0.5*rho(i);
    f2(i) = 0.5*rho(i);
end

%Collision
for k1 = 1:nstep
    for i = 1:m
        feq0 = w0 * rho(i);
        feq = w1 * rho(i);
        f0(i) = (1-omega)*f1(i) + omega*feq0;
        f1(i) = (1-omega)*f1(i) + omega*feq;
        f2(i) = (1-omega)*f2(i) + omega*feq;
    end

    %Streaming
    for i = 1:m-1
        f1(m-i+1) = f1(m-i);% f1(m) = f1(m-1) , f1(m-1) = f1(m-2) ...  f1(2)   = f1(1)
        f2(i) = f2(i+1);    % f2(1) = f2(3)   , f2(2)   = f2(3)   ...  f2(m-1) = f2(m)
    end
    %Boundary condition
    f1(1) = twall - f2(1) - f0(1);
    f1(m) = f1(m-1);
    f2(m) = f2(m-1);
    f0(m) = f0(m-1);

    for j = 1:m
        rho(j) = f1(j) + f2(j) + f0(j);
    end
end


%Flux
fluxq = zeros(m);
for k = 1:m
    flux(k) = omega*(f1(k) - f2(k))/c2;
end
for k = 1:m 
    fluxq(k) = rho(k)-rho(k+1);
end
fluxq(m) = fluxq(m-1);

figure(1)
plot(x,rho)
   title('Tempture')
   xlabel('X')
   ylabel('T')

figure(2)
plot(x,flux,'o',x,fluxq,'x')
    title('Flux, time step = 200')
    xlabel('X')
    ylabel('Flux')
