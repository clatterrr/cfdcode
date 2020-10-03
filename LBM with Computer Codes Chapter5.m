% LBM 1D diffusion equation D1Q2
clear
m = 101;
dx = 1.0;
rho = zeros(m);%这玩意就是温度T
f1 = zeros(m);
flux = zeros(m);
x = zeros(m);
x(1) = 0.0;

for i = 1:m-1
    x(i+1)=x(i)+dx;
end

alpha = 0.25;
omega = 1/(alpha + 0.5);
twall = 1.0;%左边界温度为1.0
nstep = 2;

for i = 1:m
    f1(i) = 0.5*rho(i);
    f2(i) = 0.5*rho(i);
end

%Collision
for k1 = 1:nstep
    for i = 1:m
        feq = 0.5 * rho(i);
        f1(i) = (1-omega)*f1(i) + omega*feq;
        f2(i) = (1-omega)*f2(i) + omega*feq;
    end

    %Streaming
    for i = 1:m-1
        f1(m-i+1) = f1(m-i);% f1(m) = f1(m-1) , f1(m-1) = f1(m-2) ...  f1(2)   = f1(1)
        f2(i) = f2(i+1);    % f2(1) = f2(3)   , f2(2)   = f2(3)   ...  f2(m-1) = f2(m)
    end
    %Boundary condition
    f1(1) = twall - f2(1);
    f1(m) = f1(m-1);
    f2(m) = f2(m-1);

    for j = 1:m
        rho(j) = f1(j) + f2(j);
    end
end


%Flux
for k = 1:m
    flux(k) = omega*(f1(k) - f2(k));
end

figure(1)
plot(x,rho)
   title('Tempture')
   xlabel('X')
   ylabel('T')

figure(2)
plot(x,flux,'o')
    title('Flux, time step = 200')
    xlabel('X')
    ylabel('Flux')
