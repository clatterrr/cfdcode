% LBM 2D diffusion equation D2Q5
clear
m = 51;
n = 51;
x1 = 1.0;
y1 = 1.0;
dx = x1/(m-1.0);
dy = y1/(n-1.0);
w0 = 2./6.; %中心点的权重
w = 1./6.;   %上下左右的权重

%初始化
f0 = zeros(m,n);f1 = zeros(m,n);
f2 = zeros(m,n);f3 = zeros(m,n);f4 = zeros(m,n);
rho = zeros(m,n);x = zeros(m);y = zeros(n);
Tm = zeros(m); Z = zeros(n,m);

for i = 1:m-1
    x(i+1)=x(i)+dx;
end
for j = 1:n-1
    y(j+1) = y(j)+dy;
end

alpha = 0.25;
omega = 1/(3.*alpha+0.5);
twall = 1.0;%左边界温度为1.0
nstep = 400;

for i = 1:m
    for j = 1:n
    f0(i,j) = w0*rho(i,j);
    f1(i,j) = w*rho(i,j);
    f2(i,j) = w*rho(i,j);
    f3(i,j) = w*rho(i,j);
    f4(i,j) = w*rho(i,j);
    end
end

%Collision
for k1 = 1:nstep
    for j = 1:n
      for i = 1:m
        feq0 = w0 * rho(i,j);
        feq = w * rho(i,j);
        f0(i,j) = (1.-omega)*f0(i,j) + omega*feq0;
        f1(i,j) = (1.-omega)*f1(i,j) + omega*feq;
        f2(i,j) = (1.-omega)*f2(i,j) + omega*feq;
        f3(i,j) = (1.-omega)*f3(i,j) + omega*feq;
        f4(i,j) = (1.-omega)*f4(i,j) + omega*feq;
       end
    end

    %Streaming
    for j = 1:n
      for i = 1:m-1
        f1(m-i+1,j) = f1(m-i,j);% f1(m) = f1(m-1) , f1(m-1) = f1(m-2) ...  f1(2)   = f1(1)
        f2(i,j)     = f2(i+1,j);    % f2(1) = f2(3)   , f2(2)   = f2(3)   ...  f2(m-1) = f2(m)
      end
    end

    for i = 1:m
        for j = 1:n-1
          f3(i,n-j+1) = f3(i,n-j);% f1(m) = f1(m-1) , f1(m-1) = f1(m-2) ...  f1(2)   = f1(1)
          f4(i,j)     = f4(i,j+1);    % f2(1) = f2(3)   , f2(2)   = f2(3)   ...  f2(m-1) = f2(m)
        end
      end
    %Boundary condition
    for j = 1:n
        f1(1,j) = twall - f2(1,j) - f0(1,j) - f3(1,j) - f4(1,j);
        f3(m,j) = -f1(m,j) - f0(m,j) - f2(m,j) - f4(m,j);
    end

    for i = 1:m
        f3(i,1) = f3(i,2);
        f4(i,n) = -f0(i,n) - f3(i,n) - f2(i,n) - f1(i,n);
    end


    for j = 1:n
       for i = 1:m
        rho(i,j) = f1(i,j) + f2(i,j) + f0(i,j) + f3(i,j) + f4(i,j);
       end
    end
end


%rotating matrix for contour plotting
for j = 1:n
    for i = 1:m
        Z(j,i) = rho(i,j);
    end
end

for i = 1:n
    Tm(i) = rho(i,(n-1)/2);
end

figure(1)
plot(x,Tm)
   xlabel('X')
   ylabel('T')
figure(2)
contour(Z)
