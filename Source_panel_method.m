
clc;
clear all;
close all;
data  = load('airfoil.txt'); 

v = input('Value of freestream velocity = ');
alpha = input('Value of angle of attack = ');
Xt = flip(data(:,1)); Yt = flip(data(:,2));
n = size(Xt,1); r = zeros(n-1,n-1); 
X = zeros(1,n);  Y = zeros(1,n); 
xm = zeros(1,n-1);  ym = zeros(1,n-1); 
S = zeros(1,n-1);  beta = zeros(1,n-1); 
V_1 = zeros(1,n-1);  V_2 = zeros(1,n-1);
l1 = zeros(n-1,n-1);  In = zeros(n-1,n-1);  It = zeros(n-1,n-1); 
Cp = zeros(1,n-1);  l2 = zeros(n-1,n-1);
x_x = zeros(1,n-1);  y_y = zeros(1,n-1);


for i = 1:n
    X(i) = Xt(i);
    Y(i) = Yt(i);
end

for j = 1:(n-1)
    xm(j) = 0.5*(X(j) + X(j+1));
    ym(j) = 0.5*(Y(j) + Y(j+1));
    S(j) = ((X(j)-X(j+1))^2 + (Y(j)-Y(j+1))^2)^0.5;
    beta(j) = atand((Y(j)-Y(j+1)) / (X(j)-X(j+1))) + 90 - alpha;
    for m = 1:200
           x_x(j,m) = ((m)/200)*(X(j) + X(j+1)); 
           y_y(j,m) = ((m)/200)*(Y(j) + Y(j+1));
    end
end
for t = 1:(n-1)
    V_1(t) = -2*pi*v*cosd(beta(t));
    V_2(t) = v*sind(beta(t));
end
for k = 1:(n-1)
    for m = 1:(n-1)
        for d = 1:200
            l1(k,m) = l1(k,m) + (((log(((xm(k)-x_x(m,d))^2 + (ym(k)-y_y(m,d))^2)^0.5))/...
                (v*(X(m+1) - X(m)/Y(m+1) - Y(m)))));
            l2(k,m) = l2(k,m) + (((log(((xm(k)-x_x(m,d))^2 + (ym(k)-y_y(m,d))^2)^0.5))/...
                (-v/(X(m+1) - X(m)/Y(m+1) - Y(m)))));
        end
    end
end


for p = 1:(n-1)
    for q = 1:(n-1)
        if p==q
            In(p,q) = pi;
            It(p,q) = 0;
        else
            In(p,q) = (l1(p,q))*(S(q)/200);
            It(p,q) = (l2(p,q))*(S(q)/200);
        end
    end
end
lambda = V_1/In;
Vf = V_2 + (lambda*It)./(2*pi);
total = 0;
for f = 1:(n-1)
  Cp(f) = 1 - (Vf(f)/v)^2;
  total = total + lambda(f)*S(f);
end

t = zeros(1,n-1);

for i = 1:n-1
    t(i) = 0.5*( X(n-i) + X(n+1-i) ); 
end


plot(t,Cp);
hold on
set(gca,'Ydir','reverse');
set(gca,'Xdir','reverse');
xlabel('x/c');
ylabel('C_p');