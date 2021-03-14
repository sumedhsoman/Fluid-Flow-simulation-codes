% SIMPLE Algorithm solver. This code is designed to solve Steady-State
% Navier Stokes equations and Continuity equations, in 2D, grid size 10x10
%% Author: Sumedh Soman
%% V1.2
%% Variable Definition and preallocation, grid generation
mu = 1;
rho = 1;
dx = 0.1;
dy = 0.1;
Nx = 10;
Ny = 10;
%Pressure Underrelaxation = 0.7, Velocity underrelaxation = 0.3
alpha_u = 0.9;
alpha_v = 0.9;
alpha_p = 0.5;
u_northwall = 5;
u_southwall = 5;
u_westwall = 5;
u_eastwall = 5;
v_northwall = 0;
v_southwall = 0;
v_westwall = 0;
v_eastwall = 0;
ConvergenceCondition = 0.05;
% Preallocation process
uold = zeros(Nx,Ny);
vold = zeros(Nx,Ny);
uold(1,2:2:10) = 5;
vold(10,1:2:9) = 0;
vold(2:2:8,1)= 0;
uold(3:2:9,10) = 0;
u = zeros(Nx,Ny);
v = zeros(Nx,Ny);
p = zeros(Nx,Ny);
u(1,2:2:10) = 5;
v(10,1:2:9) = 0;
v(2:2:8,1)= 0;
u(3:2:9,10) = 0;
p_prime = zeros(Nx,Ny);
a_u = ones(10,10);
a_v = ones(10,10);
a_p_prime = ones(10,10);
an_p_prime = zeros(10,10);
as_p_prime = zeros(10,10);
aw_p_prime = zeros(10,10);
ae_p_prime = zeros(10,10);
F_u = zeros(10,10);
F_v = zeros(10,10); 
ae_u = zeros(10,10);
aw_u = zeros(10,10);
an_u = zeros(10,10);
as_u = zeros(10,10);
ae_v = zeros(10,10);
as_v = zeros(10,10);
an_v = zeros(10,10);
aw_v = zeros(10,10);
bdash = zeros(10,10);
uold = u;
vold = v;
%% Solver implementation
% Coefficient Definition; 
for var = 1:2000
A = dx*dy;
for i = 3:2:7
    for j = 4:2:8
      F_u(i,j-1) = 0.5*rho*(uold(i,j)+uold(i,j-2));
      F_u(i,j+1) = 0.5*rho*(uold(i,j)+uold(i,j+2));
      F_u(i+1,j) = 0.5*rho*(vold(i-1,j+1)+vold(i-1,j-1));
      F_u(i-1,j) = 0.5*rho*(vold(i+1,j+1)+vold(i+1,j-1));
      F_u(i,1) = 0.5*rho*(uold(i,2));
      F_u(i+1,2) = 0.5*rho*(vold(i-1,3)+vold(i-1,1));
      F_u(i-1,2) = 0.5*rho*(vold(i+1,3)+vold(i+1,1));
      F_u(10,j) = 0.5*rho*(vold(8,j+1)+vold(8,j-1));
      F_u(9,j-1) = 0.5*rho*(uold(9,j)+uold(9,j-2));
      F_u(9,j+1) = 0.5*rho*(uold(9,j)+uold(9,j+2));
      F_u(9,1) = 0.5*rho*(uold(9,2));
      F_u(10,2) = 0.5*rho*(vold(8,3)+vold(8,1));
    end
      
end

for j = 3:2:7
    for i = 4:2:8
      F_v(i,j-1) = 0.5*rho*(uold(i-1,j-1)+uold(i+1,j-1));
      F_v(i,j+1) = 0.5*rho*(uold(i-1,j+1)+uold(i+1,j+1));
      F_v(i-1,j) = 0.5*rho*(vold(i,j)+vold(i-2,j));
      F_v(i+1,j) = 0.5*rho*(vold(i,j)+vold(i+2,j));
      F_v(i,10) = 0.5*rho*(uold(i-1,10)+uold(i+1,10));
      F_v(i-1,9) = 0.5*rho*(vold(i,9)+vold(i-2,9));
      F_v(i+1,9) = 0.5*rho*(vold(i,9)+vold(i+2,9));
      F_v(2,j-1) = 0.5*rho*(uold(1,j-1)+uold(3,j-1));
      F_v(2,j+1) = 0.5*rho*(uold(1,j+1)+uold(3,j+1));
      F_v(1,j) = 0.5*rho*(vold(2,j)+5);
      F_v(2,10) = 0.5*rho*(uold(1,10)+uold(3,10));
      F_v(1,9) = 0.5*rho*(vold(2,9)+5);
      
    end

end
D = mu/dx;
% Defining coefficients by upwind
for i = 3:2:9
    for j = 2:2:8
        aw_u(i,j)= max(F_u(i,j-1),0)+ mu/dx;
        ae_u(i,j)= max(-F_u(i,j+1), 0)+mu/dx;
        an_u(i,j)= max(-F_u(i-1,j), 0)+mu/dx;
        as_u(i,j)= max(F_u(i+1,j), 0)+mu/dx;
        
        a_u(i,j) =  as_u(i,j)+an_u(i,j)+ae_u(i,j)+aw_u(i,j)/alpha_u;
        
     end
end

for i = 2:2:8
    for j = 3:2:9
        aw_v(i,j)= max(F_v(i,j-1), 0)+mu/dx;
        ae_v(i,j)= max(-F_v(i,j+1), 0)+mu/dx;
        an_v(i,j)= max(-F_v(i-1,j), 0)+mu/dx;
        as_v(i,j)= max(F_v(i+1,j), 0)+mu/dx;
                
        a_v(i,j) = (as_v(i,j)+an_v(i,j)+ae_v(i,j)+aw_v(i,j))/alpha_v;
    end
end

% Solving u momentum equations
for viter  = 1:100
for i = 3:2:7
    for j = 4:2:8
        u(i,j) = (1-alpha_u)*uold(i,j)+(1/a_u(i,j))*(an_u(i,j)*u(i-2,j)+as_u(i,j)*u(i+2,j)+ae_u(i,j)*u(i,j+2)+aw_u(i,j)*u(i,j-2)+...
            (p(i,j-1)-p(i,j+1))*A);
      u(i,2) = (1-alpha_u)*uold(i,2)+(1/a_u(i,2))*(an_u(i,2)*u(i-2,2)+as_u(i,2)*u(i+2,2)+ae_u(i,2)*u(i,4)+aw_u(i,2)*0+...
            (p(i,1)-p(i,3)*A));
        u(9,j) = (1-alpha_u)*uold(9,j)+(1/a_u(9,j))*(an_u(9,j)*u(7,j)+as_u(9,j)*0+ae_u(9,j)*u(9,j+2)+aw_u(9,j)*u(9,j-2)+...
            (p(9,j-1)-p(9,j+1))*A);
        u(9,2) = (1-alpha_u)*uold(9,2)+(1/a_u(9,2))*(an_u(9,2)*u(7,2)+as_u(9,2)*0+ae_u(9,2)*u(9,4)+aw_u(9,2)*0+...
            (p(9,1)-p(9,3))*A);

    end
end
end

%Solving the V momentum equation
for uiter = 1:100
for i = 4:2:8
    for j = 3:2:7
        v(i,j) = (1-alpha_v)*vold(i,j)+(1/a_v(i,j))*(as_v(i,j)*v(i-2,j)+an_v(i,j)*v(i+2,j)+ae_v(i,j)*v(i,j+2)+aw_v(i,j)*v(i,j-2)+...
            (p(i-1,j)-p(i+1,j))*A);
        v(2,j) = (1-alpha_v)*vold(2,j)+(1/a_v(2,j))*(as_v(2,j)*0+an_v(2,j)*v(4,j)+ae_v(2,j)*v(2,j+2)+aw_v(2,j)*v(2,j-2)+...
            (p(1,j)-p(3,j))*A);
        v(i,9) = (1-alpha_v)*vold(i,9)+(1/a_v(i,9))*(as_v(i,9)*v(i-2,9)+an_v(i,9)*v(i+2,9)+ae_v(i,9)*0+aw_v(i,9)*v(i,7)+...
            (p(i-1,9)-p(i+1,9))*A);
        v(2,9) = (1-alpha_v)*vold(2,9)+(1/a_v(2,9))*(as_v(2,9)*0+an_v(2,9)*v(4,9)+ae_v(2,9)*0+aw_v(2,9)*v(2,7)+...
            (p(1,9)-p(3,9))*A);
        
    end
end
end
%Pressure Correction Equation

for i = 3:2:9
    for j= 3:2:9
        as_p_prime(i,j) = A/a_v(i-1,j);
        an_p_prime(i,j) = A/a_v(i+1,j);
        aw_p_prime(i,j) = A/a_u(i,j-1);
        ae_p_prime(i,j) = A/a_u(i,j+1);
        a_p_prime(i,j) = as_p_prime(i,j)+ aw_p_prime(i,j)+ae_p_prime(i,j)+an_p_prime(i,j);
        bdash(i,j) = A*(u(i,j-1)-u(i,j+1))+(v(i-1,j)-v(i+1,j));
        a_p_prime(1,1) = 1e30;
    end
end



for piter = 1:1100
for i = 3:2:7
    for j= 3:2:7
          p_prime(i,j) = p_prime(i,j) + (1.7/(a_p_prime(i,j)))*(ae_p_prime(i,j)*p_prime(i,j+2)+an_p_prime(i,j)*p_prime(i-2,j)+aw_p_prime(i,j)*p_prime(i,j-2)+as_p_prime(i,j)*p_prime(i+2,j)-...
        a_p_prime(i,j)*p_prime(i,j)+bdash(i,j));
          
          p_prime(9,j) = p_prime(9,j) + (1.7/(a_p_prime(9,j)))*(ae_p_prime(9,j)*p_prime(9,j+2)+an_p_prime(9,j)*p_prime(7,j)+aw_p_prime(9,j)*p_prime(9,j-2)+as_p_prime(9,j)*0-...
        a_p_prime(9,j)*p_prime(9,j)+bdash(9,j));
        
        
          p_prime(i,9) = p_prime(i,9) + (1.7/(a_p_prime(i,9)))*(ae_p_prime(i,9)*0+an_p_prime(i,9)*p_prime(i-2,9)+aw_p_prime(i,9)*p_prime(i,7)+as_p_prime(i,9)*p_prime(i+2,9)-...
        a_p_prime(i,9)*p_prime(i,9)+bdash(i,9)); 
        
        
          p_prime(9,9) = p_prime(9,9) + (1.7/(a_p_prime(9,9)))*(ae_p_prime(9,9)*0+an_p_prime(9,9)*p_prime(7,9)+aw_p_prime(9,9)*p_prime(9,7)+as_p_prime(9,9)*0-...
        a_p_prime(9,9)*p_prime(9,9)+bdash(9,9));
                
    end
end
end
        
  

for i = 3:2:9
    for j = 2:2:8
        u(i,j) = u(i,j)+ (A/a_u(i,j))*(p_prime(i,j-1)-p_prime(i,j+1));
    end
end


for i = 2:2:8
    for j = 3:2:9
        v(i,j) = v(i,j)+ (A/a_v(i,j))*(p_prime(i,j-1)-p_prime(i,j+1));
    end
end


for i = 3:2:9
    for j= 3:2:9
        p(i,j) = p(i,j) + alpha_p*p_prime(i,j);
    end
end

uold = u;
vold = v;
% Continuity Condition:
for i = 3:2:9
    for j= 3:2:9
     Cont =   A*((-u(i,j-1)+u(i,j+1))+(-v(i-1,j)+v(i+1,j)));
     fprintf('%f\n',Cont)
     if Cont < ConvergenceCondition && Cont > 0
         break
     end
    end
end
end
figure
x = 1:2:9;
y = 1:2:9;
Z = p(1:2:9,1:2:9);
contourf(x,y,Z);
figure
y1 = 1:2:9;
x1 = 2:2:10;
Z2 = u(y1,x1);
contourf(x1,y1,Z2);
figure 
Z3 = v(x1,y1);
contourf(x1,y1,Z3);
colormap jet;
