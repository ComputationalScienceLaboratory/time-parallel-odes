clear;
clc;
global ntracer;
global nx;
global ny;
global nz;
%set test values
ntracer = 1;
nx = 100;
ny = 100;
nz = 1;

% % test calls(build bump and blow it sideways)
% u =  zeros(nx, ny, nz); %wind
% u(:, 1, 1) = 1*ones(nx, 1);
% % u(1, :, 1) = ones(1, ny);
% % %u = ones(nx, ny, nz);
% kh = zeros(nx, ny, nz); %k
% % % %kh = 1000*ones(nx, ny, nz);
% concentration = zeros(nx, ny, nz, ntracer);
% x = zeros(100, 1);
% x(25:50) = 1;
% concentration(:, 1, 1, 1) = x;
% concentration_bc = zeros(2, ny, nz, ntracer); %first component is w second e
% % % %concentration_bc(:, 1, 1, 1) = zeros(2, 1);
% cdot = zeros(nx, ntracer);
% % %test calls
% transport_X_fwd(10000000000, u, kh, concentration, concentration_bc, cdot);
% %in y direction
% % concetrantion2 = zeros(nx, ny, nz, ntracer);
% % y = zeros(ny, 1);
% % y(25:50) = 1;
% % concentration(1, :, 1, 1) = y;
% % %transport_Y_fwd(100000000000, u, kh, concentration, concentration_bc, cdot);


% %%%%%%%%transport_Z_fwd test block
% w = zeros(nx, ny, nz);
% w(1, 1, :) = .1*ones(nz, 1);
% %w(1, 1, 25:50) = 1;
% kv = zeros(nx, ny, nz);
% %kv(1, 1, :) = 10*ones(nz, 1);
% concentration = zeros(nx, ny, nz, ntracer);
% %put a bump in the initial concentration
% x = zeros(nz, 1);
% x(25:50) = 1;
% concentration(1, 1, :, 1) = x;
% concentration_bc = zeros(nx, ny, ntracer);
% surfaceem = zeros(nx, ny, ntracer);
% %surfaceem(:, :, 1) = ones(nx, ny);
% volumeem = zeros(nx, ny, nz, ntracer);
% %volumeem(1, 1, :, 1) = ones(nz, 1);
% depositionvel = zeros(nx, ny, ntracer);
% transport_Z_fwd(100000, w, kv, concentration, concentration_bc, surfaceem, volumeem, depositionvel);

% %%%%%%%adj_advdiff_fdh test block
% %%%%%%%the purpose of this test block is to test that adj_advdiff_fdh
% %%%%%%%matches the numerical adjoint of advdiff_fdh
% %%%%%%%set test values
% %%%%%%%verification was successful at time of writing
% ntracer = 1;
% u =  zeros(nx, ny, nz); %wind
% %u(:, 1, 1) = ones(nx, 1);
% kh = zeros(nx, ny, nz); %k
% %kh = 1000*ones(nx, ny, nz);
% concentration = zeros(nx, ny, nz, ntracer);
% x = zeros(nx, 1);
% %x(floor(nx/4):floor(nx/2)) = 1;
% concentration(:, 1, 1, 1) = x;
% concentration_bc = zeros(2, ny, nz, ntracer); %first component is w second e
% sb(1, :) = squeeze(concentration_bc(1, 1, 1, 1:ntracer));
% sb(2, :) = squeeze(concentration_bc(2, 1, 1, 1:ntracer));
% cdot = zeros(nx, ntracer);
% dt = 15;
% nstep = floor(dt/300) + 1;
% dx = 1000;
% conc = squeeze(concentration(1:nx, 1, 1, 1:ntracer));
% lam = ones(size(conc));
% %%%%%%%first test that adj_advdiff_fdh runs without error
% adj_advdiff_fdh(dt, nstep, nx, ntracer, dx, u, kh, lam);
% %%%%%%%verify with finite differences
% %%%%%%%get numerical jacobian
% n = nx;
% A = eye(n); %grab standard basis vectors from here
% advdiff_fdh(dt, nstep, nx, ntracer, dx, u, kh, conc, sb);
% h = 1e-6;
% numjac = zeros(n, n);
% f = @(x) advdiff_fdh(dt, nstep, nx, ntracer, dx, u, kh, x, sb);
% x0 = conc;
% i = 1;
% for j = 1:n
%     e = A(:,j);
%     stepf = x0 + h*e;
%     stepb = x0 - h*e;
%     fd = (f(stepf) - f(stepb))/(2*h);
%     numjac(:, j) = fd;
% end
% adj = zeros(n, n); %full adjoint matrix as recovered from adjv products 
% adjv = @(x) adj_advdiff_fdh(dt, nstep, nx, ntracer, dx, u, kh, x);
% for j = 1:n    
%     e = A(:,j);
%     adj(:, j) = adjv(e);
% end
% disp("Relative error (Frobenius)");
% disp(norm(numjac' - adj, 'fro')/norm(numjac', 'fro'));

% %%%%%%%transport_X_fwd test block
% %%%%%%%the purpose of this test block is to test that transport_X_dadj
% %%%%%%%matches the numerical adjoint of transport_X_fwd
% %%%%%%%this block is still in progress
% %%%%%%%set test values
% ntracer = 1;
% u =  zeros(nx, ny, nz); %wind
% %u(:, 1, 1) = ones(nx, 1);
% kh = zeros(nx, ny, nz); %k
% %kh = 1000*ones(nx, ny, nz);
% concentration = zeros(nx, ny, nz, ntracer);
% x = zeros(nx, 1);
% x(floor(nx/4):floor(nx/2)) = 1;
% concentration(:, 1, 1, 1) = x;
% concentration_bc = zeros(2, ny, nz, ntracer); %first component is w second e
% sb(1, :) = squeeze(concentration_bc(1, 1, 1, 1:ntracer));
% sb(2, :) = squeeze(concentration_bc(2, 1, 1, 1:ntracer));
% cdot = zeros(nx, ntracer);
% dt = 15;
% nstep = floor(dt/300) + 1;
% dx = 1000;
% lam = ones(size(concentration));
% %%%%%%%first test that Transport_X_adj runs without error
% %adj_advdiff_fdh(dt, nstep, nx, ntracer, dx, u, kh, lam);
% %%%%%%%verify with finite differences
% %%%%%%%get numerical jacobian
% n = nx;
% A = eye(n); %grab standard basis vectors from here
% transport_X_fwd(5, u, kh, concentration, concentration_bc, cdot);
% h = 1e-6;
% numjac = zeros(n, n);
% f = @(x) transport_X_fwd(5, u, kh, x, concentration_bc, cdot);
% x0 = concentration;
% i = 1;
% for j = 1:n
%     e = A(:,j);
%     stepf = x0 + h*e;
%     stepb = x0 - h*e;
%     fd = (f(stepf) - f(stepb))/(2*h);
%     numjac(:, j) = fd;
% end
% adj = zeros(n, n); %full adjoint matrix as recovered from adjv products 
% adjv = @(x) adj_advdiff_fdh(dt, nstep, nx, ntracer, dx, u, kh, x);
% for j = 1:n    
%     e = A(:,j);
%     adj(:, j) = adjv(e);
% end
% disp("Relative error (Frobenius)");
% disp(norm(numjac' - adj, 'fro')/norm(numjac', 'fro'));



% %%%%%%%%%%%test that I can call transport_X_dadj properly
% dt = 5;
% u = zeros(nx, ny, nz);
% u(:, 1, 1) = 1*ones(nx, 1);sq
% % u(1, :, 1) = ones(1, ny);
% % %u = ones(nx, ny, nz);
% kh = zeros(nx, ny, nz); %k
% % % %kh = 1000*ones(nx, ny, nz);
% concentration = zeros(nx, ny, nz, ntracer);
% x = zeros(100, 1);
% x(25:50) = 1;
% concentration(:, 1, 1, 1) = x;
% concentration_bc = zeros(2, ny, nz, ntracer); %first component is w second e
% % % %concentration_bc(:, 1, 1, 1) = zeros(2, 1);
% cdot = zeros(nx, ntracer);
% % %test calls
% transport_X_fwd(dt, u, kh, concentration, concentration_bc, cdot);
% lambda = ones(size(concentration));
% trasport_X_dadj(dt, u, kh, concentration, lambda);

%%%%%My goal here is to successfully use the transport_fwd functions to
%%%%%perform the 2d test problem
%%%%%I want to build a circle of some radius where everything is 1
%%%%%everything 0 elsewhere
concentration = zeros(nx, ny, nz, ntracer);
cntr = [20, 20];
rad = 3;
for i = 1:nx
    for j = 1:ny
        pt = [i, j];
        if(norm(pt - cntr) <= rad)
           concentration(i, j, 1, 1) = 2; 
        end
    end
end
clf;
%%%%%I want to build the wind profiles now
omega = 2;
windcntr = [50, 50];
u = zeros(nx, ny, nz, ntracer);
v = zeros(nx, ny, nz, ntracer);
for i = 1:nx
    for j = 1:ny
        %u(i, j, 1, 1) = -omega*(j - windcntr(2));
        %v(i, j, 1, 1) = omega*(i - windcntr(1));
        %ihat = (i - windcntr(1))/nx;
        %jhat = (j - windcntr(2))/ny;
        ihat = i;
        jhat = j;
        u(i, j, 1, 1) = 1;
        v(i, j, 1, 1) = 0;
    end
end

quiver(u,v)

dt = 100;
kh = zeros(nx, ny, nz); %k
% % %kh = 1000*ones(nx, ny, nz);
concentration_bc = zeros(2, ny, nz, ntracer); %first component is w second e
% % %concentration_bc(:, 1, 1, 1) = zeros(2, 1);
cdot = zeros(nx, ntracer);
i = 0;
% %test calls
contour(concentration);
pause
while true
quiver(u, v)
hold on
contour(concentration);
pause;
hold off
concentration1 = transport_X_fwd(dt, u, kh, concentration, concentration_bc, cdot);
concentration2 = transport_Y_fwd(dt, v, kh, concentration, concentration_bc, cdot);
concentration = (concentration1 + concentration2)/2;
end

function concentration = transport_Z_fwd(dt, w, kv, concentration, concentration_bc, surfem, volem, depositionvel)
global ntracer;
global nx;
global ny;
global nz;

nsteps = floor(dt/300)+1;
%nsteps = 1000;
deltat = dt/nsteps;

deltaz = 100;
%generate uniformly spaced grid
%the vertical grid
z = linspace(0, deltaz*nz, nz);

for i = 1:nx
    for j = 1:ny
        conc(1:nz, 1:ntracer) = concentration(i, j, 1:nz, 1:ntracer);
        wind(1:nz) = w(i, j, 1:nz);
        dif(1:nz) = kv(i, j, 1:nz);
        sb(1, 1:ntracer) = concentration(i, j, 1, 1:ntracer);
        sb(2, 1:ntracer) = concentration_bc(i, j, 1:ntracer);
        surfaceem(1:ntracer) = surfem(i, j, 1:ntracer);
        volumeem(1:nz, 1:ntracer) = 0;
        depvel(1:ntracer) = depositionvel(i, j, 1:ntracer);
        
        c = advdiff_fdv(deltat, nsteps, nz, ntracer, z, wind, dif, conc, sb, surfaceem, volumeem, depvel);
        
        concentration(i, j, 1:nz, 1:ntracer) = subplus(conc); %positive part
        
    end
end

end

function concentration = transport_Y_fwd(dt, v, kh, concentration, concentration_bc, cdot)
global ntracer;
global nx;
global ny;
global nz;

dy = 100;
deltay = dy;
nsteps = floor(dt/300) + 1;
%nsteps = 1e10;
deltat = dt/nsteps;

for i = 1:nx
    for k = 1:nz
        %disp("next");
        conc = squeeze(concentration(i, 1:ny, k, 1:ntracer));
        wind = v(1, 1:ny, k)';
        dif = kh(i, 1:ny, k)';
        sb(1, :) = squeeze(concentration_bc(1, i, k, 1:ntracer));
        sb(2, :) = squeeze(concentration_bc(2, i, k, 1:ntracer));
        conc = advdiff_fdh(deltat, nsteps, ny, ntracer, deltay, ...
            wind, dif, conc, sb);
        
        concentration(i, 1:ny, k, 1:ntracer) = subplus(conc); %positive part
    end
end

end


function concentration = transport_X_fwd(dt, u, kh, concentration, concentration_bc,  cdot)
global ntracer;
global nx;
global ny;
global nz;

dx = 100;
deltax = dx;
nsteps = floor(dt/300)+1;
deltat = dt/nsteps;

for j = 1:ny
    for k = 1:nz
        conc = squeeze(concentration(1:nx, j, k, 1:ntracer))';
        wind = u(1:nx, j, k);
        dif = kh(1:nx, j, k);
        sb(1, :) = squeeze(concentration_bc(1, j, k, 1:ntracer));
        sb(2, :) = squeeze(concentration_bc(2, j, k, 1:ntracer));
        for l = 1:nsteps
            %plot(1:100, conc(:, 1)) %first tracer
            %axis([0, nx, 0, 5])
            %pause
            %cdot = advdiff_fun_fdh(nx, deltax, wind, dif, sb, conc, cdot);
            %c1 = conc + deltat*cdot;
            %cdot = advdiff_fun_fdh(nx, deltax, wind, dif, sb, c1, cdot);
            %conc = conc/2 + (c1 + deltat*cdot)/2;
            %wrap around to right
            conc = advdiff_fdh(deltat, nsteps, nx, ntracer, deltax, ...
            wind, dif, conc, sb);

        end
        concentration(1:nx, j, k, 1:ntracer) = subplus(conc)'; %positive part
        %disp("next y, z")
    end
end

end


function dc = advdiff_fun_fdh(n, dx, u, k, bdry, c, dc)
global ntracer;
global nx;
global ny;
global nz;

ap = -1/6;
bp = 1;
cp = - 1/2;
dp = -1/3;
an = 1/3;
bn = 1/2;
cn = -1;
dn = 1/6;

for j = 1:ntracer
    
    F = c(:, j).*u;
    F = [bdry(1, j)*u(1); F];
    F = [F; bdry(2, j)*u(nx)];
    
    %advection discretization
    if u(1) >= 0 %inflow
        dc(1, j) = (F(1) - F(2))/dx;
    else %outflow
        dc(1, j) = (F(2) - F(3))/dx;
    end
    
    for i = 3:n %array indices for F are shifted by 1 in MATLAB
        if u(i) >= 0 %inflow
            dc(i - 1, j) = (ap*F(i-2) + bp*F(i - 1) + cp*F(i) + dp*F(i+1))/dx;
        else
            dc(i - 1, j) = (an*F(i - 1) + bn*F(i) + cn*F(i+1)+dn*F(i+2))/dx;
        end
    end
    
    if u(n) >= 0 %outflow
        dc(n, j) = (F(n) - F(n + 1))/dx;
    else
        dc(n, j) = ((F(n + 1) - F(n + 2)))/dx;
    end
    
    %the diffusion part
    if u(1) >= 0 %inflow
        dc(1, j) = dc(1, j) + ( (k(1) + k(2)) * (c(2, j) - c(1, j)) ...
            - 2*k(1)*(c(1, j) - bdry(1, j)))/(2*dx^2);
    else %outflow
        dc(1, j) = dc(1, j) + ( (k(1) + k(2))*(c(2, j) - c(1, j)))/(2*dx^2);
    end
    
    for i=2:n-1
        dc(i, j) = dc(i, j) + ((k(i+1) + k(i))*(c(i+1, j) - c(i, j)) ...
            - (k(i) + k(i-1))*(c(i, j) - c(i-1, j)))/(2*dx^2);
    end
    
    if u(n) >= 0 %outflow
        dc(n, j) = dc(n, j) + ( ...
            - (k(n) + k(n-1))*(c(n, j)-c(n-1, j)) ...
            )/(2*dx^2);
    else %inflow
        dc(n, j) = dc(n, j) + ( 2*k(n)*(bdry(2, j) - c(n, j)) ...
            - (k(n) - k(n-1))*(c(n, j) - c(n - 1, j)) ...
            )/(2*dx^2);
    end
    
end %j (1:ntracer)
return;
end


function jac = advdiff_jac_fdh(n, dx, u, k)
ku = 2;
kl = 2;

ap = -1/6;
bp = 1;
cp = -1/2;
dp = -1/3;
an = 1/3;
bn = 1/2;
cn = -1;
dn = 1/6;

%5 rows of banded Jacobian
%used to build sparse pentadiagonal jacobian below
jac = zeros(5, n);

%advection discretization
if u(1) >= 0
    jac(jrow(1, 1), 1) = -u(1)/dx;
else
    jac(jrow(1, 1), 1) = u(1)/dx;
    jac(jrow(1, 2), 2) = -u(2)/dx;
end

for i=2:n-1
    if u(i) >= 0
        if i > 2
            jac(jrow(i, i-2), i-2) = ap*u(i-2)/dx;
        end
        jac(jrow(i, i-1), i-1) = bp*u(i-1)/dx;
        jac(jrow(i, i), i) = cp*u(i)/dx;
        jac(jrow(i, i+1), i+1) = dp*u(i+1)/dx;
    else
        jac(jrow(i, i-1), i-1) = an*u(i-1)/dx;
        jac(jrow(i, i), i) = bn*u(i)/dx;
        jac(jrow(i, i+1), i+1) = cn*u(i+1)/dx;
        if i < n - 1
            jac(jrow(i, i+2), i+2) = dn*u(i+2)/dx;
        end
    end
end

if u(n) >= 0
    jac(jrow(n, n-1), n-1) = u(n-1)/dx;
    jac(jrow(n, n), n) = -u(n)/dx;
else
    jac(jrow(n, n), n) = u(n)/dx;
end


%the diffusion part
if u(1) >= 0 %inflow
    jac(jrow(1, 1), 1) = jac(jrow(1, 1), 1) ...
        - (k(2) + 3*k(1))/(2*dx^2);
    jac(jrow(1, 2), 2) = jac(jrow(1, 2), 2) ...
        + (k(1) + k(2))/(2*dx^2);
else %outflow
    jac(jrow(1, 1), 1) = jac(jrow(1, 1), 1) - (k(1) + k(2))/(2*dx^2);
    jac(jrow(1, 2), 2) = jac(jrow(1, 2), 2) + (k(1) + k(2))/(2*dx^2);
end

for i = 2:n-1
    jac(jrow(i, i-1), i-1) = jac(jrow(i, i-1), i-1) ...
        + (k(i) + k(i - 1))/(2*dx^2);
    jac(jrow(i, i), i) = jac(jrow(i, i), i) ...
        - (k(i+1) + 2*k(i) + k(i-1))/(2*dx^2);
    jac(jrow(i, i+1), i+1) = jac(jrow(i, i+1), i+1) ...
        + (k(i+1) + k(i))/(2*dx^2);
end

if u(n) >= 0 %outflow
    jac(jrow(n, n-1), n-1) = jac(jrow(n, n-1), n-1) ...
        + (k(n) + k(n-1))/(2*dx^2);
    jac(jrow(n, n), n) = jac(jrow(n, n), n) ...
        - (k(n) + k(n-1))/(2*dx^2);
else %inflow
    jac(jrow(n, n-1), n-1) = jac(jrow(n, n-1), n-1)...
        + (k(n) + k(n-1))/(2*dx^2);
    jac(jrow(n, n), n) = jac(jrow(n, n), n) ...
        - (3*k(n) + k(n-1))/(2*dx^2);
end

jac = spdiags(jac', 2:-1:-2, n, n);
return;
end

function c = jrow(i, j)
kl = 2;
ku = 2;
if i <= 0 || j <= 0
    disp("error in advdiff_jac_fdh")
    return
end
c = ku + 1 + i - j;
end

function c = advdiff_fdh(dt, nstep, n, nspec, dx, u, k, c, bdry)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Performs Nstep timesteps of length DT
%      to solve the adv_diff equation using a
%      linearly implicit Crank-Nicholson time discretization
%
%  n   = no. of grid points
%  nspec = no. of chemical species
%  nstep = no of time steps
%  x(1:N) = grid point coordinates
%  u(1:N) = wind speeds
%  k(1:N) = diffusion coefficients
%  surfaceem  = Surface Emission intensity
%  volumeem   = Elevated Emission intensity
%  vd    = deposition velocity
%  c     = concentration of each species
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kl = 2;
ku = 2;
ldjac = kl+ku+1;
lda = 2*kl+ku+1;
jac = zeros(n, n);
A = zeros(n, n);
ctmp = zeros(n, 1);
b = zeros(n, 1);
d = zeros(2, 1);

jac = advdiff_jac_fdh(n, dx, u, k);
A = -(dt/2)*jac;
A = A + eye(n);

for ispec=1:nspec
    for istep=1:nstep
        ctmp = c(ispec, :)';
        %plot(ctmp)
        %pause;
        d = bdry(:, ispec);
        %the free term
        b = advdiff_free_fdh(n, dx, u, k, d);
        ctmp = ctmp + dt*b;
        alpha = (dt/2);
        ctmp = alpha*jac*c(ispec, :)' + ctmp;
        ctmp = A\ctmp;
        c(ispec, :) = ctmp;
        %wrap around
        %c(2, :) = c(n, :);
        %c(1, :) = c(n-1,:);
    end
    %disp("next step");
end
end

function b = advdiff_free_fdh(n, dx, u, k, bdry)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Advection-diffusion derivative function
% Free term B such that: c' = Fun_fdh = Jac_fdh*c + B
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ap = -1/6;
bp = 1;
cp = -1/2;
dp = -1/3;
an = 1/3;
bn = 1/2;
cn = -1;
dn = 1/6;

b = zeros(n, 1);

%the advection discretization
if u(1) >= 0
    b(1) = bdry(1)*u(1)/dx + 2*k(1)*bdry(1)/(2*dx^2);
end
if u(2) >= 0
    b(2) = ap*bdry(1)*u(1)/dx;
end
if u(n-1) < 0
    b(n-1) = dn*bdry(2)*u(n)/dx;
end
if u(n) < 0
    b(n) = -bdry(2)*u(n)/dx + 2*k(n)* bdry(2)/(2*dx^2);
end

end

function jac = advdiff_jac_fdv(n, z, w, k)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Advection-diffusion derivative by finite volumes
% Advection is discretized by simple upwind
% Jac is in Blas-banded-format before building
% with spdiags: Jac(1:3,j) = A(j-1:j+1,j)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jac = zeros(n, n);
% difflux/advflux = diffusive/advective fluxes through i-1/2
% flux = total flux through i-1/2
dflux = zeros(n+1, 1);
lflux = zeros(n+1, 1);

dflux(1) = 0; % VD - W(1)
lflux(1) = 0;

%Intermediate boundaries
for i = 2:n
    if w(i) >= 0
        dflux(i) = (k(i) + k(i-1))/2/(z(i) - z(i-1));
        lflux(i) = -(k(i) - k(i-1))/2/(z(i) - z(i-1)) - w(i-1);
    else
        dflux(i) = (k(i) + k(i-1))/2/(z(i) - z(i-1)) - w(i);
        lflux(i) = -(k(i) + k(i-1))/2/(z(i) - z(i-1));
    end
end

%Top of the domain
if w(n) < 0 %inflow
    dflux(n+1) = 0;
    lflux(n+1) = w(n);
else %outflow
    dflux(n+1) = 0;
    lflux(n+1) = -w(n);
end

jac(jrowfdv(1,1),1) = (lflux(2) - dflux(1))/(z(2) - z(1));
jac(jrowfdv(1,2),2) = dflux(2)/(z(2) - z(1));
for i = 2:n-1
    jac(jrowfdv(i, i-1), i-1) = (-lflux(i))/(z(i+1)-z(i-1))*2;
    jac(jrowfdv(i, i), i) = (lflux(i+1) - dflux(i))/(z(i+1) - z(i-1))*2;
    jac(jrowfdv(i, i+1), i+1) = dflux(i+1)/(z(i+1) - z(i-1))*2;
end
jac(jrowfdv(n, n-1), n-1) = (-lflux(n))/(z(n) - z(n-1));
jac(jrowfdv(n, n), n) = (lflux(n+1) - dflux(n))/(z(n) - z(n-1));

jac = spdiags(jac', 1:-1:-1, n, n);
end


function jr = jrowfdv(i, j)
%gives the row of the BLAS banded format for pentadiagonal Jacobian
kl = 1;
ku = 1;
if ((i<=0) || (j<=0))
    print("error")
end
jr = ku + 1 + i - j;

end


function dc = advdiff_fun_fdv(n, z, w, k, bdry, surfem, vd, c)
%time derivative of concentration
dc = zeros(n);

%difflux/advflux = diffusive/advective fluxes through i-1/2
flux = zeros(n+1);

%ground level boundary
flux(1) = vd*c(1) - surfem;

%intermediate boundaries
for i = 2:n
    if w(i) >= 0
        advflux = w(i-1)*c(i-1);
    else
        advflux = w(i)*c(i);
    end
    difflux = (k(i) + k(i-1))/2*(c(i)-c(i-1))/(z(i)-z(i-1));
    flux(i) = difflux - advflux;
end

%top of the domain boundary
if w(n) < 0 %inflow
    advflux = w(n)*bdry(2);
    difflux = k(n)*(bdry(2)-c(n))/(z(n) - z(n-1));
else %outflow
    advflux = w(n)*c(n);
    difflux = 0;
end
flux(n+1) = difflux - advflux;

%time derivatives
dc(1) = (flux(2) - flux(1))/(z(2) - z(1));
for i = 2:n
    dc(i) = (flux(i+1)-flux(i))/(z(i+1)-z(i-1))*2;
end
dc(n) = (flux(n+1)-flux(n))/(z(n)-z(n-1));
end

function c = advdiff_fdv(dt, nstep, n, nspec, z, w, k, c, bdry, surfaceem, volumeem, vd)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Performs Nstep timesteps of length DT
%      to solve the adv_diff equation in vertical direction
%      using finite volume method and Crank-Nicholson
%
%  N     = no. of grid points
%  Nspec = no. of chemical species
%  Nstep = no of time steps
%  Z(1:NGP) = grid point coordinates
%  W(1:NGP) = wind speeds
%  K(1:NGP) = diffusion coefficients
%  SurfaceEm  = Surface Emission intensity
%  VolumeEm   = Elevated Emission intensity
%  Vd    = deposition velocity
%  C     = concentration of each species
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jac = advdiff_jac_fdv(n, z, w, k);

A = -dt/2*jac;
A = A + eye(n);

for ispec=1:nspec
    for istep = 1:nstep
        ctmp = c(:, ispec);
        plot(ctmp);
        pause;
        d = bdry(:, ispec);
        %the free term
        b = advdiff_free_fdv(n, z, w, k, d, surfaceem(ispec), ...
            vd(ispec), ctmp);
        %ctmp := c + dt*b + dt*volumeemission
        ctmp = ctmp + dt*b' + dt*volumeem(:, ispec);
        %ctmp := dt/2*jac*c + ctmp
        %      = (I + dt/2*jac)*c + dt*b + dt*volumeemission
        ctmp = ((dt/2)*jac)*c(:, ispec) + ctmp;
        %ctmp := A\ctmp
        %      = [I - dt/2*jac]\[ (I + dt/2*jac)*c + dt*b +...
        %      dt*volumeemission]
        ctmp = A\ctmp;
        c(:, ispec) = ctmp;
        
        %wrap around
        c(2, ispec) = c(n, ispec);
        c(1, ispec) = c(n-1,ispec);
        
    end
    disp("next species");
end
return;
end

function b = advdiff_free_fdv(n, z, w, k, bdry, sem, vd, c)
%ground level
b(1) = (-vd*c(1) + sem)/(z(2) - z(1));

%top of the domain
if w(n) < 0
    advflux = w(n)*bdry(2);
    difflux = k(n)*bdry(2)/(z(n) - z(n-1));
else
    advflux = 0;
    difflux = 0;
end
b(n) = (difflux - advflux)/(z(n) - z(n-1));
return;
end












%%%%%%%%%%%%%%%%%%%%%%%adjoints
%Discrete adjoint model of X transport
function lambda = trasport_X_dadj(dt, u, kh, concentration, lambda)
global ntracer;
global nx;
global ny;
global nz;

dx = 1000;
deltax = dx;
nsteps = floor(dt/300)+1;
deltat = dt/nsteps;

for j = 1:ny
    for k = 1:nz
        lam = squeeze(lambda(1:nx, j, k, 1:ntracer));
        wind = u(1:nx, j, k);
        dif = kh(1:nx, j, k);
        
        lambda(:, j, k, :) = adj_advdiff_fdh(dt, nsteps, nx, ntracer, deltax, wind, dif, lam);
    end
end
end



%verified with finite difference
function lam = adj_advdiff_fdh(dt, nstep, n, nspec, dx, u, k, lam)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  The adjoint of advdiff_fdh
%
%  N   = no. of grid points
%  Nspec = no. of chemical species
%  Nstep = no of time steps
%  X(1:N) = grid point coordinates
%  U(1:N) = wind speeds
%  K(1:N) = diffusion coefficients
%  SurfaceEm  = Surface Emission intensity
%  VolumeEm   = Elevated Emission intensity
%  Vd      = deposition velocity
%  Lam     = adjoint of concentration of each species
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%the jacobian
jac = advdiff_jac_fdh(n, dx, u, k);
%a = I - dt/2*jac
A = (-dt/2)*jac;
A = eye(n) + A; %add 1 to diagonal terms

for ispec=1:nspec
    for istep=1:nstep
        %Lam := (I - DT/2*JAC^T)\Lam
        lam(:, ispec) = (A')\lam(:, ispec);
        %temp = lam(:, ispec);
        %Lam := (dt/2)*Jac^T*Lam + Lam
        %     = (I + dt/2)*Jac)^T * Lam
        lam(:, ispec) = (dt/2)*(jac')*lam(:, ispec) + lam(:, ispec);
    end
end
end