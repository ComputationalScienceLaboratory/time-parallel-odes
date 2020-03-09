%Try to build the model object and call the advdiff_fdh function without
%crashing.
%Start with only x
%Take first coordinate of everything
model = initAdvDiffModel();
c = squeeze(model.x0(:, 1, 1, :));
sb(1, :) = model.bc.w(1, 1, :);
sb(2, :) = model.bc.e(1, 1, :);
dif = squeeze(model.kh(:, 1, 1));
wind = model.wind.u(:, 1, 1);

c = advdiff_fdh(model.dt, model.M, model.nx, model.ntracer, model.dx, wind, dif, c, sb);





% %%%Test block for optimization with gradients only, using MATLAB
% %%%optimizer
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% model.stateestimate = [model.x0, x]; 
% cf = @(x) costfcn(x, model);
% gradfunc = @(x) gradfun(x, model);
% %findiffverify(xf, cf, gradfunc);
% optoptions = optimoptions('fminunc', 'MaxIterations', 1250, 'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'Display', 'iter');
% xg = fminunc(cf, x, optoptions);
% xg = [model.x0, xg];
% norm(xg - xt)


% % %%GN Hessian Test Block using Linv, Linvtrans applications and line
% % search
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% model.stateestimate = [model.x0, x]; 
% while true
%     [c, g] = costfcn(x, model);
%     g = reshape(gradfun(x, model), model.n, m);
%     dx = Linv(x, Linvtrans(x, -g, model), model);
%     a = backtrack(x(:), dx(:), model);
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
%     c
%     norm(x - xt(:, 2:end))
% end
% 

% 
% %%%GN Hessian Test Block using coarse Linv, coarse Linvtrans applications and line
% %%%search
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% model.stateestimate = [model.x0, x]; 
% while true
%     [c, g] = costfcn(x, model);
%     g = reshape(gradfun(x, model), model.n, m);
%     dx = cLinv(x, cLinvtrans(x, -g, model), model);
%     a = backtrack(x(:), dx(:), model);
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
%     c
%     norm(x - xt(:, 2:end))
% end


% % % 
% % % 
% %GN Hessian Test Block with full Hessian and line search
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% model.stateestimate = [model.x0, x]; 
% while true
%     [c, g] = costfcn(x, model);
%     H = fullHessian(x, model);
%     dx = H\(-g);
%     a = backtrack(x(:), dx, model);
%     dx = reshape(dx, model.n, m);
%     c
%     x = x + a*dx;
%    	norm(x - xt(:, 2:end))
% end
% 
% % %%% Test block comparing Linv, Linvtrans, Hessian solves on random
% % %%% vectors
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% model.stateestimate = [model.x0, x]; 
% while true
%      mine = 100*randn(size(x));
%      g = reshape(mine, model.n, m);
%      dx = Linv(x, Linvtrans(x, g, model), model);
%      H = fullHessian(x, model);
%      dx2 = H\mine(:);
%      norm(dx2(:) - dx(:))
% end
% % 

%Return cost function and gradient. Depends on model.stateestimate to form
%scaling matrices.
function [c, G] = costfcn(x, model)
M = numel(model.times) - 1; %number of points from trajectory minus 1 (initial value already known).
x = reshape(x, model.n, M);

x = [model.x0, x]; %prepend known initial value
u = zeros(model.n, M); %store cost function integrations. 
u = [model.x0, x]; %prepend known initial value to make indices line up
diffs = zeros(size(u)); %to store differences between integration values
scaleddiffs = zeros(size(diffs)); %differences with scaling matrices applied
adj = zeros(model.n, M + 1); %adjoint applied to appropriately scaled difference
G = zeros(model.n, M + 1); %to store gradient vectors. will be flattened before return.
%cost function integrations (integrate estimated system states forward one
%step each)
parfor i = 2:M + 1
    u(:, i) = fwdmodel([model.times(i - 1), model.times(i)], x(:, i - 1), model);
    diffs(:, i) = x(:, i) - u(:, i);
    d = rMat(model.stateestimate(:, i), model); %diagonal elements of scaling matrix
    scaleddiffs(:, i) = (1./d).*diffs(:, i); %inverse of diagonal matrix
end
%Sum cost function terms
c = 0;
for i = 1:M + 1
   c = c + diffs(:, i).'*scaleddiffs(:, i);
end
c = c/2;

%Calculate gradients
%Adjoint model runs
parfor i = 2:M
   v = diffs(:, i + 1);
   d = rMat(model.stateestimate(:, i + 1), model); %diagonal elements of scaling matrix
   v = (1./d).*diffs(:, i+1);  %inverse of diagonal matrix
   t = [model.times(i), model.times(i + 1)];
   v = adjmodel(t, x(:, i), v, model);
   adj(:, i) = v;
end

for i = 2:M + 1
    G(:, i) = scaleddiffs(:, i) - adj(:, i);
end
G = G(:, 2:end); %first column is known initial value (not fed to cost function).
G = G(:);       %throw out
return;

end

%Given u, return the diagonal of the diagonal scaling matrix given in Vishwas 
%& Sandu.
%return ones to give identity matrix (for testing)
function d = rMat(u, model)
d = (model.rtol*abs(u) + model.atol).^2;
%d = ones(size(u));
%R = diag(d)
end

%takes  timespan [t_0, t_1], initial state y, vector v that adjoint model over given
%parameters should be applied to
function v = adjmodel(tspan, y, v, model)
%adjoint run
%Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', model.atol, 'ADJ_RelTol', model.rtol);
Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', model.atol, 'ADJ_RelTol', model.rtol);
[~, ~, v] = MATLODE_SDIRK_ADJ_Integrator(model.rhs, tspan, y, Options);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
end


%takes  timespan [t_0, t_1], initial state y, vector v that TLM over given
%parameters should be applied to
function v = tlmmodel(tspan, y, v, model)
%tlm run
Options = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac, 'Y_TLM', eye(model.n), 'TLM_AbsTol', model.atol, 'TLM_RelTol',model.rtol);
%SDIRK TLM seems to break on application to a single input  vector for some
%reason, so return full sensitivity matrix and apply to v
[~, ~, sens] = MATLODE_SDIRK_TLM_Integrator(model.rhs, tspan, y, Options);
v = sens*v;
end

%given timespan [t_0, t_1], initial value y, return system state at t_1
%starting from state y at t_0
function y = fwdmodel(tspan, y, model)
fwdoptions = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac);

[~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, tspan, y, fwdoptions);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
y = y(end, :).';
end

%same as above but return a full trajectory instead of the final state
function y = fullfwd(times, y, model)
fwdoptions = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac);

[~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, times, y, fwdoptions);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
y = y.';
end


%takes  timespan [t_0, t_1], initial state y, vector v that adjoint model over given
%parameters should be applied to
function y = coarseadj(tspan, y, v, model)
%adjoint run
%Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', model.atol, 'ADJ_RelTol', model.rtol);
Options  = MATLODE_OPTIONS('AbsTol',model.catol, 'RelTol', model.crtol, 'Lambda', @(t, y) v, 'ADJ_AbsTol', model.catol, 'ADJ_RelTol', model.crtol, 'Jacobian', model.jac);
[~, ~, v] = MATLODE_ERK_ADJ_Integrator(model.rhs, tspan, y, Options);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
end

%takes  timespan [t_0, t_1], initial state y, vector v that TLM over given
%parameters should be applied to
function v = coarsetlm(tspan, y, v, model)
%tlm run
Options = MATLODE_OPTIONS('Jacobian', model.jac, 'AbsTol', model.catol, 'RelTol', model.crtol,'Y_TLM', eye(model.n), 'TLM_AbsTol', model.catol, 'TLM_RelTol',model.crtol);
%SDIRK TLM seems to break on application to a single input  vector for some
%reason, so return full sensitivity matrix and apply to v
[~, ~, sens] = MATLODE_ERK_TLM_Integrator(model.rhs, tspan, y, Options);
v = sens*v;
end

%given timespan [t_0, t_1], initial value y, return system state at t_1
%starting from state y at t_0
function y = coarsefwd(tspan, y, model)
fwdoptions = MATLODE_OPTIONS('AbsTol', model.catol, 'RelTol', model.crtol);

[~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, tspan, y, fwdoptions);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
y = y(end, :).';
end


%to set up the system resulting from the cholesky decomposition applied to
% the GN system, we need the applications of the resulting L^-1 and L^-T
% matrices. This function returns the application of L^-T applied to v. 
function Lntv = Linvtrans(x, v, model)
m = size(v, 2); 

Lntv = v(:, m);
times = model.times;
for j = (m - 1):-1:1 %outer loop represents row in L^-T matrix. Scale after integrations
    lambda = v(:, j);
    for i = j:m - 1
        tspan = [times(i), times(i+1)];
        lambda = adjmodel(tspan, x(:, i), lambda, model);
    end
    Lntv = [lambda, Lntv];
end
Lntv = reshape(Lntv, model.n, m);

for i = 1:m
    d = rMat(model.stateestimate(:, i), model);
    v(:, i) = sqrt(d).*v(:, i);
end

end

%to set up the system resulting from the cholesky decomposition applied to
% the GN system, we need the applications of the resulting L^-1 and L^-T
% matrices. This function returns the application of L^-T applied to v.
function Lnv = Linv(x, v, model)
m = size(v, 2); 

%apply scaling components before TLM runs 
for i = 1:m
    d = rMat(model.stateestimate(:, i), model);
    v(:, i) = sqrt(d).*v(:, i);
end



Lnv = v(:, 1);
times = model.times;
%outer loop represents row in L^-1 matrix. Scale after integrations.
%since the first row applied to V is simply v(:, 1), we start at 2.
for j = 2:m  
    xi = v(:, 1);
    for i = 1:j-1
        tspan = [times(i), times(i+1)];
        xi = tlmmodel(tspan, x(:, i), xi, model);
    end
    xi = xi + v(:, j);
    Lnv = [Lnv, xi];
end
Lnv = reshape(Lnv, model.n, m);

end



%to set up the system resulting from the cholesky decomposition applied to
% the GN system, we need the applications of the resulting L^-1 and L^-T
% matrices. This function returns the application of L^-T applied to v.
function Lnv = cLinv(x, v, model)
m = size(v, 2); 

%apply scaling components before TLM runs 
parfor i = 1:m
    d = rMat(model.stateestimate(:, i), model);
    v(:, i) = sqrt(d).*v(:, i);
end



Lnv = v(:, 1);
times = model.times;
%outer loop represents row in L^-1 matrix. Scale after integrations.
%since the first row applied to V is simply v(:, 1), we start at 2.
parfor j = 2:m  
    xi = v(:, 1);
    for i = 1:j-1
        tspan = [times(i), times(i+1)];
        xi = coarsetlm(tspan, x(:, i), xi, model);
    end
    xi = xi +  v(:, j);
    Lnv = [Lnv, xi];
end
Lnv = reshape(Lnv, model.n, m);

end


%to set up the system resulting from the cholesky decomposition applied to
% the GN system, we need the applications of the resulting L^-1 and L^-T
% matrices. This function returns the application of L^-T applied to v. 
function Lntv = cLinvtrans(x, v, model)
m = size(v, 2); 

Lntv = v(:, m);
times = model.times;
%for j = (m - 1):-1:1 %outer loop represents row in L^-T matrix. Scale after integrations
for j = 1:(m - 1)
    lambda = v(:, j);
    for i = j:m - 1
        tspan = [times(i), times(i+1)];
        lambda = coarseadj(tspan, x(:, i), lambda, model);
    end
    %Lntv = [lambda, Lntv];
    Lntv = [Lntv, lambda];
end
Lntv = reshape(Lntv, model.n, m);

parfor i = 1:m
    d = rMat(model.stateestimate(:, i), model);
    v(:, i) = sqrt(d).*v(:, i);
end

end

function Hes = fullHessian(x, model)
diags = {}; %diagonal elements of hessian matrix
bdiags = {}; %below diagonal elements
adiags = {}; %above diagonal elements
m = size(x, 2); %m is the number of free (vectors) along the trajectory
n = size(x, 1); %n is the dimension of state vectors
H = sparse(n*m, n*m);
%build diagonal blocks
for i = 1:m
    tspan = [model.times(i), model.times(i+1)];
    mat = tlmmodel(tspan, x(:, i), eye(n), model);
    d = rMat(model.stateestimate(:, i), model);
    mat = (1./d).*mat;
    mat = adjmodel(tspan, x(:,i), mat, model);
    mat = diag(1./d) + mat;
    diags{end+1} = mat;
end
%build below diagonal blocks
for i = 2:m
   d = rMat(model.stateestimate(:, i), model);
   tspan = [model.times(i), model.times(i+1)];
   mat = tlmmodel(tspan, x(:, i), eye(n), model);
   mat = diag(1./d) * mat;
   bdiags{end + 1} = mat;
end
%build above diagonal blocks
for i = 2:m
   d = rMat(model.stateestimate(:, i), model);
   tspan = [model.times(i), model.times(i+1)];
   mat = adjmodel(tspan, x(:, i), eye(n), model);
   mat = mat * diag(1./d);
   adiags{end + 1} = mat;
end
%build full GN Hessian from blocks
%build top row outside loop
H(1:n,1:n) = diags{1};
%If the Hessian has an above diagonal block
try
    H(1:n,n+1:2*n) = adiags{1};
end
%build bottom row outside loop
H((m - 1)*n + 1:m*n,(m-1)*n + 1:m*n) = diags{1};
%If the Hessian has a below diagonal block
try
    H((m - 1)*n + 1:m*n,(m-2)*n + 1:(m-1)*n) = bdiags{1};
end
%main loop builds all rows except first and last (i.e. those containing all
%3 blocks)
for i = 2:m-1
    H((i-1)*n + 1:i*n,(i-1)*n + 1:i*n) = diags{1};%diagonal blocks
    H((i-1)*n + 1:i*n,(i-2)*n + 1:(i-1)*n) = bdiags{1}; %below diagonal blocks
    H((i-1)*n + 1:i*n,i*n + 1:(i+1)*n) = adiags{1};
end
Hes = H;
end

%for use with finite difference verification routine
function G = gradfun(x, model)
    [c, G] = costfcn(x, model);
    return;
end



function a = backtrack(x, p, model)
r = .9;
c = 1e-4;
a = 1; %initial step size
%af = 0;
while true
[f, g] = costfcn(x, model);
try
[af, ag] = costfcn(x +  a*p, model);

if imag(af) == 0 && af <= (f + c*a*(g.'*p))
    break
end
end
if a < numel(x)*eps
    break
end

a = r*a;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below here are the advection/diffusion code functions.
%They are used in the forward, tlm, adjoint model functions above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        ctmp = c(:, ispec);
        %plot(ctmp)
        %pause;
        d = bdry(:, ispec);
        %the free term
        b = advdiff_free_fdh(n, dx, u, k, d);
        ctmp = ctmp + dt*b;
        alpha = (dt/2);
        ctmp = alpha*jac*c(:, ispec) + ctmp;
        ctmp = A\ctmp;
        c(:, ispec) = ctmp;
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

