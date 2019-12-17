model = initIntModel();

%Calculate true trajectory (for testing and performance examination)
fwdoptions = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac);
[~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, model.times, model.x0, fwdoptions);
xt = y.'; %MATLODE returns integration states as row vectors.
m = numel(model.times) - 1;
x = 8*randn(model.n, m);
model.stateestimate = [model.x0, x]; %store initial model estimate for cost function calculations.
xf = x(:);                           %this needs to be updated before each optimization 

% % %%%Apply Linv, Linvtrans to full block Hessian, and see if they match
% % %%%By comparing norm to identity.
% H = fullHessian(x, model);
% prod = zeros(size(H));
% for i = 1:size(H, 1)
%     g = reshape(H(:, i), model.n, m);
%     myprod = Linv(x, Linvtrans(x, g, model), model);
%     prod(:, i) = reshape(myprod, size(H, 1), 1);
% end
% disp(norm(prod - eye(size(H,1))));
% %%%Test block
% cf = @(x) costfcn(x, model);
% gradfunc = @(x) gradfun(x, model);
% %findiffverify(xf, cf, gradfunc);
% optoptions = optimoptions('fminunc', 'MaxIterations', 1250, 'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'Display', 'iter');
% xg = fminunc(cf, x, optoptions);
% xg = [model.x0, xg];
% norm(xg - xt)
% 
% 
% %%GN Hessian Test Block using Linv, Linvtrans applications
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% while true
%     g = reshape(gradfun(x, model), model.n, m);
%     dx = Linv(x, Linvtrans(x, -g, model), model);
%     x = x + dx;
%     norm(x - xt(:, 2:end))
% end
% 
% 
% %GN Hessian Test Block with full Hessian
% xtest = xt(:, 2:end);
% x = xtest + 1*randn(size(xtest));
% while true
%     g = reshape(gradfun(x, model), model.n, m);
%     g = gradfun(x, model);
%     H = fullHessian(x, model);
%     dx = H\(-g);
%     dx = reshape(dx, model.n, m);
%     x = x + dx;
%    	norm(x - xt(:, 2:end))
% end
% 
% 
% 
%The gradient appears to  verify properly with my own finite differences
%implementation, or the MATLAB option for fminunc.
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
for i = 2:M + 1
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
for i = 2:M
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
%SDIRK TLM doesn't seems to break on application to a single input  vector for some
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

%my own implementation of a finite differences verification routine for the
%gradient
function maxdiff = findiffverify(x, fun, grad)
    h = 1e-8;
    n = numel(x);
    maxdiff = 0;
    A = eye(n); %grab standard basis vectors from here
    gradeval = grad(x);
    diffs = zeros(n, 1);
    parfor i = 1:n %test in all coordinate directions
        stepf = x + h*A(:, i);
        stepb = x - h*A(:, i);
        fd = (fun(stepf) - fun(stepb))/(2*h);
        reldiff = abs(fd - gradeval(i))/abs(fd);
        %disp(i)
        disp(reldiff)
        diffs(i) = reldiff;
    end
    disp("Max relative difference between finite difference approximation and grad function:");
    disp(max(diffs));
    disp("Mean relative difference:");
    disp(mean(diffs));
end

function G = gradfun(x, model)
    [c, G] = costfcn(x, model);
    return;
end

%to set up the system resulting from the cholesky decomposition applied to
% the GN system, we need the applications of the resulting L^-1 and L^-T
% matrices. This function returns the application of L^-T applied to v. 
function Lntv = Linvtrans(x, v, model)
m = size(v, 2); 
for i = 1:m
    d = rMat(model.stateestimate(:, i), model);
    v(:, i) = sqrt(d).*v(:, i);
end

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
