model = initIntModel();

%Calculate true trajectory (for testing and performance examination)
fwdoptions = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac);
[~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, model.times, model.x0, fwdoptions);
xt = y.'; %MATLODE returns integration states as row vectors.
m = numel(model.times) - 1;
x = 8*randn(40, m);
model.stateestimate = [model.x0, x]; %store initial model estimate for cost function calculations.
xf = x(:);                           %this needs to be updated before each optimization 


% %See if L^-1 and L^-T run without issues, return correct sizes
% Linvtrans(x, x, model)
% Linv(x, x, model)

% %%%Test block
% cf = @(x) costfcn(x, model);
% gradfunc = @(x) gradfun(x, model);
% %findiffverify(xf, cf, gradfunc);
% optoptions = optimoptions('fminunc', 'MaxIterations', 1250, 'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'Display', 'iter');
% xg = fminunc(cf, x, optoptions);
% xg = [model.x0, xg];
% norm(xg - xt)


%%%GN Hessian Test Block
while true
    g = reshape(gradfun(x, model), 40, m);
    dx = Linv(x, Linvtrans(x, -g, model), model)
    x = x + dx;
    norm(x - xt(:, 2:end))
end


%The gradient appears to  verify properly with my own finite differences
%implementation, or the MATLAB option for fminunc.
function [c, G] = costfcn(x, model)
M = numel(model.times) - 1; %number of points from trajectory minus 1 (initial value already known).
x = reshape(x, 40, M);
x = [model.x0, x]; %prepend known initial value
u = zeros(40, M); %store cost function integrations. 
u = [model.x0, x]; %prepend known initial value to make indices line up
diffs = zeros(size(u)); %to store differences between integration values
scaleddiffs = zeros(size(diffs)); %differences with scaling matrices applied
adj = zeros(40, M + 1); %adjoint applied to appropriately scaled difference
G = zeros(40, M + 1); %to store gradient vectors. will be flattened before return.
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
%Right now I have it set to return all 1s (identity scaling matrix) so that
%I don't have to debug these parts of the L^{-1}, L^{-T} applications yet.
function d = rMat(u, model)
%d = (model.rtol*abs(u) + model.atol).^2;
d = ones(size(u));
%R = diag(d)
end

%takes  timespan [t_0, t_1], initial state y, vector v that adjoint model over given
%parameters should be applied to
function v = adjmodel(tspan, y, v, model)
%adjoint run
%Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', model.atol, 'ADJ_RelTol', model.rtol);
Options  = MATLODE_OPTIONS('AbsTol',1e-10, 'RelTol', 1e-10, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', 1e-12, 'ADJ_RelTol', 1e-12);
[~, ~, v] = MATLODE_SDIRK_ADJ_Integrator(model.rhs, tspan, y, Options);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
end


%takes  timespan [t_0, t_1], initial state y, vector v that TLM over given
%parameters should be applied to
function v = tlmmodel(tspan, y, v, model)
%tlm run
Options = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac, 'Y_TLM', eye(40));
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
m = size(x, 2); 
Lntv = v(:, m);

times = model.times;
for j = 1:m - 1 %outer loop represents row in L^-T matrix. Scale after integrations
    lambda = v(:, j);
    for i = j:m - 1
        tspan = [times(i), times(i+1)];
        lambda = adjmodel(tspan, x(:, i), lambda, model);
    end
    Lntv = [Lntv, lambda];
end
Lntv = reshape(Lntv, 40, m);
%apply scaling components after adjoint runs
% disp(size(Lntv))
%for i = 1:m
%    d = rMat(model.stateestimate(:, i), model);
%    Lntv(:, i) = sqrt(d).*Lntv(:, i);
%end

end

%to set up the system resulting from the cholesky decomposition applied to
% the GN system, we need the applications of the resulting L^-1 and L^-T
% matrices. This function returns the application of L^-T applied to v. I
% suspect that this is currently working (not including the application of 
%scaling matrices) but that Linvt is not.
function Lnv = Linv(x, v, model)
m = size(v, 2); 
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
%apply scaling components after TLM runs 
%for i = 1:m
%    d = rMat(model.stateestimate(:, i), model);
%    Lnv(:, i) = sqrt(d).*Lnv(:, i);
%end

end
