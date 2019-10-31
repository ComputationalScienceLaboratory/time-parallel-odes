model = initIntModel();

%Calculate true trajectory (for testing and performance examination)
fwdoptions = MATLODE_OPTIONS('AbsTol', model.atol, 'RelTol', model.rtol, 'Jacobian', model.jac);
[~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, model.times, model.x0, fwdoptions);
xt = y.'; %MATLODE returns integration states as row vectors.
m = numel(model.times) - 1;
x = randn(40, m);
% %%%%Test block
% cf = @(x) costfcn(x, model);
% optoptions = optimoptions('fminunc', 'MaxIterations', 1250, 'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'Display', 'iter');
% xg = fminunc(cf, x, optoptions);
% xg = [model.x0, xg];
% norm(xg - xt)
xt = xt(:, 2:end);
xt = xt(:);
x = x(:);
for i = 1:100
    H = @(v) hvp(x, v, model);
    [~, G] = costfcn(x, model);
    dx = gmres(H, -G);
    norm(dx)
    x = x + dx;
    norm(x - xt)/norm(xt)
end


%With the scaling matrix applied, this doesn't seem to verify with the
%MATLAB finite differences check in optimoptions. Actually applying the
%optimization routine with the gradient seems to work, so I assume that
%this is working.
function [c, G] = costfcn(x, model)
M = numel(model.times) - 1; %number of points from trajectory minus 1. initial value already known
x = reshape(x, 40, M);
x = [model.x0, x]; %prepend known initial value
u = zeros(40, M); %store cost function integrations. 
u = [model.x0, x]; %prepend known initial value to make indices line up
diffs = zeros(size(u)); %to store differences between integration values
scaleddiffs = zeros(size(diffs)); %differences with scaling matrices applied
adj = zeros(40, M + 1); %adjoint applied to appropriately scaled difference
G = zeros(40, M); %to store gradient vectors. will be flattened before return.
%cost function integrations (integrate estimated system states forward one
%step each)
parfor i = 2:M + 1
    u(:, i) = fwdmodel([model.times(i - 1), model.times(i)], x(:, i - 1), model);
    diffs(:, i) = x(:, i) - u(:, i);
    d = rMat(x(:, i), model);
    scaleddiffs(:, i) = (1./d).*diffs(:, i);
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
   %v = diffs(:, i + 1);
   d = rMat(x(:, i), model);
   v = (1./d).*diffs(:, i+1);
   t = [model.times(i), model.times(i + 1)];
   v = adjmodel(t, x(:, i), v, model);
   adj(:, i) = v;
end

parfor i = 2:M + 1
    G(:, i - 1) = scaleddiffs(:, i) - adj(:, i);
end
G = G(:);
return;

end

% (Approximate) Hessian vector product (Hessian approximation at x applied
% to vector v).
% function Hv = hvp(x, v, model)
% M = numel(model.times) - 1; %number of points from trajectory minus 1. initial value already known
% x = reshape(x, 40, M);
% x = [model.x0, x]; %prepend known initial value
% v = reshape(v, 40, M); %reshape vector transformation will be applied to
% v = [zeros(40, 1), v]; %prepend initial value so indices match
% tlms = zeros(size(v)); %store tlm applications to elements of v (diagonal and off diagonal)
% adjs = zeros(size(v)); %adjoint applications to v (off diagonal)
% Vts = zeros(size(v)); %diagonal blocks of H matrix
% Uts = zeros(size(v)); %below diagonal blocks of H matrix
% Wts = zeros(size(v)); %above diagonal blocks of H matrix
% Hv = zeros(size(v)); %store Hessian vector product
% TLM and adjoint runs (applied to TLM outputs)
% parfor i = 2:M
%     tspan = [model.times(i), model.times(i+1)];
%     tlms(:, i) = tlmmodel(tspan, x(:, i), v(:, i), model);
%     d = rMat(x(:,i), model);
%     scaledtlm = (1./d).*tlms(:, i);
%     adjtlm = adjmodel(tspan, x(:, i), scaledtlm, model);
%     Vts(:, i) = (1./d).*v(:, i) + adjtlm;
% end
% parfor i = 3:M
%    tspan = [model.times(i - 1), model.times(i)];
%    d = rMat(x(:, i), model);
%    Uts(:, i) = (1./d).*tlms(:, i - 1);
% end
% parfor i = 2:M
%    tspan = [model.times(i), model.times(i+1)];
%    d = rMat(x(:, i), model);
%    scaledVt = (1./d).*Vts(:, i + 1);
%    Wts(:, i) = adjmodel(tspan, x(:, i), scaledVt, model);
% end
% for i = 2:M + 1
%     Hv(:, i) = -Uts(:, i) + Vts(:, i) - Wts(:, i);
% end
% Hv = Hv(:, 2:end);
% Hv = Hv(:);
% end
%HVP attempt 2 (from scratch)
function H = hvp(x, v, model)
M = numel(model.times) - 1; %number of points from trajectory minus 1. initial value already known
x = reshape(x, 40, M);
x = [model.x0, x]; %prepend known initial value
v = reshape(v, 40, M); %reshape vector transformation will be applied to
v = [zeros(40, 1), v]; %prepend initial value so indices match
tlms = zeros(size(v)); %store tlm applications to elements of v (diagonal and off diagonal)
adjs = zeros(size(v)); %adjoint applications to v (off diagonal)
adjtlms = zeros(size(v));
Hv = zeros(size(v)); %store Hessian vector product

%all TLM runs
parfor i = 2:M
    tspan = [model.times(i), model.times(i+1)];
    tlms(:, i) = tlmmodel(tspan, x(:, i), v(:, i), model);
end
%scale TLM outputs and apply adjoint model
parfor i = 2:M
   tspan = [model.times(i), model.times(i+1)];
   d = rMat(x(:, i), model);
   scaledtlm = (1./d).*tlms(:, i);
   adjtlm(:, i) = adjmodel(tspan, x(:, i), scaledtlm, model);
end
%apply adjoint model to vectors of v
parfor i = 2:M-1
  tspan = [model.times(i), model.times(i+1)];
  d = rMat(x(:, i), model);
  sv = (1./d).*v(:, i + 1);
  adjs(:, i) = adjmodel(tspan, x(:, i - 1), sv, model);
end
%apply diagonal blocks
parfor i = 2:M
   d = rMat(x(:, i), model);
   Hv(:, i) = (1./d).*v(:, i) + adjtlm(:, i);
end
%apply above diagonal blocks
parfor i = 2:M-1
   d = rMat(x(:, i + 1), model);
   Hv(:, i+1) = Hv(:, i+1) - (1./d).*adjs(:, i); 
end
%apply below diagonal blocks
parfor i = 2:M-1
   d = rMat(x(:, i + 1), model);
   Hv(:, i) = Hv(:, i) - (1./d).*tlms(:, i);
end
for i = 1:M
    row = zeros(40, 40*M);
    row(:, 1:40 = 
end

end

%Given u, return the diagonal of the diagonal scaling matrix given in Vishwas 
%& Sandu
function d = rMat(u, model)
d = (model.rtol*abs(u) + model.atol).^2;
%d = ones(size(u));
%R = diag(d);
end


%takes  timespan [t_0, t_1], initial state y, vector v that adjoint model over given
%parameters should be applied to
function v = adjmodel(tspan, y, v, model)
%adjoint run
%Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', model.atol, 'ADJ_RelTol', model.rtol);
Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', 1e-12, 'ADJ_RelTol', 1e-12);
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
