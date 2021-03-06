model = initAdvDiffModel();
%%%%%%%%%%%%%%Moving to 2D test problem
%build initial concentration
concentration = zeros(model.nx, model.ny, model.nz, model.ntracer);
cntr = [6, 6];
rad = .25;
for i = 1:model.nx
    for j = 1:model.ny
        pt = [i, j];
        if(norm(pt - cntr) <= rad)
            concentration(i, j, 1, 1) = 2;
        end
    end
end
model.x0 = concentration;
% surf(concentration(:, :, 1, 1));
% pause;
%clf;

% %build diffusion coefficients
% model.kh = 1000*ones(size(model.kh));
% model.kv = 1000*ones(size(model.kv));

%%%%%I want to build the wind profiles now
omega = 2;
windcntr = [5, 5];
u = zeros(model.nx, model.ny, model.nz, model.ntracer);
v = zeros(model.nx, model.ny, model.nz, model.ntracer);
for i = 1:model.nx
    for j = 1:model.ny
        u(i, j, 1, 1) = -omega*(j - windcntr(2));
        v(i, j, 1, 1) = omega*(i - windcntr(1));
    end
end
model.u = u;
model.v = v;
quiver(v(:, :, 1, 1),u(:, :, 1, 1))
hold on
contour(concentration(:, :, 1, 1));
hold off

movieit = 0;
%Integration loop (normal time integrators)
while true
    movieit = movieit + 1;
    
    concentration = par_transport_X_fwd(concentration, model);
    %concentration = transport_X_fwd(concentration, model);
    quiver(v(:, :, 1, 1),u(:, :, 1, 1))
    hold on
    contour(concentration(:, :, 1, 1));
    pause(5/1000);
    hold off
    pause(5/1000);
    
    concentration = par_transport_Y_fwd(concentration, model);
    %concentration = transport_Y_fwd(concentration, model);
    quiver(v(:, :, 1, 1),u(:, :, 1, 1))
    hold on
    contour(concentration(:, :, 1, 1));
    pause(5/1000);
    hold off
    
    M(movieit) = getframe;
end

%%%%%%%%end of 2d problem test code
%%%%%%%%2d problem 1d transport functions (our GN algo, coarsened H blcks)
function concentration = par_transport_Y_fwd(concentration, model)
% dx = model.dx;
% deltax = dx;
%nsteps = floor(dt/300)+1;
%deltat = dt/nsteps;
model.n = model.ny;

for i = 1:model.nx
    for k = 1:model.nz
        conc = squeeze(concentration(i, :, k, :))';
        wind = model.v(i, :, k);
        dif = model.kh(i, :, k);
        sb(1, :) = squeeze(model.bc.s(i, k, :));
        sb(2, :) = squeeze(model.bc.n(i, k, :));
        %for l = 1:nsteps
        %specify RHS, Jacobian
        jac = advdiff_jac_fdh(model.ny, model.dy, wind, dif);
        model.jac = @(t, x) jac;
        model.rhs = @(t, x) jac*x + advdiff_free_fdh(model.ny, model.dy, wind, dif, sb);
        model.x0 = conc;
        %pass to our algo
        %Do forward Euler for initial forward integration
        yc = zeros(model.nx, model.M);
        yc(:, 1) = conc;
        for f = 2:model.M
            yc(:, f) = yc(:, f- 1) + (model.times(f) - model.times(f-1))*model.rhs(model.times(f - 1), yc(:, f - 1));
        end
        model.stateestimate = yc;
        
        x = model.stateestimate(:, 2:end);
        optiter = 0;
        [c, g] = costfcn(x, model);
        
        while c > 1
            %disp(optiter);
            %disp(c);
            optiter = optiter + 1;
            Li = fullLinv(x, model);
            Lit = fullLinvT(x, model);
            dx = Li*Lit*(-g);
            a = backtrack(x(:), dx, model);
            dx = reshape(dx, model.n, model.M-1);
            x = x + a*dx;
            [c, g] = costfcn(x, model);
        end
        
        
        
        xt = x(:, end);
        concentration(i, :, k, :) = subplus(xt)'; %positive part
        %disp("next y, z")
        %pause;
    end
end

return;
end

function concentration = par_transport_X_fwd(concentration, model)
% dx = model.dx;
% deltax = dx;
%nsteps = floor(dt/300)+1;
%deltat = dt/nsteps;
model.n = model.nx;

for j = 1:model.ny
    for k = 1:model.nz
        conc = squeeze(concentration(:, j, k, :));
        wind = model.u(:, j, k);
        dif = model.kh(:, j, k);
        sb(1, :) = squeeze(model.bc.w(j, k, :));
        sb(2, :) = squeeze(model.bc.e(j, k, :));
        %for l = 1:nsteps
        %specify RHS, Jacobian
        jac = advdiff_jac_fdh(model.nx, model.dx, wind, dif);
        model.jac = @(t, x) jac;
        model.rhs = @(t, x) jac*x + advdiff_free_fdh(model.nx, model.dx, wind, dif, sb);
        %pass to our algo
        model.x0 = conc;
        
        %Do forward Euler for initial forward integration
        yc = zeros(model.nx, model.M);
        yc(:, 1) = conc;
        for f = 2:model.M %all the other indices are taken
            yc(:, f) = yc(:, f- 1) + (model.times(f) - model.times(f-1))*model.rhs(model.times(f - 1), yc(:, f - 1));
        end
        model.stateestimate = yc;
        
        x = model.stateestimate(:, 2:end);
        optiter = 0;
        [c, g] = costfcn(x, model);
        
        while c > 1
            %disp(optiter);
            %disp(c);
            optiter = optiter + 1;
            Li = fullLinv(x, model);
            Lit = fullLinvT(x, model);
            dx = Li*Lit*(-g);
            a = backtrack(x(:), dx, model);
            dx = reshape(dx, model.n, model.M-1);
            x = x + a*dx;
            [c, g] = costfcn(x, model);
        end
        
        
        
        xt = x(:, end);
        concentration(:, j, k, :) = subplus(xt)'; %positive part
        %disp("next y, z")
        %pause;
    end
end

return;
end



%%%%%%%%2d problem 1d transport functions (MATLODE time integrators)
function concentration = transport_Y_fwd(concentration, model)
% dx = model.dx;
% deltax = dx;
%nsteps = floor(dt/300)+1;
%deltat = dt/nsteps;

for i = 1:model.nx
    for k = 1:model.nz
        conc = squeeze(concentration(i, :, k, :))';
        wind = model.v(i, :, k);
        dif = model.kh(i, :, k);
        sb(1, :) = squeeze(model.bc.s(i, k, :));
        sb(2, :) = squeeze(model.bc.n(i, k, :));
        %for l = 1:nsteps
        %specify RHS, Jacobian
        jac = advdiff_jac_fdh(model.ny, model.dy, wind, dif);
        model.jac = @(t, x) jac;
        model.rhs = @(t, x) jac*x + advdiff_free_fdh(model.ny, model.dy, wind, dif, sb);
        %pass to MATLODE
        fwdoptions = MATLODE_OPTIONS('AbsTol', 1e-3, 'RelTol', 1e-3);
        [~, y] = MATLODE_ERK_FWD_Integrator(model.rhs, model.times, conc, fwdoptions);
        xt = y(end, :).'; %MATLODE returns integration states as row vectors.
        
        
        %end
        concentration(i, :, k, :) = subplus(xt)'; %positive part
        %disp("next y, z")
    end
end

return;
end


function concentration = transport_X_fwd(concentration, model)
dx = model.dx;
deltax = dx;
%nsteps = floor(dt/300)+1;
%deltat = dt/nsteps;

for j = 1:model.ny
    for k = 1:model.nz
        conc = squeeze(concentration(:, j, k, :));
        wind = model.u(:, j, k);
        dif = model.kh(:, j, k);
        sb(1, :) = squeeze(model.bc.w(j, k, :));
        sb(2, :) = squeeze(model.bc.e(j, k, :));
        %for l = 1:nsteps
        %specify RHS, Jacobian
        jac = advdiff_jac_fdh(model.nx, model.dx, wind, dif);
        model.jac = @(t, x) jac;
        model.rhs = @(t, x) jac*x + advdiff_free_fdh(model.nx, model.dx, wind, dif, sb);
        %pass to MATLODE
        fwdoptions = MATLODE_OPTIONS('AbsTol', 1e-3, 'RelTol', 1e-3);
        [~, y] = MATLODE_ERK_FWD_Integrator(model.rhs, model.times, conc, fwdoptions);
        xt = y(end, :).'; %MATLODE returns integration states as row vectors.
        
        
        %end
        concentration(:, j, k, :) = subplus(xt)'; %positive part
        %disp("next y, z")
    end
end

return;
end
%%%%%%%%%%end of 2d problem 1d transport functions

% %%%%Test that I can do the integrations in one direction with MATLODE
% %%%%Throw out everything but the first x direction
% u = model.wind.u;
% kh = model.kh;
% concentration = model.x0full;
% x = squeeze(model.x0full(:, 1, 1, 1));
% concentration_bc(1, :, :, :) = model.bc.w;
% concentration_bc(2, :, :, :) = model.bc.e;
% %%%%%%%%%%%%%%%%%%%%
% u = squeeze(u(:, 1, 1));
% u = ones(size(u));
% k = squeeze(kh(:, 1, 1));
% %jac = advdiff_jac_fdh(model.nx, model.dx, u, k);
% model.jac = @(t, x) advdiff_jac_fdh(model.nx, model.dx, u, k);
% jac = advdiff_jac_fdh(model.nx, model.dx, u, k);
% bdry(1) = concentration_bc(1, 1, 1, 1);
% bdry(2) = concentration_bc(2, 1, 1, 1);
% model.rhs = @(t, x) jac*x + advdiff_free_fdh(model.nx, model.dx, u, k, bdry);
% model.x0 = zeros(model.nx, 1);
% model.x0(25:50, 1) = 1;
%
% model.n = model.nx;
% %Calculate reference trajectory
% fwdoptions = MATLODE_OPTIONS('AbsTol', 1e-12, 'RelTol', 1e-12, 'Jacobian', model.jac);
% [~, y] = MATLODE_SDIRK_FWD_Integrator(model.rhs, model.times, model.x0, fwdoptions);
% xt = y.'; %MATLODE returns integration states as row vectors.
% plot(y)
% pause;
% m = numel(model.times) - 1;
% %x = norm(model.x0)*randn(model.n, m);
% %model.stateestimate = [model.x0, x]; %store initial model estimate for cost function calculations.
% % xf = x(:);                           %this needs to be updated before each optimization
%
%
% %Do forward Euler for initial forward integration
% yc = zeros(size(y'));
% yc(:, 1) = model.x0;
% for i = 2:model.M
%     yc(:, i) = yc(:, i- 1) + (model.times(i) - model.times(i-1))*model.rhs(model.times(i - 1), yc(:, i - 1));
% end
% model.stateestimate = yc;

%%%%NOTE TO VISHWAS
%%%%The chunks below here and be block uncommented to run each optimization
%%%%test.

% %%%Test block for optimization with gradients only, using MATLAB
% %%%optimizer
% x = model.stateestimate(:, 2:end);
% cf = @(x) costfcn(x, model);
% gradfunc = @(x) gradfun(x, model);
% optoptions = optimoptions('fminunc', 'MaxIterations', 1250, 'SpecifyObjectiveGradient', true, 'CheckGradients', false, 'Display', 'iter', 'StepTolerance', 1e-100);
% xg = fminunc(cf, x, optoptions);

% %As of 3/6/20, this is beating LBFGS pretty convincingly by iteration
% %GN Hessian Test Block with full Hessian and line search
% %(solves H\-g to get search direction)
% x = model.stateestimate(:, 2:end);
% i = 0;
% while true
%     [c, g] = costfcn(x, model);
%     disp(i);
%     disp(c);
%     i = i + 1;
%     H = fullHessian(x, model);
%     dx = H\(-g);
%     a = backtrack(x(:), dx, model);
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
% end

% %GN Hessian Test Block with full block L^{-T}, L^{-1} and line search
% %(Search direction by p = L^{-1}L^{-T}-g)
% %As of 3/9/20, this seems to be handily beating LBFGS and explicit H
% %solves
% x = model.stateestimate(:, 2:end);
% i = 0;
% while true
%     [c, g] = costfcn(x, model);
%     disp(i);
%     disp(c);
%     i = i + 1;
%     Li = fullLinv(x, model);
%     Lit = fullLinvT(x, model);
%     dx = Li*Lit*(-g);
%     a = backtrack(x(:), dx, model);
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
% end

% %GN Hessian Test Block with full block coarest block L^{-T}, L^{-1} and line search
% %(Search direction by p = L^{-1}L^{-T}-g)
% %As of 3/9/20, this seems to be handily beating LBFGS and explicit H
% %solves
% x = model.stateestimate(:, 2:end);
% i = 0;
% while true
%     [c, g] = costfcn(x, model);
%     disp(i);
%     disp(c);
%     i = i + 1;
%     Li = cfullLinv(x, model);
%     Lit = cfullLinvT(x, model);
%     dx = Li*Lit*(-g);
%     a = backtrack(x(:), dx, model);
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
% end

% %As of 3/6/20, this is beating LBFGS pretty convincingly by iteration
% %GN Hessian Test Block with full Hessian and line search
% %(solves H\-g to get search direction)
% %coarsest TLM, adj approximations (identity)
% x = model.stateestimate(:, 2:end);
% i = 0;
% while true
%     [c, g] = costfcn(x, model);
%     disp(i);
%     disp(c);
%     i = i + 1;
%     H = cfullHessian(x, model);
%     dx = H\(-g);
%     a = backtrack(x(:), dx, model);
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
% end

%%%%NOTE TO VISHWAS
%%%%This is the end of the tests we discussed today (3/9/20)

% %%%Verify the Hessian
% %%%As of 3/9/20, seems to check out
% A = eye(model.n*m);
% %xtest = xt + .001*randn(size(xt));
% %model.stateestimate = xtest;
% e = 1e-12;
% maxrelerr = 0;
% for i = 1:model.n*m
% j = randi([1, model.n*m], 1, 1);
% x = reshape(model.stateestimate(:, 2:end), model.n, m);
% v = A(:, j);
% v = reshape(v, size(x));
% g = reshape(gradfun(x, model), model.n, m);
% x = model.stateestimate(:, 2:end);
% [~, stepf] = costfcn(x + e * v, model);
% [~, stepb] = costfcn(x - e * v, model);
% fd = (stepf - stepb)./(2*e);
% %s = Linv(x, Linvtrans(x, v, model), model);
% H = fullHessian(x, model);
% s = H*v(:);
% relerr = norm(fd - s(:))/norm(fd);
% disp(j)
% disp(relerr)
% if relerr > maxrelerr
%    maxrelerr = relerr;
% end
% end
% disp("Max rel err")
% disp(maxrelerr);


% % %%GN Hessian Test Block using Linv, Linvtrans applications and line
% % search
% %xtest = xt(:, 2:end);
% x = model.stateestimate(:, 2:end);
% i = 0;
% while true
%     [c, g] = costfcn(x, model);
%     disp(i)
%     c
%     i = i + 1;
%     g = reshape(gradfun(x, model), model.n, m);
%     dx = Linv(x, Linvtrans(x, -g, model), model);
%     a = backtrack(x(:), dx(:), model);
%     %a = 1;
%     dx = reshape(dx, model.n, m);
%     x = x + a*dx;
%     %norm(x - xt(:, 2:end))
% end
% % % %

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


%
% % %%% Test block comparing Linv, Linvtrans, Hessian solves on random
% % %%% vectors
% %xtest = xt(:, 2:end);
% %x = xtest + 1*randn(size(xtest));
% %model.stateestimate = [model.x0, x];
% x = model.stateestimate(:, 2:end);
% while true
%      mine = 100*randn(size(x));
%      g = reshape(mine, model.n, m);
%      dx = Linv(x, Linvtrans(x, g, model), model);
%      H = fullHessian(x, model);
%      dx2 = H\mine(:);
%      norm(dx2(:) - dx(:))/norm(dx(:))
% end
% %

% x = model.stateestimate(:, 2:end);
% H = fullHessian(x, model);
% Li = fullLinv(x, model);
% %spy(L)
% Lit = fullLinvT(x, model);
% pause;
% mydim = size(H, 2);
% prod = zeros(size(H));
% %should be identity
% for i = 1:mydim
% h = reshape(H(:, i), size(x));
% col = Linv(x, Linvtrans(x, h, model), model);
% prod(:, i) = col(:);
% col(:)
% pause;
% end
% norm(prod)




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
function v = coarseadj(tspan, y, v, model)
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


%takes  timespan [t_0, t_1], initial state y, vector v that adjoint model over given
%parameters should be applied to
function v = coarsestadj(tspan, y, v, model)
%adjoint run
%Options  = MATLODE_OPTIONS('AbsTol',model.atol, 'RelTol', model.rtol, 'Lambda', @(t, y) v, 'Jacobian', model.jac, 'ADJ_AbsTol', model.atol, 'ADJ_RelTol', model.rtol);
%Options  = MATLODE_OPTIONS('AbsTol',model.catol, 'RelTol', model.crtol, 'Lambda', @(t, y) v, 'ADJ_AbsTol', model.catol, 'ADJ_RelTol', model.crtol, 'Jacobian', model.jac);
%[~, ~, v] = MATLODE_ERK_ADJ_Integrator(model.rhs, tspan, y, Options);
%MATLODE integration returns the initial value (y) and the integrated value
%as row vectors stacked on top of each other. So, we remove the top row and
%transpose before returning.
return;
end

%takes  timespan [t_0, t_1], initial state y, vector v that TLM over given
%parameters should be applied to
function v = coarsesttlm(tspan, y, v, model)
%tlm run
%Options = MATLODE_OPTIONS('Jacobian', model.jac, 'AbsTol', model.catol, 'RelTol', model.crtol,'Y_TLM', eye(model.n), 'TLM_AbsTol', model.catol, 'TLM_RelTol',model.crtol);
%SDIRK TLM seems to break on application to a single input  vector for some
%reason, so return full sensitivity matrix and apply to v
%[~, ~, sens] = MATLODE_ERK_TLM_Integrator(model.rhs, tspan, y, Options);
%v = sens*v;
return;
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
    %disp(j);
    lambda = v(:, j);
    for i = m-1:-1:j
        %disp(i);
        tspan = [times(i), times(i+1)];
        lambda = adjmodel(tspan, x(:, i), lambda + v(:, i), model);
    end
    lambda = lambda + v(:, j);
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
for j = 1:m-1
    xi = v(:, 1);
    for i = 1:j
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
    tspan = [model.times(i - 1), model.times(i)];
    mat = tlmmodel(tspan, x(:, i-1), eye(n), model);
    mat = -diag(1./d) * mat;
    bdiags{end + 1} = mat;
end
%build above diagonal blocks
for i = 2:m
    d = rMat(model.stateestimate(:, i), model);
    tspan = [model.times(i - 1), model.times(i)];
    mat = adjmodel(tspan, x(:, i-1), eye(n), model);
    mat = -mat * diag(1./d);
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
H((m - 1)*n + 1:m*n,(m-1)*n + 1:m*n) = diags{end};
%If the Hessian has a below diagonal block
try
    H((m - 1)*n + 1:m*n,(m-2)*n + 1:(m-1)*n) = bdiags{end};
end
%main loop builds all rows except first and last (i.e. those containing all
%3 blocks)
for i = 2:m-1
    H((i-1)*n + 1:i*n,(i-1)*n + 1:i*n) = diags{i};%diagonal blocks
    H((i-1)*n + 1:i*n,(i-2)*n + 1:(i-1)*n) = bdiags{i}; %below diagonal blocks
    H((i-1)*n + 1:i*n,i*n + 1:(i+1)*n) = adiags{i};
end
Hes = H;
end

%for use with finite difference verification routine
function G = gradfun(x, model)
[c, G] = costfcn(x, model);
return;
end



function a = backtrack(x, p, model)
r = .1;
c = 1e-8;
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
    %disp(a)
    
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


%%%%Return the full L^{-1} matrix (for testing)
%%%%Right multiply by sqrt(R) matrices at the end for scaling
function Lnv = fullLinv(x, model)
[n, m] = size(x); %model dim, free vectors
%each column cell will contain a cell array of that column's blocks
cols = [];
col = [];
for i = 1:m %columns of matrix
    col = [];
    for j = 1:i - 1
        col = [zeros(n, n); col];
    end
    col = [col; eye(n)];
    for j = i + 1:m %rows of column
        tspan = [model.times(i), model.times(j)];
        mat = tlmmodel(tspan, x(:, i), eye(n), model);
        col = [col; mat];
    end
    cols = [cols, col];
    col = [];
end

Lnv = cols;


%Multiply by scaling matrices on right
scale = [];
for i = 1:m %block rows
    d =  rMat(model.stateestimate(:, i), model);
    scale = [scale; sqrt(d)];
end
scale = diag(scale);
%spy(scale);

Lnv = cols*scale;
end


%%%%Return the full L^{-T} matrix (for testing)
%%%%Left multiply by sqrt(R) matrices for scaling
function Lnvt = fullLinvT(x, model)
[n, m] = size(x); %model dim, free vectors
rows = [];
row = [];
for i = 1:m %rows of matrix
    row = [];
    for j = 1:i - 1
        row = [zeros(n, n), row];
    end
    row = [row, eye(n)];
    for j = i + 1:m %rows of column
        tspan = [model.times(i), model.times(j)];
        mat = adjmodel(tspan, x(:, i), eye(n), model);
        row = [row, mat];
    end
    rows = [rows; row];
    row = [];
end

Lnvt = rows;


%Multiply by scaling matrices on right
scale = [];
for i = 1:m %block rows
    d =  rMat(model.stateestimate(:, i), model);
    scale = [scale; sqrt(d)];
end
scale = diag(scale);

Lnvt = scale*rows;
end

%%%%Return the full coarse L^{-1} matrix (for testing)
%%%%Right multiply by sqrt(R) matrices at the end for scaling
function Lnv = cfullLinv(x, model)
[n, m] = size(x); %model dim, free vectors
%each column cell will contain a cell array of that column's blocks
cols = [];
col = [];
for i = 1:m %columns of matrix
    col = [];
    for j = 1:i - 1
        col = [zeros(n, n); col];
    end
    col = [col; eye(n)];
    for j = i + 1:m %rows of column
        tspan = [model.times(i), model.times(j)];
        mat = coarsesttlm(tspan, x(:, i), eye(n), model);
        col = [col; mat];
    end
    cols = [cols, col];
    col = [];
end

Lnv = cols;


%Multiply by scaling matrices on right
scale = [];
for i = 1:m %block rows
    d =  rMat(model.stateestimate(:, i), model);
    scale = [scale; sqrt(d)];
end
scale = diag(scale);
%spy(scale);

Lnv = cols*scale;
end


%%%%Return the full coarse L^{-T} matrix (for testing)
%%%%Left multiply by sqrt(R) matrices for scaling
function Lnvt = cfullLinvT(x, model)
[n, m] = size(x); %model dim, free vectors
rows = [];
row = [];
for i = 1:m %rows of matrix
    row = [];
    for j = 1:i - 1
        row = [zeros(n, n), row];
    end
    row = [row, eye(n)];
    for j = i + 1:m %rows of column
        tspan = [model.times(i), model.times(j)];
        mat = coarsestadj(tspan, x(:, i), eye(n), model);
        row = [row, mat];
    end
    rows = [rows; row];
    row = [];
end

Lnvt = rows;


%Multiply by scaling matrices on right
scale = [];
for i = 1:m %block rows
    d =  rMat(model.stateestimate(:, i), model);
    scale = [scale; sqrt(d)];
end
scale = diag(scale);

Lnvt = scale*rows;
end

%full GN hessian using coarsest TLM, ADJ approximations (identity)
function Hes = cfullHessian(x, model)
diags = {}; %diagonal elements of hessian matrix
bdiags = {}; %below diagonal elements
adiags = {}; %above diagonal elements
m = size(x, 2); %m is the number of free (vectors) along the trajectory
n = size(x, 1); %n is the dimension of state vectors
H = sparse(n*m, n*m);
%build diagonal blocks
for i = 1:m
    tspan = [model.times(i), model.times(i+1)];
    mat = coarsesttlm(tspan, x(:, i), eye(n), model);
    d = rMat(model.stateestimate(:, i), model);
    mat = (1./d).*mat;
    mat = coarsestadj(tspan, x(:,i), mat, model);
    mat = diag(1./d) + mat;
    diags{end+1} = mat;
end
%build below diagonal blocks
for i = 2:m
    d = rMat(model.stateestimate(:, i), model);
    tspan = [model.times(i - 1), model.times(i)];
    mat = coarsesttlm(tspan, x(:, i-1), eye(n), model);
    mat = -diag(1./d) * mat;
    bdiags{end + 1} = mat;
end
%build above diagonal blocks
for i = 2:m
    d = rMat(model.stateestimate(:, i), model);
    tspan = [model.times(i - 1), model.times(i)];
    mat = coarsestadj(tspan, x(:, i-1), eye(n), model);
    mat = -mat * diag(1./d);
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
H((m - 1)*n + 1:m*n,(m-1)*n + 1:m*n) = diags{end};
%If the Hessian has a below diagonal block
try
    H((m - 1)*n + 1:m*n,(m-2)*n + 1:(m-1)*n) = bdiags{end};
end
%main loop builds all rows except first and last (i.e. those containing all
%3 blocks)
for i = 2:m-1
    H((i-1)*n + 1:i*n,(i-1)*n + 1:i*n) = diags{i};%diagonal blocks
    H((i-1)*n + 1:i*n,(i-2)*n + 1:(i-1)*n) = bdiags{i}; %below diagonal blocks
    H((i-1)*n + 1:i*n,i*n + 1:(i+1)*n) = adiags{i};
end
Hes = H;
end

