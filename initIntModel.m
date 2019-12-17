function [model] = initIntModel()
model.n = 2; %model dimension
model.jac = @(t, y) fulljacobian(1, y, 5, 5); %RHS Jacobian
model.rhs = @(t, y) fullf(1, y, 5, 5); %ODE RHS
tspan = [0, 1]; %full integration timespan
model.tspan = tspan;
M = 3; %subintervals to divide tspan into. If you set M=1 it will crash.
model.times = linspace(tspan(1), tspan(2), M);
model.atol = 1e-6; %integration tolerances
model.rtol = 1e-6; 
model.x0 = randn(model.n, 1); %specify initial values
model.stateestimate = zeros(model.n, M); %store estimate of states before optimization
                                    %(used to build R_i scaling matrices,
                                    %fixed before optimization).
end