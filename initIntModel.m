function [model] = initIntModel()
model.jac = @(t, y) lorenz96_jacobian(t, y); %RHS Jacobian
model.rhs = @(t, y) lorenz96_rhs(t, y); %ODE RHS
tspan = [0, .1]; %full integration timespan
model.tspan = tspan;
M = 3; %subintervals to divide tspan into
model.times = linspace(tspan(1), tspan(2), M);
model.atol = 1e-3; %integration tolerances
model.rtol = 1e-3; 
model.x0 = randn(40, 1); %specify initial values
model.stateestimate = zeros(40, M); %store estimate of states before optimization
                                    %(used to build R_i scaling matrices,
                                    %fixed before optimization).
end