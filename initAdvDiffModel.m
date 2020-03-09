function [model] = initAdvDiffModel()
model.nx = 100; %model dimension
model.ny = 100;
model.nz = 100;
model.ntracer = 1;

tspan = [0, 1]; %full integration timespan
model.tspan = tspan;
model.dt = model.tspan(2) - model.tspan(1);

M = 5;
model.M = M; %subintervals to divide tspan into. If you set M=1 it will crash.
model.times = linspace(tspan(1), tspan(2), M);
model.dts = diff(model.times);

model.x0full = randn(model.nx, model.ny, model.nz, model.ntracer); %specify initial values
model.stateestimate = zeros(model.nx, model.ny, model.nz, model.ntracer, M);
%store estimate of states before optimization
%(used to build R_i scaling matrices,
%fixed before optimization).

%wind tensors
model.wind.u = randn(model.nx, model.ny, model.nz);
model.wind.v = randn(model.nx, model.ny, model.nz);
model.wind.w = randn(model.nx, model.ny, model.nz);

%diffusion coefficients
model.kh = randn(model.nx, model.ny, model.nz);
model.kv = randn(model.nx, model.ny, model.nz);

%boundary conditions
model.bc.w = zeros(model.ny, model.nz, model.ntracer);
model.bc.e = zeros(model.ny, model.nz, model.ntracer);
model.bc.s = zeros(model.nx, model.nz, model.ntracer);
model.bc.n = zeros(model.nx, model.nz, model.ntracer);
model.bc.top = zeros(model.nx, model.ny, model.ntracer);

%spatial step sizes
model.dx = 1000;
model.dy = 1000;
model.dz = 1000;

%integration tolerances
model.atol = 1e-6;
model.rtol = 1e-6;
end