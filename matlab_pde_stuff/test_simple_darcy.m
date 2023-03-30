% Darcy flow with realistic parameters
%
%   u(1,:) = vel_x
%   u(2,:) = vel_y
%   u(3,:) = p
%   u(4,:) = C
%
%   u1:     (ux) = - K * grad(p,x)
%   u2:     (uy) = - K * grad(p,y)
%   u3:     grad( K * grad(p) ) = 0
%   u4:     eps*d(u)/dt - grad * (eps * D * grad(u)) = -eps*ux*grad(u,x) - eps*uy*grad(u,y)
%                                              -eps*u*(grad(ux) + grad(uy))

% Constants
p_exit = 101350;   % Pa
ux_inlet = 2;      % m/s
C_input = 0.001;   % mol/m^3
n_steps = 40;
end_time = 2;      % s

% System params
t_span = linspace(0,end_time,n_steps);  % s
N = 4; % num PDEs 

% Object and calculated constants
obj = porous_flow_adsorption_chemistry();
mu = obj.ViscosityWater();  % kg/m/s
D = obj.EffectiveDispersionWater(1, 298, ux_inlet);  % m^2/s
eps = obj.bulk_porosity;
K = obj.KozenyCarmannConst();        % m^3*s/kg


% PDE coefficients
dtsqrd = 0;
dtcoeff = zeros(N,1);

% NOTE: Only the tracer has a time coeff
dtcoeff(4,1) = eps;

%              u(1)    u(2)  u(3)  u(4)
diff_tensor = [mu;mu; mu;mu; K;K; eps*D;eps*D];

rxn_coeff = [1;1;0;0];

nonlin_fun= @(location,state) [ -K*state.ux(3,:);
                               -K*state.uy(3,:);
                               0*state.u(3,:);
                       -eps*state.u(1,:).*state.ux(4,:) - eps*state.u(2,:).*state.uy(4,:) - eps*state.u(4,:).*(state.ux(1,:) + state.uy(2,:))];

% Info for specifying functions
% -----------------------------
% location is a structure with these fields:
% 
%     location.x
% 
%     location.y
% 
%     location.z
% 
%     location.subdomain
% 
% The fields x, y, and z represent the x-, y-, and z- coordinates of points 
% for which your function calculates coefficient values. The subdomain field 
% represents the subdomain numbers, which currently apply only to 2-D models. 
% The location fields are row vectors.
% 
% state is a structure with these fields:
% 
%     state.u
% 
%     state.ux
% 
%     state.uy
% 
%     state.uz
% 
%     state.time
% 
% The state.u field represents the current value of the solution u. The 
% state.ux, state.uy, and state.uz fields are estimates of the solution’s 
% partial derivatives (∂u/∂x, ∂u/∂y, and ∂u/∂z) at the corresponding points 
% of the location structure. The solution and gradient estimates are N-by-Nr 
% matrices. The state.time field is a scalar representing time for 
% time-dependent models.

% NOTE: There is no acess to derivatives in BCs


% --------- Define the geometry ------------

% Create geometry for the pde
%   2D geometry 

% Rectangle is code 3, 4 sides,
% followed by x-coordinates and then y-coordinates
R1 = [3,4,  -1,1,1,-1,  -.4,-.4,.4,.4]';
% Circle is code 1, center (.5,0), radius .2
C1 = [1,.5,0,.2]';
% Pad C1 with zeros to enable concatenation with R1
C1 = [C1;zeros(length(R1)-length(C1),1)];
geom = [R1,C1];

% Names for the two geometric objects
ns = (char('R1','C1'))';

% Set formula
sf = 'R1 - C1';

% Create geometry
geo = decsg(geom,sf,ns);


% ------------- PDE Model Objects ---------------
% Create the PDE model
model = createpde(N);

% Include the geometry in the model
geometryFromEdges(model,geo);

% generate the mesh 
generateMesh(model,"Hmax",0.05);

% Specify what the model coefficients are
specifyCoefficients(model,"m",dtsqrd,"d",dtcoeff,"c",diff_tensor,"a",rxn_coeff,"f",nonlin_fun);


% -------------- Set BCs ---------------------
% BC Formats
%
%       Dirichlet:  h*u=r
%
%       Neumann:    n * (c * grad(u)) = g - q*u

% For "Mixed" BCs with multiple coupled variables, you have to supply BC
%   parameters in a Matrix format (https://www.mathworks.com/help/pde/ug/steps-to-specify-a-boundary-conditions-object.html)

% Edge 3 = outlet
% Edge 1 = inlet
% Edge 2, 4:8 = walls

% Apply BCs
% Exit --> u(1) and u(2) are neumann with zero slope (g=0, q=0)
%          u(3) is dirichlet with value pinned to p_exit (h=1, r=0)
%          u(4) is neumann with zero slope (g=0, q=0)
he = [0 0 0 0; 
     0 0 0 0; 
     0 0 1 0;
     0 0 0 0];

re = [0;0;p_exit;0];

qe = [0 0 0 0; 
     0 0 0 0; 
     0 0 0 0;
     0 0 0 0];

ge = [0;0;0;0];

applyBoundaryCondition(model,"mixed", ...
                             "Edge",3, ...
                             "h",he,"r",re,"q",qe,"g",ge);


% Walls
%       u(1) and u(2) are dirichlet with value pinned to 0 (h=1, r=0)
%       u(3) is neumann with zero slope (g=0, q=0)
%       u(4) is neumann with zero slope (g=0, q=0)
hw = [1 0 0 0; 
     0 1 0 0; 
     0 0 0 0;
     0 0 0 0];

rw = [0;0;0;0];

qw = @(region,state) [0 0 0 0; 
                     0 0 0 0; 
                     0 0 0 0;
                     0 0 0 0];

gw =  [0; 0; 0;0];

applyBoundaryCondition(model,"mixed", ...
                             "Edge",[2,4:8], ...
                             "h",hw,"r",rw,"q",qw,"g",gw);

% Inlet BC
%       u(1) and u(2) are neumann with zero slope (g=0, q=0)
%       u(3) is neumann with slope equal to ux_inlet
%       u(4) is dirichlet with value r=C_input
%               -- or --
%       u(4) is neumann/cauchy with g = -(n*v)*(u_input-u)
hi = [0 0 0 0; 
     0 0 0 0; 
     0 0 0 0;
     0 0 0 1];

%hi = [0 0 0 0; 
%     0 0 0 0; 
%     0 0 0 0;
%     0 0 0 0];

ri = [0;0;0;C_input];

%ri = [0;0;0;0];

qi = [0 0 0 0; 
     0 0 0 0; 
     0 0 0 0;
     0 0 0 0];

% Use the 1.255 as a correction factor for the pressure at inlet 
gi = @(region,state) [0;0; ux_inlet*1.255; 0];

% NOTE: Cauchy BC possible, but much slower
%gi = @(region,state) [0;0; ux_inlet*1.255; ux_inlet*(C_input-state.u(4,:))];

applyBoundaryCondition(model,"mixed", ...
                             "Edge",1, ...
                             "h",hi,"r",ri,"q",qi,"g",gi);


% Set initial conditions 
u0 = [0;0;p_exit;0];
setInitialConditions(model,u0);


% Solve the model
model.SolverOptions.ReportStatistics = 'on';
model.SolverOptions.AbsoluteTolerance = 1e-4; % ODE opt
model.SolverOptions.RelativeTolerance = 1e-4; % ODE opt
model.SolverOptions.ResidualTolerance = 1e-6; % Nonlinear opt
model.SolverOptions.MaxIterations = 30;       % Nonlinear opt
model.SolverOptions.MinStep = 0.001;          % Min step size 
model.SolverOptions.ResidualNorm = 2;         % L-2 norm
model.SolverOptions.MaxShift = 500;           % Lanczos solver shift
model.SolverOptions.BlockSize = 50;           % Block size for Lanczos recurrence


results = solvepde(model,t_span);
u = results.NodalSolution;


% plot single state flow pressure and flow
f1 = figure;
pdeplot(model,"XYData",u(:,1,end))
f2 = figure;
pdeplot(model,"XYData",u(:,2,end))
f3 = figure;
pdeplot(model,"XYData",u(:,3,end),FlowData=[u(:,1,end),u(:,2,end)], ColorMap="jet", FaceAlpha=0.5)
f4 = figure;
pdeplot(model,"XYData", sqrt( u(:,1,end).^2 + u(:,2,end).^2),Mesh="on", ColorMap="jet")
%f5 = figure;
%pdeplot(model,FlowData=[u(:,1,end),u(:,2,end)])

% Plot gif for tracer results
umax = max(max(u(:,4,:)));
umin = min(min(u(:,4,:)));
delete("tracer_flow.gif")
for i = 1:n_steps
    p = pdeplot(model,"XYData",u(:,4,i),"Mesh","off");
    axis([-1 1 -1 1 umin umax]); 
    clim([umin umax]);
    xlabel x
    ylabel y
    zlabel u
    M(i) = getframe;
    exportgraphics(gca, "tracer_flow.gif", "Append",true);
end


% Interpolate solution on a line
%
%   line from point (-1,-.4) to (-1, .4)  [exit]
%                   (x1, y1) to (x2, y2)
xq = -1*ones(1,10);
yq = linspace(-.4,.4,10);

C_intrp = interpolateSolution(results,xq,yq,4,1:length(t_span));
ux_intrp_exit = interpolateSolution(results,xq,yq,1,length(t_span));
uy_intrp_exit = interpolateSolution(results,xq,yq,2,length(t_span));

%   line from point (-1,-.4) to (-1, .4)  [exit]
%                   (x1, y1) to (x2, y2)
xq = 1*ones(1,10);
yq = linspace(-.4,.4,10);

ux_intrp_enter = interpolateSolution(results,xq,yq,1,length(t_span));
uy_intrp_enter = interpolateSolution(results,xq,yq,2,length(t_span));

% C_intrp will be a 10x1x50 (SpacePoints x Vars x TimePoints) structure 
breakthrough = zeros(1,1,length(t_span));
for i=1:length(t_span)
    breakthrough(1,1,i) = mean(C_intrp(:,1,i));
end
times = results.SolutionTimes;
C_out = breakthrough(1,1,:);
C_out = reshape(C_out,[1,length(t_span)]);

f5 = figure;
plot(times,C_out)
xlabel('Time (s)')
ylabel('Exit Concentration (mol/m^3)')

% Check input and output fluxes
avg_ux_input = mean(ux_intrp_enter);
avg_ux_output = mean(ux_intrp_exit);

avg_uy_input = mean(uy_intrp_enter);
avg_uy_output = mean(uy_intrp_exit);
