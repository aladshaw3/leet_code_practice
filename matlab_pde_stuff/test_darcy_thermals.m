% Darcy flow with thermal balances
%
%   u(1,:) = vel_x
%   u(2,:) = vel_y
%   u(3,:) = p
%   u(4,:) = T
%   u(5,:) = Tc
%
%   u1:     (ux) = - K * grad(p,x)
%   u2:     (uy) = - K * grad(p,y)
%   u3:     grad( K * grad(p) ) = 0
%   u4:     εb ρ cpg ∂T/∂t + εb ρ cpg v ∂T/∂z = -(1-εb) Ga hc (T - Tc) [ - εb α hwg (T - Tw) ]
% 
%   u5:     (1-εb) ρc cpc ∂Tc/∂t = (1-εb) Kc ∂2Tc/∂z2 + (1-εb) Ga hc (T - Tc) [ - (1-εb) α hwc (Tc - Tw) ] 

% Create properties object 
prop = porous_flow_adsorption_chemistry();

% Constants
p_exit = 101350;   % Pa
ux_inlet = 0.25;      % m/s
T_input = 298;     % K
T_init = 298;      % K
n_steps = 50;
end_time = 5;      % s

% System params
t_span = linspace(0,end_time,n_steps);  % s
N = 5; % num PDEs 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rxn_coeff = [1;1;0;0;0];
dtcoeff = @(location,state) dcoeffunction(location,state,prop);
diff_tensor = @(location,state) ccoeffunction(location,state,prop);
nonlin_fun = @(location,state) fcoeffunction(location,state,prop);

% ------------- PDE Model Objects ---------------
% Create the PDE model
model = createpde(N);

% Include the geometry in the model
geometryFromEdges(model,geo);

% generate the mesh 
generateMesh(model,"Hmax",0.05,"GeometricOrder","linear");

% Specify what the model coefficients are
specifyCoefficients(model,"m",0,"d",dtcoeff,"c",diff_tensor,"a",rxn_coeff,"f",nonlin_fun);


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
%          u(4) and u(5)is neumann with zero slope (g=0, q=0)
he = [0 0 0 0 0; 
     0 0 0 0 0; 
     0 0 1 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

re = [0;0;p_exit;0;0];

qe = [0 0 0 0 0; 
     0 0 0 0 0; 
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

ge = [0;0;0;0;0];

applyBoundaryCondition(model,"mixed", ...
                             "Edge",3, ...
                             "h",he,"r",re,"q",qe,"g",ge);


% Walls
%       u(1) and u(2) are dirichlet with value pinned to 0 (h=1, r=0)
%       u(3) is neumann with zero slope (g=0, q=0)
%           ---Temp--
%       u(4) and u(5) is neumann with zero slope (g=0, q=0)
hw = [1 0 0 0 0; 
     0 1 0 0 0; 
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

rw = [0;0;0;0;0];

qw = @(region,state) [0 0 0 0 0; 
                     0 0 0 0 0; 
                     0 0 0 0 0;
                     0 0 0 0 0;
                     0 0 0 0 0];

Twall = 470;
scale = 1;
hwall = prop.WallHeatTransferWater(p_exit, T_init, ux_inlet, 0, 0)/1e6*scale;
area_bulk = prop.bulk_porosity;
area_solid = (1-prop.bulk_porosity);
gw =  @(region,state) [0; 0; 0; area_bulk*hwall*(Twall-state.u(4,:)); area_solid*hwall*(Twall-state.u(5,:))];

applyBoundaryCondition(model,"mixed", ...
                             "Edge",[2,4:8], ...
                             "h",hw,"r",rw,"q",qw,"g",gw);

% Inlet BC
%       u(1) and u(2) are neumann with zero slope (g=0, q=0)
%       u(3) is neumann with slope equal to ux_inlet
%
%           --- Temp ---
%       u(4) is dirichlet with value r=C_input
%               -- or --
%       u(4) is neumann/cauchy with g = -(n*v)*(u_input-u)
%
%       u(5) is neumann with zeros
hi = [0 0 0 0 0; 
     0 0 0 0 0; 
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

%hi = [0 0 0 0 0; 
%     0 0 0 0 0; 
%     0 0 0 0 0;
%     0 0 0 0 0;
%     0 0 0 0 0];

ri = [0;0;0;0;0];

%ri = [0;0;0;0;0];

qi = [0 0 0 0 0; 
     0 0 0 0 0; 
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

% Use the 1.255 as a correction factor for the pressure at inlet 
%gi = @(region,state) [0;0; ux_inlet*1.255; 0; 0];

% NOTE: Cauchy BC possible, but much slower
%gi = @(region,state) [0;0; ux_inlet*1.255; εb ρ cpg v*(T_input-state.u(4,:)); 0];
gi = @(region,state) [0;0; ux_inlet*1.255; prop.bulk_porosity*prop.DensityWater().*prop.SpecificHeatWater()/ 1e6*ux_inlet*(T_input-state.u(4,:)); 0];

applyBoundaryCondition(model,"mixed", ...
                             "Edge",1, ...
                             "h",hi,"r",ri,"q",qi,"g",gi);


% Set initial conditions 
u0 = [0;0;p_exit;T_init;T_init];
setInitialConditions(model,u0);


% Solve the model
model.SolverOptions.ReportStatistics = 'on';
model.SolverOptions.AbsoluteTolerance = 1e-4; % ODE opt
model.SolverOptions.RelativeTolerance = 1e-4; % ODE opt
model.SolverOptions.ResidualTolerance = 1e-6; % Nonlinear opt
model.SolverOptions.MaxIterations = 30;       % Nonlinear opt
model.SolverOptions.MinStep = 0.0001;          % Min step size 
model.SolverOptions.ResidualNorm = 2;         % L-2 norm
model.SolverOptions.MaxShift = 500;           % Lanczos solver shift
model.SolverOptions.BlockSize = 50;           % Block size for Lanczos recurrence


results = solvepde(model,t_span);
u = results.NodalSolution;



% plot single state flow pressure and flow
f1 = figure;
pdeplot(model,"XYData",u(:,1,end))
axis([-1 1 -1 1]); 
f2 = figure;
pdeplot(model,"XYData",u(:,2,end))
axis([-1 1 -1 1]); 
f3 = figure;
pdeplot(model,"XYData",u(:,3,end),FlowData=[u(:,1,end),u(:,2,end)], ColorMap="jet", FaceAlpha=0.5)
axis([-1 1 -1 1]); 

% Plot gif for tracer results
f4 = figure;
umax = max(max(u(:,4,:)));
umin = min(min(u(:,4,:)));
delete("fluid_heat_flow.gif")
for i = 1:n_steps
    p = pdeplot(model,"XYData",u(:,4,i),"Mesh","off");
    axis([-1 1 -1 1]); 
    clim([umin umax]);
    xlabel x
    ylabel y
    zlabel u
    Ma(i) = getframe;
    exportgraphics(gca, "fluid_heat_flow.gif", "Append",true);
end

% % Plot gif for tracer results
f5 = figure;
umax = max(max(u(:,5,:)));
umin = min(min(u(:,5,:)));
delete("solid_heat_flow.gif")
for i = 1:n_steps
    p = pdeplot(model,"XYData",u(:,5,i),"Mesh","off");
    axis([-1 1 -1 1]); 
    clim([umin umax]);
    xlabel x
    ylabel y
    zlabel u
    Mb(i) = getframe;
    exportgraphics(gca, "solid_heat_flow.gif", "Append",true);
end




function cmatrix = ccoeffunction(location,state,prop)
    % The C matrix is an N1xNr matrix where N1 is the number of
    % coefficients and Nr is the number of locations.
    %
    %   Nr always equals numel(location.x)
    %
    %   N1 can vary depending on the form of the C matrix and
    %   the number of pdes being solved.
    %
    %   For 5 pdes in 2D, we can do just 5 for N1 (representing isotropic
    %   behaviors). 
    n1 = 5;
    nr = numel(location.x);
    cmatrix = zeros(n1,nr);
    % mu 
    cmatrix(1,:) = prop.ViscosityWater(state.u(4,:));
    % mu
    cmatrix(2,:) = cmatrix(1,:);
    % K
    cmatrix(3,:) = prop.KozenyCarmannDarcyCoeffient(state.u(4,:));
    % eb*Kw
    cmatrix(4,:) = prop.bulk_porosity*prop.ThermalConductivityWater(state.u(4,:))/50;
    % (1-eb)*Kc
    cmatrix(5,:) = (1-prop.bulk_porosity)*prop.bulk_solids_cond/50;
end

function fmatrix = fcoeffunction(location,state,prop)
    n1 = 5;
    nr = numel(location.x);
    fmatrix = zeros(n1,nr);

    % pressure grad x
    fmatrix(1,:) = -prop.KozenyCarmannDarcyCoeffient(state.u(4,:)).*state.ux(3,:);

    % pressure grad y
    fmatrix(2,:) = -prop.KozenyCarmannDarcyCoeffient(state.u(4,:)).*state.uy(3,:);

    % Zeros

    % -(1-εb) Ga hc (T - Tc) - εb ρ cpg v ∂T/∂z
    transfer = -(1-prop.bulk_porosity)*prop.SphericalSurfaceVolumeRatio()*...
        prop.SolidsHeatTransferWater(state.u(3,:), state.u(4,:), 2, 0, 0).*...
            (state.u(4,:)-state.u(5,:)) / 1e6;
    advcoeff = prop.bulk_porosity*prop.DensityWater(state.u(3,:), state.u(4,:)).*prop.SpecificHeatWater(state.u(4,:)) / 1e6;
    fmatrix(4,:) = transfer ...
        -advcoeff.*state.u(1,:).*state.ux(4,:) - advcoeff.*state.u(2,:).*state.uy(4,:) ...
            - advcoeff.*state.u(4,:).*(state.ux(1,:) + state.uy(2,:));

    % (1-εb) Ga hc (T - Tc)
    fmatrix(5,:) = -transfer;

end

function dmatrix = dcoeffunction(location,state,prop)
    n1 = 5;
    nr = numel(location.x);
    dmatrix = zeros(n1,nr);

    % 4: εb ρ cpg
    dmatrix(4,:) = prop.bulk_porosity*prop.DensityWater(state.u(3,:), state.u(4,:)).*prop.SpecificHeatWater(state.u(4,:)) / 1e6;

    % 5: (1-εb) ρc cpc
    dmatrix(5,:) = (1-prop.bulk_porosity)*prop.bulk_solids_dens*prop.bulk_solids_Cp / 1e6;
end