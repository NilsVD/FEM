function [x,elem,nel,nnp,ndm,ndf,nen,matparam,drlt,robin,loadcurve,Q]=input_Rocket_Engine()

% import mesh:
x = dlmread('Rocket_Engine.nodes'); % 1st colum containes node number
x = x(:,2:end);                     % remove 1st column

conn = dlmread('Rocket_Engine.conn'); % 1st colum containes element number
conn = conn(:,2:end);                 % remove 1st column
[nel,nen] = size(conn);

% number of node points and number of spatial dimensions
[nnp, ndm] = size(x);

% number of degrees of freedom per node
ndf = 1;  % temperature is a scalar!

% allocate memory for the element structure 
% 1. organize connectivity list as a struct
% 2. determine edof, element degrees of freedom
% 3. reserve memory for state variables
%    gradT  [gradT_x    gradT_y] and 
%    q      [    q_x        q_y]
% number of state variables
nsv = 2 + 2;

elem = repmat(struct('cn',zeros(nen,1), ...
                     'edof', zeros(nen*ndf,1), ...
                     'stateVar0', zeros(1,nsv), ...
                     'stateVar', zeros(1,nsv) ...
                    ),  nel, 1 );

for e = 1:nel
    elem(e).cn = conn(e,:);
    elem(e).edof = elem(e).cn; 
end

% material parameter
% heat conductivity
matparam(1) = 22;

% thickness of the plate (m)
matparam(2) = 0.001;

% density of titanium
matparam(3) = 4540;

% mass speci?c heat capacity of titanium 
matparam(4) = 520;

% Initial boundary temperature
matparam(5) = 20;


%-----------------------------------------------------------------------
% load curves
%-----------------------------------------------------------------------
id = 1;
loadcurve(id).id = id;
loadcurve(id).scalefactor =  1.0;
loadcurve(id).time =  [0.0   1.0];
loadcurve(id).value = [1.0   1.0];

%Retrieves all boundary nodes for the different boundaries.
[NodesAN,NodesAE,NodesEF,NodesFI,NodesIJ,NodesJN] = BCNodes();

% boundary conditions
% node: Node where condition is applied
% loadcurveID: id of load curve
% ldof: load degree of freedom (1 for heat, 1 or 2 for mechanical)
% scale: magnitude of bc

% Dirichlet boundary condition on boundaries A-N
drlt = repmat(struct('nodes',zeros(1,length(NodesAN)),'loadcurveID',0,'ldof',0,'scale',0),1,1);
drlt.nodes = NodesAN;
drlt.loadcurveID = 1;
drlt.ldof = 1;
drlt.scale = -150;

%Robin boundary conditions on all remaining boundaries
robin = repmat(struct('nodes',zeros(1,length(NodesFI)),'alpha',0,'Tinf',0),5,1);

%Fills in all different conditions:
%A-E boundary
robin(1).nodes = NodesAE;
robin(1).alpha = 20;     
robin(1).Tinf = 20; 

%E-F boundary
robin(2).nodes = NodesEF;
robin(2).alpha = 20;     
robin(2).Tinf = 20; 

%F-I boundary
robin(3).nodes = NodesFI;
robin(3).alpha = 120;     
robin(3).Tinf = 2000;     

%I-J boundary
robin(4).nodes = NodesIJ;
robin(4).alpha = 20;     
robin(4).Tinf = 20; 

%J-N boundary
robin(5).nodes = NodesJN;
robin(5).alpha = 20;     
robin(5).Tinf = 20; 

% constant heat source 
Q = 0;

end