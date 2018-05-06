%%  FEM assignment 1

clear all
close all

[x,elem,nel,nnp,ndm,ndf,nen,matparam, drlt, neum, robin_1, robin_2, loadcurve, Q] = input_Rocket_Engine();

plotTriaMesh(elem, x);

% separate dof types
allDofs = (1:1:nnp*ndf)';

% dofs with Dirichlet boundary conditions
numDrltDofs = length(drlt.nodes);
drltDofs = zeros(numDrltDofs,1);
for i=1:numDrltDofs
    node = drlt.nodes(i);
    ldof = drlt.ldof;
    drltDofs(i) = (node-1)*ndf + ldof;
end

%     for i=1:length(drlt_1.nodes)
%         node = drlt_1.nodes(i);
%         ldof = drlt_1.ldof;
%         drltDofs(j) = (node-1)*ndf + ldof;
%         j = j + 1;
%     end
% 
%     for i=1:length(drlt_2.nodes)
%         node = drlt_2.nodes(i);
%         ldof = drlt_2.ldof;
%         drltDofs(j) = (node-1)*ndf + ldof;
%         j = j + 1;
%     end

% free dofs
freeDofs = setdiff(allDofs, drltDofs);

% dofs with Neumann boundary conditions
numNeumDofs = length(neum.nodes);
neumDofs = zeros(numNeumDofs,1);
for i=1:numNeumDofs
    node = neum.nodes(i);
    ldof = neum.ldof;
    neumDofs(i) = (node-1)*ndf + ldof;
end

% dofs with Robin boundary conditions
numRobinDofs = length(robin_1.nodes)+length([robin_2(1).nodes,robin_2(2).nodes,robin_2(3).nodes,robin_2(4).nodes]);
RobinDofs = zeros(numRobinDofs,1);
j = 1;
for i=1:length(robin_1.nodes)
    node = robin_1.nodes(i);
    ldof = 1;

    RobinDofs(j) = (node-1)*ndf + ldof;
    j = j + 1;
end
for k=1:4
    for i=1:length(robin_2(k).nodes)
        node = robin_2(k).nodes(i);
        ldof = 1;

        RobinDofs(j) = (node-1)*ndf + ldof;
        j = j + 1;
    end
end

%=========================================================================
% fe-analysis
%=========================================================================

% unknown node temperatures
a = zeros(ndf*nnp,1);

for i=1:length(a)
    a(i) = 20;
end
%a = [20;0;0;0];

% simulation time and time/load step
stopTime = 50.0;
dt = 0.1;
% numerical tolerance to detect equilibrium
TOL = 1e-10;
%TOL = 2;

% loop over all time/load steps
time = 0;
step = 0;
while (time < stopTime)
    time = time + dt;
    step = step + 1;
    
    fprintf(1,'step= %2d time= %8.4e dt= %8.4e\n', step, time, dt);
    % Neumann loads for this timestep
    fsur = zeros(ndf*nnp,1);
    numNeumDofs = size(neum,1);
    neumValue = zeros(numNeumDofs,1);
    for i=1:numNeumDofs
        lcID = neum.loadcurveID;
        scale = neum.scale * loadcurve(lcID).scalefactor;
        neumValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    %fsur(neumDofs) = neumValue;
    
    RobinValue = zeros(numRobinDofs,1);
    T_0 = matparam(5);
    j = 1;
    for i=1:length(robin_1.nodes)
        q_0 = (-1)*robin_1.alpha*(a(robin_1.nodes(i))-robin_1.Tinf);
        %q_0 = (-1)*robin_1.alpha*(robin_1.Tinf);
        RobinValue(j) = q_0*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
        j = j + 1;
    end
    
    for k=1:4
        for i=1:length(robin_2(k).nodes)
            q_0 = (-1)*robin_2(k).alpha*(a(robin_2(k).nodes(i))-robin_2(k).Tinf);
            %q_0 = (-1)*robin_2(k).alpha*(robin_2(k).Tinf);
            RobinValue(j) = q_0*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
            j = j + 1;
        end
    end
    fsur(RobinDofs) = RobinValue;
    
    % Dirichlet boundary conditions
    drltValue = zeros(numDrltDofs,1);
    for i=1:numDrltDofs
        lcID = drlt.loadcurveID;
        scale = drlt.scale * loadcurve(lcID).scalefactor;
        drltValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    drltValueIncr = drltValue - a(drltDofs);
    
    % always run through the while loop once to check equilibrium
    rsn = 1.0; % residuum norm
    
    iter = 0;
    while (rsn > TOL && iter < 2)
        
        % initialise global internal force vector
        fint = zeros(ndf*nnp,1);
        
        % initialise global volume force vector
        fvol = zeros(ndf*nnp,1);
        
        % initialise global tangent stiffness matrix
        K = zeros(ndf*nnp,ndf*nnp);
        
        % initialise global residuum
        rsd = zeros(ndf*nnp,1);
        
        % element loop
        for e=1:nel
            
            % dofs belonging to element e
            edof = elem(e).edof;
            
            % coordinates of the element nodes
            ex = x(elem(e).cn,1);   % x-coord
            ey = x(elem(e).cn,2);   % y-coord
            
            % temperature values at the element nodes
            ae = a(edof);
                        
            % element volume force vector
            %fvole = zeros(size(edof));
            
            % area of the element 
            Ae=1/2*det([ones(3,1) ex ey]);
  
            % B-matri8ex for a triangular element with 3 nodes
            %B = [X(elem(e).cn(2),2)-X(elem(e).cn(3),2)  X(elem(e).cn(3),2)-X(elem(e).cn(1),2)   X(elem(e).cn(1),2)-X(elem(e).cn(2),2); 
            %     X(elem(e).cn(3),1)-X(elem(e).cn(2),1)  X(elem(e).cn(1),1)-X(elem(e).cn(3),1)   X(elem(e).cn(2),1)-X(elem(e).cn(1),1)] * 1/(2*Ae);
            B = [ey(2)-ey(3)  ey(3)-ey(1)   ey(1)-ey(2); 
                 ex(3)-ex(2)  ex(1)-ex(3)   ex(2)-ex(1)] * 1/(2*Ae);
            
            ex_mean = (ex(1)+ex(2)+ex(3))/3;
            ey_mean = (ey(1)+ey(2)+ey(3))/3;
            
            Ne = 1/(2*Ae)*[ex(2)*ey(3)-ex(3)*ey(2)+(ey(2)-ey(3))*ex_mean+(ex(3)-ex(2))*ey_mean, ...
                ex(3)*ey(1)-ex(1)*ey(3)+(ey(3)-ey(1))*ex_mean+(ex(1)-ex(3))*ey_mean, ...
                ex(1)*ey(2)-ex(2)*ey(1)+(ey(1)-ey(2))*ex_mean+(ex(2)-ex(1))*ey_mean];
            % temperature gradient
            gradT = B*ae;
            
            % temperature
            T = Ne*ae;
                
            % Fourier's law
            [q, D, elem(e).stateVar(1,:)] = fourierLaw2d(gradT, matparam, elem(e).stateVar0(1,:));
            
            % thickness of the plate 
            t = matparam(2);
            
            % element internal force vector 
            finte = B'* q * Ae * t;
            %finte = [0,0,0]';
            
            % element volume force vector 
            % fvole = 
            
            % --- thermal conductivity ----
            k = matparam(1);
            
            % element stiffness matrix 
            Ke = k*(B')*B*t*Ae*matparam(3);
            
            alpha = 0;
            T_inf = 0;
            BoolRobin_1 = ismember(elem(e).cn,robin_1.nodes);
            BoolRobin_2 = 0;
            %BoolRobin_2 = ismember(elem(e).cn, [robin_2(1).nodes,robin_2(2).nodes,robin_2(3).nodes,robin_2(4).nodes]);
            for i=1:4
                if(sum(ismember(elem(e).cn,robin_2(i).nodes) == 2))
                    BoolRobin_2 = ismember(elem(e).cn,robin_2(i).nodes);
                end
            end
            
            if (sum(BoolRobin_1) == 2)
                alpha = robin_1.alpha;
                T_inf = robin_1.Tinf;
            end
            % && ~isequal(elem(e).cn,[74,26,83]) && ~isequal(elem(e).cn,[53,11,55])
            if (sum(BoolRobin_2) == 2)
                alpha = robin_2(1).alpha;
                T_inf = robin_2(1).Tinf;
            end
            if (sum(BoolRobin_1) == 2 && sum(BoolRobin_2) == 2)
                alpha = (robin_1.alpha + robin_2(1).alpha)/2;
                T_inf = (robin_1.Tinf + robin_2(1).Tinf)/2;
            end

            disp('Temp at odes 79: ')
            a(79)
            Kce = zeros(3,3);
            if (alpha ~= 0)
                Kce = plantml(ex', ey', t*alpha)*matparam(3);
                fsur(edof) = fsur(edof) - alpha*(ae-T_inf)*Ae;
            end
            
            % assemble finte into global internal force vector fint
            fint(edof) = fint(edof) + finte;
            
            % assemble fvole into global volume force vector fvol
            % fvol(gdof) = fvol(gdof) + fvole;
            
            % assemble Ke into the global stiffness matrix K
            K(edof,edof) = K(edof,edof) + Ke + Kce;
            
        end % element loop

        % calculate residuum and its norm
        rsd(freeDofs) = -(fint(freeDofs) - fvol(freeDofs) - fsur(freeDofs));
        rsn = norm(rsd) + norm(drltValueIncr);
        fprintf(1, ' %2d. residuum norm= %e\n', iter, rsn);
        
        if (rsn >= TOL && iter < 2)
            disp('     solve and update');
            % modified right hand side
            rsd(freeDofs) = rsd(freeDofs) - K(freeDofs,drltDofs) * drltValueIncr;
            % calculate increment of node displacements
            da = zeros(size(a));
            da(freeDofs) = K(freeDofs,freeDofs)\rsd(freeDofs);
            %da(freeDofs) = rsd(freeDofs)\K(freeDofs,freeDofs);
            a(freeDofs) = a(freeDofs) + da(freeDofs);
            a(drltDofs) = a(drltDofs) + drltValueIncr;
            drltValueIncr = zeros(size(drltValueIncr));
        end
         
        iter = iter + 1;     
    end
    
    % save state variables for next time/load step
    for e = 1:nel
        elem(e).stateVar0 = elem(e).stateVar;
    end;
    
    postprocessing(step, time, ndm, elem, x, drltValue, drltDofs, fsur, RobinDofs, a);
    
  
                      
end;