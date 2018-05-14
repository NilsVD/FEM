% ---- Mechanical problem ----

clear all
close all

[ndm, ndf, nnp, nel, nen, x, elem, matparam, drlt, neum, loadcurve, b] = input_Rocket_Engine_Mechanical();

plotTriaMesh(elem,x);

% separate dof types
allDofs = (1:1:nnp*ndf)';

% dofs with Dirichlet boundary conditions
numDrltDofs = ndf*length(drlt.nodes);
drltDofs = zeros(numDrltDofs,1);
for i=1:numDrltDofs/2
    node = drlt.nodes(i);
    drltDofs(2*i-1) = node*ndf - 1;
    drltDofs(2*i) = node*ndf;
end

% free dofs
freeDofs = setdiff(allDofs, drltDofs);

% dofs with Neumann boundary conditions
numNeumDofs = ndf*length(neum.nodes);
neumDofs = zeros(numNeumDofs,1);
for i=1:numNeumDofs/2
    node = neum.nodes(i);
    neumDofs(2*i-1) = node*ndf - 1;
    neumDofs(2*i) = node*ndf;
end

%-----calculate normal vectors along the boundary F-I for each element----- 
% This assumes that the deformations are small and therefor that the
% normalvectors along the boundary barely will change direction with the
% passage of time.

%Initialise matrix containing normal coordinates for each element along FI.
%[x-coord, y-coord]
normalVectors_FI = zeros(length(neum.elements), ndm+1);

for i=1:length(neum.elements)
    x1=0;
    x2=0;
    y1=0;
    y2=0;
    
    % coordinates of the element nodes
    ex = x(elem(neum.elements(i)).cn,1)*0.001;   % x-coord
    ey = x(elem(neum.elements(i)).cn,2)*0.001;   % y-coord
    
    %wnob = which nodes on the boundary
    wnob = ismember(elem(neum.elements(i)).cn, neum.nodes);
    
      if(wnob == [1,1,0])
         x1=ex(1);
         x2=ex(2);
         y1=ey(1);
         y2=ey(2);               
      elseif(wnob == [1,0,1])
         x1=ex(1);
         x2=ex(3);
         y1=ey(1);
         y2=ey(3);              
      elseif(wnob == [0,1,1])
         x1=ex(2);
         x2=ex(3);
         y1=ey(2);
         y2=ey(3);               
      end 
      
      %vector along the element boundary
      vec = [x2-x1; y2-y1];
      vec = 1/norm(vec) * vec;
      
      %Define matrix that transforms vector into its normal vector
      ntm = [0,-1;1,0];
      
      %Create normalvector to vector vec above.
      norm_vec = ntm*vec;
      
      %Fill boundaryVectors_FI with coordinates
      normalVectors_FI(i,1) = norm_vec(1);
      normalVectors_FI(i,2) = norm_vec(2);
      normalVectors_FI(i,3) = neum.elements(i);
end
%-------Done calculating normalvectors-------------------------------------

%=========================================================================
% fe-analysis
%=========================================================================

%code to check boundary and normalvectors along 
%ghgh = [boundaryVectors_FI, normalVectors_FI];

% unknown node displacements
a = zeros(ndf*nnp,1);

%Initialise global Von-Mises vector
Vm = zeros(nel,1);

% simulation time and time/load step
stopTime = 1.0;
dt = 0.1;
% numerical tolerance to detect equilibrium
TOL = 1e-10;

% loop over all time/load steps
time = 0;
step = 0;

figure
while (time < stopTime)
    time = time + dt;
    step = step + 1;
    
    %Initialise global fsur (traction forces)
    fsur = zeros(ndf*nnp,1);
    
    %Initialise global fload (load forces, t.ex gravity)
    fload = zeros(ndf*nnp,1);
    
    % initialise global tangent stiffness matrix
    K = zeros(ndf*nnp,ndf*nnp);
    
    % initialise global residuum
    rsd = zeros(ndf*nnp,1);
    
    fprintf(1,'step= %2d time= %8.4e dt= %8.4e\n', step, time, dt);
    % Neumann loads for this timestep
    %numNeumDofs = size(neum,1);
    neumValue = zeros(numNeumDofs,1);
    for i=1:numNeumDofs
        lcID = neum.loadcurveID;
        scale = neum.pmax * loadcurve(lcID).scalefactor;
        neumValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    %fsur(neumDofs) = neumValue;
    
    % Dirichlet boundary conditions
    drltValue = zeros(numDrltDofs,1);
    for i=1:numDrltDofs
        lcID = drlt.loadcurveID;
        scale = drlt.scale * loadcurve(lcID).scalefactor;
        drltValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    drltValueIncr = drltValue - a(drltDofs);
 
    % element loop
    for e=1:nel
        
        % dofs belonging to element e
        edof = elem(e).edof;
        
        % coordinates of the element nodes
        ex = x((elem(e).cn),1)*0.001;   % x-coord
        ey = x((elem(e).cn),2)*0.001;   % y-coord
        
        % displacement values at the element nodes
        ae = a(edof);
        
        % area of the element
        Ae=1/2*det([ones(3,1) ex ey]);
              
        % B-matrix for a triangular element with 3 nodes
        B = [ ey(2)-ey(3)  0  ey(3)-ey(1)  0  ey(1)-ey(2)  0;
              0  ex(3)-ex(2)  0  ex(1)-ex(3)  0  ex(2)-ex(1);
              ex(3)-ex(2) ey(2)-ey(3) ex(1)-ex(3) ey(3)-ey(1) ex(2)-ex(1) ey(1)-ey(2) ]* 1/(2*Ae);
            
        % strain
        et = B*ae;
        
        % linear elasticity
        % ep(1) = 1;    % plane stress
         ep(1) = 2;    % plane strain
        [es, D, stateVar] = linElast_planeStress(et,ep,matparam, elem(e).stateVar0(1,:));
        sig_zz = stateVar(7);
        
        % thickness of the plate
        t = matparam(4);
        
        %Calculate Von_mises fro element e
        Vm(e) = sqrt( es(1)^2 + es(2)^2 + sig_zz^2 - es(1)*es(2) - es(1)*sig_zz - es(2)*sig_zz + 3*es(3)^2);
        
        %------Skapa och fyll fsure beroende på neum villkoren--------------
        
        %Initialise local fsur
        fsure = zeros(ndf*nen,1);
        
        %check if element is on boundary F-I
        Boolneum = ismember(elem(e).cn, neum.nodes);
        if sum(Boolneum) == 2
            
            normalvector_e = [0,0];
            k=1;
            while k<=length(neum.elements)
                if e == neum.elements(k)
                    normalvector_e = normalVectors_FI(k,1:2);
                end
                k = k+1;
            end
            
            %Instansiate coordinates to create length to integrate over
            x1=0;
            x2=0;
            y1=0;
            y2=0;
            
            if Boolneum(1)==0
                        x1 = ex(2);
                        x2 = ex(3);
                        y1 = ey(2);
                        y2 = ey(3);
            end
            if Boolneum(2)==0
                        x1 = ex(1);
                        x2 = ex(3);
                        y1 = ey(1);
                        y2 = ey(3);
            end
            if Boolneum(3)==0
                        x1 = ex(1);
                        x2 = ex(2);
                        y1 = ey(1);
                        y2 = ey(2);
            end
            
            line([x1 x2],[y1 y2])
            hold on
            %boundary length to be integrated over
            L = sqrt((x2-x1)^2 + (y2-y1)^2);
            
               fsureVec = zeros(2,6);
                    if(Boolneum == [1,1,0])
                        fsureVec = [1,0,1,0,0,0;
                                    0,1,0,1,0,0];
                    elseif(Boolneum == [1,0,1])
                        fsureVec = [1,0,0,0,1,0;
                                    0,1,0,0,0,1];
                    elseif(Boolneum == [0,1,1])
                        fsureVec = [0,0,1,0,1,0;
                                    0,0,0,1,0,1];
                    end
                    
                    fsure = t*L/2*fsureVec'*neum.pmax*normalvector_e';
        end
        
        
        % element stiffness matrix
        Ke = B'*D*B*Ae*t;
        
        % element load vector
        floade = b(2)*t*Ae/3*[0;1;0;1;0;1];
        
%----Assemble vectors och matrices into global variants--------------------
        
        % assemble fsure into global surface force vector fsur
        fsur(edof) = fsur(edof) + fsure;
        
        % assemble floade into global load vector fload
        fload(edof) = fload(edof) + floade;
        
        % assemble Ke into the global stiffness matrix K
        K(edof,edof) = K(edof,edof) + Ke;
        
    end % element loop
%--------------------------------------------------------------------------
if step==1
    fsur(neumDofs)
end
    
    % calculate residuum and its norm
    rsd(freeDofs) =  fsur(freeDofs) + fload(freeDofs) - K(freeDofs,freeDofs)*a(freeDofs); 
    
    % modified right hand side
    rsd(freeDofs) = rsd(freeDofs) - K(freeDofs,drltDofs) * drltValueIncr;
    rsn = norm(rsd(freeDofs));
    fprintf(1,'residuum norm %e \n',rsn)
    if rsn > 1e-8	
    % calculate increment of node displacements
       da = zeros(size(a)); 
       da(freeDofs) = K(freeDofs,freeDofs) \ rsd(freeDofs);
       a(freeDofs) = a(freeDofs) + da(freeDofs);
       a(drltDofs) = a(drltDofs) + drltValueIncr;
       drltValueIncr = zeros(size(drltValueIncr));
    end
    % save state variables for next time/load step
    for e = 1:nel
        elem(e).stateVar0 = elem(e).stateVar;
    end;
    VmNodes = zeros(nnp,1);
    Conn = zeros(nel,3);
    for i=1:nel
        Conn(i, 1:3) = elem(i).cn;
    end
    for i=1:size(x,1)
        [c0, c1] = find(Conn == i);
        VmNodes(i) = sum(Vm(c0))/size(c0,1);
    end
    
    reaction_forces = K(drltDofs,freeDofs)*a(freeDofs);
    reaction_forces_x = 0;
    reaction_forces_y = 0;
    for i=1:length(drltDofs)
        if mod(i,2) == 1
            reaction_forces_x = reaction_forces_x + reaction_forces(i);
        else
            reaction_forces_y = reaction_forces_y + reaction_forces(i);
        end
    end
        
    
    postprocessingMovie(step, time, ndm, elem, x, drltValue, drltDofs, fsur, neumDofs, a, normalVectors_FI, neum, VmNodes);
                      
end;