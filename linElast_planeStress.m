function  [es, D, stateVar] = linElast_planeStress(et, ep, matparam, stateVar)
%
% linear elasticity
%
    %----- material parameter
    E =    matparam(1);        % Young's modulus
    nu =   matparam(2);        % Poisson ratio

    %----- analysis type
    % ptype == 1 plane stress
    % ptype == 2 plane strain
    ptype = ep(1);
    
    if (ptype == 1)
        % plane stress
        % constitutive matrix
        kappa = E/(1-nu*nu); 
        D =  kappa * [ 1     nu           0; ...
                      nu      1           0; ...
                       0      0    (1-nu)/2 ];
    
        % stress [sig_xx sig_yy tau_xy] 
        es = D * et;
        % strain and stress in z-direction
        eps_zz = -nu/E*(es(1) + es(2));
        sig_zz = 0;
    
    elseif (ptype == 2)
        %Plane strain
        et = [et(1);et(2);et(3)];
        
        %constitutive matrix
        kappa = E/(1+nu)/(1-2*nu); 
        D =  kappa * [ 1-nu    nu        0; ...
                       nu      1-nu      0; ...
                       0       0         1-2*nu];
                   
        es = D*et;
        
        % strain and stress in z-direction
        eps_zz = 0;
        sig_zz = kappa*[nu,nu]*[et(1);et(2)];
    end;
    
    % save strain for postprocessing
    stateVar(1:4) = [et(1)  et(2)  eps_zz  et(3)];
    % save stress for postprocessing
    stateVar(5:8) = [es(1)  es(2)  sig_zz  es(3)]; 
    
return;
    
