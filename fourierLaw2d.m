function  [q, D, stateVar] = fourierLaw2d(gradT, matparam, stateVar)
%
% linear elasticity
%
    %----- material parameter
    k = matparam(1);        % heat condutivity
    
    % constitutive matrix for isotropic material
    D = k*eye(2);
 
    % Fourier Law
    q = -D*gradT;
    
   % save temperature gradient gradT for postprocessing
    stateVar(1:2) = gradT;
   % save heat flux for postprocessing
    stateVar(3:4) = q; 
    
return;
    
