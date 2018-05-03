function  postprocessingMovie(step, time, ndm, elem, X, drltValue, drltDofs, fsur, neumDofs, a)

% number of elements in the mesh
nel = length(elem);

%=========================================================================
% post-processing - visualising results
%=========================================================================

scale_factor = 0.005;

figure(2);
clf;
for e=1:nel
    nen = length(elem(e).cn);
    Xe = zeros(ndm, nen);
    for idm = 1:ndm
        Xe(idm,:) = X(elem(e).cn,idm);
    end
    
    % dofs belonging to element e
    edof = elem(e).edof;
    Te = a(edof);
    
    % extract stress and strain tensor
    gradT = elem(e).stateVar(1:2);
    q = elem(e).stateVar(3:4);
    
    % plot deformed configuration  
    patch(Xe(1,:), Xe(2,:), Te); hold on;
    Xecenter = 1/3 * sum(Xe(1,:));
    Yecenter = 1/3 * sum(Xe(2,:));
    %quiver(Xecenter, Yecenter, q(1), q(2), scale_factor);
    hold on;
end

title('temperature and heat flux vector');
caxis('auto');
h = colorbar;
xlabel(h, 'temperature T');

axis equal;
hold off
drawnow;

%fname = sprintf('pictures/1tri3_%03d.png', step);
%print(gcf, '-dpng', fname); 
