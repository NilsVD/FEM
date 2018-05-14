function  plotTriaMesh(elem, X)

% nel number of elements in the mesh
nel = length(elem);
% nnp number of node points
% ndm number of spatial dimensions
[nnp, ndm] = size(X);

%=========================================================================
% plot triangular mesh
%=========================================================================

figure(1);
clf;
%draw the elements with element numbers
for e=1:nel
    nen = length(elem(e).cn);
    Xe = zeros(ndm, nen);
    for idm = 1:ndm
        Xe(idm,:) = X(elem(e).cn,idm);
    end
    
    % plot undeformed configuration  
    patch(Xe(1,:), Xe(2,:), e);
    hold on
    Xcenter = 1/3 * sum(Xe(1,:));
    Ycenter = 1/3 * sum(Xe(2,:));
    elemNum = num2str(e, '%d');
%     h = text(Xcenter, Ycenter, elemNum);
%     set(h, 'Color', 'r'); 
    hold on;
end

for in=1:nnp
    nodeNum = num2str(in, '%d');
    h = text(X(in,1), X(in,2), nodeNum);
    set(h, 'Color', 'b');
end;



Xmax = max(X(:,1)); 
Xmin = min(X(:,1)) - 0.5;
Ymax = max(X(:,2)); 
Ymin = min(X(:,2));

title('mesh');
caxis('auto');
h = colorbar;

axis equal;
hold off
drawnow;

%fname = sprintf('pictures/1tri3_%03d.png', step);
%print(gcf, '-dpng', fname); 
