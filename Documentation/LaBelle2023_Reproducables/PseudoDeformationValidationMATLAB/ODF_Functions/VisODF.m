function [ out ] = VisODF( ODF, Ncart, ax_tit)

% in order to better visualize the ODF, we will display the ODF such that
% the sum of the point-wise probability values times their representative
% area adds up to 1. This is in contrast to their representation after
% min-max normalization, where the points represent the probability
% averaged over their areas. That is convenienvt for measuring FR distance.

% This way of visualization makes the resultant ODFs indpendent of the
% number of points used to sample the sphere.

% Assume the unit sphere is sampled such that every point has equal weighting
% (almost true with the dodecahedron sampling we use)
PtArea = 4*pi/size(Ncart, 1); % area on the unit sphere each point represents

VUC = sum(PtArea*ODF); % Volume under the curve of the probability function
ODF = ODF./VUC; % normalized probability function (pdf)

% check 
% if abs(1.0 - sum(ODF*PtArea)) < 0.001
%     disp('ODF point probabilities times surface area sums to 1');
% end

g = gca; hold on; g.LineWidth = 2;
colormap('jet');
scatter3(Ncart(:,1), Ncart(:,2), Ncart(:,3), 110, ODF, '.');
view([-60, 30]); hold off
xlim([-1,1]); xticks([]); 
ylim([-1,1]); yticks([]);
zlim([-1,1]); zticks([]);
%X and Y flipped because matlab messes with the eigenvector order...
xlabel('Y','FontSize', 24,'FontWeight','bold','Position',[0 -1 -1]);
ylabel('X','FontSize', 24,'FontWeight','bold','Position',[-1 0 -1]);
zlabel('Z','FontSize', 24,'FontWeight','bold','Position',[-1.5 1 0],'Rotation',0);
title(ax_tit,'FontSize',24,'FontWeight','bold');
axis image;  colorbar; 
c = colorbar;
c.LineWidth=1; caxis([0.02 0.15]); c.Ticks=[0.02 0.085 0.15];
%caxis([0.02 0.15]);
% try
%     caxis([min(ODF), max(ODF)]);
% catch
%     caxis([0,max(ODF)])
% end
grid off; set(gca, 'FontSize',20); box on;
out = false;

end
