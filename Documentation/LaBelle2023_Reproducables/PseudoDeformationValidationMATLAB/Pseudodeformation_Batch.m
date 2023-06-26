% Steven LaBelle
% MRL - angiogenesis project
%
% Code used to verify Pseudodeformation by using a MonteCarlo sampling approach. Set
% the SPD by changing the principal axes and rotation angles.

% Methods:
% 1) Set up the SPD by applying a rotation matrix.
% 2) Generate samples of position and reject those that are exterior to the
% SPD.

%% Init

clc; clear all; close all;
addpath('ODF_Functions','UniformSphereSampling');
warning('off');

% misc variables used. 
% TR_v: integer. Icosahedron mesh discretization increases with this.
% TR: Icosahedron mesh generated using the external function provided by Anton Semechko (https://www.mathworks.com/matlabcentral/fileexchange/37004-suite-of-functions-to-perform-uniform-sampling-of-a-sphere)
% rS: Random stream using mersenne twister.
% q: Angle in degrees used to search for local points on an icosahedron mesh.
TR_v = 6;
TR = IcosahedronMesh; TR = SubdivideSphericalMesh(TR,TR_v);
rS = RandStream('mt19937ar');
q = 60;
Icos = TR.Points;

% set variables to call the default matlab colors
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.9290, 0.6940, 0.1250];
purple = [0.4940, 0.1840, 0.5560];
green = [0.4660, 0.6740, 0.1880];
sky = [0.3010, 0.7450, 0.9330];
maroon = [0.6350, 0.0780, 0.1840];
black = [0.25, 0.25, 0.25];


% set the semiprinciapl axes of the SPD
a = [1 3 3]; b = [1 1 3]; c = [1 1 1];
lam = 0.5:0.1:1.5;
kap = 0:0.05:0.5;
ODF_True_uniaxial = zeros(length(Icos),3); 
ODF_Pseudo_uniaxial = zeros(length(Icos),3);
GFA_True_uniaxial=zeros(length(lam),3);
GFA_Pseudo_uniaxial=zeros(length(lam),3);
GFA_Diff_uniaxial=zeros(length(lam),3);
distDeg_uniaxial=zeros(length(lam),3);

ODF_True_biaxial = zeros(length(Icos),3); 
ODF_Pseudo_biaxial = zeros(length(Icos),3);
GFA_True_biaxial=zeros(length(lam),3);
GFA_Pseudo_biaxial=zeros(length(lam),3);
GFA_Diff_biaxial=zeros(length(lam),3);
distDeg_biaxial=zeros(length(lam),3);

ODF_True_ss = zeros(length(Icos),3); 
ODF_Pseudo_ss = zeros(length(Icos),3);
GFA_True_ss=zeros(length(kap),3);
GFA_Pseudo_ss=zeros(length(kap),3);
GFA_Diff_ss=zeros(length(kap),3);
distDeg_ss=zeros(length(kap),3);

ODF_True_ps = zeros(length(Icos),3); 
ODF_Pseudo_ps = zeros(length(Icos),3);
GFA_True_ps=zeros(length(kap),3);
GFA_Pseudo_ps=zeros(length(kap),3);
GFA_Diff_ps=zeros(length(kap),3);
distDeg_ps=zeros(length(kap),3);

for i=1:3
    v = eye(3); d = diag([a(i),b(i),c(i)]);
    P = v*d*transpose(v);

    % provide the rotation angles then rotate the SPD
    t1 = 45;
    t2 = 45;
    t3 = 45;
    Rotx = [1 0 0; 0 cosd(t1) -sind(t1); 0 sind(t1) cosd(t1)];
    Roty = [cosd(t2) 0 sind(t2); 0 1 0; -sind(t2) 0 cosd(t2)];
    Rotz = [cosd(t3) -sind(t3) 0; sind(t3) cosd(t3) 0; 0 0 1];
    Q = Rotz*Roty*Rotx;
    P = Q*P*transpose(Q);

    % determine the eigensolution and get the indices of descending order along the diagonal of d. 
    [v, d] = eig(P);
    [R, ~, ~] = poldecomp(v);

    %% Uniaxial Tension/Compression


    %Perform Uniaxial Tension/Compression
    for j=1:length(lam)
        F=[lam(j) 0                 0              ; 
           0      1/sqrt(lam(j))    0              ;
           0      0                 1/sqrt(lam(j))];
       % Get the rotational component of F
       [RF, ~, ~] = poldecomp(F);
       ODF_True = getODFKnown(Icos,P,F);
       ODF_Pseudo = getODFPseudo(Icos, P, F);
       GFA_True_uniaxial(j,i) = std(ODF_True)/rms(ODF_True);
       GFA_Pseudo_uniaxial(j,i) = std(ODF_Pseudo)/rms(ODF_Pseudo);
       GFA_Diff_uniaxial(j,i) = GFA_True_uniaxial(j,i)-GFA_Pseudo_uniaxial(j,i);
       distDeg_uniaxial(j,i) = computeFisherRao(ODF_True, ODF_Pseudo);   
       if j==length(lam)
           ODF_True_uniaxial(:,i) = ODF_True;
           ODF_Pseudo_uniaxial(:,i) = ODF_Pseudo;
       end
    end
    %% Strip Biaxial Tension/Compression

    %Perform Biaxial Tension/Compression
    for j=1:length(lam)
        F=[lam(j) 0                 0        ; 
           0      1                 0        ;
           0      0                 1/lam(j)];
       % Get the rotational component of F
       [RF, ~, ~] = poldecomp(F);
       ODF_True = getODFKnown(Icos,P,F); 
       ODF_Pseudo = getODFPseudo(Icos, P, F);
       GFA_True_biaxial(j,i) = std(ODF_True)/rms(ODF_True);
       GFA_Pseudo_biaxial(j,i) = std(ODF_Pseudo)/rms(ODF_Pseudo);
       GFA_Diff_biaxial(j,i) = GFA_True_biaxial(j,i)-GFA_Pseudo_biaxial(j,i);
       distDeg_biaxial(j,i) = computeFisherRao(ODF_True, ODF_Pseudo);   
       if j==length(lam)
           ODF_True_biaxial(:,i) = ODF_True;
           ODF_Pseudo_biaxial(:,i) = ODF_Pseudo;
       end
    end
    %% Simple Shear
    maxhit = false;
    %Perform Simple Shear
    for j=1:length(kap)
        F=[1      kap(j)            0 ; 
           0      1                 0 ;
           0      0                 1];
       % Get the rotational component of F
       [RF, ~, ~] = poldecomp(F);
       ODF_True = getODFKnown(Icos,P,F); 
       ODF_Pseudo = getODFPseudo(Icos, P, F);
       GFA_True_ss(j,i) = std(ODF_True)/rms(ODF_True);
       GFA_Pseudo_ss(j,i) = std(ODF_Pseudo)/rms(ODF_Pseudo);
       GFA_Diff_ss(j,i) = GFA_True_ss(j,i)-GFA_Pseudo_ss(j,i);
       distDeg_ss(j,i) = computeFisherRao(ODF_True, ODF_Pseudo);   
       if (distDeg_ss(j,i) > 6 && ~maxhit)
           ODF_True_ss(:,i) = ODF_True;
           ODF_Pseudo_ss(:,i) = ODF_Pseudo;
           maxhit = true;
       end
       if (j==length(lam) && ~maxhit)
           ODF_True_ss(:,i) = ODF_True;
           ODF_Pseudo_ss(:,i) = ODF_Pseudo;
       end       
    end

    %% Pure Shear
    maxhit = false;
    %Perform Pure Shear
    for j=1:length(kap)
        F=[1        kap(j)  0              ; 
           kap(j)   1       0              ;
           0        0       1/(1-kap(j)^2)];
       % Get the rotational component of F
       [RF, ~, ~] = poldecomp(F);
       ODF_True = getODFKnown(Icos,P,F); 
       ODF_Pseudo = getODFPseudo(Icos, P, F);
       GFA_True_ps(j,i) = std(ODF_True)/rms(ODF_True);
       GFA_Pseudo_ps(j,i) = std(ODF_Pseudo)/rms(ODF_Pseudo);
       GFA_Diff_ps(j,i) = GFA_True_ps(j,i)-GFA_Pseudo_ps(j,i);
       distDeg_ps(j,i) = computeFisherRao(ODF_True, ODF_Pseudo);   
       if (distDeg_ps(j,i) > 6 && ~maxhit)
           ODF_True_ps(:,i) = ODF_True;
           ODF_Pseudo_ps(:,i) = ODF_Pseudo;
           maxhit = true;
       end
       if (j==length(lam) && ~maxhit)
           ODF_True_ps(:,i) = ODF_True;
           ODF_Pseudo_ps(:,i) = ODF_Pseudo;
       end       
    end    
end

figure; 
for i=1:3 
    %% Plotting Uniaxial
    subplot(3,4,(i-1)*4+1);
    plot(lam,abs(GFA_Diff_uniaxial(:,i)),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('GFA Difference','FontSize', 24);
    end
    xlabel('\lambda','FontSize',24,'Position',[1,-0.0001,1]); ylabel('Difference (GFA)','FontSize', 24); xlim([0.5 1.5]); ylim([-1e-5 1e-3]);
    g = gca; g.LineWidth=3; g.XTick=[0.5 1.5]; g.YTick=[0 1e-3]; g.YTickLabel={["0"],["1e-3"]}; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';

    subplot(3,4,(i-1)*4+2);
    plot(lam,distDeg_uniaxial(:,i),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('Fisher-Rao Distance','FontSize', 24); 
    end
    xlabel('\lambda','FontSize', 24,'Position',[1,-0.0001,1]); ylabel('Distance (\circ)','FontSize', 24); xlim([0.5 1.5]); ylim([-0.1 6]);
    g = gca; g.LineWidth=3; g.XTick=[0.5 1.5]; g.YTick=[0 6]; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
    
    subplot(3,4,(i-1)*4+3);
    if(i==1)
        VisODF(ODF_True_uniaxial(:,i),Icos,'True Deformation');
    else
        VisODF(ODF_True_uniaxial(:,i),Icos,'');
    end
    
    subplot(3,4,(i-1)*4+4);
    if(i==1)
        VisODF(ODF_Pseudo_uniaxial(:,i),Icos,'Pseudo-Deformation');
    else
        VisODF(ODF_Pseudo_uniaxial(:,i),Icos,'');
    end
end

figure; 
for i=1:3 
    %% Plotting Biaxial
    subplot(3,4,(i-1)*4+1); 
    plot(lam,abs(GFA_Diff_biaxial(:,i)),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('GFA Difference','FontSize', 24);
    end
    xlabel('\lambda','FontSize', 24,'Position',[1,-0.0001,1]); ylabel('Difference (GFA)','FontSize', 24); xlim([0.5 1.5]); ylim([-1e-5 1e-3]);
    g = gca; g.LineWidth=3; g.XTick=[0.5 1.5]; g.YTick=[0 1e-3]; g.YTickLabel={["0"],["1e-3"]}; g.YAxis.Exponent = 0; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
    
    subplot(3,4,(i-1)*4+2);
    plot(lam,distDeg_biaxial(:,i),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('Fisher-Rao Distance','FontSize', 24);
    end
    xlabel('\lambda','FontSize', 24,'Position',[1,-0.0001,1]); ylabel('Distance (\circ)','FontSize', 24); xlim([0.5 1.5]); ylim([-0.1 6]);
    g = gca; g.LineWidth=3; g.XTick=[0.5 1.5]; g.YTick=[0 6]; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
    
    subplot(3,4,(i-1)*4+3);
    if(i==1)
        VisODF(ODF_True_biaxial(:,i),Icos,'True Deformation');
    else
        VisODF(ODF_True_biaxial(:,i),Icos,'');
    end
    
    subplot(3,4,(i-1)*4+4);
    if(i==1)
        VisODF(ODF_Pseudo_biaxial(:,i),Icos,'Pseudo-Deformation');
    else
        VisODF(ODF_Pseudo_biaxial(:,i),Icos,'');
    end
end


figure;
for i=1:3 
    %% Plotting Simple Shear
    subplot(3,4,(i-1)*4+1); 
    plot(kap,abs(GFA_Diff_ss(:,i)),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('GFA Difference:','FontSize', 24); 
    end
    xlabel('\kappa','FontSize', 24,'Position',[0.25,-0.00015,1]); ylabel('Difference (GFA)','FontSize', 24); xlim([0 0.5]); ylim([-1e-5 1e-3]);
    g = gca; g.LineWidth=3; g.YTick=[0 1e-3]; g.YTickLabel={["0"],["1e-3"]}; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
        
    subplot(3,4,(i-1)*4+2);
    plot(kap,distDeg_ss(:,i),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('Fisher-Rao Distance','FontSize', 24); 
    end
    xlabel('\kappa','FontSize', 24,'Position',[0.25,-0.00015,1]); ylabel('Distance (\circ)','FontSize', 24); xlim([0 0.5]); ylim([-0.1 6]);
    g = gca; g.LineWidth=3; g.YTick=[0 6]; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
    
    subplot(3,4,(i-1)*4+3);
    if(i==1)
        VisODF(ODF_True_ss(:,i),Icos,'True Deformation');
    else
        VisODF(ODF_True_ss(:,i),Icos,'');
    end
    
    subplot(3,4,(i-1)*4+4);
    if(i==1)
        VisODF(ODF_Pseudo_ss(:,i),Icos,'Pseudo-Deformation');
    else
        VisODF(ODF_Pseudo_ss(:,i),Icos,'');
    end
end

figure;
for i=1:3 
    %% Plotting Pure Shear
    subplot(3,4,(i-1)*4+1); 
    plot(kap,abs(GFA_Diff_ps(:,i)),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('GFA Comparison','FontSize', 24); 
    end
    xlabel('\kappa','FontSize', 24,'Position',[0.25,-0.00015,1]); ylabel('Difference (GFA)','FontSize', 24); xlim([0 0.5]); ylim([-1e-5 1e-3]);
    g = gca; g.LineWidth=3; g.YTick=[0 1e-3]; g.YTickLabel={["0"],["1e-3"]}; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
    
    subplot(3,4,(i-1)*4+2); 
    plot(kap,distDeg_ps(:,i),'-','Color',blue,'Linewidth',6);
    if(i==1)
        title('Fisher-Rao Distance','FontSize', 24); 
    end
    xlabel('\kappa','FontSize', 24,'Position',[0.25,-0.00015,1]); ylabel('Distance (\circ)','FontSize', 24); xlim([0 0.5]); ylim([-0.1 6]);
    g = gca; g.LineWidth=3; g.YTick=[0 6]; g.FontWeight='bold'; g.FontSize=20; g.Layer='top'; g.PlotBoxAspectRatio=[1,1,1]; g.Clipping='on';
    
    subplot(3,4,(i-1)*4+3); 
    if(i==1)
        VisODF(ODF_True_ps(:,i),Icos,'True Deformation');
    else
        VisODF(ODF_True_ps(:,i),Icos,'');
    end
    
    subplot(3,4,(i-1)*4+4); 
    if(i==1)
        VisODF(ODF_Pseudo_ps(:,i),Icos,'Pseudo-Deformation');
    else
        VisODF(ODF_Pseudo_ps(:,i),Icos,'');
    end
end