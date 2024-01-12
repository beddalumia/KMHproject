clear,clc

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/RDMFT/for_paper';
cd(DATA)

%% Main body

figure('Name','Potential [flakes]')

Nflake = [24,54,96,150];

for i = 2:5

   dE = load(sprintf("%dN_ediff.txt",i)); Ne = length(dE);
   Kz = load(sprintf("%dN_kinZ.txt",i));  Nz = length(Kz);
   Kx = load(sprintf("%dN_kinX.txt",i));  Nx = length(Kx);
   Nk = min([Ne,Nz,Nx]); U = (0:0.1:0.1*(Nk-1))';
   dE = dE(1:Nk);  Kz = Kz(1:Nk);  Kx = Kx(1:Nk); 

   s = plot(U,dE+Kz*Nflake(i-1)-Kx*Nflake(i-1),'d:','LineWidth',1); hold on

   c = get_palette('matter',16);

   switch i 
      case(2)
         s.Color = str2rgb('matlab3');
         s.Marker = '^';
         s.MarkerEdgeColor = str2rgb('matlab3');%c(3,:);%
         s.MarkerFaceColor = [1.00000,1.00000,0.06667];
      case(3)
         s.Color = str2rgb('matlab2');
         s.Marker = 's'; 
         s.MarkerEdgeColor = str2rgb('matlab2');%c(5,:);%
         s.MarkerFaceColor = [1.00000,0.41176,0.16078];
      case(4)
         s.Color = str2rgb('matlab4');
         s.Marker = 'd'; 
         s.MarkerEdgeColor = str2rgb('matlab4');%c(10,:);%
         s.MarkerEdgeColor = [0.70196,0.56078,0.72941];
      case(5)
         s.Color = str2rgb('matlab1');
         s.Marker = 'v'; 
         s.MarkerEdgeColor = str2rgb('matlab1');%c(12,:);%
   end

   labels(i-1) = string(i)+"N-flake";

end

legend(labels,'Location','southwest','Interpreter','latex','box','off');
box on;
xlim([0,5]); 
xlabel('$U/t$','Interpreter','latex');
ylim([-0.75,0.75]);
ylabel('$\displaystyle\sum_i \Bigl( U_{i,\parallel} - U_{i,\perp} \Bigr)$','Interpreter','latex');

cd(CODE)

%% Export to TikZ
matlab2tikz('filename','nanopotential.tex','width','4cm','height','5cm');

%% Reset path
rmpath ../lib/m2tex/src





