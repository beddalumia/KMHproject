clear,clc

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/RDMFT/for_paper';
cd(DATA)

%% Main body

figure('Name','AFMx-AFMz [flakes]')

for i = 2:5

   dE = load(sprintf("%dN_ediff.txt",i));
   Ne = length(dE); U = (0:0.1:0.1*(Ne-1))';

   s = scatter(U,dE,'d','LineWidth',1); hold on

   c = get_palette('matter',16);

   switch i 
      case(2)
         s.Marker = '^'; 
         s.MarkerEdgeColor = str2rgb('matlab3');%c(3,:);%
         s.MarkerFaceColor = [1.00000,1.00000,0.06667];
      case(3)
         s.Marker = 's'; 
         s.MarkerEdgeColor = str2rgb('matlab2');%c(5,:);%
         s.MarkerFaceColor = [1.00000,0.41176,0.16078];
      case(4)
         s.Marker = 'd'; 
         s.MarkerEdgeColor = str2rgb('matlab4');%c(10,:);%
         s.MarkerEdgeColor = [0.70196,0.56078,0.72941];
      case(5)
         s.Marker = 'v'; 
         s.MarkerEdgeColor = str2rgb('matlab1');%c(12,:);%
   end

   labels(i-1) = string(i)+"N-flake";

end

legend(labels,'Location','southwest','Interpreter','latex','box','off');
box on;
xlim([0,7.5]); xlabel('$U/t$','Interpreter','latex');
ylim([-1.25,0.75]); ylabel('$\displaystyle\sum_i \Bigl( E_{i,\parallel} - E_{i,\perp} \Bigr)$',...
      'Interpreter','latex');

cd(CODE)

%% Export to TikZ
matlab2tikz('filename','nanoenergy.tex','width','4cm','height','5cm');

%% Reset path
rmpath ../lib/m2tex/src





