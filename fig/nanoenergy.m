clear,clc

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/RDMFT/for_paper';
cd(DATA)

%% Main body
for i = 2:5

   dE = load(sprintf("%dN_ediff.txt",i));
   Ne = length(dE); U = (0:0.1:0.1*(Ne-1))';

   s = scatter(U,dE,'d','LineWidth',2); hold on

   c = get_palette('matter',16);

   switch i 
      case(2)
         s.Marker = '^'; s.MarkerEdgeColor = c(3,:);%str2rgb('matlab3')%
      case(3)
         s.Marker = 's'; s.MarkerEdgeColor = c(5,:);%str2rgb('matlab2')%
      case(4)
         s.Marker = 'd'; s.MarkerEdgeColor = c(10,:);%str2rgb('matlab2')%
      case(5)
         s.Marker = 'v'; s.MarkerEdgeColor = c(12,:);%str2rgb('matlab1')%
   end

end

cd(CODE)

xlim([0,7.2])
ylim([-1.2,1])
