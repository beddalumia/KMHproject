clear,clc

do_post = false;

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/DMFT/AFMx/';
cd(DATA)

SOIvals = [(0:9)/30,0.4,0.5,0.6];
Uvals = 0:0.1:10;

plane = zeros(length(Uvals),length(SOIvals));
for i=1:13
   soi = SOIvals(i);
   swd = sprintf('SOI=%f',soi);
   cd(swd); disp(swd)
   U = QcmP.post.get_list('U');
   if do_post
      QcmP.post.observables_line;
   end
   Mx = load('magX_1.txt');
   if do_post
      Z2 = QcmP.post.custom_line('Z2_inv.txt');
   else
      Z2 = load('Z2_inv.txt');
   end
   if(length(U)~=length(Uvals))
      Usubset = Uvals;
      if do_post
         QcmP.post.observables_line('U',[],Usubset);
      end
      Mxsubset = load('magX_1.txt');
      if do_post
         Z2subset = QcmP.post.custom_line('Z2_inv.txt','U',Usubset);
      else
         Z2 = load('Z2_inv.txt');
      end
   else
      Usubset = U;
      Mxsubset = Mx;
      Z2subset = Z2;
   end
   if soi==0
      Z2 = zeros(size(Z2));
      Z2subset=Z2;
   end
   cd('..')
   %threshold = 0.002*(length(SOIvals)-i/2)
   threshold = 0.02;
   if SOIvals(i)==0
      threshold = 0.05;
   end
   % if SOIvals(i)==4/30
   %    disp here
   %    threshold = 0.005;
   % end
   for j = 1:length(Usubset)
      if Z2subset(j) == 0 && abs(Mxsubset(j))>threshold
      plane(j,i) = 0.15; % Trivial AFMx solution
      elseif Z2subset(j) == 1 && abs(Mxsubset(j))>threshold
         plane(j,i) = 0.3-abs(Mxsubset(j)); % Topological AFMx solution
      else
         plane(j,i) = +0.36; % Topological nonmagnetized
      end
   end
end
%[X,Y] = meshgrid(SOIvals,Uvals);
[X,Y] = meshgrid(Uvals,SOIvals);
%imagesc(Uvals,SOIvals,plane')
surf(X,Y,plane'); view(2); hold on
shading flat
set(gca,'YDir','normal')
xlabel('$U/t$','Interpreter','latex')
ylabel('$\lambda_\mathrm{so}/t$','Interpreter','latex')
QcmP.plot.import_colorlab
cmap = diverging_cmap('ocean blue','peach',101);
cmap = [cmap(1:40,:);cmap(101,:)];
colormap(cmap)
ylim([0,0.6])
xlim([0,10]) 
grid off
box on
% for i = 1:length(SOIvals)
%    yline(SOIvals(i));
% end

cd(CODE)

matlab2tikz('filename','phase_diagram.tex','width','7cm','heigth','4cm');

%% Clean path
rmpath ../lib/m2tex/src