clear,clc

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/';
cd(DATA)

%% OUT OF PLANE AFM [DMFT]

AFMZ = 'DMFT/AFMz_rotateMag/';
AFMZ = 'DMFT/AFMz_normalSOI/';

cd(AFMZ)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);
Umat = repmat(0:0.1:15,Nso,1);
Mmat = zeros(size(Umat));

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   M{i} = get_order_parameter();
   try
      Z{i} = load('Z2_inv.txt');
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   [~,Mtemp] = extra_points(U{i},M{i},"grid");
   %Umat(i,:) = [Umat(i,:),Umat(i,end):0.1:15];
   Mmat(i,:) = [Mtemp',repmat(Mtemp(end),1,length(Umat(i,:))-length(Mtemp))];
   % > prepare smoothed data for later surface plot...

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      fill([U{i};flipud(U{i})],[Z{i};0.*Z{i}],str2rgb("peach"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},M{i},'o:','MarkerSize',13,'MarkerEdgeColor',str2rgb('ocean blue'),...
         'MarkerFaceColor',str2rgb('aqua blue'),'Color',str2rgb('ocean blue'),'Linewidth',2);
      xlim([0,10]); ylim([0,1]); box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end

% figure('Name','AFMz intraorbital correlations')
% [X,Y] = meshgrid(0:0.1:10,so_vals);
% surf(X,Y,Mmat,'EdgeColor','none','FaceColor','interp');
% caxis([0,1]); set_palette('viridis')
% zlim([0.,1]); view(-15,60); grid off 
% xlabel('$U/t$','Interpreter','latex');
% ylabel('$\lambda_\mathrm{so}/t$','Interpreter','latex');
% zlabel('$I(\,\uparrow\,:\,\downarrow\,)$ [bit]','Interpreter','latex');

% Export to TikZ
%matlab2tikz('filename',[CODE,'/zCorrSurf.tex'],'width','5cm','heigth','6cm');

%close all
cd(CODE);
cd(DATA);

%% OUT OF PLANE AFM [HF]

AFMZ = 'MF/AFMz/';

cd(AFMZ)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);
Umat = repmat(0:0.1:15,Nso,1);
Mmat = zeros(size(Umat));

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   M{i} = get_order_parameter();
   Ztmp = zeros(size(U{i}));
   try
      Zraw = load('Z2_inv.txt');
      Ztmp(1:length(Zraw)) = Zraw;
      Z{i} = Ztmp;
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   [~,Mtemp] = extra_points(U{i},M{i},"grid");
   %Umat(i,:) = [Umat(i,:),Umat(i,end):0.1:15];
   Mmat(i,:) = [Mtemp',repmat(Mtemp(end),1,length(Umat(i,:))-length(Mtemp))];
   % > prepare smoothed data for later surface plot...

   if so_vals(i)==0.3
      fill([U{i};flipud(U{i})],[Z{i};0.*Z{i}],str2rgb("salmon"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},M{i},'d:','MarkerSize',7,'MarkerEdgeColor',str2rgb("ocean blue"),...
         'MarkerFaceColor',str2rgb('light cyan'),'Color',str2rgb("ocean blue"),'Linewidth',2);
      xlim([0,10]); ylim([0,1]); box on;
      xlabel('$U/t$','Interpreter','latex');      
   end

   cd('..')

end

elements = get(gca,'Children'); uistack(elements(2),'down',1);

legend(["$\mathbf{Z}_2$ ~(DMFT)",
        "$\mathbf{Z}_2$ ~(HF)",
        "$\mathcal{M}_i$ (DMFT)",
        "$\mathcal{M}_i$ (HF)"],...
        'Location','Northwest','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/afmz_mag.tex'],'width','8cm','height','6cm');

%close all
cd(CODE);
cd(DATA);

%% IN-PLANE AFM

AFMX = 'DMFT/AFMx/';

cd(AFMX)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);
Umat = repmat(0:0.1:15,Nso,1);
Mmat = zeros(size(Umat));

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   M{i} = get_order_parameter();
   try
      Z{i} = load('Z2_inv.txt');
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   [~,Mtemp] = extra_points(U{i},M{i},"grid");
   %Umat(i,:) = [Umat(i,:),Umat(i,end):0.1:15];
   Mmat(i,:) = [Mtemp',repmat(Mtemp(end),1,length(Umat(i,:))-length(Mtemp))];
   % > prepare smoothed data for later surface plot...

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      fill([U{i};flipud(U{i})],[Z{i};0.*Z{i}],str2rgb("peach"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},M{i},'o:','MarkerSize',13,'MarkerEdgeColor',str2rgb('ocean blue'),...
         'MarkerFaceColor',str2rgb('aqua blue'),'Color',str2rgb('ocean blue'),'Linewidth',2);
      xlim([0,10]); ylim([0,1]); box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end

% figure('Name','AFMx intraorbital correlations')
% [X,Y] = meshgrid(0:0.1:10,so_vals);
% surf(X,Y,Mmat,'EdgeColor','none','FaceColor','interp');
% caxis([0,0.3]); set_palette('viridis')
% zlim([0.,0.3]); view(-15,60);  grid off
% xlabel('$U/t$','Interpreter','latex');
% ylabel('$\lambda_\mathrm{so}/t$','Interpreter','latex');
% zlabel('$I(\rightarrow:\leftarrow)$ [bit]','Interpreter','latex');

% % Export to TikZ
% matlab2tikz('filename',[CODE,'/xCorrSurf.tex'],'width','5cm','heigth','6cm');

%close all
cd(CODE);
cd(DATA);

%% IN-PLANE AFM [HF]

AFMX = 'MF/AFMxy/';

cd(AFMX)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);
Umat = repmat(0:0.1:15,Nso,1);
Mmat = zeros(size(Umat));

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   M{i} = get_order_parameter();
   Ztmp = zeros(size(U{i}));
   try
      Zraw = load('Z2_inv.txt');
      Ztmp(1:length(Zraw)) = Zraw;
      Z{i} = Ztmp;
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   [~,Mtemp] = extra_points(U{i},M{i},"grid");
   %Umat(i,:) = [Umat(i,:),Umat(i,end):0.1:15];
   Mmat(i,:) = [Mtemp',repmat(Mtemp(end),1,length(Umat(i,:))-length(Mtemp))];
   % > prepare smoothed data for later surface plot...

   if so_vals(i)==0.3
      fill([U{i};flipud(U{i})],[Z{i};0.*Z{i}],str2rgb("salmon"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},M{i},'d:','MarkerSize',7,'MarkerEdgeColor',str2rgb("ocean blue"),...
         'MarkerFaceColor',str2rgb('light cyan'),'Color',str2rgb("ocean blue"),'Linewidth',2);
      xlim([0,10]); ylim([0,1]); box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end

elements = get(gca,'Children'); uistack(elements(2),'down',1);

legend(["$\mathbf{Z}_2$ ~(DMFT)",
        "$\mathbf{Z}_2$ ~(HF)",
        "$\mathcal{M}_i$ (DMFT)",
        "$\mathcal{M}_i$ (HF)"],...
        'Location','Northwest','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/afmx_mag.tex'],'width','8cm','height','6cm');

%close all
cd(CODE);

%% Reset path
rmpath ../lib/m2tex/src

%% contains

function [Mi] = get_order_parameter()

   try % DMFT
      try   % NORMAL ED calculations
         Mi = load('mag_1.txt');
      catch % NONSU2 ED calculations
         Mx = load('magX_1.txt');
         Mz = load('magZ_1.txt');
         Mi = max(abs(Mx),abs(Mz));
      end
   catch
      % Hartree Fock
      try
         Mi = sqrt(load('Rx.txt').^2 + load('Ry.txt').^2) / 2;
      catch
         Mi = load('Rz.txt') / 2;
      end
   end

end


function [x_,y_] = extra_points(x,y,mode)
   switch mode
   case "tail"
      model = fit(x(end-6:end),y(end-6:end),'smoothingspline');
      x_ = x(end):0.1:(1.1*x(end)); 
   case "grid"
      model = fit(x,y,'smoothingspline');
      x_ = x(1):0.1:x(end); 
   end
   x_ = x_';
   y_ = model(x_);
end


