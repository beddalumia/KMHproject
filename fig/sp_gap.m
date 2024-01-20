clear,clc

global do_post
do_post = false;

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

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   G{i} = get_sp_fap(U{i});
   try
      Z{i} = QcmP.post.custom_line('Z2_inv.txt','U',U{i});
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      fill([U{i};flipud(U{i})],[Z{i}*10;0.*Z{i}],str2rgb("peach"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},G{i},'o:','MarkerSize',13,'MarkerEdgeColor',str2rgb('matlab4'),...
         'MarkerFaceColor',str2rgb('pale purple'),'Color',str2rgb('matlab4'),'Linewidth',2);
      xlim([0,10]); ylim([0,max(G{i})]);  box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end


cd(CODE);
cd(DATA);

%% OUT OF PLANE AFM [HF]

AFMZ = 'MF/AFMz/';

cd(AFMZ)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   G{i} = get_sp_fap(U{i});
   Ztmp = zeros(size(U{i}));
   try
      Zraw = QcmP.post.custom_line('Z2_inv.dat','U',U{i});
      Zraw(isnan(Zraw)) = 0;
      Z{i} = Zraw;
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   if so_vals(i)==0.3
      fill([U{i};flipud(U{i})],[Z{i}*10;0.*Z{i}],str2rgb("salmon"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},G{i},'d:','MarkerSize',7,'MarkerEdgeColor',str2rgb("matlab4"),...
         'MarkerFaceColor',str2rgb('pale lilac'),'Color',str2rgb("matlab4"),'Linewidth',2);
      xlim([0,10]); ylim([0,max(G{i})]);  box on;
      xlabel('$U/t$','Interpreter','latex');      
   end

   cd('..')

end

elements = get(gca,'Children'); uistack(elements(2),'down',1);

legend(["$\mathbf{Z}_2$ \enspace\,(DMFT)",
        "$\mathbf{Z}_2$ \enspace\,(HF)",
        "$\Delta/t$ (DMFT)",
        "$\Delta/t$ (HF)"],...
        'Location','Northwest','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/zgap.tex'],'width','6cm','height','6cm');

%close all
cd(CODE);
cd(DATA);

%% IN-PLANE AFM

AFMX = 'DMFT/AFMx/';

cd(AFMX)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   G{i} = get_sp_fap(U{i});
   try
      if do_post
         Z{i} = QcmP.post.custom_line('Z2_inv.txt')
      else
         Z{i} = load('Z2_inv.txt');
      end
   catch
      Z{i} = NaN*ones(size(U{i}));
   end


   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      fill([U{i};flipud(U{i})],[Z{i}*10;0.*Z{i}],str2rgb("peach"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},G{i},'o:','MarkerSize',13,'MarkerEdgeColor',str2rgb('matlab4'),...
         'MarkerFaceColor',str2rgb('pale purple'),'Color',str2rgb('matlab4'),'Linewidth',2);
      xlim([0,10]); ylim([0,max(G{i})]);  box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end

cd(CODE);
cd(DATA);

%% IN-PLANE AFM [HF]

AFMX = 'MF/AFMxy/';

cd(AFMX)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   G{i} = get_sp_fap(U{i});
   Ztmp = zeros(size(U{i}));
   try
      Zraw = load('Z2_inv.txt');
      Ztmp(1:length(Zraw)) = Zraw;
      Z{i} = Ztmp;
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   if so_vals(i)==0.3
      fill([U{i};flipud(U{i})],[Z{i}*10;0.*Z{i}],str2rgb("salmon"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},G{i},'d:','MarkerSize',7,'MarkerEdgeColor',str2rgb("matlab4"),...
         'MarkerFaceColor',str2rgb('pale lilac'),'Color',str2rgb("matlab4"),'Linewidth',2);
      xlim([0,10]); ylim([0,max(G{i})]);  box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end

elements = get(gca,'Children'); uistack(elements(2),'down',1);

legend(["$\mathbf{Z}_2$ \enspace\,(DMFT)",
        "$\mathbf{Z}_2$ \enspace\,(HF)",
        "$\Delta/t$ (DMFT)",
        "$\Delta/t$ (HF)"],...
        'Location','Northwest','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/xgap.tex'],'width','6cm','height','6cm');

%close all
cd(CODE);
cd(DATA);

%% PARAMAGNETIC DMFT

PARA = 'DMFT/Para/';

cd(PARA)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   U{i} = QcmP.post.get_list('U'); 
   G{i} = get_sp_fap(U{i});
   Ztmp = zeros(size(U{i}));
   try
      Zraw = load('Z2_inv.txt');
      Ztmp(1:length(Zraw)) = Zraw;
      Z{i} = Ztmp;
   catch
      Z{i} = NaN*ones(size(U{i}));
   end

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      fill([U{i};flipud(U{i})],[Z{i}*10;0.*Z{i}],str2rgb("salmon"),...
         'EdgeColor','none','FaceAlpha',0.7)
      plot(U{i},G{i},'o:','MarkerSize',13,'MarkerEdgeColor',str2rgb("matlab4"),...
         'MarkerFaceColor',str2rgb('pale purple'),'Color',str2rgb("matlab4"),'Linewidth',2);
      xlim([0,15]); box on;
      xlabel('$U/t$','Interpreter','latex');
   end

   cd('..')

end

legend(["$\Delta/t$ (DMFT)"],...
        'Location','Northwest','Interpreter','latex');

% Export to TikZ
matlab2tikz('filename',[CODE,'/paragap.tex'],'width','6cm','height','6cm');

%close all
cd(CODE);

%% Reset path
rmpath ../lib/m2tex/src

%% contains

function gap = get_sp_fap(Uvec)

   global do_post

   if do_post
      gap = zeros(length(Uvec),1);
      for i=1:length(Uvec)
         cd(sprintf("U=%f",Uvec(i)))
         try
            G = QcmP.plot.spectral_load('Greal_l11_s1_realw_indx000001.dat');
         catch
            G = QcmP.plot.spectral_load('Greal_l11_s1_realw__indx000001.dat');
         end
         w = G.zeta; R = -G.real; 
         wh = w(w>=0); Rh = abs(R(w>=0));
         wp = w(w<=0); Rp = abs(R(w<=0));
         Ro = 0;
         for j = length(wp):-1:1
            if Rp(j) <= Ro && Rp(j) > 0.05
               break
            else
               Ro = Rp(j);
            end
         end
         HOMO = wp(j);
         Ro = 0;
         for j = 1:+1:length(wh)
            if Rh(j) <= Ro && Rh(j) > 0.05
               break
            else
               Ro = Rh(j); 
            end
         end
         LUMO = wh(j);
         gap(i) = 2*max(-HOMO,LUMO);
         cd('..')
      end
      writematrix(gap);
   else
      gap = load('gap.txt');
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


