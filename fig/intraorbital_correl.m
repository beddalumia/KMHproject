clear,clc

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/';
cd(DATA)

%% Paramagnetic

PARA = 'DMFT/Para/';

cd(PARA)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   [S{i},Sup{i},Sdw{i}] = build_local_entropies();
   U{i} = QcmP.post.get_list('U'); 
   I{i} = intraorbital(S{i},Sup{i},Sdw{i});

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      [x,y] = extra_points(U{i},S{i},"grid");
      [t,r] = extra_points(U{i},S{i},"tail");
      plot([x;t],[y;r],':','Color',str2rgb('mauve'),'Linewidth',1.2);
      plot(U{i},S{i},'s','MarkerSize',5,'MarkerEdgeColor',str2rgb('mauve'),...
         'MarkerFaceColor',str2rgb('light rose'),'Linewidth',0.5);
      annotation('textbox',[0.18,0.88,0.1,0.05],'string','$s_i$',...
         'Interpreter','latex','EdgeColor','none')
      [x,y] = extra_points(U{i},Sup{i},"tail");
      plot([U{i};x],[Sup{i};y],'-','Color',str2rgb('watermelon'),'Linewidth',2);
      [x,y] = extra_points(U{i},Sdw{i},"tail");
      plot([U{i};x],[Sdw{i};y],'--','Color',str2rgb('azure'),'Linewidth',2);
      annotation('textbox',[0.18,0.57,0.1,0.05],'string','$s_{i,\uparrow}=s_{i,\downarrow}$',...
         'Interpreter','latex','EdgeColor','none')
      [x,y] = extra_points(U{i},I{i},"grid");
      [t,r] = extra_points(U{i},I{i},"tail");
      fill([x;t;flipud(t);flipud(x)],[y;r;0.*y;0.*r],str2rgb("light khaki"),...
         'EdgeColor','none')
      plot([x;t],[y;r],':','Color',str2rgb('grass'),'Linewidth',1.2);
      plot(U{i},I{i},'d','MarkerSize',3,'MarkerEdgeColor',str2rgb('grass'),...
         'MarkerFaceColor',str2rgb('yellowish green'),'Linewidth',0.5);
      annotation('textbox',[0.18,0.19,0.1,0.05],'string','$I(\,\uparrow\,:\,\downarrow\,)$',...
         'Interpreter','latex','EdgeColor','none')
      xlim([0,15]); ylim([0,2]); box on;
      xlabel('$U/t$','Interpreter','latex');
      ylabel('[bit]','Interpreter','latex');

      % Export to TikZ
      matlab2tikz('filename',[CODE,'/para_correlation.tex'],'width','6cm','height','8cm');

   end

   cd('..')

end



%figure('Name','Paramagnetic intraorbital correlations')


close all
cd(CODE);
cd(DATA);

%% OUT OF PLANE AFM

AFMZ = 'DMFT/AFMz_normalSOI/';

cd(AFMZ)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   [S{i},Sup{i},Sdw{i}] = build_local_entropies();
   U{i} = QcmP.post.get_list('U'); 
   I{i} = intraorbital(S{i},Sup{i},Sdw{i});

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      plot(U{i},S{i},':','Color',str2rgb('mauve'),'Linewidth',1.2);
      plot(U{i},S{i},'s:','MarkerSize',5,'MarkerEdgeColor',str2rgb('mauve'),...
         'MarkerFaceColor',str2rgb('light rose'),'Linewidth',0.5);
      annotation('textbox',[0.18,0.88,0.1,0.05],'string','$s_i$',...
         'Interpreter','latex','EdgeColor','none')
      plot(U{i},Sup{i},'-','Color',str2rgb('watermelon'),'Linewidth',2);
      plot(U{i},Sdw{i},'--','Color',str2rgb('azure'),'Linewidth',2);
      annotation('textbox',[0.18,0.57,0.1,0.05],'string','$s_{i,\uparrow}=s_{i,\downarrow}$',...
         'Interpreter','latex','EdgeColor','none')
      fill([U{i};flipud(U{i})],[I{i};0.*I{i}],str2rgb("light khaki"),...
         'EdgeColor','none')
      plot(U{i},I{i},':','Color',str2rgb('grass'),'Linewidth',1.2);
      plot(U{i},I{i},'d','MarkerSize',3,'MarkerEdgeColor',str2rgb('grass'),...
         'MarkerFaceColor',str2rgb('yellowish green'),'Linewidth',0.5);
      annotation('textbox',[0.18,0.19,0.1,0.05],'string','$I(\,\uparrow\,:\,\downarrow\,)$',...
         'Interpreter','latex','EdgeColor','none')
      xlim([0,10]); ylim([0,2]); box on;
      xlabel('$U/t$','Interpreter','latex');
      ylabel('[bit]','Interpreter','latex');
      % Export to TikZ
      matlab2tikz('filename',[CODE,'/afmz_correlation.tex'],'width','6cm','height','8cm');
   end

   cd('..')

end

close all
cd(CODE);
cd(DATA);

%% EASY PLANE AFM

AFMX = 'DMFT/AFMx/';

cd(AFMX)

[so_vals, so_dirs] = QcmP.post.get_list('SOI');

Nso = length(so_vals);

for i = 1:Nso

   cd(so_dirs(i));

   [S{i},Sup{i},Sdw{i}] = build_local_entropies();
   U{i} = QcmP.post.get_list('U'); 
   I{i} = intraorbital(S{i},Sup{i},Sdw{i});

   if so_vals(i)==0.3
      figure("Name",so_dirs(i))
      hold on
      plot(U{i},S{i},':','Color',str2rgb('mauve'),'Linewidth',1.2);
      plot(U{i},S{i},'s:','MarkerSize',5,'MarkerEdgeColor',str2rgb('mauve'),...
         'MarkerFaceColor',str2rgb('light rose'),'Linewidth',0.5);
      annotation('textbox',[0.18,0.88,0.1,0.05],'string','$s_i$',...
         'Interpreter','latex','EdgeColor','none')
      plot(U{i},Sup{i},'-','Color',str2rgb('watermelon'),'Linewidth',2);
      plot(U{i},Sdw{i},'--','Color',str2rgb('azure'),'Linewidth',2);
      annotation('textbox',[0.18,0.57,0.1,0.05],'string','$s_{i,\rightarrow}=s_{i,\leftarrow}$',...
         'Interpreter','latex','EdgeColor','none')
      fill([U{i};flipud(U{i})],[I{i};0.*I{i}],str2rgb("light khaki"),...
         'EdgeColor','none')
      plot(U{i},I{i},':','Color',str2rgb('grass'),'Linewidth',1.2);
      plot(U{i},I{i},'d','MarkerSize',3,'MarkerEdgeColor',str2rgb('grass'),...
         'MarkerFaceColor',str2rgb('yellowish green'),'Linewidth',0.5);
      annotation('textbox',[0.18,0.19,0.1,0.05],'string','$I(\rightarrow:\leftarrow)$',...
         'Interpreter','latex','EdgeColor','none')
      xlim([0,10]); ylim([0,2]); box on;
      xlabel('$U/t$','Interpreter','latex');
      ylabel('[bit]','Interpreter','latex');
      % Export to TikZ
      matlab2tikz('filename',[CODE,'/afmx_correlation.tex'],'width','6cm','height','8cm');
   end

   cd('..')

end

close all
cd(CODE);

%% Reset path
rmpath ../lib/m2tex/src

%% contains

function [Si,Sup,Sdw] = build_local_entropies()

   try   % if there's no docc file then
      Di = load('docc_1.txt');
   catch % we need to retrieve all obs.
      QcmP.post.observables_line('U');
      Di = load('docc_1.txt');
   end

   try   % NORMAL ED calculations
      Mi = load('mag_1.txt');
   catch % NONSU2 ED calculations
      Mx = load('magX_1.txt');
      Mz = load('magZ_1.txt');
      Mi = max(abs(Mx),abs(Mz));
   end

   Si = (Di + (Mi-1)/2).*log2((1-Mi)/2-Di) ...
      + (Di - (Mi+1)/2).*log2((1+Mi)/2-Di) ...
      - 2*Di.*log2(Di);
   
   assert(all(Si>0));

   Sup_plus_Sdw = (Mi-1).*log2((1-Mi)/2) ...
                - (Mi+1).*log2((1+Mi)/2);

   assert(all(Sup_plus_Sdw>0))

   Sup = Sup_plus_Sdw/2;
   Sdw = Sup_plus_Sdw/2;

end

function i = intraorbital(s,sup,sdw)
   i = sup + sdw - s; 
end

function [x_,y_] = extra_points(x,y,mode)
   switch mode
   case "tail"
      model = fit(x(end-6:end),y(end-6:end),'smoothingspline');
      x_ = x(end):0.01:(1.1*x(end)); 
   case "grid"
      model = fit(x,y,'smoothingspline');
      x_ = x(1):0.01:x(end); 
   end
   x_ = x_';
   y_ = model(x_);
end


