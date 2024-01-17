clear,clc

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/m2tex/src

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/RDMFT/for_paper';
cd(DATA)

%% Main body

figure('Name','Nanocorrelation')

t = tiledlayout(4, 2);

Nflake = [24,54,96,150];

for i = 2:5

   Ix = load(sprintf("corr_x_R%d.txt",i-1));
   Nx = size(Ix,2);
   Ux = 0:0.1:(Nx-1)/10;
   Iz = load(sprintf("corr_z_R%d.txt",i-1));
   Nz = size(Iz,2);
   Uz = 0:0.1:(Nz-1)/10;
   
   ruler = Ix(:,45);
   [~,k] = sort(ruler,'descend');
   Ix = Ix(k,:);
   Iz = Iz(k,:);
      
   Nedge = Nflake(i-1)/i;
   
   ax(i-1,1) = nexttile;
   for j = 1:Nedge
      patchline(Ux,Ix(j,:),'LineWidth',2,...
         'EdgeColor',str2rgb('lemon green'),'EdgeAlpha',5/Nedge);
      hold on
   end
   for j = Nedge+1:Nflake(i-1)
      patchline(Ux,Ix(j,:),'LineWidth',2,...
         'EdgeColor',str2rgb('jade'),'EdgeAlpha',3/(Nflake(i-1)-Nedge));
      hold on
   end; box on
   if i==2
      title("$I(\rightarrow:\leftarrow)$",'Interpreter','latex');
   elseif i==5
      xlabel("$U/t$",'Interpreter','latex')
   end
   xlim([0,round(max(Ux))])
   ylabel(sprintf('%dN-flake',i),'Interpreter','latex');

   ax(i-1,2) = nexttile;
   for j = 1:Nedge
      patchline(Uz,Iz(j,:),'LineWidth',2,...
         'EdgeColor',str2rgb('lemon green'),'EdgeAlpha',5/Nedge);
      hold on
   end
   for j = Nedge+1:Nflake(i-1)
      patchline(Uz,Iz(j,:),'LineWidth',2,...
         'EdgeColor',str2rgb('jade'),'EdgeAlpha',3/(Nflake(i-1)-Nedge));
      hold on
   end; box on
   if i==2
      title("$I(\,\uparrow\,:\,\downarrow\,)$",'Interpreter','latex');
   elseif i==5
      xlabel("$U/t$",'Interpreter','latex')
   end
   xlim([0,round(max(Uz))])

end

linkaxes(ax,'y'); ylim([0,0.35])

t.TileSpacing = 'compact';
t.Padding = 'compact';

cd(CODE)

%% Export to TikZ
matlab2tikz('filename','nanocorrelation.tex','width','6cm','height','15cm');

%% Reset path
rmpath ../lib/m2tex/src




%% contains

function p = patchline(xs,ys,varargin)
   % Plot lines as patches (efficiently)
   %
   % SYNTAX:
   %     patchline(xs,ys)
   %     patchline(xs,ys,zs,...)
   %     patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
   %     p = patchline(...)
   %
   % PROPERTIES: 
   %     Accepts all parameter-values accepted by PATCH.
   % 
   % DESCRIPTION:
   %     p = patchline(xs,ys,zs,'PropertyName',propertyvalue,...)
   %         Takes a vector of x-values (xs) and a same-sized
   %         vector of y-values (ys). z-values (zs) are
   %         supported, but optional; if specified, zs must
   %         occupy the third input position. Takes all P-V
   %         pairs supported by PATCH. Returns in p the handle
   %         to the resulting patch object.
   %         
   % NOTES:
   %     Note that we are drawing 0-thickness patches here,
   %     represented only by their edges. FACE PROPERTIES WILL
   %     NOT NOTICEABLY AFFECT THESE OBJECTS! (Modify the
   %     properties of the edges instead.)
   %
   %     LINUX (UNIX) USERS: One test-user found that this code
   %     worked well on his Windows machine, but crashed his
   %     Linux box. We traced the problem to an openGL issue;
   %     the problem can be fixed by calling 'opengl software'
   %     in your <http://www.mathworks.com/help/techdoc/ref/startup.html startup.m>.
   %     (That command is valid at startup, but not at runtime,
   %     on a unix machine.)
   %
   % EXAMPLES:
   %%% Example 1:
   %
   % n = 10;
   % xs = rand(n,1);
   % ys = rand(n,1);
   % zs = rand(n,1)*3;
   % plot3(xs,ys,zs,'r.')
   % xlabel('x');ylabel('y');zlabel('z');
   % p  = patchline(xs,ys,zs,'linestyle','--','edgecolor','g',...
   %     'linewidth',3,'edgealpha',0.2);
   %
   %%% Example 2: (Note "hold on" not necessary here!)
   %
   % t = 0:pi/64:4*pi;
   % p(1) = patchline(t,sin(t),'edgecolor','b','linewidth',2,'edgealpha',0.5);
   % p(2) = patchline(t,cos(t),'edgecolor','r','linewidth',2,'edgealpha',0.5);
   % l = legend('sine(t)','cosine(t)');
   % tmp = sort(findobj(l,'type','patch'));
   % for ii = 1:numel(tmp)
   %     set(tmp(ii),'facecolor',get(p(ii),'edgecolor'),'facealpha',get(p(ii),'edgealpha'),'edgecolor','none')
   % end
   %
   %%% Example 3 (requires Image Processing Toolbox):
   %%%   (NOTE that this is NOT the same as showing a transparent image on 
   %%%         of the existing image. (That functionality is
   %%%         available using showMaskAsOverlay or imoverlay).
   %%%         Instead, patchline plots transparent lines over
   %%%         the image.)
   %
   % img = imread('rice.png');
   % imshow(img)
   % img = imtophat(img,strel('disk',15));
   % grains = im2bw(img,graythresh(img));
   % grains = bwareaopen(grains,10);
   % edges = edge(grains,'canny');
   % boundaries = bwboundaries(edges,'noholes');
   % cmap = jet(numel(boundaries));
   % ind = randperm(numel(boundaries));
   % for ii = 1:numel(boundaries)
   % patchline(boundaries{ii}(:,2),boundaries{ii}(:,1),...
   %     'edgealpha',0.2,'edgecolor',cmap(ind(ii),:),'linewidth',3);
   % end
   %
   % Written by Brett Shoelson, PhD
   % brett.shoelson@mathworks.com
   % 5/31/2012
   % 
   % Revisions:
   % 6/26 Improved rice.png example, modified FEX image.
   % 11/30/30/2020 Fixed an issue on line 113 if varargin is empty.
   %
   % Copyright 2012 MathWorks, Inc.
   %
   % See also: patch, line, plot
   [zs,PVs] = parseInputs(varargin{:});
   if rem(numel(PVs),2) ~= 0
       % Odd number of inputs!
       error('patchline: Parameter-Values must be entered in valid pairs')
   end
   % Facecolor = 'k' is (essentially) ignored here, but syntactically necessary
   if isempty(zs)
       p = patch([xs(:);NaN],[ys(:);NaN],'k');
   else
       p = patch([xs(:);NaN],[ys(:);NaN],[zs(:);NaN],'k');
   end
   % Apply PV pairs
   for ii = 1:2:numel(PVs)
       set(p,PVs{ii},PVs{ii+1})
   end
   if nargout == 0
       clear p
   end
   function [zs,PVs] = parseInputs(varargin)
   if ~isempty(varargin) && isnumeric(varargin{1})
       zs = varargin{1};
       PVs = varargin(2:end);
   else
       PVs = varargin;
       zs = [];
   end

end
end