clear,clc

%% INPUT VARS
U = [4.5,6.0];
M = ["X","Z"];

%% LIBRARIES
QcmP.plot.import_colorlab
addpath ../lib/arrow
% addpath ../lib/m2tex/src
addpath ../lib/export_fig

%% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/RDMFT/for_paper';
cd(DATA)

%% Main body
for i = 5:-1:2
    
    flake = load(sprintf("%dN_flake.txt",i));

    x = flake(:,1);
    y = flake(:,2);

    for j = 1:2

        id = sprintf("%dN Flake @ U/t=%s",i,string(U(j)));
        filename = sprintf("%s/%dNflakeU%s",CODE,i,string(U(j)));
        frame = single_frame(id,x,y,i,U(j));

        % matlab2tikz('filename',[char(filename),'.tex'],'width','5cm','height','4cm');
        % >> DOES NOT WORK: 
        %     1. Rectangle does not accept a curvature in TikZ (or m2tex does not parse it...)
        %     2. Arrow colors are not reproduced well (weird stuff happening I guess)

        % exportgraphics(frame,filename+'.pdf','ContentType','vector')
        % >> HUGE FILE SIZE (also some weird artifacts on arrow heads)

        set(frame, 'Color', 'White')
        export_fig(filename+'.pdf','-pdf','-painters')
        % >> PERFECTION <3

    end

end

cd(CODE)

%% Reset path
rmpath ../lib/arrow
% rmpath ../lib/m2tex/src
rmpath ../lib/export_fig

%% Contains

function h = single_frame(id,x,y,n,u)

    % load magnetic data
    magX = load(sprintf("%dN_magX_U%s.txt",n,string(u)));
    magZ = load(sprintf("%dN_magX_U%s.txt",n,string(u)));

    % Represent as double bubblechart
    h = figure('Name',id); z = magX;
    draw_lattice(x(z>0),y(z>0),sqrt(z(z>0)/2),'flat'); hold on
    draw_lattice(x(z<0),y(z<0),sqrt(-z(z<0)/2),'flat'); hold on

    xlim([min(x)-0.5,max(x)+0.5])
    ylim([min(y)-0.5,max(y)+0.5])
    zlim([-1,1])
    axis image
    axis off
    
    % AFMz
    z = magZ;
    s = 1/max(x).*z(z>0); 
    c = palette.brewer(9,'YlOrRd');
    colormap(c(1:7,:))
    arrow3([x(z>0),y(z>0),0*z(z>0)],[x(z>0),y(z>0),0.5*z(z>0)],'|',3*s,3*s);
    s = 1/max(x).*z(z<0); 
    c = palette.brewer(9,'YlGnBu');
    colormap(c(1:7,:))
    arrow3([x(z<0),y(z<0),0*z(z<0)],[x(z<0),y(z<0),0.5*z(z<0)],'|',3*s,3*s);

    % AFMx
    z = magX;
    s = 1/max(x).*z(z>0);
    c = palette.brewer(11,'-PiYG');
    colormap(c(6:9,:))
    arrow3([x(z>0),y(z>0),0*z(z>0)],[x(z>0)+0.5*z(z>0),y(z>0),0*z(z>0)],'|',3*s,3*s);
    s = 1/max(x).*z(z<0); 
    c = palette.brewer(11,'PiYG');
    colormap(c(6:9,:))
    arrow3([x(z<0),y(z<0),0*z(z<0)],[x(z<0)+0.5*z(z<0),y(z<0),0*z(z<0)],'|',3*s,3*s);
    view(-20,45)

end

function sites = draw_lattice(x,y,z,mode,varargin)
    %% Draws a lattice, with site dimensions given in proper x,y units
    %  i.e. rescalign the figure would not change the aspect ratio, as
    %  the circles must maintain their dimensions in the rescaled x,y
    %  axes. This way we can ensure that, no matter what screensize or
    %  system specs (HPC, laptop, whatever) the sites do not overlap.
    %
    %  >> lattice = draw_site(x,y,r,varargin) 
    %
    %   x :: a 1d array of site coordinates (x axis)
    %   y :: a 1d array of site coordinates (y axis)
    %   r :: a 1d array of values for the site radii (hint: data ∝ area!)
    %
    %   lattice :: a cell array of handles to all the rectangles (to allow 
    %              OOP-style selection of properties, with granularity)
    %              > you can use the varargin to pass global settings, as
    %                keyargs or name-value pairs.
    %
    %  See also rectangle
    N = length(x);
    switch mode
    case("radius")
        r = z;
        a = repmat(1.0,N,1);
    case("color")
        r = repmat(0.5,N,1);
        a = z;
    case("flat")
        r = repmat(0.5,N,1);
        a = repmat(0.0,N,1);
    end
    sites = cell(N,1);
    for i = 1:N
        pos = [x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)];
        sites{i} = rectangle('Position',pos,'Curvature',[1 1],varargin{:});
        c = sites{i}.FaceColor; c = [c,a(i)]; sites{i}.FaceColor = c;
        hold on
    end
    hold off
end

