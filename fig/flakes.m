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
        [zframe, xframe ] = double_frame(id,x,y,i,U(j));

        % matlab2tikz('filename',[char(filename),'.tex'],'width','5cm','height','4cm');
        % >> DOES NOT WORK: 
        %     1. Rectangle does not accept a curvature in TikZ (or m2tex does not parse it...)
        %     2. Arrow colors are not reproduced well (weird stuff happening I guess)

        %exportgraphics(zframe,filename+'_z.pdf','ContentType','vector')
        %exportgraphics(xframe,filename+'_x.pdf','ContentType','vector')
        % >> HUGE FILE SIZE (also some weird artifacts on arrow heads)

        set(zframe, 'Color', 'White')
        export_fig(zframe,filename+'_z.pdf','-pdf','-painters')
        set(xframe, 'Color', 'White')
        export_fig(xframe,filename+'_x.pdf','-pdf','-painters')
        % >> PERFECTION <3

    end

end

cd(CODE)

%% Reset path
rmpath ../lib/arrow
% rmpath ../lib/m2tex/src
rmpath ../lib/export_fig

%% Contains

function [Z,X] = double_frame(id,x,y,n,u)

    % load magnetic data
    magX = -load(sprintf("%dN_magX_U%s.txt",n,string(u)));
    magZ = -load(sprintf("%dN_magX_U%s.txt",n,string(u)));

    % load correlation data
    corX = load(sprintf("corr_x_R%d.txt",n-1));
    try 
        corX = corX(:,1:61);
    catch
        warning(id+"X")
    end
    corZ = load(sprintf("corr_z_R%d.txt",n-1));
    try 
        corZ = corZ(:,1:61);
    catch
        warning(id+"Z")
    end
    if u < 5
        corX = corX(:,u*10+1);
        corZ = corZ(:,u*10+1);
    else
        corX = corX(:,end);
        corZ = corZ(:,end);
    end

    %% AFMz
    % Represent correlations as color
    Z = figure('Name',id); z = corZ;
    draw_lattice(x,y,z,'radius'); hold on
    % Set camera angle and similar
    xlim([min(x)-0.5,max(x)+0.5])
    ylim([min(y)-0.5,max(y)+0.5])
    zlim([-1,1])
    view(-20,48)
    axis image
    axis off
    % Represent magnetization with arrows
    z = magZ;
    s = 1/max(x).*z(z>0); 
    c(1:7,:)=repmat(str2rgb("watermelon"),7,1);
    colormap(c(1:7,:))
    arrow3([x(z>0),y(z>0),0*z(z>0)],[x(z>0),y(z>0),0.5*z(z>0)],'|',3*s,3*s);
    s = 1/max(x).*z(z<0); 
    c(1:7,:)=repmat(str2rgb("azure"),7,1);
    colormap(c(1:7,:))
    arrow3([x(z<0),y(z<0),0*z(z<0)],[x(z<0),y(z<0),0.5*z(z<0)],'|',3*s,3*s);

    %% AFMx
    % Represent correlations as color
    X = figure('Name',id); z = corX;
    draw_lattice(x,y,z,'radius'); hold on
    % Set camera angle and similar
    xlim([min(x)-0.5,max(x)+0.5])
    ylim([min(y)-0.5,max(y)+0.5])
    zlim([-1,1])
    view(-20,48)
    axis image
    axis off
    % Represent magnetization with arrows
    z = magX;
    s = 1/max(x).*z(z>0);
    c(6:9,:)=repmat(str2rgb("Watermelon"),4,1);
    colormap(c(6:9,:))
    arrow3([x(z>0),y(z>0),0*z(z>0)],[x(z>0)+0.5*z(z>0),y(z>0),0*z(z>0)],'|',3*s,3*s);
    s = 1/max(x).*z(z<0); 
    c(6:9,:)=repmat(str2rgb("azure"),4,1);
    colormap(c(6:9,:))
    arrow3([x(z<0),y(z<0),0*z(z<0)],[x(z<0)+0.5*z(z<0),y(z<0),0*z(z<0)],'|',3*s,3*s);


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
        r = sqrt(abs(z))*.8;
        c = repmat(0.0,N,1);
    case("color")
        r = repmat(0.5,N,1);
        c = z;
    case("flat")
        r = repmat(0.5,N,1);
        c = repmat(0.0,N,1);
    end
    sites = cell(N,1);
    ruler = linspace(0,max(c),25);
    color = palette.crameri('vanimo',100);
    color = color(51:75,:);
    for i = 1:N
        pos = [x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)];
        sites{i} = rectangle('Position',pos,'Curvature',[1 1],varargin{:});
        [~,minloc] = min(abs(ruler-c(i)));
        sites{i}.EdgeColor = [color(minloc,:)];
        hold on
    end
    hold off
end

