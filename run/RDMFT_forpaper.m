clear,clc

%% INPUT
U = [4.5,6.0];
M = ["X","Z"];

%% LIBRARIES
QcmP.plot.import_colorlab;
addpath ../lib/arrow

% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/RDMFT/for_paper';
cd(DATA)

    for i = 2:5%2:5
        
        flake = load(sprintf("%dN_flake.txt",i));

        x = flake(:,1);
        y = flake(:,2);

        for j = 2:2%2

            for k = 1:2%1:2

                mag = load(sprintf("%dN_mag%s_U%s.txt",i,M(k),string(U(j))));
                id = sprintf("AFM$_%s$ %dN flake @ $U=%s$",M(k),i,string(U(j)));

                frame = single_frame(id,x,y,mag,M(k));
                
            end

        end

    end

cd(CODE)
rmpath ../lib/arrow

function h = single_frame(id,x,y,z,ax)

    % Represent as double bubblechart
    h = figure('Name',id);
    %draw_lattice(x(z>0),y(z>0),sqrt(z(z>0)/2),'radius''FaceColor',rgb.xkcd('red'));
    %draw_lattice(x(z>0),y(z>0),sqrt(z(z>0)/2),'color','FaceColor',rgb.xkcd('red'));
    draw_lattice(x(z>0),y(z>0),sqrt(z(z>0)/2),'flat');
    hold on
    %draw_lattice(x(z<0),y(z<0),sqrt(-z(z<0)/2),'radius','FaceColor',rgb.xkcd('blue'));
    %draw_lattice(x(z<0),y(z<0),sqrt(-z(z<0)/2),'color','FaceColor',rgb.xkcd('blue'));
    draw_lattice(x(z<0),y(z<0),sqrt(-z(z<0)/2),'flat');
    title(id,'Interpreter','latex');

    xlim([min(x)-0.5,max(x)+0.5])
    ylim([min(y)-0.5,max(y)+0.5])
    hold on
    switch ax
        case('Z')
            s = 3/max(x).*z(z>0);
            arrow3([x(z>0),y(z>0),0*z(z>0)],[x(z>0),y(z>0),z(z>0)],'r',s,s);
            s = 3/max(x).*z(z<0);
            arrow3([x(z<0),y(z<0),0*z(z<0)],[x(z<0),y(z<0),z(z<0)],'b',s,s);
        case('X')
            s = 3/max(x).*z(z>0);
            arrow3([x(z>0),y(z>0),0*z(z>0)],[x(z>0)+0.5*z(z>0),y(z>0),0*z(z>0)],'r',s,s);
            s = 3/max(x).*z(z<0);
            arrow3([x(z<0),y(z<0),0*z(z<0)],[x(z<0)+0.5*z(z<0),y(z<0),0*z(z<0)],'b',s,s);
    end
    zlim([-1,1])
    axis off

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

