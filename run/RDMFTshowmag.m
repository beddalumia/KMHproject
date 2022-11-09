Ulist = postDMFT.get_list('U');
lattice_file = 'flake.txt';
gif_name = 'animation.gif';
plotDMFT.import_colorlab;

if not(isempty(Ulist))
    for iU = 1:length(Ulist)
        
        U = Ulist(iU);
        UDIR= sprintf('U=%f',U);
        if ~isfolder(UDIR)
            errstr = 'U_list appears to be inconsistent: ';
            errstr = [errstr,UDIR];
            errstr = [errstr,' folder has not been found.'];
            error(errstr);
        end
        cd(UDIR);
        
        frame = single_frame(U,lattice_file);
        plotDMFT.push_frame(gif_name,iU,length(Ulist),0.1);
        close(frame)
        
        cd ..
        
    end
else
    % If Ulist is empty then we are inside a Udir already and want just a frame
    [~,Udir,decimals] = fileparts(pwd); % hate this
    Udir = [Udir,decimals];
    Uval = sscanf(Udir,'U=%f');
    frame = single_frame(Uval,lattice_file);
end

function h = single_frame(U,lattice_file)

    % POSIX dependent, sorry
    system("awk '{print $5}' observables_last_ineq0* > magX.txt");
    z = load('magX.txt');
    delete('magX.txt');

    % Load coordinates
    lattice = load(lattice_file);
    x = lattice(:,1);
    y = lattice(:,2);

    % Represent as double bubblechart
    h = figure('Name',"U="+string(U));
    %bubblechart(x(z>0),y(z>0),sqrt(z(z>0)),rgb.xkcd('red'));
    draw_lattice(x(z>0),y(z>0),sqrt(z(z>0)/2),FaceColor=rgb.xkcd('red'));
    hold on
    %bubblechart(x(z<0),y(z<0),sqrt(-z(z<0)),rgb.xkcd('blue'));
    draw_lattice(x(z<0),y(z<0),sqrt(-z(z<0)/2),FaceColor=rgb.xkcd('blue'));
    axis equal
    axis off
    title("U="+string(U));

end

function sites = draw_lattice(x,y,r,varargin)
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
    N = length(r);
    sites = cell(N,1);
    for i = 1:N
        pos = [x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)];
        sites{i} = rectangle('Position',pos,'Curvature',[1 1],varargin{:});
        hold on
    end
    hold off
end
