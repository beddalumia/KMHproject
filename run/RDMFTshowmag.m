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
        plotDMFT.push_frame(gif_name,iU,length(Ulist),0.1,frame);
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
    bubblechart(x(z>0),y(z>0),sqrt(z(z>0)),rgb.xkcd('red'));
    hold on
    bubblechart(x(z<0),y(z<0),sqrt(-z(z<0)),rgb.xkcd('blue'));
    axis equal
    axis off
    title("U="+string(U));
    
    % What to do about xlim and ylim??

end
