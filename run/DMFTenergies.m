clear,clc

%% INPUT

doPost   = 0    %  true  |  false  ->  Forces fresh postDMFT run
doSmooth = 1    %  true  |  false  ->  Smoothing on energy differences
ignConv  = 0    %  true  |  false  ->  Ignores converged-U lists

mode = 'map'    % 'line' | 'map' (unimplemented: only maps)

whichDIR = {'AFMz','AFMx'}  %  where to look for DMFT-lines

%% MAIN

% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/KMH-DMFT_Data/';
cd(DATA)

for iDIR = 1:length(whichDIR)

    cd(whichDIR{iDIR});
    disp(whichDIR{iDIR});

    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    SOI_list = SOI_list(1:end);
    SOI_names = SOI_names(1:end);
    Nlines = length(SOI_list);
    tots = cell(Nlines,1);

    for iSOI = 1:Nlines

        SOI = SOI_list(iSOI);
        lineID = SOI_names(iSOI);

        fprintf('••••••••••••\n');
        fprintf(lineID);
        fprintf('\n••••••••••••\n');

        cd(lineID);

        if isfile('U_conv.txt') && not(ignConv)
            U_list = postDMFT.get_list('U');
            U_conv = load('U_conv.txt');
            U_shit = setdiff(U_list,U_conv);
        else
            U_list = postDMFT.get_list('U');
            U_conv = U_list;
            U_shit = [];
        end

        fprintf('> E_kin.');
        if isfile('kinetic_energy.txt') && not(doPost)
            kins = load('kinetic_energy.txt');
        else
            kins = postDMFT.kinetic_line(U_list); 
        end
        fprintf('..DONE\n');

        fprintf('> E_pot.');
        if isfile('Hi.txt') && not(doPost)
            pots = load('Hi.txt');
        else
            [~,ens] = postDMFT.energy_line(U_list); 
            pots = ens{1};
        end
        fprintf('..DONE\n');

        tots{iSOI} = kins + pots; % Total <GS|H|GS>
        writematrix(tots{iSOI},'groundstate_energy.txt')

        %••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
        %% 2D-FIGURES
        %••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
        if iSOI==1
            figure("Name","DMFT Energy: lines "+whichDIR{iDIR});
            ax = axes;
        end
        % Get the map data
        U = U_list;
        SOI = SOI_list(iSOI)*ones(length(U),1);
        z = tots{iSOI};
        %z(ismember(U,U_shit)) = NaN;
        %z = repnan(z,'spline');
        % Plot the map
        Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
        Plt = plot3(ax,SOI,U,z,'g','LineWidth',2); hold on
        % •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

        Etot(iSOI,:) = z;
        
        cd('..');

    end

    % Title, legend, all of that
    title(ax,"DMFT Energy ["+whichDIR{iDIR}+"]");
    xlabel(ax,'\lambda_{SO} / t');
    ylabel(ax,'U / t');
    ax.Box = 'on';
    colormap(ax,'copper');
    cb = colorbar(ax);
    view(ax,-70,52);
    fig = gcf;
    fig.Renderer='Painters';
    %set(gca, 'ZScale', 'log');
    
    %••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
    %% 2D-SURFACES 
    %••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
    % Get the map data
    [X,Y] = meshgrid(SOI_list,U_list); % Assaad's convention
    Z = Etot'; 
    if doSmooth, Z = smoothdata(Z,'SmoothingFactor',0.05); end
    % Plot the map data
    plotDMFT.colorlab_importall
    figure("Name","DMFT Energy: surfaces "+whichDIR{iDIR});
    switch whichDIR{iDIR}
        case 'AFMx'
            AFMx = Z;
            surf(X,Y,AFMx,'FaceAlpha',0.5);
            colormap(palette.cmocean('-matter'))
        case 'AFMz'
            AFMz  = Z;
            surf(X,Y,AFMz,'FaceAlpha',0.5);
            colormap(palette.cmocean('-algae'))
        otherwise
            AFM   = Z;
            surf(X,Y,AFM,'FaceAlpha',0.5);
            colormap(palette.cmocean('tempo'))
    end
    view(136,30)
    % •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

    cd ..

end

if size(AFMx) == size(AFMz)
    figure('Name','GSenergy difference @ DMFT')
    title('Groundstate Energy Difference [AFMx - AFMz]','Interpreter','latex')
    ylabel('$U/t$','Interpreter','latex')
    xlabel('$\lambda_{SO}/t$','Interpreter','latex')
    zlabel('$E_0$(AFMx) - $E_0$(AFMz)','Interpreter','latex')
    meshz(X,Y,AFMx-AFMz,'EdgeColor','k','FaceColor',str2rgb('powder blue'));
    view(136,30)
%     figure('Name','GSenergy difference @ DMFT')
%     title('Groundstate Energy Difference [AFMx - AFMz]','Interpreter','latex')
%     xlabel('$U/t$','Interpreter','latex')
%     ylabel('$\lambda_{SO}/t$','Interpreter','latex')
%     imagescn(SOI_list,U_list,AFMx-AFMz);
%     palette.cmocean('matter');
%     colorbar
%     caxis([-0.1,0])
end

% Dirty path reset
cd(CODE);

