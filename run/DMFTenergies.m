clear,clc

%% INPUT

doPost   = 0    %  true  |  false  ->  Forces fresh QcmP.post run
doSmooth = 1    %  true  |  false  ->  Smoothing on energy differences
ignConv  = 1    %  true  |  false  ->  Ignores converged-U lists

mode = 'map'    % 'line' | 'map' (unimplemented: only maps)

whichDIR = {'AFMx','AFMz_rotateMag'}  %  where to look for DMFT-lines

%% MAIN
cd
% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/DMFT/';
cd(DATA)

for iDIR = 1:length(whichDIR)

    cd(whichDIR{iDIR});
    disp(whichDIR{iDIR});

    [SOI_list, SOI_names] = QcmP.post.get_list('SOI');
    SOI_list = SOI_list(1:end);
    SOI_names = SOI_names(1:end);
    Nlines = length(SOI_list);
    tots = cell(Nlines,1);
    U_list_old = [];
    Etot = [];

    for iSOI = 1:Nlines

        SOI = SOI_list(iSOI);
        lineID = SOI_names(iSOI);

        fprintf('••••••••••••\n');
        fprintf(lineID);
        fprintf('\n••••••••••••\n');

        cd(lineID);

        if isfile('U_conv.txt') && not(ignConv)
            U_list = QcmP.post.get_list('U');
            U_conv = load('U_conv.txt');
            U_shit = setdiff(U_list,U_conv);
        else
            U_list = QcmP.post.get_list('U');
            U_conv = U_list;
            U_shit = [];
        end
        if iSOI > 1
            U_shit = [U_shit,setdiff(U_list,U_list_old)];
        end
        U_list_old = U_list;

        fprintf('> E_kin.');
        if isfile('kinetic_energy.txt') && not(doPost)
            kins = load('kinetic_energy.txt');
        else
            kins = QcmP.post.kinetic_line('U',U_list); 
        end
        fprintf('..DONE\n');

        fprintf('> E_pot.');
        if isfile('Hi.txt') && not(doPost)
            pots = load('Hi.txt');
        else
            [~,ens] = QcmP.post.energy_line('U',[],U_list); 
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

        Etot(iSOI,:) = z(~ismember(U,U_shit));
        
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
    close
    
    %••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
    %% 2D-SURFACES 
    %••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
    % Get the map data
    [X,Y] = meshgrid(SOI_list,U_list); % Assaad's convention
    Z = Etot'; 
    if doSmooth, Z = smoothdata(Z,'SmoothingFactor',0.1); end
    % Plot the map data
    QcmP.plot.import_colorlab
    figure("Name","DMFT Energy: surfaces "+whichDIR{iDIR});
    switch whichDIR{iDIR}
        case 'AFMx'
            AFMx = Z;
            surf(X,Y,AFMx,'FaceAlpha',0.5);
            colormap(palette.cmocean('-matter'))
        case 'AFMz_rotateMag'
            AFMz  = Z;
            surf(X,Y,AFMz,'FaceAlpha',0.5);
            colormap(palette.cmocean('-algae'))
        otherwise
            AFM   = Z;
            surf(X,Y,AFM,'FaceAlpha',0.5);
            colormap(palette.cmocean('tempo'))
    end
    view(136,30)
    close
    % •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

    cd ..

end

figure('Name','GSenergy difference @ DMFT')
title('Groundstate Energy Difference [AFMx - AFMz]','Interpreter','latex')
ylabel('$U/t$','Interpreter','latex')
xlabel('$\lambda_{SO}/t$','Interpreter','latex')
zlabel('$E_0$(AFMx) - $E_0$(AFMz)','Interpreter','latex')
meshz(X,Y,AFMx(:,1:10)-AFMz,'EdgeColor','none','FaceColor','interp');%,str2rgb('powder blue'));
palette.cmocean('matter');
view(110,30)
xlim([0,0.3]); 
xlabel('$\lambda_\mathrm{so}/t$','Interpreter','latex');
ylim([0,10]);
ylabel('$U/t$','Interpreter','latex');
zlim([-0.0565,0.0020]);
zlabel('$E_\parallel - E_\perp$','Interpreter','latex');
caxis([-0.0565,0.0020]);
grid off

% figure('Name','GSenergy difference @ DMFT')
% title('Groundstate Energy Difference [AFMx - AFMz]','Interpreter','latex')
% xlabel('$U/t$','Interpreter','latex')
% ylabel('$\lambda_{SO}/t$','Interpreter','latex')
% imagescn(SOI_list,U_list,AFMx(:,1:10)-AFMz);
% palette.cmocean('matter');
% colorbar
%caxis([-0.1,0])


% Dirty path reset
cd(CODE);

% Export to LaTeX
addpath ../lib/m2tex/src
matlab2tikz('filename','../fig/DMFT_energy.tex','width','\textwidth');
rmpath ../lib/m2tex/src

