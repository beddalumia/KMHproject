%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear,clc

whichMF = 'AFMz';  % 'AFMz' | 'AFMxy' | 'AFMxyz'

% Select order parameter (varID==0 means everything)
varID = 2; 

% Select mode ('line' | 'map')
mode = 'map';

%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validate input
switch whichMF

    case 'AFMz'

        maxID = 2;

    case 'AFMxy'

        maxID = 4;

    case 'AFMxyz'

        maxID = 6;

otherwise

    error('Nonvalid MF-scheme!');

end

if (varID > maxID) || (varID < 0)

    error('Nonvalid varID!');

else
    
    mustBeInteger(varID);

end

% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/KMH-MF_Data/';
cd([DATA,whichMF]);

% Plotting
switch mode

    case 'line'

        phase_line(varID); 

    case 'map'

        switch whichMF
            
            case 'AFMz'
                
                phase_map_1(varID); 
                
            case 'AFMxy'
                
                phase_map_2(varID,varID+1);
                
            otherwise
                
                error('AFMxyz still not implemented!')
        end

otherwise

    error('Nonvalid mode!'); 

end

% Dirty path reset
cd(CODE);

%% Contains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Single Phase-Lines
function phase_line(varID)
    % Select ordpmservable 
    iOBS = varID; % (iOBS==0 means everything)
    % Get SOI value list
    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    Nlines = length(SOI_list);
    for iSOI = 1:Nlines
        lineID = SOI_names(iSOI);
        cd(lineID); fprintf(lineID);
        clear('ids','ordpms','U_list');
        load('order_parameter_line.mat','ids','ordpms','U_list');
        if iOBS~=0
            figName = strcat(lineID,' | ',ids{iOBS})';
            figure("Name",figName);
            plot(U_list,ordpms{iOBS},'LineWidth',2); drawnow
        else
            for iOBS = 1:length(ids)
                figName = strcat(lineID,' | ',ids{iOBS})';
                figure("Name",figName);
                plot(U_list,ordpms{iOBS},'LineWidth',2); drawnow
            end
        end
        cd('..'); fprintf('..DONE!\n\n');
    end  
end
    
%% Full Phase Diagram | Just one channel
function [ax,cb] = phase_map_1(varID)
    if varID == 0
       error('All ordpmservables option not allowed for phase maps!') 
    end
    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    Nlines = length(SOI_list);
    phaseVAR = cell(Nlines,1);
    transLine = zeros(2,Nlines);
    for iSOI = 1:Nlines
        lineID = SOI_names(iSOI);
        cd(lineID);
        clear('ids','ordpms','U_list');
        load('order_parameter_line.mat','ids','ordpms','U_list');
        if iSOI==1
            figure("Name",ids{varID});
            ax = axes;
        end
        % Get the map data
        U = U_list;
        SOI = SOI_list(iSOI)*ones(length(U),1);
        phaseVAR{iSOI} = ordpms{varID};
        z = phaseVAR{iSOI};
        % Plot the map
        %Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
        Z(iSOI,:) = z;
        Plt = plot3(ax,SOI,U,z,'.','LineWidth',2,'Color',hex2rgb('E91D63')); hold on
        cd('..');
    end
    % Imagesc plot
    X = SOI_list;
    Y = U;
    imagescn(X,Y,smoothdata(Z','SmoothingFactor',0.1));
    % Title, legend, all of that
    zlabel(ax,'AFMz order parameter','Interpreter','latex');
    zticks([]); caxis([0,1]);
    xlabel(ax,'$\lambda_{SO} / t$','Interpreter','latex');
    ylabel(ax,'$U / t$','Interpreter','latex');
    ax.Box = 'on';
    plotDMFT.colorlab_importall
    %colormap(ax,palette.brewer([],'PuOr'));
    palette.cmocean('matter')
    cb = colorbar(ax,'eastoutside');
    view(ax,-70,52); set(gca,'Color','none')
    axis tight    
end 

%% Full Phase Diagram | Two channels
function [ax,cb] = phase_map_2(varID1,varID2)
    if varID1 == 0 && varID2 == 0
       error('All ordpmservables option not allowed for phase maps!') 
    end
    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    Nlines = length(SOI_list);
    phaseVAR = cell(Nlines,1);
    for iSOI = 1:Nlines
        lineID = SOI_names(iSOI);
        cd(lineID);
        clear('ids','ordpms','U_list');
        load('order_parameter_line.mat','ids','ordpms','U_list');
        if iSOI==1
            figure("Name",sprintf('%s and %s',ids{varID1},ids{varID2}));
            ax = axes;
        end
        % Get the map data
        U = U_list;
        SOI = SOI_list(iSOI)*ones(length(U),1);
        phaseVAR1{iSOI} = ordpms{varID1}./2;
        phaseVAR2{iSOI} = ordpms{varID2}./2;
        z = sqrt((phaseVAR1{iSOI}).^2+(phaseVAR2{iSOI}).^2);
        % Plot the map
        %Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
        Z(iSOI,:) = z;
        Plt = plot3(ax,SOI,U,z,'.','LineWidth',2,'Color',hex2rgb('92D050')); hold on
        cd('..');
    end
    % Imagesc plot
    X = SOI_list;
    Y = U;
    imagescn(X,Y,Z');
    % Title, legend, all of that
    zlabel(ax,'AFMxy order parameter','Interpreter','latex');
    zticks([]); caxis([0,1]);
    xlabel(ax,'$\lambda_{SO} / t$','Interpreter','latex');
    ylabel(ax,'$U / t$','Interpreter','latex');
    ax.Box = 'on';
    plotDMFT.colorlab_importall
    %colormap(ax,palette.brewer([],'BrBG'));
    palette.cmocean('speed')
    cb = colorbar(ax,'eastoutside');
    view(ax,-70,52); set(gca,'Color','none')
    axis tight   
end 




