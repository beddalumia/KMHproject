%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global doPost
plotDMFT.colorlab_importall

    % Select directory
    whichDIR = 'AFMx'; %'AFMx', 'AFMz_normal', 'AFMz_replica'
    doPost   = 0    %  true  |  false  ->  Forces fresh postDMFT run

    % Select observable (varID==0 means everything)
    varID = 12;
    
    % Select mode ('line' | 'map')
    mode = 'map';
    % Do you want a transition line ( true | false )
    doTransLine = false;

    % Dirty path selector
    CODE = fileparts(mfilename('fullpath'));
    DATA = '../../Data/KMH-DMFT_Data/';
    cd(DATA)
    cd(whichDIR)
    
    % Plotting
    if strcmp(mode,'line')
       var = phase_line(varID); 
    elseif strcmp(mode,'map')
       var = phase_map(varID); 
    else
       error('Nonvalid mode!'); 
    end
    
    disp(cell2table(var))
    
    % Dirty path reset
    cd(CODE);
    
%% Single Phase-Lines
function ids = phase_line(varID)
                                                              global doPost
    % Select observable 
    iOBS = varID; % (iOBS==0 means everything)
    % Get SOI value list
    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    Nlines = length(SOI_list);
    for iSOI = 1:Nlines
        lineID = SOI_names(iSOI);
        cd(lineID); fprintf(lineID);
        fprintf('> Read observables.');
        if isfile('observables_line.mat') && not(doPost)
            load('observables_line.mat','ids','obs','U_list');
        else
            [ids,obs,U_list] = postDMFT.observables_line;
            save('observables_line.mat','ids','obs','U_list');
        end
        fprintf('..DONE\n');
        if iOBS~=0
            figName = strcat(lineID,' | ',ids{iOBS})';
            figure("Name",figName);
            plot(U_list,obs{iOBS},'LineWidth',2); drawnow
        else
            for iOBS = 1:length(ids)
                figName = strcat(lineID,' | ',ids{iOBS})';
                figure("Name",figName);
                plot(U_list,obs{iOBS},'LineWidth',2); drawnow
            end
        end
        cd('..'); fprintf('..DONE!\n');
    end  
end
    
%% Full Phase Diagram
function ids = phase_map(varID)
                                                              global doPost
% varID \in [1,15]
if varID == 0
   error('All observables option not allowed for phase maps!') 
end
[SOI_list, SOI_names] = postDMFT.get_list('SOI');
Nlines = length(SOI_list);
phaseVAR = cell(Nlines,1);
for iSOI = 1:Nlines
    lineID = SOI_names(iSOI);
    cd(lineID);
    fprintf('> Read observables.');
    if isfile('observables_line.mat') && not(doPost)
        load('observables_line.mat','ids','obs','U_list');
    else
        try %#ok
            [ids,obs,U_list] = postDMFT.observables_line;
        end
            save('observables_line.mat','ids','obs','U_list');
    end
    fprintf('..DONE\n');
    if iSOI==1
        figure("Name",ids{varID});
        ax = axes;
    end
    % Get the map data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    phaseVAR{iSOI} = obs{varID};
    z = abs(phaseVAR{iSOI});
    Z(iSOI,:) = z;
    % Plot the map
    %Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
    Plt = plot3(ax,SOI,U,z,'.','LineWidth',2,'Color',hex2rgb('92D050')); hold on
    %Plt = plot3(ax,SOI,U,z,'.','LineWidth',2,'Color',hex2rgb('E91D63')); hold on
    cd('..');
end
% Imagesc plot
X = SOI_list;
Y = U;
imagescn(X,Y,Z');
% Title, legend, all of that
%zlabel(ax,'AFMz order parameter','Interpreter','latex');
%zticks([]); 
zlim([0,1]);
zlabel(ax,ids{varID},'Interpreter','latex')
xlabel(ax,'$\lambda_{SO} / t$','Interpreter','latex');
ylabel(ax,'$U / t$','Interpreter','latex');
ax.Box = 'on';
palette.cmocean('speed')
%palette.cmocean('matter')
cb = colorbar(ax,'eastoutside'); caxis([0,1]);
view(ax,-70,52); set(gca,'Color','none')
axis tight 
end 


