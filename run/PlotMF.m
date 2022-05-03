%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

whichMF = 'AFMxy';  % 'AFMz' | 'AFMxy' | 'AFMxyz'

% Select order parameter (varID==0 means everything)
varID = 3; 

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
DATA = '../../Data/KMH-MF_Data/';
cd([DATA,whichMF]);

% Plotting
switch mode

    case 'line'

        phase_line(varID); 

    case 'map'

        phase_map(varID); 

otherwise

    error('Nonvalid mode!'); 

end

% Dirty path reset
CODE = '../../../KMproj[git]/run/';
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
    
%% Full Phase Diagram | Just one channel (spin resolved in testing/)
function [ax,cb] = phase_map(varID)
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
    % Get the line data
    try
        ztrim = z(z<1e-2);
        ztrans = max(ztrim);
        transID = find(z==ztrans);
        transLine(2,iSOI) = U(transID);
        transLine(1,iSOI) = SOI(transID);
    catch
        transLine(2,iSOI) = NaN;
        transLine(1,iSOI) = NaN;     
    end
    % Plot the map
    Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
    Plt = plot3(ax,SOI,U,z,'k','LineWidth',2); hold on
    cd('..');
end
% Plot transition line
plot3(transLine(1,:),transLine(2,:),min(z)*ones(2,Nlines),'r','LineWidth',2.5);
% Title, legend, all of that
title(ax,ids{varID});
xlabel(ax,'$\lambda_{SO} / t$','Interpreter','latex');
ylabel(ax,'$U / t$','Interpreter','latex');
ax.Box = 'on';
colormap(ax,plotDMFT.colorlab.brewer.cmap([],'BrBG'));
cb = colorbar(ax);
view(ax,-70,52);
fig = gcf;
fig.Renderer='Painters';
clc    
end 


