%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

    % Select observable (varID==0 means everything)
    varID = 5;
    
    % Select mode ('line' | 'map')
    mode = 'line';
    % Do you want a transition line ( true | false )
    doTransLine = false;
    
    %%% Structure: obs{varID} (for ed_kane_mele) %%%%%%%%%%%%%%%%
    % obs{1}    Density (at half filling = 1)                   %
    % obs{2}    Docc (Double Occupancy \in [0,1])               %
    % obs{3}    Nup (intensive)                                 %
    % obs{4}    Ndown (intensive)                               %
    % obs{5}    Magn(Nup-Ndown)                                 %
    % obs{6}    S2 (Impurity magnetic dipole: <S^2>)            %
    % obs{7}    Egs (Ground-State Energy)                       %
    % obs{8}    Sz2_11 (Impurity magnetic dipole:<Sz^2>)        %
    % obs{9}    N2_11 (??????)                                  %
    % obs{10}   Z: Quasiparticle Weight for spin-up             %
    % obs{11}   Z: Quasiparticle Weight for spin-down           %
    % obs{12}   Sig1_s1(\Im\Sigma(iw=0) for spin-up)            %
    % obs{13}   Sig1_s1(\Im\Sigma(iw=0) for spin-down)          %
    % obs{14}   Nph (??????)                                    %
    % obs{15}   Wph (??????)                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting
    
    if strcmp(mode,'line')
       phase_line(varID); 
    elseif strcmp(mode,'map')
       phase_map(varID,doTransLine); 
    else
       error('Nonvalid mode!'); 
    end
    
%% Single Phase-Lines
function phase_line(varID)
    % Select observable 
    iOBS = varID; % (iOBS==0 means everything)
    % Get SOI value list
    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    Nlines = length(SOI_list);
    for iSOI = 1:Nlines
        lineID = SOI_names(iSOI);
        cd(lineID); fprintf(lineID);
        clear('ids','obs','U_list');
        load('observables_line.mat','ids','obs','U_list');
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
    
%% Full Phase Diagram | Just one channel (spin resolved in testing/)
function phase_map(varID,doTransLine)
% varID \in [1,15]
if varID == 0
   error('All observables option not allowed for phase maps!') 
end
[SOI_list, SOI_names] = postDMFT.get_list('SOI');
Nlines = length(SOI_list);
phaseVAR = cell(Nlines,1);
transLine = zeros(2,Nlines);
for iSOI = 1:Nlines
    lineID = SOI_names(iSOI);
    cd(lineID);
    clear('ids','obs','U_list');
    load('observables_line.mat','ids','obs','U_list');
    if iSOI==1
        figure("Name",ids{varID});
        ax = axes;
    end
    % Get the map data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    phaseVAR{iSOI} = obs{varID};
    z = phaseVAR{iSOI};
    % Get the line data
    if doTransLine
    ztrim = z(z<1e-2);
    ztrans = max(ztrim);
    transID = find(z==ztrans);
    transLine(2,iSOI) = U(transID);
    transLine(1,iSOI) = SOI(transID);
    end
    % Plot the map
    Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
    Plt = plot3(ax,SOI,U,z,'g','LineWidth',2); hold on
    cd('..');
end
% Plot transition line
if doTransLine
plot3(transLine(1,:),transLine(2,:),min(z)*ones(2,Nlines),'r','LineWidth',2.5);
end
% Title, legend, all of that
title(ax,ids{varID});
xlabel(ax,'\lambda_{SO} / t');
ylabel(ax,'U / t');
ax.Box = 'on';
colormap(ax,'copper');
cb = colorbar(ax);
view(ax,-70,52);
fig = gcf;
fig.Renderer='Painters';
%set(gca, 'ZScale', 'log');
clc    
end 


