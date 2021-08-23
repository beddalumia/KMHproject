%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

whichMF = 'AFMxy';  % 'AFMz' | 'AFMxy' | 'AFMxyz'
doBands = false;    %  true  |  false
doRaster = false;   %  true  |  false

cd(whichMF);

[SOI_list, SOI_names] = get_list('SOI');
Nlines = length(SOI_list);
phaseVAR = cell(Nlines,1);
transLine = zeros(2,Nlines);
for iSOI = 1:Nlines
    SOI = SOI_list(iSOI);
    lineID = SOI_names(iSOI);
    fprintf('************\n');
    fprintf(lineID);
    fprintf('\n************\n');
    cd(lineID);
    clear('ids','ordpms','U_list');
    load('order_parameter_line.mat','ids','ordpms','U_list');
    Npar = length(ordpms);
    Npoints = length(U_list);
    for iHubb = 1:Npoints
        U = U_list(iHubb);
        UDIR= sprintf('U=%f',U);
        cd(UDIR);
        %
        Emb(iSOI,iHubb) = 0;
        %
        Sx = 0;
        Sy = 0;
        Sz = 0;
        %
        Rx = 0;
        Ry = 0;
        Rz = 0;
        %
        Pz = 0;
        Nel = 2;
        %
        if      Npar == 2
            Sz = ordpms{1}(iHubb);
            Rz = ordpms{2}(iHubb);
        elseif  Npar == 4
            Sx = ordpms{1}(iHubb);
            Sy = ordpms{2}(iHubb);
            Rx = ordpms{3}(iHubb);
            Ry = ordpms{4}(iHubb);
        elseif  Npar == 6
            Sx = ordpms{1}(iHubb);
            Sy = ordpms{2}(iHubb);
            Sz = ordpms{3}(iHubb);
            Rx = ordpms{4}(iHubb);
            Ry = ordpms{5}(iHubb);
            Rz = ordpms{6}(iHubb);
        end
        %
        Emb(iSOI,iHubb) = Emb(iSOI,iHubb) + (Sz+Rz)^2 - (Pz+Nel)^2 + (Sz-Rz)^2 - (Pz-Nel)^2;
        Emb(iSOI,iHubb) = Emb(iSOI,iHubb) - (Sx+Rx)^2 - (Sy+Ry)^2 - (Sx-Rx)^2 - (Sy-Ry)^2;
        Emb(iSOI,iHubb) = Emb(iSOI,iHubb) * U/4;
        %
        Eigenbands = load('Eigenbands.nint');
        cd('..');
        Ncell = length(Eigenbands);
        Ev = Eigenbands(1:round(Ncell/2),:);
        fprintf('Computing GS energy for U=%f..',U);
        Egs(iSOI,iHubb) = sum(Ev(:,2))/Ncell;
        fprintf('.DONE!\n');
        %% Plotting Bands
        %------------------------------------------------------------------
        if doBands
            Ec = Eigenbands(round(Ncell/2):end,:);
            scatter(Ev(:,1),Ev(:,2),'r'); hold on
            scatter(Ec(:,1),Ec(:,2),'b');
            % Title, legend, all of that
            title(sprintf('Mean-Field Bands [SOI=%f U=%f]',SOI,U));
            xticks([0 2.418399 4.836798 7.245524])
            xlim([0,7.245524]);
            xticklabels({'\Gamma','K','K`','\Gamma'})
            ylabel('\epsilon / t');
            ax = gca; ax.Box = 'on';
        
           % [These two lines ensure filling of the fig]
            InSet = get(ax, 'TightInset');
            set(ax, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3),...
                1-InSet(2)-InSet(4)]);
            if doRaster
               filename = sprintf('MF_bands_SOI=%f_U=%f.png',SOI,U);
            else
               filename = sprintf('MF_bands_SOI=%f_U=%f.pdf',SOI,U);            
            end
            fprintf('Printing %s..',filename);
            if doRaster
               print(gcf,filename,'-dpng','-r600');
            else
               print(gcf,filename,'-dpdf','-fillpage');
            end
            fprintf('.DONE!\n');
        end
        %------------------------------------------------------------------
    end
    cd('..');
end
 
cd('..');

% Get the map data
[X,Y] = meshgrid(SOI_list,U_list);
Z = Egs' + Emb';
% Plot the map data
surf(X,Y,Z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [flist, strlist] = get_list(VARNAME)
%% Getting a list of variable values, from directories.
%  VARNAME: a string, identifying the listed variable (e.g. 'U')
%  flist: an array of float_values (e.g. U=[:] )
%  strlist: an array of dir_name strings (e.g. ['U=%f'] )
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subthings = dir('.'); % Retrieves every subdir and file inside pwd
    subfolders = subthings([subthings(:).isdir]); % Keeps only subfolders
    subfolders = subfolders(~ismember({subfolders(:).name},{'.','..'}));
    N = length(subfolders); flist = zeros(N,1); strlist = strings(N,1);
    for i = 1:N
        DIR = subfolders(i).name; % Let's get the indexed string...
        flist(i) = sscanf(DIR, [VARNAME,'=%f']); %...and extract the value!
        strlist(i) = DIR;
    end
    % We need to sort the lists by floats (not strings, as it is now)
    [flist, sortedIDX] = sort(flist); strlist = strlist(sortedIDX);
end