%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global ignUlist
ignUlist = false; %  true  |  false  
% >>>>>>>>>>>> ACHTUNG <<<<<<<<<<<<
% ignoring the U_list will give
% problems for every nonconverged
% DMFT point (Ekin is evaluated 
% at convergence!)
% >>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<

[SOI_list, SOI_names] = get_list('SOI');
Nlines = length(SOI_list);
totalEnergies = cell(Nlines,1);
for iSOI = 1:Nlines
    SOI = SOI_list(iSOI);
    lineID = SOI_names(iSOI);
    fprintf('************\n');
    fprintf(lineID);
    fprintf('\n************\n');
    cd(lineID);
    if isfile('U_list.txt')
        U_list = load('U_list.txt');
    else
        U_list = [];
    end
    fprintf('> E_kin.');
    [kins,U_list] = kinetic_line(U_list); fprintf('..DONE\n');
    fprintf('> E_pot.');
    [ids,ens,U_list] = energy_line(U_list); fprintf('..DONE\n');
    pots = ens{1}; % The total potential energy is always the 1st value!
    totalEnergies{iSOI} = kins + pots;
    save('energy_line.mat','ids','ens','kins','U_list');
    % FIGURE ..............................................................
    if iSOI==1
        figure("Name",'DMFT Energy');
        ax = axes;
    end
    % Get the map data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    z = totalEnergies{iSOI};
    % Plot the map
    %Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
    Plt = plot3(ax,SOI,U,z,'g','LineWidth',2); hold on
    % .....................................................................
    cd('..');
end

% Title, legend, all of that
title(ax,'DMFT Energy');
xlabel(ax,'\lambda_{SO} / t');
ylabel(ax,'U / t');
ax.Box = 'on';
%colormap(ax,'copper');
%cb = colorbar(ax);
view(ax,-70,52);
fig = gcf;
fig.Renderer='Painters';
%set(gca, 'ZScale', 'log');



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

function [kins,U_list]  = kinetic_line(U_LIST)
%% Getting a list of energy values, from directories.
%  U_LIST: an array of values for Hubbard interaction U (could be empty!) 
%  kins: an array of values for Kinetic energies, forall U
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ignUlist
    if isempty(U_LIST) || ignUlist == true
       [U_LIST, ~] = get_list('U'); 
    else
       U_LIST = sort(U_LIST);
    end
    % Then we can proceed spanning all the U-values
    Nu = length(U_LIST);
    kins = zeros(Nu,1);
    for iU = 1:length(U_LIST)
        U = U_LIST(iU);
        UDIR= sprintf('U=%f',U);
        if ~isfolder(UDIR)
           errstr = 'U_list file appears to be inconsistent: ';
           errstr = [errstr,UDIR];
           errstr = [errstr,' folder has not been found.'];
           error(errstr);
        end
        cd(UDIR); 
        
        % The dmft_kinetic_energy.dat file has a weird structure...
        strCell = readcell('dmft_kinetic_energy.dat'); % cell of strings
        tempStr = strCell{1};                          % single string
        tempVec = sscanf(tempStr,'%f');                % extract all the %f
        kins(iU) = tempVec(1); % For sure the 1st value is the total K.E.
        
        cd('..');
    end
    
    U_list = U_LIST;
end

function [ids,ens,U_list]  = energy_line(U_LIST)
%% Getting a list of energy values, from directories.
%  U_LIST: an array of values for Hubbard interaction U (could be empty!)
%  ids: a cell of strings, the QcmPlab names for the pot-energy terms 
%  ens: a cell of float-arrays, corresponding to the names above
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    global ignUlist
    if isempty(U_LIST) | ignUlist == true
       [U_LIST, ~] = get_list('U'); 
    else
       U_LIST = sort(U_LIST);
    end
    % Then we can proceed spanning all the U-values
    Nu = length(U_LIST);
    cellEn = cell(Nu,1);
    for iU = 1:length(U_LIST)
        U = U_LIST(iU);
        UDIR= sprintf('U=%f',U);
        if ~isfolder(UDIR)
           errstr = 'U_list file appears to be inconsistent: ';
           errstr = [errstr,UDIR];
           errstr = [errstr,' folder has not been found.'];
           error(errstr);
        end
        cd(UDIR); 
        [ids, cellEn{iU}] = get_energies();
        cd('..');
    end
    % We need some proper reshaping
    Nids = length(ids);
    ens = cell(1,Nids);
    for jEn = 1:Nids
        ens{jEn} = zeros(Nu,1);
        for iU = 1:Nu
           ens{jEn}(iU) = cellEn{iU}(jEn);
        end
    end
    U_list = U_LIST;
end

function [names, energies] = get_energies()
%% Getting all information from energy_last.ed and energy_info.ed
%  names: a cell of strings, the QcmPlab names for the pot-energy terms
%  energies: an array of float values, corresponding to the names above
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    names = readcell('energy_info.ed','FileType','fixedwidth');
    names(strcmp(names,'#'))=[];
    for i = 1:length(names)
        tempstr = names{i};                 % Temporary string variable
        head = sscanf(tempstr,'%d');        % Extracts the initial integer
        head = int2str(head);               % Int to Str conversion
        names{i} = erase(tempstr,head);     % Proper beheading ;D
    end
    energies = load('energy_last.ed');
end