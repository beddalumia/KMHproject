clear all
clc

global ignUlist
ignUlist = false;

whichMF = 'AFMxy';  % 'AFMz' | 'AFMxy' | 'AFMxyz'

% Dirty path selector
DATA = '../../../Data/KMH-MF_Data/';
cd([DATA,whichMF]);

% We don't have a SOI-values list, but we can obtain that by just
% inspecting the subdirectories...
[SOI_list, SOI_names] = postDMFT.get_list('SOI');
Nlines = length(SOI_list);
for iSOI = 1:Nlines
    SOIDIR = SOI_names(iSOI);
    fprintf(strcat(SOIDIR));
    cd(SOIDIR);
    if isfile('U_list.txt')
        U_list = load('U_list.txt');
    else
        U_list = [];
    end
    [ids,ordpms,U_list] = postDMFT.order_parameter_line(U_list); fprintf('...DONE\n');
    save('order_parameter_line.mat','ids','ordpms','U_list');
    mfens = mf_energy_line(U_list);
    cd('..');
end

% Dirty path reset
CODE = '../../../KMproj[git]/KMH-MF/KMH-MF_mat/';
cd(CODE);

%% contains

function [mfens,U_list]  = mf_energy_line(U_LIST)
    %% Getting a list of energy values, from directories.
    %
    %     [mfens,U_list]  = mf_energy_line(U_LIST)
    %
    %  U_LIST: an array of values for Hubbard interaction U (could be empty!)
    %  mfens:  an array of values for mean-field energies, forall U
    %  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ~exist('U_LIST','var') || isempty(U_LIST)
           [U_LIST, ~] = postDMFT.get_list('U'); 
        else
           U_LIST = sort(U_LIST);
        end
        % Then we can proceed spanning all the U-values
        Nu = length(U_LIST);
        mfens = zeros(Nu,1);
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
            
            % Loading the relevant file
            mfens(iU) = load('mean_field_gs.dat');
            
            cd('..');
        end
        filename = 'mf_energy.txt';
        writematrix(mfens,filename,'Delimiter','tab');
        U_list = U_LIST;
    end


