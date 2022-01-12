%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global ignUlist
ignUlist = false;

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
    [ids,ordpms,U_list] = postDMFT.order_parameter_line(U_list); fprintf('..DONE\n');
    save('order_parameter_line.mat','ids','ordpms','U_list');
    cd('..');
end




