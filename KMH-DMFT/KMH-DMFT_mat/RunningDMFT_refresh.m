%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOW TO USE?
%
%  - Put this script in the home path of a pre-existent calculation
%    > A path already containing the folders for all the computed U values.
%  - Set-up the name of your *updated* driver program (no .f90 extension)
%    > e.g. driver = 'ed_kane_mele_SomeNewOption';
%  - Select doMPI (true.or.false) to run with openMPI or not
%  - Select ignConv (true.or.false) to ignore unvconverged points or not
%  - Define how much loops the refresh will have at most: Nnew int variable
%  - Run everything with $ matlab -batch RunningDMFT_refresh
%  - At the end you will find some additional output in the U=%f folders
%    > LOG_dmft_refr.txt giving just the stdout, for the new DMFT loops
%    > LOG_time_refr.txt giving the wall-clock time, for the new DMFT loops
%  - Also you may have some update to
%    > U_list.txt, if some points that were unconverged now are converged!
%      i.e. only the *newly* converged calculations will update the list :D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver = 'ed_kane_mele';    doMPI = true;    ignConv = true;    Nnew = 1;

% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Open file to write/update converged U-values
fileID_Ulist = fopen('U_list.txt','a');

%% Retrieve the list of the *converged* U-values
U_converged = unique(sort(load('U_list.txt')));

%% Build the list of all U-values or use the converged ones only
if isempty(U_converged) || ignConv == true
   [U_list, ~] = get_list('U');  
else
   U_list = U_converged;
end

%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for U = U_list'                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Define the U-folder name.
if isfolder(UDIR)
   cd(UDIR);                   % Enter the U-folder (if exists)
else
   errstr = 'U_list file appears to be inconsistent: ';
   errstr = [errstr,UDIR];
   errstr = [errstr,' folder has not been found.'];
   error(errstr);
end

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
mpi = 'mpirun ';                            % Control of MPI
else                                        % boolean flag...
mpi = [];
end
HUBBARD = sprintf(' uloc=%f',U);            % OVERRIDE of U-VALUE
NLOOP = sprintf(' nloop=%d',Nnew);          % OVERRIDE of #{loops}
outLOG = ' > LOG_dmft_refr.txt';            % STDOUT destination
dmft_ed_call = [mpi,driver,HUBBARD,NLOOP,outLOG];
tic
system(dmft_ed_call);                       % Fortran-call
chrono = toc;
file_id = fopen('LOG_time_refr.txt','w');   % Write on time-log
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A *NEWLY* CONVERGED DMFT POINT
if ~isfile('ERROR.README') && isempty(find(U_converged == U,1))
    fprintf(fileID_Ulist,'%f\n', U);        % Update U-list, only
end                                         % if *newly* converged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fileID_Ulist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional Functions
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
