% !!! > DO NOT USE, STILL BROKEN! (don't know yet how to get those slurms)

% ALSO I HAVE NO IDEA ON HOW TO STOP SENDING JOBS IF DMFT DOES NOT CONVERGE!
% -> most probably you *must* control it upstream from slurm, bleah D:

% Let MATLAB see the goddamn environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
arr_id    = str2double(getenv('SLURM_ARRAY_TASK_ID')); 
arr_count = str2double(getenv('SLURM_ARRAY_TASK_COUNT')); 
arr_max   = str2double(getenv('SLURM_ARRAY_TASK_MAX'));
arr_min   = str2double(getenv('SLURM_ARRAY_TASK_MIN'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase Line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOI  = 0.01;                   % Input Spin-Orbit 
Umin = 0; Umax = 10;           % Input Hubbard 
Ustep = (Umax-Umin)/arr_count; % Hubbard Step
Wmix = 0.3;		 	% Input Self-Mixing

Ulist = [-1,Umin:Ustep:Umax];  % Hubbard loop [controlled by slurm]~~~~~~>
U = Ulist(arr_id+1);		% <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Uold=Ulist(arr_id);	

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
if isfolder(oldDIR)            % If it exist a "previous" folder: 
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../input*             % Copy inside the **external** input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
mpi = 'mpirun ';
driver = 'ed_kane_mele';
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE
T2 =sprintf(' t2=%f',SOI);			% of
MIXING = sprintf(' wmixing=%f',Wmix);		% PARAMETERS
outLOG = ' > LOG.dmft';
dmft_ed_call = [mpi,driver,HUBBARD,T2,MIXING,outLOG];
tic
system(dmft_ed_call);
chrono = toc;
file_id = fopen('LOG.time','w');
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE NEED TO CATCH SOMEHOW A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README') 
    ed_status = -1
    return ed_status % <--- I THINK HERE IS THE PROBLEM
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder

ed_status = 0
return ed_status % <--- I THINK HERE IS THE PROBLEM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% - It's broken;
% - I have no idea on how to stop sending concatenated jobs it DMFT does not converge a point...
%   > most probably you *must* control it upstream from slurm, bleah D:
%
% => Complete redesign may follow, either
%       > the single-point DMFT evaluation has to be a *function* to which we feed the SLURM environment variables; 
%       > no use of slurm arrays: instead a "main" matlab script has to call both a SLURM-parser-function (using $sbatch instead of
%         #SBATCH) and the DMFT executable. Appears to be ugly but should work fine. I think I prefer this second solution.
