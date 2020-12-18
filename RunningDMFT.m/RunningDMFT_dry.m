%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOW TO USE?
%
%  - Put this script, together with an input-file for your model, in a path
%    > This path will contain directories for all the U values you set-up.
%  - Set-up the name of your driver program (without .f90 extension)
%    > e.g. driver = 'ed_kane_mele';
%  - Adjust Umin, Umax and Ustep to your desire --> U = Umin:Ustep:Umax
%  - Set SOI to your desire: --> you will get a fixed-SOI linear span
%  - Adjust Uold to catch a 'restart-folder' in the path [!applies -> -1]
%  - Select doMPI (true.or.false) to run with openMPI or not
%  - Run everything with $ matlab -batch RunningDMFT_dry
%  - At the end you will find some additional output in the U=%f folders
%    > a LOG_dmft.txt which is just a mirror of the DMFT output (via tee)
%    > a LOG_time.txt which is a wall-clock-time value for the whole DMFT
%  - Also an additional output file in the main (external) path
%    > a U_list.txt that stores all the used U-values (for easier post..)
%      NB. Only the converged calculations will update the U_list :D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver = 'ed_kane_mele';	doMPI = false;

% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ulist = fopen('U_list.txt','a');

%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOI  = 0;                      % Input Spin-Orbit 
Umin = 1; Umax = 10;           % Input Hubbard 
Ustep = 0.25;			% Phase-line step
Uold = -1;			% Restart option

U = Umin; 
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
if isfolder(oldDIR)            % If it exist a "previous" folder: 
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../inputKANEMELE.conf % Copy inside the **external** input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
mpi = 'mpirun ';			        % Control of MPI
else					        % boolean flag...
mpi = [];
end
HUBBARD =sprintf(' uloc=%f',U);	        % OVERRIDE of
T2 =sprintf(' t2=%f',SOI);		        % PARAMETERS
outLOG = ' > LOG.out';
dmft_ed_call = [mpi,driver,HUBBARD,T2,outLOG];
tic
system(dmft_ed_call);				% Fortran-call
chrono = toc;
file_id = fopen('LOG_time.txt','w');		% Write on time-log
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
if ~isfile('ERROR.README') 
    fprintf(Ulist,'%f\n', U);	          % Update U-list (only if converged)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder

Uold = U;
U = U + Ustep;              	% Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(Ulist);

