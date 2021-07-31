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
%  - Run everything with $ matlab -batch DryMF_line
%  - At the end you will find some additional output in the U=%f folders
%    > a LOG_mf.txt which is just a mirror of the MF output (via tee)
%    > a LOG_time.txt which is a wall-clock-time value for the whole MF
%  - Also an additional output file in the main (external) path
%    > a U_list.txt that stores all the used U-values (for easier post..)
%      NB. Only the converged calculations will update the U_list :D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver = 'mf_km_2d';	doMPI = true;

% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ulist = fopen('U_list.txt','a');

%% SOC-Line: outer loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spin_orbit = [0,0.02,0.04,0.06,0.08,0.1,0.2,0.3];
Nsoi = length(spin_orbit);
for iSOI = 1:Nsoi		% Spin-Orbit loop =========================>

SOI  = spin_orbit(iSOI);       % Input Spin-Orbit 

soiDIR= sprintf('SOI=%f',SOI); % Make a folder named 'SOI=...', where '...'
mkdir(soiDIR);                 % is the given value for SpinOrb interaction
cd(soiDIR);                    % Enter the SOI-folder

copyfile ../input*             % Copy inside the **external** input file

%% Phase-Line: inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Umin = 0; Umax = 8;            % Input Hubbard 
Ustep = 0.1;			% Phase-line step
Uold = -1;			% Restart option

U = Umin; 
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);      % ------------------------------------
if isfolder(oldDIR)                  % If it exist a "previous" folder: 
restartpack = [oldDIR,'/*.restart']; % Copy all the restart files from the
copyfile(restartpack);               % last mf evaluation...
end                                  % ------------------------------------

copyfile ../input*             % Copy inside the **external** input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
mpi = 'mpirun ';			        % Control of MPI
else					        % boolean flag...
mpi = [];
end
HUBBARD =sprintf(' uloc=%f',U);	        % OVERRIDE of
T2 =sprintf(' t2=%f',SOI);		        % PARAMETERS
outLOG = ' > LOG.mf';
mf_call = [mpi,driver,HUBBARD,T2,outLOG];
tic
system(mf_call);				% Fortran-call
chrono = toc;
file_id = fopen('LOG_time.txt','w');		% Write on time-log
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A FAILED (unconverged) MF LOOP
if ~isfile('ERROR.README') 
    fprintf(Ulist,'%f\n', U);	          % Update U-list (only if converged)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder

Uold = U;
U = U + Ustep;              	% Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the SOI-folder

end                            % <=========================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(Ulist);

