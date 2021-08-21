%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOW TO USE?
%
%  - Put this script, together with an input-file for your driver in a path
%    > This path will contain directories for all the U values you set-up.
%  - Set-up the name of your driver program (without .f90 extension)
%    > e.g. driver = 'ed_kane_mele'; 
%  - Set SOI to your desire: --> you will get a fixed-SOI linear span
%  - Adjust Umin and Umax to your desire --> U \in [Umin, Umax]
%  - Adjust Uold to catch a 'restart-folder' in the path [!applies -> -1]
%  - Select doMPI (true.or.false) to run with openMPI or not
%  - Run everything with $ matlab -batch RunningDMFT_autostep
%  - At the end you will find some additional output in the U=%f folders
%    > a LOG_dmft.txt which is just a mirror of the DMFT output (via tee)
%    > a LOG_time.txt which is a wall-clock-time value for the whole DMFT
%  - Also a additional output files in the main (external) path
%    > a U_list.txt that stores all the used U-values (for easier post..)
%    > possibly some error-flag files in the format 'ERROR_U=%f'
%      if you see them you *may* have to discard the corresponding folder
%      > check it for convergence! (look at LOG_dmft.txt)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver = 'ed_kane_mele';	doMPI = true;

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
Umin = 0; Umax = 10;           % Input Hubbard 
Wmix = 0.3;		 	% Input Self-Mixing

Ustep = [0.5, 0.25, 0.1, 0.05, 0.01]; % To be auto-determined...hopefully!
NUstep = length(Ustep);

notConvFlag = false;		% Convergence-fail *flag*
notConvCount = 0;		% Convergence-fail *counter*
notConvThreshold = NUstep-1;   % Maximum #{times} we accept DMFT to fail

U = Umin; Uold = -1;
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);      % ------------------------------------
if isfolder(oldDIR)                  % If it exist a "previous" folder: 
restartpack = [oldDIR,'/*.restart']; % Copy all the restart files from the
copyfile(restartpack);               % last dmft evaluation...
end                                  % ------------------------------------

copyfile ../input*             % Copy inside the **external** input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
mpi = 'mpirun ';				% Control of MPI
else						% boolean flag...
mpi = [];
end
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE
T2 =sprintf(' t2=%f',SOI);			% of
MIXING = sprintf(' wmixing=%f',Wmix);		% PARAMETERS
outLOG = ' > LOG_dmft.txt';
dmft_ed_call = [mpi,driver,HUBBARD,T2,MIXING,outLOG];
tic
system(dmft_ed_call);				% Fortran-call
chrono = toc;
file_id = fopen('LOG_time.txt','w');
fprintf(file_id,'%f\n', chrono);		% Write on time-log
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README')
    notConvFlag = true;
    notConvCount = notConvCount + 1;
    movefile('ERROR.README',sprintf('../ERROR_U=%f',U));
else
    fprintf(Ulist,'%f\n', U);	               % Write on U-log
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder


if notConvCount > notConvThreshold
   error('DMFT not converged: phase-span stops now!');         
end 

if notConvFlag == true
   U = Uold; 			% if nonconverged we don't want to update 
   notConvFlag = false;	% > but we want to reset the flag(!)
else
   Uold = U; 			% if converged we update Uold and proceed to...         
end

U = U + Ustep(notConvCount+1); % ...Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(Ulist);
