% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOI  = 0;                      % Input Spin-Orbit 
Umax = 10; Ustep = 0.25;       % Input Hubbard 
Wmix = 0.3;		 	% Input Self-Mixing

notConverged = 0;		% Convergence-fail *counter*
notConvThreshold = 4;          % Maximum #{times} we accept DMFT to fail

U = 0; Uold = -1;
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

if Uold >= 0                   % If it exist a "previous" folder:
oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../inputKANEMELE.conf % Copy inside the *external* input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
mpi = 'mpirun ';
driver = 'ed_kane_mele';
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE
T2 =sprintf(' t2=%f',SOI);			% of
MIXING = sprintf(' wmixing=%f',Wmix);		% PARAMETERS
outLOG = ' | tee LOG_dmft.txt';
dmft_ed_call = [mpi,driver,HUBBARD,T2,MIXING,outLOG];
tic
system(dmft_ed_call);
chrono = toc;
file_id = fopen('LOG_time.txt','w');
fprintf(file_id,'%f\n', chrono);
fclose(file_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE NEED TO CATCH SOMEHOW A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README') 
    notConverged = notConverged + 1;
    delete ERROR.README
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder


if notConverged > notConvThreshold
   error('DMFT not converged: phase-span stops now!');
elseif notConverged > 0 && notConverged < notConvThreshold
   % You either reduce the U grain or the selfmixing parameter.
   % Or maybe both (?!) ---> TO FIGURE OUT YET
   Ustep = Ustep/5;           
   % Wmix = Wmix/10;   
   U = Uold; % if nonconverged we don't want to update!          
else
   Uold = U; % if converged we update Uold and proceed to...         
end

U = U + Ustep;                 % ...Hubbard update  

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
