HERE = pwd;

Rvals = 1:1:4
Nvals = [24,54,96,150,216,294,384];

commonUmax = 100;

for i = 1:length(Rvals)

   R = Rvals(i);
   RDIR = sprintf("R=%d",R)

   cd(RDIR)

   U_list = QcmP.post.get_list('U');

   commonUmax = min(commonUmax,max(U_list));

   Etot = 0;

   kins = QcmP.post.kinetic_line('U',U_list);

   for j = 1:Nvals(i)

      site = sprintf('ineq%.4d',j)
      
      [~,ens] = QcmP.post.energy_line('U',site,U_list);
      pots = ens{1};

      Etot = Etot + pots + kins; % kins is intensive

   end

   writematrix(Etot,sprintf('Etot_R%d.txt',R));

   cd('..')

end

writematrix(0:0.1:commonUmax,'commonUgrid.txt')

