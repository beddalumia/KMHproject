load('spectral_line.mat','DOS','SE','U_list');
L = length(U_list); 
for i = 46:L
    
    U = U_list(i);
    TitleString = sprintf('U=%s', U)

    %% G(w)

    %Gw = importdata(append(TitleString,'\Greal_l11_s1_realw.dat'));
    Gw = importdata('Greal_l11_s1_realw_indx000001.dat');
    %%% Structure %%%%%%%%%%%%%%%%%%
    % Gw(:,1) gives frequencies w  %
    % Gw(:,2) gives Im(G(w))       %
    % Gw(:,3) gives Re(G(w))       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %f = figure("Name", 'DOS','visible','off');
    w = Gw(:,1);
    Aw = -Gw(:,2);
    FilledStates = w(w<=0);
    FilledDOS = Aw(w<=0);
    EmptyStates = w(w>0);
    EmptyDOS = Aw(w>0);
    area(FilledStates, FilledDOS, 'FaceColor', [0.7 0.7 0.7]); hold on
    area(EmptyStates, EmptyDOS, 'FaceColor', [1 1 1]);
    %xlabel('\omega')
    %ylabel('\piA(\omega)')
    xlabel('Real Frequency')
    ylabel('Density of States')
    %title(TitleString);
    id = 'SOI = 0 | U = 10.50';
    title(id)
    xlim([-6,6]);
    %print(append('DOS\',TitleString,'.png'),'-dpng','-r1200')
    print(append(id,'.png'),'-dpng','-r300')
    %close(f);
    clc

    %f = figure("Name", 'smoothedDOS','visible','off');
    FilledDOS = smoothdata(Aw(w<=0),'gaussian');
    EmptyDOS = smoothdata(Aw(w>0),'gaussian');
    area(FilledStates, FilledDOS, 'FaceColor', [0.7 0.7 0.7]); hold on
    area(EmptyStates, EmptyDOS, 'FaceColor', [1 1 1]);
    xlabel('\omega')
    ylabel('\piA(\omega)')
    %title(TitleString);
    xlim([-6,6]);
    %print(append('smoothedDOS\',TitleString,'.png'),'-dpng','-r1200')
    print(append('DOS.png'),'-dpng','-r1200')
    close(f);

    %% Sigma(w)

    %Ew = importdata(append(TitleString,'\impSigma_l11_s1_realw.ed'));
    Gw = importdata('Greal_l11_s1_realw_indx000001.dat');
    Ew = importdata('impSigma_l11_s1_realw.dat');
    %%% Structure %%%%%%%%%%%%%%%%%%
    % Ew(:,1) gives frequencies w  %
    % Ew(:,2) gives Im(\Sigma(w))  %
    % Ew(:,3) gives Re(\Sigma(w))  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %f = figure("Name",'DOS_and_Self-Energy','visible','off');
    w = Ew(:,1);
    Re = Ew(:,3);
    Im = Ew(:,2);
    w = Gw(:,1);
    Aw = -Gw(:,2);
    %oldmax = max(abs(Im));
    %newmax = max(Aw);
    %Re = Re*newmax/oldmax;
    %Im = Im*newmax/oldmax;
    plot(w,Aw,'LineWidth',1.1); hold on
    plot(w,Re,'LineWidth',1.1, 'Color', [0.9290    0.6940    0.1250]); 
    plot(w,Im,'LineWidth',1.1,'Color', [0.4660    0.6740    0.1880]);
    legend('\piA(\omega)','Re(\Sigma(\omega))', 'Im(\Sigma(\omega))');
    xlabel('\omega')
    %title(TitleString);
    xlim([-6,6]);
    %yticks([]);
    %ylabel('Arb. Units')
    %print(append('All\',TitleString,'.png'),'-dpng','-r1200')
    %close(f);

%     f = figure("Name",'smoothed_DOS_and_Self-Energy','visible','off');
%     ReSmooth = smoothdata(Re,'gaussian');
%     ImSmooth = smoothdata(Im,'gaussian');
%     Aw = smoothdata(Aw,'gaussian');
%     plot(w,Aw,'LineWidth',1.1); hold on
%     plot(w,ReSmooth,'LineWidth',1.1, 'Color', [0.9290    0.6940    0.1250]); 
%     plot(w,ImSmooth,'LineWidth',1.1,'Color', [0.4660    0.6740    0.1880]);
%     legend('\piA(\omega)','Re(\Sigma(\omega))', 'Im(\Sigma(\omega))');
%     xlabel('\omega')
%     title(TitleString);
%     xlim([-6,6]);
%     yticks([]);
%     ylabel('Arb. Units')
%     print(append('smoothedALL\',TitleString,'.png'),'-dpng','-r1200')
%     close(f);

    %% G(iw)
    Giw = importdata(append(TitleString,'\Gloc_l11_s1_iw.dat'));
    %%% Structure %%%%%%%%%%%%%%%%%%%%
    % Giw(:,1) gives frequencies w   %
    % Giw(:,2) gives Im(G(iw))       %
    % Giw(:,3) gives Re(G(iw))       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = figure("Name",'smoothed_DOS_and_Self-Energy','visible','off');
    plot(Giw(:,1),Giw(:,2),'r','LineWidth', 1.1)
    xlabel('i\omega')
    ylabel('Im(G(i\omega))')
    title(TitleString);
    print(append('Giw\',TitleString,'.png'),'-dpng','-r1200')
    close(f);

end
