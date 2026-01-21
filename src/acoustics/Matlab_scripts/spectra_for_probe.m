%**********************************************************************
% Example of Matlab scrip for signal post-processing
% - read probe time history from LES solver
% - compute and plot the corresponding spectra
%
% G. A. Brès - Cascade Technologies
%**********************************************************************

clear all;
close all;

% Provide path and name of files
LEScase = strvcat('10M-300000');
jet_case = ['/Users/Guillaume/Desktop/Cascade/Projects/SBIR/SubsonicJet/2-Simulations/'];
surface = strvcat('/probe');

% loadmatlab=1 if spectra already computed and saved; 0 otherwise 
loadmatlab=0;

% choose variables
var=strvcat('U-x', 'U-y', 'U-z');
nvar=3;
var_start=1;
var_end=nvar;

xprobe=[0.05 0.5 1 2 5 10];
nprobes = 6;

color=strvcat('r','b','g','m','c');
color2=strvcat('r','b','g','--k','--m','--c');
color3=strvcat('--m','--b','--m','--b','--c');

% Mach number for velocity scaling
Ma=0.9;

DSt=0.05;
bandwidth = DSt*Ma;


for l=1:size(LEScase,1)
    
    LES_name = strtrim(LEScase(l,:))
    
    for k=1:size(surface,1)
        
        surface_name = strtrim(surface(k,:))
        
        for ivar=var_start:var_end;
            
            for i=1:nprobes;
                
                if (loadmatlab==0)
                    
                    % load time history
                    filename=[jet_case,LES_name,surface_name,'/Probe_x',num2str(xprobe(i)),'.',strtrim(var(ivar,:))];
                    [h, aux] = hdrload(filename);
                    t=aux(:,2);
                    nmic = size(aux,2)-3
                    for j=1:nmic;
                        p(j,:) = aux(:,j+3);
                    end
                    Dt = t(2)-t(1);
                    Df = 1/Dt;
                    
                    NFFT_total=length(t);
                    Df_total=1/(t(end)-t(1));
                    
                    FFToverlap =0.5;
                    
                    
                    for j = 1:nmic;
                        %=============================================
                        % LES probes
                        %=============================================
                        %%% compute narroband spectra
                        [Phat(j,:), P2(j,:), f]=Spectra(p(j,:),NFFT_total,Df);
                        
                        %%% overlap spectra
                        [PSD_o(j,:), foverlap] = SpectraOverlap(p(j,:),Df,bandwidth,FFToverlap);
                        
                        %%% DSt bin-averaged spectra
                        [df, fbin, fc, fc_l, fc_u, PSD_b(j,:)] = BinAveraging(P2(j,:), f, f(1), f(end), bandwidth, 0, 0);
                        
                    end
                    
                    %%% azimuthal averaging
                    PSD(i,:) = mean(P2,1)/Df_total;
                    PSDoverlap(i,:) = mean(PSD_o,1);
                    PSDbin(i,:) = mean(PSD_b,1);
                    
                    %%% save spectra to matlab file
                    matname = [jet_case,LES_name,surface_name,'/Probe-',LES_name,'_x',num2str(xprobe(i)),'_',strtrim(var(ivar,:)),'.mat']
                    save(matname,'Ma','DSt','bandwidth','fbin','PSDbin','f','PSD','Df_total','PSDoverlap','foverlap','Df');
                    
                else
                    %%% load existing matlab file
                    matname = [jet_case,LES_name,surface_name,'/Probe-',LES_name,'_x',num2str(xprobe(i)),'_',strtrim(var(ivar,:)),'.mat'];
                    load(matname)
                end
              
                figure(i)
                hold on;
                %%% narrowband
                %plot(f/Ma,PSD(i,:)/(Ma*Ma),'k','LineWidth',2);
                %%% bin-average
                plot(fbin/Ma,PSDbin(i,:)/(Ma*Ma),strtrim(color2((ivar-var_start+1)+(l-1)*(var_end-var_start+1),:)),'LineWidth',2,'MarkerSize', 10);
                
                clear aux
                clear p PSD PSDoverlap PSDbin P2 Phat f foverlap fbin PSD_b
            end
            
        end
    end
end


for i=1:nprobes;
    
    %=============================================
    % LES probes
    %=============================================
    savefile=['Probe-',LES_name,'_x',num2str(xprobe(i)),'.png'];
    figure(i)
    hold on;
    ylim([0.00001 0.1]);
    xlim([0.05 10]);
    legend('10M Ux','10M Uy','10M Uz','Location','SouthWest'); 
    set(gca,'XScale','Log');
    set(gca,'YScale','Log');        
    ylabel('PSD');  
    xlabel('St=fD/Uj');
    box on;
    set(gca,'layer','top')
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,savefile,'png')
end





