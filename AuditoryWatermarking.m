function gui_2()

clear all
close all
clc

global data signal Fs L y signal_watermarked ber
%////////////////////////////////////////////////////

figure_1 = figure('MenuBar','none','Units','centimeters','Name','watermarking','NumberTitle','off','Position',[2,2,27,18]);
axes_img_1 = axes('parent',figure_1,'units','centimeters','position',[8 10 7 7]);
axes_img_2 = axes('parent',figure_1,'units','centimeters','position',[17 10 7 7]);
axes_img_3 = axes('parent',figure_1,'units','centimeters','position',[8 1.5 7 7]);
axes_img_4 = axes('parent',figure_1,'units','centimeters','position',[17 1.5 7 7]);

panel_1 = uipanel('parent',figure_1,'BorderType','etchedin','Title','functions','Units','centimeters','Position',[0.5 7 5 10]);
sel_button_audio =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','Select Audio File','Position',[0.2,8.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@sel_audio_rot);
sel_button_watemark =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','Select Watermark','Position',[0.2,7.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@sel_watermark);

start_embed =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','Embed Watermark','Position',[0.2,6.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@embed_watermark);
prameter_calc1 =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','SNR','Position',[0.2,5.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@snr_rot);
prameter_calc2 =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','BER','Position',[0.2,4.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@ber_rot);

play_button1 =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','play orignal','Position',[0.2,2.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@play_orignal);
play_button2 =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','play watermarked','Position',[0.2,1.5,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@play_watermarked);

exit_button =uicontrol('parent',panel_1,'Style','Pushbutton','Units','centimeters','String','EXIT','Position',[0.2,0.1,4,0.6],...
 'FontSize',11,'HorizontalAlignment','Center','CallBack',@exit_rot);

progress_text = uicontrol('parent',figure_1,'Style','Text','String','Status','Units','centimeters','Position',[0.5,6,7.5,0.5],'HorizontalAlignment','Left','FontSize',11);
ber_text = uicontrol('parent',figure_1,'Style','Text','String','BER:','Units','centimeters','Position',[0.5,5,7.5,0.5],'HorizontalAlignment','Left','FontSize',11);
snr_text = uicontrol('parent',figure_1,'Style','Text','String','SNR:','Units','centimeters','Position',[0.5,4,7.5,0.5],'HorizontalAlignment','Left','FontSize',11);
%/////////// callback functions /////////////////////

function sel_audio_rot(varargin)
    [name,path]=uigetfile('*.wav');
    P=fullfile(path,name);
    [signal,Fs] = audioread(P);
    
    T = 1/Fs;                    
    L = length(signal);           
    t = (0:L-1)*T;                
    y = signal;  
    
    
    plot(axes_img_1,Fs*t,y)
    title(axes_img_1,'orignal signal')
    drawnow
    set(progress_text,'String','Status');
end

function sel_watermark(varargin)
    [name,path]=uigetfile('*.txt');
    P1=fullfile(path,name);
    
    embedding_data=1000;
    disp('reading data .....')
    fid=fopen(P1); 
    combine=[];
    while(1)
        tline = fgetl(fid);   
            if ~ischar(tline)
                break
            end  
            combine=[combine unicode2native(tline)];        
    end
    bytes=round(embedding_data/8);
    data=combine(1:bytes);
end

function embed_watermark(varargin)
    NFFT = 2^nextpow2(L); 
    Y = fft(y,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);

    X=f;
    Y=2*abs(Y(1:NFFT/2+1));
    
    plot(axes_img_2,X,Y) 
    title(axes_img_2,'Single-Sided Amplitude Spectrum of y(t)')

    mean_freq=mean(Y);
    drawnow
    
    %////////// applying GA to get optimized frequencies //////////

    com=[X' Y];
    chromosomes=100;
    generations=100;
    [r1,c1]=find(com(:,2)<mean_freq);
    possible_freq=com(r1,1);

    disp(strcat('detected possible frequices to watermark data:[',num2str(length(possible_freq)),']'));
    disp(' ')
    disp(strcat('applying GA to predict best:[',num2str(length(data)),']'));

    for i=1:chromosomes
        P_chromo(i,:)=possible_freq(randi([1 length(possible_freq)],length(data),1));
    end
    
    %////////// fit function ////////////////////////////////////////

        for j=1:chromosomes
                curr_chromo=P_chromo(j,:);
                signal_temp=signal;
                for i=1:length(data)
                    curr_data=(double(data(i))/255)/100;
                    pos=i:ceil(1000000/curr_chromo(i)):length(signal);
                    signal_temp(pos)=signal_temp(pos) + curr_data;
                end
                Err=sum(abs(signal - signal_temp));
                P_fit(j,1)=Err;
        end
        pbest=min(P_fit);
        gbest=min(P_fit);
        max_val=gbest;

        pbest=1-(pbest/max_val);
        gbest=gbest/max_val;

        S=(strcat('Generation:[',num2str(0),'] Gbest:[',num2str(gbest),']'));
        disp(S)
        set(progress_text,'String',S);
        %////////////////////////////////////////////////////////////////

        for k=1:generations
            for j=1:chromosomes
                curr_chromo=P_chromo(j,:);
                x=randi([1 length(curr_chromo)],10,1);
                y=randi([1 length(possible_freq)],10,1);

                mut_chromo=curr_chromo;
                for jj=1:length(x)
                    mut_chromo(x(jj))=possible_freq(y(jj));
                end
                C_chromo(j,:)=mut_chromo;
            end

            for j=1:chromosomes
                curr_chromo=C_chromo(j,:);
                signal_temp=signal;
                for i=1:length(data)
                    curr_data=(double(data(i))/255)/100;
                    pos=i:ceil(1000000/curr_chromo(i)):length(signal);
                    signal_temp(pos)=signal_temp(pos) + curr_data;
                end
                Err=sum(abs(signal - signal_temp));
                C_fit(j,1)=Err/max_val;
                p_fit=min(C_fit);
            end

            for j=1:chromosomes
                if(P_fit(j)>C_fit(j))
                    P_chromo(j,:)=C_chromo(j,:);
                    P_fit(j)=C_fit(j);
                end
            end
            gbest=min(P_fit);
            S=(strcat('Generation:[',num2str(k),'] Gbest:[',num2str(gbest),']'));
            disp(S)
            set(progress_text,'String',S);
            
            %/////////// implementing best solution to signal ////////////

            [r1,c1]=find(P_fit==min(P_fit));
            best_chromo=P_chromo(r1(1),:);

            signal_watermarked=signal;
            for i=1:length(data)
                curr_data=(double(data(i))/255)/100;
                pos=i:ceil(1000000/best_chromo(i)):length(signal);
                signal_watermarked(pos)=signal_watermarked(pos) + curr_data;
                decoded_data(1,pos)=curr_data;
            end
            decoded_data=decoded_data * 255 * 100;
            decoded_data=decoded_data(1:length(data));
            ber=(sum(abs((double(data) - double(decoded_data))))) * 10e+8;
            
            S=(strcat('Generation:[',num2str(k),'] Gbest:[',num2str(gbest),']'));
            disp(S)
            set(progress_text,'String',S);
            
            plot(axes_img_3,signal,'-b')
            title(axes_img_3,'orignal signal')
            
            plot(axes_img_4,signal_watermarked,'-r')
            title(axes_img_4,'watermarked signal')
            drawnow
            
            clear decoded_data
        end
end

function play_orignal(varargin)
        sound(signal,Fs);
end

function play_watermarked(varargin)
        sound(signal_watermarked,Fs);
end

function snr_rot(varargin)
    SNR = 100*(abs(20*log10(std(signal)/std(signal_watermarked)))*100);
    S=(strcat('SNR:[',num2str(SNR),'db]'));
    disp(S)
    set(snr_text,'String',S);
end

function ber_rot(varargin)
    S=(strcat('BER:[',num2str(ber),']'));
    disp(S)
    set(ber_text,'String',S);
end

end