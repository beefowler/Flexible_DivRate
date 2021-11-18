filelist = dir('\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2010\euk_model\dawnstart_inputs_2019\*data.mat');
filepath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2010\euk_model\dawnstart_inputs_2019\';
savepath = '\\sosiknas1\Backup\Overflow_Outputs_BLF\MVCO_Jan2010_Redo2019\'; 
ts = 0 % Decide if this is 6 for syn or 0 for euk 


%start while loop 
i = 68; 
while i <= length(filelist)
     
    %load days data 
    eval(['load ' filepath filelist(i).name])
    savename = [filelist(i).name(1:9) 'output.mat']; 
    
    
    if size(N_dist, 2) == 25
        
        %Interpolate light data
        time=0:(1/6):25;
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE, eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(isnan(Einterp)) = 0;
        

        %apply model
        [modelresults, modelfits, allstarts, simPROPS, simCONC] = OnePopMultiStart(day, Einterp, volbins, N_dist, 24, ts);
        
        
        %save results
        save([savepath savename], 'CONC', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp')
               
    else
        disp([filelist(i).name(1:9) 'not enough hours'])
    end %if we have 25 hours
    
   
   clearvars('-except', 'Writerobj1', 'i', 'filelist', 'filepath', 'savepath'); 
   
i = i + 1; 

end

