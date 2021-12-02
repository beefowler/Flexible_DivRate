
%cd('\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\EukWork\Model_1pop\')
% Model now has capacity to fit two subpopulations or jsut fit 1 population.
%for now, let's fit 2 populations as before. 

function ApplyModel_Cruise(inputspath)

filelist = dir([inputspath filesep '*input.mat']);

% set ts to 0 for euks or 6 for syn IF day starts at dawn. 
% Division is limited when 
% hr < ts(2) && hr >= ts(1)
% or if ts is length 1, when 
% hr < ts && hr >= 0;  


%start while loop 
i = 1; 
while i <= length(filelist)
     savename = [filelist(i).name(1:end-9) 'output.mat']; 

     if exist([filelist(i).folder filesep savename])
         i = i+1;
     else
    %load days data 
    eval(['load ' filelist(i).folder filesep filelist(i).name])
    
    if size(N_dist, 2) == 25
        
        %Interpolate light data
        time=0:(1/6):25;
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE, eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(isnan(Einterp)) = 0;
       
        day = floor(datenum(daystarttime)); 
        
        if exist([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat']) %get previous result if available
          load([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat'], 'modelresults')
         informedguess = modelresults(2:15); 
         clear modelresults
            [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayMultiStart(day, Einterp, synvolbins, N_dist, 24, ts, informedguess);     
        else
            [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayMultiStart(day, Einterp, synvolbins, N_dist, 24, ts);     
        end
                  
        %save results
        save([filelist(i).folder filesep savename], 'N_dist', 'Vhists', 'cellsperml', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
               
    else
        disp([filelist(i).name(1:9) 'not enough hours'])
    end %if we have 25 hours
    
   
   clearvars('-except', 'Writerobj1', 'i', 'filelist', 'filepath', 'savepath', 'ts'); 
   
i = i + 1; 
     end

end

%cd('\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\Scripts')

end