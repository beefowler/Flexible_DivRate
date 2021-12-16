
%cd('\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\EukWork\Model_1pop\')
% Model now has capacity to fit two subpopulations or jsut fit 1 population.
%for now, let's fit 2 populations as before. 

function ApplyModel_SpeedVersion(inputspath)

filelist = dir([inputspath filesep '*input.mat']);

if ~exist([inputspath filesep 'speedver'], 'dir')
    mkdir([inputspath filesep 'speedver'])
end

% set ts to 0 for euks or 6 for syn IF day starts at dawn. 
% Division is limited when 
% hr < ts(2) && hr >= ts(1)
% or if ts is length 1, when 
% hr < ts && hr >= 0;


%start while loop
i = 1;
while i <= length(filelist)
    savename = [filelist(i).name(1:end-9) 'output.mat'];
    
    if exist([filelist(i).folder filesep 'speedver' filesep savename]) || exist([filelist(i).folder filesep savename])
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
            
            if i > 1 && (exist([filelist(i).folder filesep 'speedver' filesep filelist(i-1).name(1:end-9) 'output.mat']) || exist([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat'])) %get previous result if available
                if exist([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat'])
                     load([filelist(i).folder filesep filelist(i-1).name(1:end-9) 'output.mat'], 'modelresults')
                else
                     load([filelist(i).folder filesep 'speedver' filesep filelist(i-1).name(1:end-9) 'output.mat'], 'modelresults')
                end
                
                informedguess = modelresults(2:15)
                informedguess(9) = min(informedguess(9), 1-informedguess(9));
                
                clear modelresults
                
                
                opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 3000,'MaxFunEvals',10000);
                %Define parameter bounds
                lb=-[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 5 5 .5 .5 1e-4]; %negative lower bounds
                ub=[1 15 1000 1 1 15 1000 1 0.5 50 50 15 15 1e4];
                
                a1=-1*eye(14); %set parameter bounds to be interpretted by fmincon
                a2=eye(14);
                A=zeros(28,14);
                A(1:2:27,:)=a1;
                A(2:2:28,:)=a2;
                
                B=zeros(27,1);
                B(1:2:27)=lb;
                B(2:2:28)=ub;
                
                COUNTS = N_dist;
                hr2 = 24;
                hr1 = 1;
                
                
                problem = createOptimProblem('fmincon','x0',informedguess,'objective',@(theta) negloglike_calc(Einterp,COUNTS,theta,synvolbins,hr1,hr2, ts),'Aineq',A,'bineq',B, 'options',opts);
                
                % [x, fval, exitflag, soln] = fmincon(problem);
                %we don't really need multistart, but for some reason it's working better with bounds
                ms=MultiStart('Display','off','TolX',1e-5,'UseParallel', false ,'StartPointsToRun','bounds-ineqs');
                
                tpoints = CustomStartPointSet(informedguess);
                [~,~,~,~,soln] = run(ms,problem,tpoints);
                
                trial = 1;
                while size(soln, 2) == 0 && trial< 4
                    [~,~,~,~,soln] = run(ms,problem,tpoints);
                    trial = trial+1
                end
                
                if trial == 4 && size(soln, 2) == 0 %need to do it the long way
    
                    [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayFastStart2(day, Einterp, synvolbins, N_dist, 24, ts);
                    %save to long way outputs 
                    save([filelist(i).folder filesep savename], 'N_dist', 'Vhists', 'cellsperml', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
                    disp(['saved: ' savename])
                
                else %process short way ansers
                    
                    x = soln.X;
                    exitflag = soln.Exitflag;
                    fval = soln.Fval;
                    
                    temp = zeros(1,17);
                    temp(1:14) = x;
                    temp(15) = fval;
                    [mu, mu1, mu2] = growth_rate(Einterp,synvolbins,COUNTS,temp(1:13),hr1,hr2, ts);
                    
                    if mu < 0.05 %Flag if division rate is very low 
                        [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayFastStart2(day, Einterp, synvolbins, N_dist, 24, ts);
                         %save to slow outputs
                        save([filelist(i).folder filesep savename], 'N_dist', 'Vhists', 'cellsperml', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
                        disp(['saved: ' savename])
                    else 
                        temp(16) = mu;
                        temp(17) = exitflag;
                    
                        largepopn=zeros(1,7); %large population has the larger mean starting volume bin
                        smallpopn=zeros(1,7);
                        if temp(10) > temp(11)
                            largepopn = temp([1:4 9 10 12]);
                            smallpopn= [temp(5:8) 1-temp(9) temp([11 13])];
                        else %11 > 10
                            largepopn = [temp(5:8) 1-temp(9) temp([11 13])];
                            smallpopn = temp([1:4 9 10 12]);
                        end
                    
                        %modelfits stores the proportion of population with smaller mean for theta(9),
                        modelfits=[smallpopn(1:4) largepopn(1:4) smallpopn(5) smallpopn(6) largepopn(6) smallpopn(7) largepopn(7) temp(14:end)];
                        allstarts=[informedguess];
                    
                        modelresults = [day x fval mu mu1 mu2 exitflag length(modelfits) ts];
                    
                        [simCONC, simPROPS] = simdata(Einterp,COUNTS,x,synvolbins,hr1,hr2, ts);
                        simPROPS = [NaN*ones(length(synvolbins), hr1-1) simPROPS]; %fill in null hours before model kicked in
                        simCONC = [NaN*ones(length(synvolbins), hr1-1) simCONC]; %fill in null hours before model kicked in
                    
                    
                   
                        %save to speed outputs
                        save([filelist(i).folder filesep 'speedver' filesep savename], 'N_dist', 'Vhists', 'cellsperml', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
                        disp(['saved: ' savename])
                    end %if mu was low or high
                end %end if short way didn't give an answer
                
            else %if don't have previous day so have to do long way 
                [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayFastStart2(day, Einterp, synvolbins, N_dist, 24, ts);
                %save to long way outputs 
                save([filelist(i).folder filesep savename], 'N_dist', 'Vhists', 'cellsperml', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'ts')
                disp(['saved: ' savename])
            
            end %if don't have previous day   
            
        else %didn't have enough hours
            disp([filelist(i).name(1:9) 'not enough hours'])
        end %if we have 25 hours
        
        
        clearvars('-except', 'Writerobj1', 'i', 'filelist', 'filepath', 'savepath', 'ts');
        
        i = i + 1;
    end
    
end

%cd('\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\Scripts')

end
