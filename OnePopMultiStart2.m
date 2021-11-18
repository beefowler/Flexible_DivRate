%Goal: Find best fit of model to 1 day of data, for a single subpopulation
%only. 

%Inputs: 
%day = matlab date number for the day of interest 
%Solar is a two collumn light input, with time and irradiance, in LOCAL time.
% must contain but not necessarily be limited to the day of interst. 
%volbins = designation of size classes of cells in population
%COUNTS = matrix with the number of cells in each size class (rows) as
%specificed by volbins, for each hour (collumns) of day of interest

%Outputs: 
%modelresults is 1 x 21 vector = [day xmin fmin mu mu1 mu2 exitflag length(modelfits)]
%where xmin is optimal theta value. Therefore Theta = modelresults(2:15)
%all values of theta that relate to pop 2 are designated NaN
%fmin is negative log likelihood given these optimal parameters
% mu mu1 mu2 are results of growth_rate function for these paramters
% exitflag is the exitflag for that optimal run, and length(modelfits)
% gives us a sense of how long it took the function to find the answer

%modelfits stores [ theta negloglike mu exitflag ] for each result of fmincon
%allstarts stores initial parameter value for every run
%daySimPROPS is matrix with the expected proportion of cells in each size class
%(rows) as specified by volbins, for each hour (collumns) according to the
%simulation with best fit parameters 

%Einterp is the interpolated light data for the particular day of interest 

function [modelresults, onepopfits, allstarts, simPROPS, simCONC, Einterp] = OnePopMultiStart2(day, Einterp, volbins, COUNTS, hr2, ts)

hr1 = 1; %choose start hour for model fitting
 
%% Fit Model to Data


ms=MultiStart('Display','off','TolX',1e-5,'UseParallel','always','StartPointsToRun','bounds-ineqs');
opts=optimset('Display','off','TolX',1e-8,'Algorithm','interior-point','UseParallel','always','MaxIter', 3000,'MaxFunEvals',10000);
icsTol=0.2;
tolvec=[0.01 0.01 100 0.005 0.5 0.5 10];

%Define parameter bounds 
lb=-[1e-4 1e-4 1e-4 1e-4 5 .5 1e-4]; %negative lower bounds
ub=[1 15 1000 1 50 15 1e4]; %upper bounds 

    a1=-1*eye(7); %set parameter bounds to be interpretted by fmincon
    a2=eye(7);
    A=zeros(14,7);
    A(1:2:14,:)=a1;
    A(2:2:14,:)=a2;

    B=zeros(14,1);
    B(1:2:14)=lb;
    B(2:2:14)=ub;
    
%Create random initial parameter values within bounds 
x0 = zeros(1,7);
num_starts = 40; %choose the number of starting values we want, number of trials of fmincon
x2 = zeros(num_starts, 7); %matrix to create tpoints
for i = 1:7
    x0(i) = -lb(i) + (ub(i) + lb(i))*rand(1);  %remember lb has negative lower bounds. so this is lower bound + rand * difference 
    x2(:,i) = -lb(i) + (ub(i) + lb(i))*rand(num_starts,1);
end
clear i
tpoints = CustomStartPointSet(x2); 


problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta_onepop) negloglike_calc(Einterp, COUNTS,theta_onepop,volbins,hr1,hr2, ts),'Aineq',A,'bineq',B,'options',opts);
[~,~,~,~,soln] = run(ms,problem,tpoints);
    
%%

    %open up the soln structure:
    start_points=zeros(40,7);
    temp=zeros(40,10);
    c=1;
    
    for j=1:length(soln)
        %check to see if all start points led to an individual solution or
        %not (MultiSTart will only return unique solutions)
        g=cell2mat(soln(j).X0);
        if length(g)==7 %only one start_point led to that solution
            start_points(c,:)=g;
            temp(c,1:7)=soln(j).X;
            temp(c,8)=soln(j).Fval;
            temp(c,9)=growth_rate(Einterp,volbins,COUNTS,temp(c,1:7),hr1,hr2, ts);
            temp(c,10)=soln(j).Exitflag;
            c=c+1;
        else
            num=length(g)/7;
            start_points(c:c+num-1,:)=squeeze(reshape(g',1,7,num))';
            temp(c:c+num-1,1:7)=repmat(soln(j).X,num,1);
            temp(c:c+num-1,8)=repmat(soln(j).Fval,num,1);
            temp(c:c+num-1,9)=repmat(growth_rate(Einterp,volbins,COUNTS,temp(c,1:7),hr1,hr2, ts),num,1);
            temp(c:c+num-1,10)=repmat(soln(j).Exitflag,num,1);
            c=c+num;
        end
    end
    %just in case have rows left as zeros
    qq=find(temp(:,1)~=0);
    temp=temp(qq,:);

 
    onepopfits = temp;
    start_points=start_points(qq,:);
    allstarts=start_points;
    
    %let's now ask, in the first batch run, did the solver "converge"?
    [sortlogL, ii]=sort(onepopfits(:,8));
    if abs(sortlogL(min(5,size(sortlogL,1)))-sortlogL(1)) < icsTol  %Did likelihood converge within tolerance
        flag1 = 0; %Yes it did converge
    else
        disp(num2str(sortlogL(1:min(5,size(sortlogL,1)))))
        flag1 = 1; %no it didn't 
    end

       %did parameters converge within a tolerance 
    partol=max(onepopfits(ii(1:min(5,size(sortlogL,1))),1:7))-min(onepopfits(ii(1:min(5,size(sortlogL,1))),1:7));
    if sum(abs(partol) < tolvec)==7 || sum((abs(partol./onepopfits(ii(1),1:7)) < 0.05))==7 %either the modelfits are within an absolute tolerance or within a relative tolerance
        flag2 = 0; %yes they did converge 
    else
        flag2 = 1; %no they didn't converge 
    end

    disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    k=1; %batch number
    while ((flag1 || flag2) || size(onepopfits,1) <= 20) && k <= 5 %run MultiStart 4 more times D: 

        disp(['k: ' num2str(k)])
        k=k+1;

        %Create random initial parameter values within bounds 
        x0 = zeros(1,7);
        num_starts = 40; %choose the number of starting values we want, number of trials of fmincon
        x2 = zeros(num_starts, 7); %matrix to create tpoints
        for i = 1:7
             x0(i) = -lb(i) + (ub(i) + lb(i))*rand(1);  %remember lb has negative lower bounds. so this is lower bound + rand * difference 
             x2(:,i) = -lb(i) + (ub(i) + lb(i))*rand(num_starts,1);
        end
        clear i
        tpoints = CustomStartPointSet(x2); 

        problem = createOptimProblem('fmincon','x0',x0,'objective',@(theta) negloglike_calc(Einterp,COUNTS,theta,volbins,hr1,hr2, ts),'Aineq',A,'bineq',B,'options',opts);
        [~,~,~,~,soln] = run(ms,problem,tpoints);

        %open up the soln structure:
        start_points=zeros(40,7);
        temp=zeros(40,10);
        c=1;

        for j=1:length(soln)
            %check to see if all start points led to an individual solution or
            %not (MultiSTart will only return unique solutions)
            g=cell2mat(soln(j).X0);
            if length(g)==7 %only one start_point led to that solution
                start_points(c,:)=g;
                temp(c,1:7)=soln(j).X;
                temp(c,8)=soln(j).Fval;
                temp(c,9)=growth_rate(Einterp,volbins,COUNTS,temp(c,1:7),hr1,hr2, ts);
                temp(c,10)=soln(j).Exitflag;
                c=c+1;
            else
                num=length(g)/7;
                start_points(c:c+num-1,:)=squeeze(reshape(g',1,7,num))';
                temp(c:c+num-1,1:7)=repmat(soln(j).X,num,1);
                temp(c:c+num-1,8)=repmat(soln(j).Fval,num,1);
                temp(c:c+num-1,9)=repmat(growth_rate(Einterp,volbins,COUNTS,temp(c,1:7),hr1,hr2, ts),num,1);
                temp(c:c+num-1,10)=repmat(soln(j).Exitflag,num,1);
                c=c+num;
            end
        end
        
        %just in case have rows left as zeros
    qq=find(temp(:,1)~=0);
    temp=temp(qq,:);
    start_points=start_points(qq,:);


    onepopfits=[onepopfits; temp];
    allstarts=[allstarts; start_points];


        %okay, now see after this batch run, did the solver "converge"?
        [sortlogL, ii]=sort(onepopfits(:,8));

        if abs(sortlogL(5)-sortlogL(1)) < icsTol
            flag1 = 0;
        else
            disp(num2str(sortlogL(1:5))) %should be 5, but occassionally get less than 5 solver runs returned...
            flag1 = 1;
        end

        partol=max(onepopfits(ii(1:5),1:7))-min(onepopfits(ii(1:5),1:7));
        if sum(abs(partol) < tolvec)==7 || sum((abs(partol./onepopfits(ii(1),1:7)) < 0.05))==7 %either the modelfits are within an absolute tolerance or within a relative tolerance
            flag2 = 0;
        else
            flag2 = 1;
        end
        disp(['flag1 = ' num2str(flag1) ' flag2=' num2str(flag2)])

    end  %while loop

    
    %Assign Outputs
    [~, jj]=sort(onepopfits(:,8));
    xmin=onepopfits(jj(1),1:7);
    fmin=onepopfits(jj(1),8);
    exitflag=onepopfits(jj(1),10);
    
    %calculate growth rate according to ideal theta
    [mu, mu1, mu2] =growth_rate(Einterp,volbins,COUNTS,xmin,hr1,hr2);

    modelresults=[day xmin fmin mu mu1 mu2 exitflag length(onepopfits) ts]; %this is how the optimization results are stored in model output files 

    
% %Step 5: Simulate according to results of optimization 
% 
[simCONC, simPROPS] = simdata(Einterp,COUNTS,xmin,volbins,hr1,hr2, ts);
simPROPS = [NaN*ones(length(volbins), hr1-1) simPROPS]; %fill in null hours before model kicked in 
simCONC = [NaN*ones(length(volbins), hr1-1) simCONC]; %fill in null hours before model kicked in 


end 



















