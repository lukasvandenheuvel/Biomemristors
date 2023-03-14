%% Stepfinder functions

classdef StepFinder
    
    methods

        function delim = Get_Delim(obj)
            delim = '\';
            if ismac
                delim = '/';
            end
        end

        function [SaveName,initval] = load_initval(obj,source)
            try
                old  = load('config_last_run.mat', 'initval');
                initval = old.initval;
            catch
                errorfile=['No previous run of AutoStepfinder was found, missing config_last_run.mat.'];
                msgbox(errorfile,'ERROR', 'error')
                return
            end
            delim = obj.Get_Delim();
            SplitSource = split(source, delim);
            FileName = SplitSource{end};
            SaveName = FileName(1:length(FileName)-4); % filename without .txt
        end

        function [data,SaveName,initval] = Load_Data_LH(obj,source)

            try
                old  = load('config_last_run.mat', 'initval');
                initval = old.initval;
            catch
                errorfile=['No previous run of AutoStepfinder was found, missing config_last_run.mat.'];
                msgbox(errorfile,'ERROR', 'error')
                return
            end
            delim = obj.Get_Delim();
            SplitSource = split(source, delim);
            FileName = SplitSource{end};
            SaveName = FileName(1:length(FileName)-4); % filename without .txt

            disp('Loading data ...')
            data = double(dlmread(source));

            % Check data size
            [~,cols]=size(data);
            if cols>2
                errorfile=[num2str(FileName),' contains more than two columns.'];
                msgbox(errorfile,'ERROR', 'error');
                return;
            end
            if cols>1
                data=data(:,2); 
            end %if array,it is assumed first column is time
        end

        %%
        function [FinalSteps, FinalFit] = autostepfinder_mainloop(obj,Data,SaveName,initval)
            % This is the main, multi-pass loop of the autostepfinder
            %while initval.nextfile>0 
                %[Data,SaveName,initval]=Get_Data(initval, handles);

             infcheck=isinf(Data);
             nancheck=isnan(Data);
             [LD,cols]=size(Data);
            if sum(infcheck)>0
                display('ERROR: AutoStepfinder detected Inf values, please remove Inf values from data and try again')
                return    
            end
            if sum(nancheck)>0
                display('ERROR: AutoStepfinder detected NaN values, please remove Inf values from data and try again')
                return    
            end
            if cols>2
                display('ERROR: AutoStepfinder detected more than two columns, please reformat data and try again')
                return    
            end
            LD=length(Data);
            IndexAxis=(1:LD)';     
            stepnumber_firstrun=min([ceil(LD/4) initval.fitrange]);        
            Residu=Data;  Fit=0*Data;
            S_Curves=zeros(stepnumber_firstrun+1,2); 
            N_found_steps_per_round=zeros(2,1); 
            full_split_log=[];

            %% core dual pass
            for fitround=1:2
                initval.stepnumber=stepnumber_firstrun;                         
                [FitResidu,~,S_Curve,split_indices,best_shot]=obj.StepfinderCore(Residu,initval);       
                steproundaccept=(max(S_Curve)>initval.SMaxTreshold);
                if steproundaccept
                   N_found_steps_per_round(fitround)=best_shot;         
                   full_split_log=obj.expand_split_log(full_split_log,split_indices,fitround,best_shot);           
                 end   
                S_Curves(:,fitround)=S_Curve;
                Residu=Residu-FitResidu;  %new residu  
                Fit=Fit+FitResidu ;       %new fit
            end    

            %% Final analysis: 
            % init outputs for function call
            FinalFit = zeros(size(Data)) + mean(Data);
            FinalSteps = NaN;
            if max(N_found_steps_per_round)==0, cla;
                 %nostepmessage=['No steps found in: ', SaveName];
                 %nostepbox=msgbox(nostepmessage,'Operation complete', 'help'); 
               %if initval.hand_load == 2
               %    pause(1);
               %    delete(nostepbox);
               %    disp(char(nostepmessage));
               %end
               nostepmessage=['No steps found in: ', SaveName];
               disp(char(nostepmessage))
            else
                [FinalSteps, FinalFit]=obj.BuildFinalfit(IndexAxis,Data,full_split_log,initval);                 
            % SaveAndPlot(   initval,SaveName,handles,...
            %                        IndexAxis, Data, FinalFit,...
            %                        S_Curves, FinalSteps,N_found_steps_per_round);     

            end        
            disp('done!');
            % end
        end

        %% This section contains the 'Core' function of the stepfinder; 
        %it can be cut and autorun independently (on a simple simulated curve) for demo purposes

        function [FitX,stepsX,S_fin,splitlog,best_shot]=StepfinderCore(obj,X,initval)
            %This function splits data in a quick fashion.
            %This one is a compact version of the first quick 2007 version
            %output: list of stepsizes: [index time  levelbefore levelafter step dwelltimeafter steperror]
            if nargin<2
                X=2*round(0.5+0.1*sin(2*pi*(1:20000)'/500))+rand(20000,1); 
                initval.stepnumber=300; initval.overshoot=1;   
            end
                %% 1 split, estimate best fit, repeat
                initval.stepnumber=min([ceil(length(X)/4) initval.stepnumber]);
                [~,~,S_raw,splitlog]=obj.Split_until_ready(X,initval); %run 1: full iteration

                [s_peakidx,S_fin]=obj.Eval_Scurve(S_raw);
                best_shot=round(min([(s_peakidx-1) ceil(length(X)/4)])); 
                indexlist=sort(splitlog(1:best_shot)); 

                FitX=obj.Get_FitFromStepsindices(X,indexlist,initval);   
                stepsX=obj.Get_Steps(FitX); 

            if nargin<2
                close all;
                subplot(2,1,1); plot(X,'r'),hold;plot(FitX,'k','LineWidth',2);
                title('Data and Fit');xlabel('time');ylabel('position,a.u.');
                subplot(2,1,2); semilogx(S_fin,'-o');
                title('S-curve');xlabel('Stepnumber');ylabel('S-value, a.u.');
            end
        end

        function [bestshot,S_fin]=Eval_Scurve(obj,S_raw)
            S_raw(S_raw<1)=1; S2=S_raw-1;  %remove base line
            BaseLine=linspace(0,S2(end),length(S2)); S3=S2-BaseLine';
            %S4=smooth(S3,ceil(ix/25));
            [~,i1]=max(S3);  %peak
            %sel=find(S3>0.9*pk1);i2=max(sel(sel>=i1)); bestshot=i2;
            bestshot=i1;
            S_fin=S3;
        end


        function stepsX=Get_Steps(obj,FitX)
        %list of stepsizes: [index time step levelbefore levelafter dwelltimeafter]    
            lx=length(FitX);
            T=(1:lx)';
            difX=FitX(2:lx)-FitX(1:lx-1);
            sel=find(difX~=0);    
            lsel=length(sel);
            if lsel>0
                dwellX=T(sel(2:lsel))-T(sel(1:lsel-1)); dwellX=[dwellX' T(lx)-T(sel(lsel))]';
                stepsX=[sel T(sel) FitX(sel) FitX(sel+1) difX(sel) dwellX];
            else
                stepsX=[NaN NaN NaN NaN NaN NaN];
            end
        end

        function [FitX,f,S,splitlog]=Split_until_ready(obj,X,initval)
             c=1; stop=0;
             N=length(X);    
             FitX=mean(X)*ones(N,1); 
             S=ones(initval.stepnumber,1);
             splitlog=zeros(initval.stepnumber,1);
             %Create the first plateau------------------------------------------
             istart=1; istop=length(X);
             [inxt, avl, avr,rankit]=obj.Splitfast(X(istart:istop));           
             f=[[1, 1, 1, 0, 0,0];
                [istart, istop, inxt+istart-1, avl, avr,rankit]; ...
                [N, N, N,0, 0,0,]]; 

            %build first counterfit:
              cFitX=0*FitX; i1=1; i2=f(2,3); i3=N;
              cFitX(i1:i2)=avl; cFitX(i2+1:i3)=avr;


             %parameters needed for calculating S(1):-----------------
            qx=sum(X.^2);                                   %sum of squared data
            qm=N*(mean(X))^2;                               %sum of squared averages plateaus, startvalue
            aqm=(inxt-istart+1)*avl^2+(istop-inxt)*avr^2;   %sum of squared averages anti-plateaus, startvalue
            S(c)=(qx-aqm)/(qx-qm);                          %S: ratio of variances of fit and anti-fit        
            %---------------------------------       
            minimumwindowsize=3;  %minimum plateau length to split
            wm=minimumwindowsize;
             while stop==0 %Split until ready 
                c=c+1;
                fsel=find((f(:,2)-f(:,1)>wm)&f(:,6)~=0);        %among those plateaus sensibly long..
                [~,idx2]=max(f(fsel,6)); idx=(fsel(idx2));   %...find the best candidate to split. 
                splitlog(c-1)=f(idx,3);                          %keep track of index order
                FitX=obj.Adapt_Fit(f,idx,FitX);                     %adapt fit-curve
                [f,qm,aqm]=obj.expand_f(f,qm,aqm,idx,X);            %adapt plateau-table; adapt S

                cFitX=obj.Adapt_cFit(f,idx,cFitX,X);                     
                %adapt fit-curve (note we use the updated split-table as we need
                %the new fields)

                if 1
                    S(c)=mean((X-cFitX).^2)/mean((X-FitX).^2);         %direct fit   
                else
                    S(c)=(qx-aqm)/(qx-qm);                             %Calculate new S-function 
                end
                stop=(1.0*c>initval.stepnumber);
            end   %-------------------------------------------------------------------
        end

        function [f,qm,aqm]=expand_f(obj,f,qm,aqm,idx,X)
        %this function inserts two new plateau-property rows on the location of one old one
        %....and adapts the S-function nominator and denominator

           %1) Label new levels.FLR locatess 'plateau fit  left right' etc  
            nFLR=f(idx-1,2)-f(idx-1,3);       avFLR=f(idx-1,5);      %FLR
            nFML=f(idx,3)-f(idx,1)+1;         avFML=f(idx,4);        %FML
            nFMR=f(idx,2)-f(idx,3);           avFMR=f(idx,5);        %FMR
            nFRL=f(idx+1,3)-f(idx+1,1)+1;     avFRL=f(idx+1,4);      %FRL

            %remove contribution from old plateau(s) from S-function terms
            qm=qm-1/(nFML+nFMR)*(nFML*avFML+nFMR*avFMR)^2;            %FM
            aqm=aqm-    1/(nFLR+nFML)*(nFLR*avFLR+nFML*avFML)^2-...   %CL
                        1/(nFMR+nFRL)*(nFMR*avFMR+nFRL*avFRL)^2;      %CR

            %2a) construct new first plateau entry, left
            istart=f(idx,1); istop=f(idx,3);
            [inxt, avl, avr,rankit]=obj.Splitfast(X(istart:istop));
            n1=[istart istop inxt+istart-1, avl, avr,rankit];

            %2b) construct new first plateau entry, right
            istart=f(idx,3)+1; istop=f(idx,2);
            [inxt, avl, avr,rankit]=obj.Splitfast(X(istart:istop));
            n2=[istart istop inxt+istart-1, avl, avr,rankit];

            %3) Insert these two new plateaus in place of the old one
            [lf,~]=size(f);      
            block1=f(1:idx,:);      block1(idx,:)=n1;
            block2=f(idx:lf,:); 	block2(1,:)=n2;
            f=[block1; block2];

            %Label newly defined levels.FMLR locates 'plateaufit mid/left/right' etc
            nFMLL=f(idx,3)-f(idx,1)+1;          avFMLL=f(idx,4);    %FMLL
            nFMLR=f(idx,2)-f(idx,3);            avFMLR=f(idx,5);    %FMLR
            nFMRL=f(idx+1,3)-f(idx+1,1)+1;      avFMRL=f(idx+1,4);  %FMRL
            nFMRR=f(idx+1,2)-f(idx+1,3);        avFMRR=f(idx+1,5);  %FMRR

            %4) add contribution from new plateau(s) to S-function terms
            qm=qm...
                +1/(nFMLL+nFMLR)*(nFMLL*avFMLL+nFMLR*avFMLR)^2 ...  %FML
                +1/(nFMRL+nFMRR)*(nFMRL*avFMRL+nFMRR*avFMRR)^2;     %FMR
            aqm=aqm ...
                +1/(nFLR+nFMLL)*(nFLR*avFLR+nFMLL*avFMLL)^2 ...      %CTL
                +1/(nFMLR+nFMRL)*(nFMLR*avFMLR+nFMRL*avFMRL)^2 ...   %CTM
                +1/(nFMRR+nFRL)*(nFMRR*avFMRR+nFRL*avFRL)^2;         %CTR
        end

        function FitX=Adapt_Fit(obj,f,idx,FitX)
            %This function creates step  and property curves adds new plateaus
            i1=f(idx,1); i2=f(idx,3);av1=f(idx,4);
            i3=f(idx,3)+1; i4=f(idx,2);av2=f(idx,5);
            FitX(i1:i2)=av1; FitX(i3:i4)=av2;
        end

        function cFitX=Adapt_cFit(obj,f,idx,cFitX,X)
            %This function adapts the counterfit locally
            i1=f(idx-1,3); 
            i2=f(idx,3);
            i3=f(idx+1,3); 
            i4=f(idx+2,3);
            cFitX(i1+1:i2)=mean(X(i1+1:i2)); 
            cFitX(i2+1:i3)=mean(X(i2+1:i3));
            cFitX(i3+1:i4)=mean(X(i3+1:i4));    
        end

         function [idx, avl, avr,rankit,errorcurve]=Splitfast(obj,Segment)              %
        %this function also adresses a one-dim array 'Segment'
        %and determines the best step-fit there
        %To save time, functions like 'mean' are avoided
            w=length(Segment);   
            Chisq=(1:w-1)*0;  
            if w>3
                Chisq=(1:w-1)*0;                           
                AvL=Segment(1);    AvR=sum(Segment(2:w))/(w-1); AvAll=sum(Segment)/w;  
                for t=2:w-2
                    AvL=(AvL*(t-1)+Segment(t))/t;     AvR=(AvR*(w-t+1)-Segment(t))/(w-t);
                    DStepL=AvL-AvAll;           DStepR=AvR-AvAll;
                    DeltaChisq=((DStepL.^2)*t+(DStepR.^2)*(w-t));            
                    Chisq(t)=-DeltaChisq;       
                end
                 [~,idx]=min(Chisq(2:w-2)); idx=idx+1;
                 avl=mean(Segment(1:idx));                    avr=mean(Segment(idx+1:w));
                 %rankit=(avr-avl)^2/(1/(idx-1)+1/(w-idx-1));  %quantity expresing predicted accuracy step (squared)
                 rankit=(avr-avl)^2*w;
            else                                            %low-limit cases
                 rankit=0;               
                 switch w
                     case 3 
                        a=(Segment(2)+Segment(3))/2-Segment(1);            b=Segment(3)-(Segment(1)+Segment(2))/2;
                        cL=[Segment(1) , (Segment(1)+Segment(2))/2];      cR=[(Segment(2)+Segment(3))/2 , Segment(3)]; 
                        [~,idx]=max([a b]);    avl=cL(idx);  avr=cR(idx);
                    case 2
                        idx=1;  avl=Segment(1); avr=Segment(2); 
                end
            end
            errorcurve=Chisq/(w-1);
         end

          function FitX=Get_FitFromStepsindices(obj,X,indexlist,initval) 

              % This function builds plateau data
            %list of levels: [startindex  stopindex starttime stoptime level dwell stepbefore stepafter]
                if initval.fitmean&&~initval.fitmedian, modus='mean';end
                if ~initval.fitmean&&initval.fitmedian, modus='median';end
            lx=length(X);
                lsel=length(indexlist); %note: index points to last point before step

                %Build a 'FitX' based on the median levels (ot on the averages)
                idxes=[0 ; indexlist ; lx];
                FitX=0*X;
                for ii=1:lsel+1
                    ixlo=idxes(ii)+1;  %first index of plateau
                    ixhi=idxes(ii+1);  %last index
                    switch modus
                        case 'mean', FitX(ixlo:ixhi)=nanmean(X(ixlo:ixhi));
                        case 'median', FitX(ixlo:ixhi)=nanmedian(X(ixlo:ixhi));
                    end
                end 
          end

        %% This section contains code related to the multipass steps

        function full_split_log=expand_split_log(obj,full_split_log,split_indices,fitround,best_shot)
           %expand split log; remove double entries (when a residu of a
           %step is found again, we label the location by its first
           %round occurence)
           LS=length(split_indices);
           new_split_log=[split_indices 3*ones(LS,1)];
           new_split_log(1:best_shot,2)=fitround;          
               if fitround==1
                    full_split_log=[full_split_log ; new_split_log];
               else
                    % Changed by Lukas; this section raised an error when
                    % no steps are found in the first fitround.
                    if isempty(full_split_log)                              % Added by Lukas
                        full_split_log=[full_split_log ; new_split_log];    % Added by Lukas
                    else                                                    % Added by Lukas
                        already_found=find(full_split_log(:,2)<3);  %steps already spotted
                        new_found=find(~ismember(new_split_log(:,1),full_split_log(already_found,1)));                  
                        full_split_log= [full_split_log(already_found,:);...
                                         new_split_log(new_found,:)];
                    end
               end
        end

        function StepsX=AddStep_Errors(obj,X,StepsX,initval)   
        %This function calculates step errors associated with the steps.
        % Two options: 
        %     1) 'measured' ;the standard deviation of the adjacent plateaus
        %     2) 'predicted'; using the global noise level and the length of these plateaus
        steperrorestimate='measured';                           

        if strcmp(steperrorestimate,'predicted')
            shft=2;
            globalnoise=nanstd((X(shft:end)-X(1:end-shft+1)))/sqrt(2);
        end

        [ls,col]=size(StepsX); i1=0;
        for i=1:ls
            i2=StepsX(i);
            if i<ls
                i3=StepsX(i+1);
            else 
                i3=length(X);
            end    
            Nbefore=i2-i1;    
            Nafter=i3-i2;

            if strcmp(steperrorestimate,'measured')
                rmsbefore=std(X(i1+1:i2));
                rmsafter=std(X(i2+1:i3)) ;
                StepsX(i,col+1)=2*((rmsbefore^2/Nbefore+rmsafter^2/Nafter)^0.5)/2^0.5; %plus minus 95%
            end
            if strcmp(steperrorestimate,'predicted')
                StepsX(i,col+1)=2*(globalnoise^2/Nbefore+globalnoise^2/Nafter)^0.5/2^0.5; %plus minus 95%
            end
            i1=i2;
        end
        end


        function [FinalSteps, FinalFit]=BuildFinalfit(obj,T,X,splitlog,initval)
        %build a step fit based on all retained indices. Perform step-by-step error
        %analysis to accept second-round (residual) steps or not 
            best_shot=length(find(splitlog(:,2)<3)); %all non-duplicates

            if (initval.manualoff==1 && initval.manualon==0)
                steps_to_pick=round(initval.overshoot*best_shot);
            end
            if (initval.manualoff==0 && initval.manualon==1)
                steps_to_pick=initval.setsteps;
            end
            %select indices to use
            bestlist=(splitlog(1:steps_to_pick,:));  
            [candidate_loc,ix]=sort(bestlist(:,1));
            candidateround_no=bestlist(ix,2);

            %2)Rebuild fit from all indices.
            candidate_fit=obj.Get_FitFromStepsindices(X,candidate_loc,initval); 
            candidate_steps=obj.Get_StepTableFromFit(T,candidate_fit); 
            %3 get errors
            candidate_steps=obj.AddStep_Errors(X,candidate_steps,initval); 
            candidate_relsteperror=(candidate_steps(:,8)./abs(candidate_steps(:,5))); 

            %4 Reject weird steps
            %Local Stepmerge: if larger than 0, weird steps are removed from fit. the
            %treshold of removing is based on the avarage erros of the steps found in
            %round 1.      
            localstepmerge=1;
            if (localstepmerge && ~initval.setsteps)
                % Default: Keep round 1-steps AND 'good' round 2 steps:
                % Get a   measure for the error of steps in the first round. 
                % This can be used for reference of errors from second-round steps
                [~,~,FinalErrorTreshold]=obj.Outlier_flag(candidate_relsteperror,2,0.8,'positive',0);
                sel_merge=find((candidate_relsteperror<2*FinalErrorTreshold)|candidateround_no==1);        
                final_idxes=candidate_steps(sel_merge,1); 
            else
                final_idxes=candidate_steps(:,1); 
            end
            %Re-build the fit from the selected indices (Effectively, rejected
            %indices are 'merged' in this step)
            FinalFit=obj.Get_FitFromStepsindices(X,final_idxes,initval);
            FinalSteps=obj.Get_StepTableFromFit(T,FinalFit); 
            FinalSteps=obj.AddStep_Errors(X,FinalSteps,initval);
            LF=length(FinalSteps(:,1));
            FinalRoundNo=zeros(LF,1);
            for ii=1:LF
                idxC=FinalSteps(ii,1);
                sel=find(candidate_steps(:,1)==idxC);
                FinalRoundNo(ii)=candidateround_no(sel(1));
            end        
            FinalSteps(:,9)=FinalRoundNo;  

            if initval.bootstraprepeats>0
                [err_st, err_t, ~,~ ]=obj.bootstrap_get_errors(X, final_idxes,initval.bootstraprepeats);
                FinalSteps(:,10)=err_st;
                FinalSteps(:,11)=err_t;
            end
        end      

        function [data,SaveName,initval]=Get_Data(obj,initval, handles)
        % This function loads the data, either standard or user-choice
            disp('Loading..');
            CurrentFolder=pwd;
            switch initval.hand_load     
                case 0     
                filecheck = exist('config_last_run.mat');
                if filecheck > 0
                old = load('config_last_run.mat', 'initval');
                initval.source=old.initval.source;
                source   = initval.source;
                FileName = old.initval.FileName;
                PathName = old.initval.PathName;
                initval.FileName=FileName;
                initval.PathName=PathName;
                else 
                    errorfile=['No previous run of AutoStepfinder was found, missing config_last_run.mat.'];
                    msgbox(errorfile,'ERROR', 'error')
                     return
                end
                try
                data=double(dlmread(source));
                save('config_last_run.mat','initval');
                catch
                    errorfile=[num2str(FileName),' is not formatted properly.'];
                    msgbox(errorfile,'ERROR', 'error')
                     return
                end
                SaveName=FileName(1:length(FileName)-4);
                initval.nextfile=0;   
                case 1      
                cd(initval.datapath);
                [FileName,PathName] = uigetfile('*.*','Select the signal file');
                cd(CurrentFolder);
                initval.FileName=FileName;
                initval.PathName=PathName;
                source=strcat(PathName,FileName);
                initval.source=source;
                try
                data=double(dlmread(source));
                save('config_last_run.mat','initval');
                set(handles.rerun,'enable','On');
                catch
                    errorfile=[num2str(FileName),' is not formatted properly.'];
                    msgbox(errorfile,'ERROR', 'error')
                     return
                end
                SaveName=FileName(1:length(FileName)-4);
                initval.nextfile=0;      
                case 2
                    fileindex=initval.nextfile;
                    cd(initval.datapath);
                    AllFileNames=dir('*.txt');
                    cd(CurrentFolder);
                    chckflder=length(AllFileNames);
                    if chckflder == 0,  msgbox('The provided input folder is empty.','ERROR', 'error')
                    return; end;
                    AllFiles=length(AllFileNames);
                    FileName=AllFileNames(fileindex).name;
                    try
                    if ismac
                    data=double(dlmread(strcat(initval.datapath,'/',FileName)));
                    else
                    data=double(dlmread(strcat(initval.datapath,'\',FileName)));
                    end
                    save('config_last_run.mat','initval');
                    set(handles.rerun,'enable','Off');
                    catch
                    errorfile=[num2str(FileName),' is not formatted properly.'];
                    msgbox(errorfile,'ERROR', 'error')
                     return
                   end
                    SaveName=FileName(1:length(FileName)-4);
                    if initval.nextfile==AllFiles
                        initval.nextfile=0;
                    else
                        initval.nextfile=fileindex+1;
                    end           
            end
            [~,cols]=size(data);
            if cols>2, errorfile=[num2str(FileName),' contains more than two columns.'];
                        msgbox(errorfile,'ERROR', 'error');
                        return;
                        end
            if cols>1, data=data(:,2); end %if array,it is assumed first column is time


        %    check for NaN + Inf values
             infcheck=isinf(data);
             nancheck=isnan(data);

             if sum(infcheck)>0
                 errorfile=[num2str(FileName),' contains infinite values.'];
                 msgbox(errorfile,'ERROR', 'error')
                 return;
             end

             if sum(nancheck)>0
                 errorfile=[num2str(FileName),' contains NaN values.'];
                 msgbox(errorfile,'ERROR', 'error')
                 return;
             end 

             disp('Analyzing:'), disp(SaveName);
        end

        function [StepsX,levelX, histX]=Get_StepTableFromFit(obj,T,FitX)
        %This function builds tables of steps or levels properties from a step fit
        %Values are based on Averages of plateaus
            lx=length(FitX);
            difX=FitX(2:lx)-FitX(1:lx-1);
            sel=find(difX~=0);  %note: index points to last point before step
            lsel=length(sel);

            %dwell time after
            dwellT=T(sel(2:lsel))-T(sel(1:lsel-1)); 
            dwellTafter=[dwellT' T(lx)-T(sel(lsel))]';
            dwellTbefore=[T(sel(1))-T(1) dwellT']'; 
            StepsX=[sel T(sel) FitX(sel) FitX(sel+1) difX(sel) dwellTbefore dwellTafter]; 
            %list of stepsizes: [index time levelbefore levelafter step dwelltime before dwelltimeafter]

            Levstartidx=[1; sel+1];  %first index of level
            Levstopidx=[sel; lx];    %last index of level
            LevstartT=T(Levstartidx);
            LevstopT=T(Levstopidx);
            LevLevel=[FitX(sel); FitX(lx)];
            LevDwell=Levstopidx-Levstartidx+1;
            LevStepBefore=[0; difX(sel)];
            LevStepAfter=[difX(sel); 0];
            levelX=[Levstartidx Levstopidx LevstartT LevstopT LevLevel LevDwell LevStepBefore LevStepAfter];
            %list of levels: [startindex starttime stopindex stoptime level dwell stepbefore stepafter]   

            %histogram
            stepsizes=StepsX(:,4)-StepsX(:,3);
            mx=max(stepsizes); mn=min(stepsizes);
            stp=(mx-mn)/20;
            hx=(mn:stp:mx);
            if ~isempty(hx)
                histX=hist(stepsizes,hx)';
                histX=[hx' histX];
            else
                histX=[1 1];
            end
        end

         function [flag,cleandata,treshold]=Outlier_flag(obj,data,tolerance,sigchange,how,sho)
        %this function is meant to find a representative value for a standard
        %deviation in a heavily skewed distribution (typically, flat data with
        % %peaks). It calculates the standard deviation and average the data;
        % Based on these, outliers are determined and excluded for a new calculation
        % of average and SD; this is repeated until sigma does not change anymore too much
        % . This is repeated until the new sigma does not change much
        % %anymore
        %output: positions of outliers
        %Jacob Kers 2013 and before---------------------------------------------
        binz=50;

        sigma=1E20;            %at start, use a total-upper-limit 
        ratio=0;
        ld=length(data);
        flag=ones(ld,1);  %at start, all points are selected
        while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
            sigma_old=sigma;
            selc=find(flag==1);
            data(flag==1); 
            av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
            sigma=nanstd(data(selc));
            ratio=sigma/sigma_old;
            treshold=tolerance*sigma+av;
            switch how
                case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
                case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
            end
            %plot menu------------------  
            if sho==1
                cleandata=data(selc); 
                hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
                sthst=hist(cleandata,hx);
                figure;
                bar(hx,sthst);
                title('Histogram');
                pause(0.5);  
                close(gcf);
            end
            %---------------------------- 
        end
        cleandata=data(selc); 
         end

        function [error_st_boot, error_t_boot, t_refit,error_t_refit]=bootstrap_get_errors(obj,data, indices,bootstraprepeats)
        % bootstrap 
        %Aim: bootstrap to get location error (plus a time estimate)

        %use: [error_t_boot,error_t_refit,error_st_boot]=get_errors_by_bootstrap(data, indices)

        %Input: 
        %data: original data, sigle column, 
        %indices: locations of stepfits, in pts
        %

        %Output:
        % error_t_boot: 95% confidence range of error in time (pts) by bootstrapping
        % t_refit: location obtained by re-fitting
        % error_t_refit: difference of input step location and location obtained by re-fitting
        % error_st_boot:  95% confidence range of error in sep size by bootstrapping

        %References: 
        %following method:
        %[1] Received 18 Jul 2015 | Accepted 16 Nov 2015 | Published 17 Dec 2015 
        % ATP hydrolysis assists phosphate release and promotes reaction ordering in F1-ATPase 
        % Chun-Biu Li1, Hiroshi Ueno2, Rikiya Watanabe2,3,4, Hiroyuki Noji2,4 & Tamiki Komatsuzaki1


        %% bootstrap main cycle
        %loop all segments
        N_steps=length(indices-2);
        error_t_refit=zeros(N_steps,1);
        error_t_boot=zeros(N_steps,1);
        error_st_boot=zeros(N_steps,1);

        t_refit=zeros(N_steps,1);

        indices_ext=[0; indices; length(data)];

        for ii=1:length(indices)
            %% 1) our usual step fit to get step location (plus extra export of error curve)
            Segment=data(indices_ext(ii)+1:indices_ext(ii+2));
            idx_old=indices(ii);  %input location
            [idx, ~, ~,~, error_curve]=obj.Splitfast(Segment) ;

            %note that we also obtain a new estimate for the best fit in this
            %segment:
            t_refit(ii)=indices_ext(ii)+idx;
            error_t_refit(ii,1)=abs(t_refit(ii)-idx_old);

            %% 2 bootstrapping
            % repeat plateau fits left and right by bootstrapping 
            %(location is kept constant)
            tic
             %repeat many times:
            leftpart=Segment(1:idx);        L_left=length(leftpart);
            rightpart=Segment(idx+1:end);   L_right=length(rightpart);

            if ((L_left>20)&(L_right>20))    
                [bootstat_av_left,bootsam_left] = bootstrp(bootstraprepeats,@mean,leftpart);  %resample left, get many left averages
                [bootstat_av_right,bootsam_right] = bootstrp(bootstraprepeats,@mean,rightpart); %resample right, get many right averages
                stepsize_boot=bootstat_av_right-bootstat_av_left;

                %% 3 Get the value Chi-square(idx) for all resamplings.
                %To gain time, work matrix-wise :
                %expand the plateau results in a block, such that every column is a new step-fit
                %build with the resampled left and right averages
                newstep_indices=[bootsam_left; idx+bootsam_right];  %matrix of indices
                NewSegments=Segment(newstep_indices);
                NewFits=[repmat(bootstat_av_left,1,L_left) repmat(bootstat_av_right,1,L_right)]';

                all_chi_squares=(mean((NewSegments-NewFits).^2)).^0.5;
                bootstrap_error_of_minimum_value=1.96*std(all_chi_squares);



                %% 4 get estimate of 95% confidence range of value of Errcurv at minimum
                %to increase precision, we interpolate 10-fold near the minimum
                LL=length(error_curve);
                nearmin_lox=max([1 idx-50]);
                nearmin_hix=min([LL idx+50]);
                axz=nearmin_lox:nearmin_hix;
                axz_ip=nearmin_lox:0.1:nearmin_hix;
                Errcurv_blowup=interp1(axz,error_curve(axz),axz_ip);

                sel=find(abs(Errcurv_blowup-error_curve(idx))<=bootstrap_error_of_minimum_value);
                lox_ip=min(sel);  
                hix_ip=max(sel);

                lox=axz_ip(lox_ip); 
                hix=axz_ip(hix_ip);


                error_t_boot(ii,1)=(hix-lox)/2;
                error_st_boot(ii,1)=1.96*std(stepsize_boot)/(2^0.5);
            else %too short plateaus, do not bootstrap
                error_t_boot(ii,1)=NaN;
                error_st_boot(ii,1)=NaN;
            end
        end
        end
    end
end