% rsRefineVesicleFits3
% Revised version with acceleration options.
% Given an info structure mi, find and subtract vesicles
% and return the updated mi containing the vesicle coordinates.
% A vesicle model file is stored in Temp/ as <basename>v.mrc.
%
fittingMode=3;
% modes are:  1 - fit overall amplitude
%             2 - fit amplitude and shift only, don't change radius
%             3 - fit amplitude, shift and radius
forceNewModel=0;   % Always ask the user to select a new refined model
% on the first micrograph (can be from the same micrograph)
useOkField=1;      % refine every vesicle for ok(:,1) is true.
doDownsampling=1;  % Downsample for speed
maxVesiclesToFit=inf;
disA=800;          % size of displayed/fitted window in angstroms
writeMiFile=1;     % Save the updated mi file
writeSubFile=0;    % Write models and subtracted images into Temp
setBasePath=1;     % Replace the mi.basePath with the path from the file selector.
vm=struct;
vmGood=0;

% if nargin<1  % put up a file selector
[fname pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
if isnumeric(fname)  % Cancel
    return
end;

[rootPath infoPath]=ParsePath(pa);

if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);

for fileIndex=1:numel(fname)
    %%
    disp(['Reading ' infoPath fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    
    %     check if there is something to do.
    if numel(mi.vesicle.x)>0
        if ~isfield(mi,'log')
            mi.log=cell(0);
        end;
        
        %     Check to see if we have a generic vesicle model (central 5 points are
        %     essentially equal)
        nvm=numel(mi.vesicleModel);
        vmCtr=ceil((nvm+1)/2);
        if nvm < 3 || std(mi.vesicleModel(vmCtr-1:vmCtr+1))<1e-6...
                || forceNewModel
            if fileIndex==1  % Check the first one.
                [mname, pa]=uigetfile('*mi.mat','Select an mi file for membrane model');
                if ischar(mname) % user selected something
                    disp(['Loading the membrane model from ' pa mname]);
                    vm=load([pa mname]);
                    vm=vm.mi;
                    vmGood=1;
                end;
            end;
            if vmGood  % replace the model with the one in vm
                disp(['Using the membrane model from ' mname]);
                if vm.pixA == mi.pixA  % Check if the same pixel size
                    mi.vesicleModel=vm.vesicleModel;  % copy the model
                else                   % resample the model
                    disp('Resampling the model');
                    mi.vesicleModel=meDownsampleVesicleModel(...
                        vm.vesicleModel,mi.pixA/vm.pixA);
                end;
%                 Replace the ampFactors too.
                for i=1:numel(vm.ctf)
                    mi.ctf(i).ampFactor=vm.ctf(i).ampFactor;
                end;
            else
                disp('Using the existing generic membrane model');
            end;
        else  % the previous model was a good one.
            disp('Using the existing refined membrane model');
            mname=fname{fileIndex};
            vmGood=1;
            vm=mi;
        end;
        
        % Read the image and normalize to fractional contrast
        [m, mergePath]=meReadMergedImage(mi);
        
        %     Check that we have a temp directory
        if ~isfield(mi,'tempPath');
            mi.tempPath='Temp/';
        end;
        if ~exist(mi.tempPath,'dir')
            mkdir('Temp');
        end;
        if setBasePath
            mi.basePath=rootPath;
            disp(['Resetting mi.basePath to: ' mi.basePath]);
        end;
        
        %%  ------All the work is done here-------

        mi1=rsRefineVesicleFitsSub(mi,m,fittingMode);
        
        
        %%  Outputting
        mi=mi1;
        if writeMiFile
            %%
            outName=[infoPath mi.baseFilename 'mi.mat'];  % replace original
            if ~isfield(mi,'log')
                mi.log=cell(0);
            end;
            mi.log{numel(mi.log)+1}=['meRefineVesicleFits ' TimeStamp];
            save(outName,'mi');
            disp([outName ' saved']);
            
        end;
        %%             Compute and store model vesicles
        figure(2); clf; SetGrayscale;
        imacs(GaussFilt(m,.1));
        title('Original image');
        drawnow;
        
        disp('Making the final vesicle models');
        vs1=meMakeModelVesicles(mi,size(m),find(mi.vesicle.ok(:,3)));
        
        imacs(GaussFilt(m-vs1,.1));
        title('Subtracted');
        drawnow;
        
        if writeSubFile
            outVesName=[mi.tempPath mi.baseFilename 'v'];
            WriteMRC(vs1,pixA0,[outVesName '.mrc']);
            WriteJpeg(vs1,outVesName);
            imwrite(uint8(imscale(rot90(vs1),256,0)),[outVesName '.jpg']);
            disp([outVesName ' saved']);
            %               Store the subtracted micrograph
            outSubName=[mi.tempPath mi.baseFilename 'mv'];
            WriteMRC(m-vs1,pixA0,[outSubName '.mrc']);
            WriteJpeg(m-vs1,outSubName);
            %             imwrite(uint8(imscale(rot90(m-vs1),256,1e-3)),[outSubName '.jpg']);
            disp([outSubName ' saved']);
        end;
        
    else  % No vesicles have been found to refine
        if numel(mi.vesicle.x)<1
            disp('  ...no vesicles found');
        end;
        if isfield(mi.vesicle,'refined') && mi.vesicle.refined>0
            disp('  ...already refined.');
        end;
    end;
    disp(' ');
end; % for fileIndex
disp('Done.');
disp(' ');
