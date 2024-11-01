%%
clc; close all; clearvars -except dataset_name     

cd(userpath); cd ..; cd ..;

homepath = pwd;
scripts_path = fullfile(homepath,'Box\CMR_Group\Users\Fikunwa Kolawole\Scripts');
addpath(genpath(scripts_path));

%% folders

dataset_name                = input("data set name? eg 'HEALTHY_006'\n");


filePath    = matlab.desktop.editor.getActiveFilename;
parts       = regexp(filePath, '[\\/]', 'split');
currentFolder = strjoin(parts(1:end-1), '/');
cd(currentFolder); cd('..');
dataAnalysis_master_folder  = pwd;
dataAnalysis_folder         = fullfile(dataAnalysis_master_folder,dataset_name);


trackedPoints_folder        = fullfile(dataAnalysis_folder, '02_grid/06_trackedPoints/');
remesh_folder               = fullfile(dataAnalysis_folder, '01_bSSFP/04_LV_model/02_remesh'); 

febioModel_folder           = fullfile(dataAnalysis_folder, '04_PassiveMechanics/');

grid_LAX_folder             = fullfile(dataAnalysis_folder, '02_grid/03_reg2cardiac/02_sliceAligned_data/v2/03_lax_align/03_data_aligned_nrrd_sorted/LAX');
if ~exist(grid_LAX_folder, 'dir')
    grid_LAX_folder = fullfile(dataAnalysis_folder, '02_grid/03_reg2cardiac/01_raw_data_sorted_nrrd/LAX');
end

grid_SAX_folder             = fullfile(dataAnalysis_folder, '02_grid/03_reg2cardiac/02_sliceAligned_data/v2/02_sax_align_w_bssfp/03_data_aligned_nrrd_sorted/SAX/');
if ~exist(grid_SAX_folder, 'dir')
    grid_SAX_folder = fullfile(dataAnalysis_folder, '02_grid/03_reg2cardiac/01_raw_data_sorted_nrrd/SAX');
end

optimization_folder         = fullfile(dataAnalysis_folder, '04_PassiveMechanics/01_FEM/02_optimization'); 

if ~exist(optimization_folder, 'dir')    
    optimization_folder = fullfile(optimization_folder, 'v1');
    mkdir(optimization_folder)
end
%% Plot settings

fontSize    = 30;
markerSize  = 10;
position    = [200 200 2000 1500];
lineWidth   = 2; 

%% user input DS and ES phases

switch (dataset_name)
    case {'HEALTHY_006'}
        ES_phase_number = 10;
        DS_phase_number = 21;
        % slices to skip
        noTrackSlices.base     = 1; 
        noTrackSlices.apex     = 1;
    otherwise
        assert(false,'Unknown dataset. Need to identify end systole and diastasis phase number')
end

%% write a note about the optimization

cd(optimization_folder)
diary notes.txt
fprintf([num2str(noTrackSlices.apex) ' less apical and ' num2str(noTrackSlices.base) ' less basal slice. Anisotropy ratios from Nasopoulou \n']);
diary off

%% copy the model file from the febio Model folder to the optimization folder

copyfile(fullfile(febioModel_folder, [dataset_name, '_opt.feb']), fullfile(optimization_folder, [dataset_name '_opt.feb']))

%% EXPERIMENTAL DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LV mesh with labeled faces                             %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[lvMyo_mesh.facets,lvMyo_mesh.nodes,lvMyo_mesh.facet_labels,~,lvMyo_mesh.elements,~,~,lvMyo_mesh.phys_names,~] = gmsh_read_mesh(fullfile(remesh_folder, '02_tet/DS/lvMyo_DS_3p00_v2.msh'));

% define febio model file location
febioFebFileNamePart    = [dataset_name '_opt.feb'];
febioFebFileName        = fullfile(optimization_folder,febioFebFileNamePart); %FEB file name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% import tracked points                                  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tracked points
trackedPoints_folder_info = dir(fullfile(trackedPoints_folder, '02_reformatted/*.mat'));
load(fullfile(trackedPoints_folder_info.folder, trackedPoints_folder_info.name)) % trackedPoints

% visualize tracked points relative to surface model
figure;

subplot(1,2,1)
H = patch('vertices', lvMyo_mesh.nodes, 'faces', lvMyo_mesh.facets,...
        'FaceColor',       [.8 .8 .8], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15,...
        'FaceAlpha', .2, ...
        'EdgeAlpha', .2);
hold on
pcshow(trackedPoints.DS.all_QPt(:,:,DS_phase_number), 'w')
title("DS", 'color', 'white')
axis on
xlim([-40 70])
ylim([-40 40])
ylim([-40 40])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% find phase with max displacement from diastasis phase %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for phase_idx = DS_phase_number:size(trackedPoints.DS.all_QPt,3) 
    displacement(:,:,phase_idx)     = trackedPoints.DS.all_QPt(:,:,phase_idx) - trackedPoints.DS.all_QPt(:,:,DS_phase_number);
    displacement_mag(:,phase_idx)   = vecnorm(displacement(:,:,phase_idx),2,2);
end

% find max displacment phase
[~,I]               = max(displacement_mag');
max_disp_phase_no   = mode(I);

% visualize tracked points at max disp
subplot(1,2,2)
pcshow(trackedPoints.DS.all_QPt(:,:,max_disp_phase_no), 'r')
title('ED', 'color', 'red')
axis on
xlim([-40 70])
ylim([-40 40])
ylim([-40 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get the optimization criteria for only the nodal points closest to the %%%% 
%%%% tracked slices but bounded by the most apical and basal slices %%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % load grid nrrd images
% grid_SAX_image_info = dir(fullfile(grid_SAX_folder, ['*phase' num2str(DS_phase_number, '%03.0f') '.nrrd']));
% grid_SAX_image      = nhdr_nrrd_read(fullfile(grid_SAX_image_info.folder, grid_SAX_image_info.name), true);
% 
% long_vector     = grid_SAX_image.spacedirections_matrix(:,3); % vector pointing in longitudinal dir
% no_slices       = size(grid_SAX_image.data, 3)-2;

long_vector     = [-9.9444 0.026546 -0.55673]'; % for HEALTHY_006
no_slices       = 9; % for HEALTHY_006
meshNodes       = febioExtractNodalPositions(febioFebFileName);

% create interpolant of positions at end inflation phase
interpolant_pos_maxInf_x = scatteredInterpolant(trackedPoints.DS.all_QPt(:,:,DS_phase_number),...
                                                trackedPoints.DS.all_QPt(:,1,max_disp_phase_no), 'natural');
interpolant_pos_maxInf_y = scatteredInterpolant(trackedPoints.DS.all_QPt(:,:,DS_phase_number),...
                                                trackedPoints.DS.all_QPt(:,2,max_disp_phase_no), 'natural');
interpolant_pos_maxInf_z = scatteredInterpolant(trackedPoints.DS.all_QPt(:,:,DS_phase_number),...
                                                trackedPoints.DS.all_QPt(:,3,max_disp_phase_no), 'natural');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interpolate to find the mesh nodal positions at max Inflation %%%%%%%%%
%%% (only at the mesh nodes close to the tracked slices) %%%%%%%%%%%%%%%%%%
%%% (use all slices except most basal and most apical slice) %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes_exp.nodePos_referential = []; % mesh nodes that would be used in optimization
clear nodes_exp.nodeID

figure; % visualize nodal "slices" 
pcshow(meshNodes, 'w'); 
axis on
xlim([-40 70])
ylim([-40 40])
ylim([-40 40])
for slice_idx = 1+noTrackSlices.apex : no_slices-noTrackSlices.base 
    % point_on_plane = grid_SAX_image.spaceorigin + grid_SAX_image.spacedirections_matrix(:,3)*slice_idx;

    point_on_plane = [65.244 203.85 55.652]' + [-9.9444 0.026546 -0.55673]'*slice_idx;
    points_dist_from_plane  = abs(dot(long_vector'.*ones(length(meshNodes),3),meshNodes-point_on_plane',2)); % point distance from plane is abs(dot(N,point-point_on_plane))
    point_dist_min_idx      = points_dist_from_plane < 1;
    slice_meshNodes         = meshNodes(point_dist_min_idx,:);

    hold on; pcshow(slice_meshNodes, 'r')

    nodes_exp.nodePos_referential = [nodes_exp.nodePos_referential; slice_meshNodes];

end

% get the febio corresponding nodeID of the meshNodes for opt criteria
[~,~,row_idx]       = intersect(nodes_exp.nodePos_referential, meshNodes, 'rows','stable');
nodes_exp.nodeID    = row_idx;

nodes_exp.nodePos_maxInf(:,1) = interpolant_pos_maxInf_x(nodes_exp.nodePos_referential);
nodes_exp.nodePos_maxInf(:,2) = interpolant_pos_maxInf_y(nodes_exp.nodePos_referential);
nodes_exp.nodePos_maxInf(:,3) = interpolant_pos_maxInf_z(nodes_exp.nodePos_referential);

% export experimental nodes
save(fullfile(optimization_folder, 'nodes_exp'), 'nodes_exp');

%% Defining file names

% export the check points node set (ie points used for optimization criteria)
febioFebFileNamePart    = febioExportNodeSetData('position',nodes_exp.nodeID,optimization_folder, febioFebFileNamePart); % write into febio model file, the nodeset

% define full file names
febioFebFileName        = fullfile(optimization_folder,febioFebFileNamePart); %FEB file name
febioLogFileName        = fullfile(optimization_folder,[febioFebFileNamePart(1:end-4),'.log']); %FEBio log file name
febioNodePosFileName    = fullfile(optimization_folder,'nodePos.log'); %Log file name for exporting LV endocardium positions

%% Running one forward simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% parameter initial guess %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set initial guess to params from Nasopoulou improved... fix alpha, rf and rt
c           = 1.7e-3; % in MPa
alpha       = 15;
rf          = .55;
rt          = .25;

param_values_transIso   = [c alpha rf rt];
param_names_transIso    = {'c', 'alpha', 'rf', 'rt'};

[param_values_ortho, param_names_ortho] = transIso_2_ortho(param_values_transIso,param_names_transIso,1);


febioReplaceMatValue(febioFebFileName,param_names_ortho,param_values_ortho); % replace material properties in file

febioAnalysis.run_filename  = febioFebFileName; %The input file name
febioAnalysis.run_logname   = febioLogFileName; %The name for the log file
febioAnalysis.disp_on       = 1; %Display information on the command window
febioAnalysis.runMode       = 'external';%'internal';%'external'
febioAnalysis.maxLogCheckTime = 120;

[runFlag]   = runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Import FEBio results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if runFlag==1 %i.e. a succesful run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Import simulated nodal positions from log file %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [nodes_sim.nodePos_maxInf,nodes_sim.nodeID] = febioExtractNodalPositions(febioNodePosFileName);

    % Visualize forward simulated nodal positions curve
    figure; hold on;
    title('Nodal positions optimisation','FontSize',fontSize);
    
    Hn(1)=plot3(nodes_sim.nodePos_maxInf(:,1),...
                nodes_sim.nodePos_maxInf(:,2),...
                nodes_sim.nodePos_maxInf(:,3),'.r', MarkerSize=markerSize);    
    
    Hn(2)=plot3(nodes_exp.nodePos_maxInf(:,1),...
                nodes_exp.nodePos_maxInf(:,2),...
                nodes_exp.nodePos_maxInf(:,3),'.b', MarkerSize=markerSize);    
    ax = gca;
    ax.View = [5 120];
    axis off; 
    axis tight
    axis equal
    legend('Simulation','Experiment','Location','northwest','TextColor','black');
    set(gca,'FontSize',fontSize);
    set(gcf, 'Position',position);
    drawnow;

    % fix and set x,y and z limits
    xlim(xlim) 
    ylim(ylim)
    zlim(zlim)
    
    % clear images folder so can start fresh everytimes
    cd(optimization_folder)
    cmd_rmdir('images\png');   
    diary off
    if exist('opt.log', 'file') ; delete('opt.log'); end
end

%% Run series of optimizations fixing alpha to find optimal c

opt_param_combo.c       = [];
opt_param_combo.alpha   = [];

alpha_range = logspace(log10(10),log10(200), 3);

% alpha_range = param_values_transIso(2);

for alpha_idx = 1:length(alpha_range) 

    
    alpha   = alpha_range(alpha_idx);

    param_values_transIso = [param_values_transIso(1) alpha param_values_transIso(3:4)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Create structures for optimization  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Material structure %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mat_struct.id           = 1; %Material id
    
    % reformulated params. Can reconstruct actual parameters using these
    mat_struct.param_names_transIso          = param_names_transIso; %Parameter names
    mat_struct.param_values_transIso         = num2cell(param_values_transIso);%Parameter values for optimization
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% What should be known to the objective function: %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    objectiveStruct.h               = Hn(1); % to update the plot of simulated nodal positions
    if isequal(nodes_sim.nodeID, nodes_exp.nodeID )
        objectiveStruct.nodes_exp   = nodes_exp; % 
    else
        fprintf(2,'\n Warning: simulated and experimental nodeIDs do not match \n')
    end
    febioAnalysis.disp_on       = 0; 
    febioAnalysis.disp_log_on   = 0;

    objectiveStruct.febioAnalysis       = febioAnalysis;
    objectiveStruct.febioFebFileName    = febioFebFileName;
    objectiveStruct.febioNodePosFileName= febioNodePosFileName;
    objectiveStruct.mat_struct          = mat_struct;
    objectiveStruct.parNormFactors      = cell2mat(mat_struct.param_values_transIso); % This will normalize the parameters to ones(size(P))
    objectiveStruct.Pb_struct.xx_c      = cell2mat(mat_struct.param_values_transIso); %Parameter constraining centre
    objectiveStruct.Pb_struct.xxlim     = [[0.001*cell2mat(mat_struct.param_values_transIso(1)) 0 0 0]' 1000*[cell2mat(mat_struct.param_values_transIso(1)) 1 1 1]']; %Parameter bounds [[lowbound] [upperbound]]
    objectiveStruct.Pb_struct.xxlim_n   = objectiveStruct.Pb_struct.xxlim./ objectiveStruct.parNormFactors'; % normalize the xxlim to use as upper and lower bound in lsqnonlin
    objectiveStruct.method              = 2; 
    objectiveStruct.febioModel_folder   = optimization_folder;
    objectiveStruct.ax                  = ax;
    
    %Optimisation settings
    maxNumberIterations             = 20; %10; 100; %Maximum number of optimization iterations
    maxNumberFunctionEvaluations    = maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
    functionTolerance               = 1e-6; %1e-3; %Tolerance on objective function value
    parameterTolerance              = 1e-6; %1e-3; %Tolerance on parameter variation
    displayTypeIterations           = 'iter';
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% start optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params_transIso_normalized   = param_values_transIso./objectiveStruct.parNormFactors;
    
    switch objectiveStruct.method
        case 1 %fminsearch and Nelder-Mead
            OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
            OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                'MaxIter',maxNumberIterations,...
                'TolFun',functionTolerance,...
                'TolX',parameterTolerance,...
                'Display',displayTypeIterations,...
                'FinDiffRelStep',1e-2,...
                'DiffMaxChange',0.5);
            [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),params_transIso_normalized,OPT_options);
        case 2 %lsqnonlin and Levenberg-Marquardt
            OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
            OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                'MaxIter',maxNumberIterations,...
                'TolFun',functionTolerance,...
                'TolX',parameterTolerance,...
                'Display',displayTypeIterations,...
                'FinDiffRelStep',1e-2,...
                'DiffMaxChange',0.5);
            [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),params_transIso_normalized(1),[],[],OPT_options);   
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% pull out the optimization results %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn_opt,objectiveStruct); 
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Unnormalize and constrain parameters %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_opt = [Pn_opt 1 1 1].*objectiveStruct.parNormFactors; %Scale back, undo normalization
    
    % Constraining parameters
    for q = 1:1:numel(P_opt)
        [P_opt(q)] = boxconstrain(P_opt(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));  
    end
    
    disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
    
    % write to diary, the distance betweem spatial sim and opt locations and the optimzed parameters
    cd(fullfile(optimization_folder))
    diary opt.log
    
    fprintf('*******************************************************************************************************\n');
    fprintf(['P_opt   = ',disp_text, '\n']);
    fprintf(['alpha   = ',sprintf('%6.16e ',P_opt(2)), '\n']);
    fprintf(['c (opt) = ',sprintf('%6.16e ',1000*P_opt(1)), 'kPa \n']);
    fprintf(['final residual norm = ',sprintf(repmat('%6.16e ',[1,numel(OPT_out.resnorm)]),OPT_out.resnorm) '\n']);
    % objective values
    fprintf(['objective sum = ' num2str(sum(Fopt),'%6.16e') '\n']);
    fprintf(['objective sum normalized = ' num2str(sum(Fopt)/numel(Fopt),'%6.16e') '\n']);
    fprintf('*******************************************************************************************************\n\n');
    diary off
    
    % Plot and store results
    
    figure; hold on;
    
    title(strcat('Nodal positions final ',...
        [' c=' num2str(1000*P_opt(1), '%0.3f') 'kPa'], [', alpha=' num2str(P_opt(2), '%0.3f') ],...
        [', rf=' num2str(P_opt(3), '%0.3f')],         [', rt=' num2str(P_opt(4), '%0.3f')]),...
        FontSize=fontSize)   

    plot3(OPT_stats_out.nodePos_spatial_sim(:,1),...
                OPT_stats_out.nodePos_spatial_sim(:,2),...
                OPT_stats_out.nodePos_spatial_sim(:,3),'.r', MarkerSize=markerSize);   
    
    plot3(nodes_exp.nodePos_maxInf(:,1),...
                nodes_exp.nodePos_maxInf(:,2),...
                nodes_exp.nodePos_maxInf(:,3),'.b', MarkerSize=markerSize);    
    
    ax = gca;
    ax.View = [5 120];
    axis off; 
    axis tight
    axis equal
    legend({'Simulation','Experiment'},'Location','northwest','TextColor','black');
    set(gca,'FontSize',fontSize);
    set(gcf, 'Position',position);
    drawnow;

    % export final fig
    cd(fullfile(optimization_folder,['images\png\alpha=' num2str(P_opt(2))]))
    export_fig('final', '-png', '-nocrop', '-trans')

    % store optimized c and alpha parameters
    opt_param_combo.c       = [opt_param_combo.c; P_opt(1)];
    opt_param_combo.alpha   = [opt_param_combo.alpha; P_opt(2)];

end

%%% plot optimized material parameters combination

figure;
plot(opt_param_combo.c*1000, opt_param_combo.alpha, 'or',...
    MarkerSize=markerSize, MarkerFaceColor='r')

xlabel('c [kPa]'); ylabel('alpha')
title({'optimal c and alpha'; ...
        'anisotropy ratios: rf=0.55, rt=0.25'})
set(gca, "FontSize", fontSize)
set(gcf, 'Position',position/2);

export_fig(fullfile(optimization_folder, 'result optimal c and alpha'), '-png', '-nocrop')

save(fullfile(optimization_folder, 'opt_param_combo'), 'opt_param_combo')

%% In script functions 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Objective Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for bulk stiffness params
function [Fopt,OPT_stats_out]=objectiveFunctionIFEA(c_params_transIso_normalized,objectiveStruct)
    c_normalized        = c_params_transIso_normalized(1);
    alpha_normalized    = 1;
    rf_normalized       = 1;
    rt_normalized       = 1;

    params_transIso_normalized = [c_normalized,alpha_normalized,rf_normalized,rt_normalized];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% get from objectiveStruct, the needed information for the optimization %%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    febioFebFileName    = objectiveStruct.febioFebFileName; 
    febioNodePosFileName= objectiveStruct.febioNodePosFileName;
    
    nodePos_spatial_exp = objectiveStruct.nodes_exp.nodePos_maxInf;
    nodeID_exp          = objectiveStruct.nodes_exp.nodeID;
    
    febioModel_folder   = objectiveStruct.febioModel_folder;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Unnormalize and constrain parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param_transIso          = params_transIso_normalized.*objectiveStruct.parNormFactors; %Scale back, undo normalization
    
    % Constraining parameters
    for param_idx=1:1:numel(param_transIso)
        [param_transIso_constrained(param_idx)]=boxconstrain(param_transIso(param_idx),objectiveStruct.Pb_struct.xxlim(param_idx,1),objectiveStruct.Pb_struct.xxlim(param_idx,2),objectiveStruct.Pb_struct.xx_c(param_idx));    
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Setting material parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Acces material parameters
    mat_struct                               = objectiveStruct.mat_struct;
    mat_struct.param_transIso_constrained    = param_transIso_constrained;
    
    % transverse Isotropic material params
    c_constrained       = (mat_struct.param_transIso_constrained(1));
    alpha_constrained   = (mat_struct.param_transIso_constrained(2));
    rf_constrained      = (mat_struct.param_transIso_constrained(3));
    rt_constrained      = (mat_struct.param_transIso_constrained(4));

    % convert material params from trans Iso reformulated to orthropic
    [mat_struct.param_values_ortho, mat_struct.param_names_ortho] = transIso_2_ortho(mat_struct.param_transIso_constrained ,mat_struct.param_names_transIso,1);
    
    % replace material parameters in febio file accordingly
    febioReplaceMatValue(febioFebFileName,mat_struct.param_names_ortho,mat_struct.param_values_ortho); 
    
    fprintf('Done setting material parameters \n')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% START FEBio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [runFlag]=runMonitorFEBio(objectiveStruct.febioAnalysis);
    
    %pause(0.1); 
    
    if runFlag==1    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Import simulated nodal positions from log file %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        [~,nodeID_sim] = febioExtractNodalPositions(febioNodePosFileName);
        if isequal(nodeID_exp,nodeID_sim) % compare if nodeID extracted from function is the same as order we have going now
            [nodePos_spatial_sim,~] = febioExtractNodalPositions(febioNodePosFileName);
        else 
            fprintf(2,'\n Warning: extracted (from node pos log file) and model nodeIDs do not match \n')
        end
    
        % update the plot
        if ~isempty(objectiveStruct.h)
            objectiveStruct.h.XData = nodePos_spatial_sim(:,1)';
            objectiveStruct.h.YData = nodePos_spatial_sim(:,2)';
            objectiveStruct.h.ZData = nodePos_spatial_sim(:,3)';
            
            drawnow;        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Derive Objective value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Derive Fopt (objective value)
        nodePos_spatial_Dev     = nodePos_spatial_exp-nodePos_spatial_sim;  
           
        switch objectiveStruct.method
            case 1
                Fopt=sum((nodePos_spatial_Dev).^2); %Sum of squared differences
            case 2
                % Fopt=nodePos_spatial_Dev(:); %Squared differences
                Fopt=vecnorm(nodePos_spatial_Dev,2,2); 
        end
        
        OPT_stats_out.nodePos_spatial_sim   = nodePos_spatial_sim;
        OPT_stats_out.Fopt                  = Fopt;

    else %Output NaN
        switch objectiveStruct.method
            case 1
                Fopt=NaN; 
            case 2
                Fopt=NaN(size(nodePos_spatial_exp(:))); %Squared differences
        end
        OPT_stats_out=[];
    end

    % export the plot
    cd(febioModel_folder)
    if not(isfolder(['images\png\alpha=' num2str(alpha_constrained)]))
        mkdir(['images\png\alpha=' num2str(alpha_constrained)])
    end
    cd(fullfile(febioModel_folder,['images\png\alpha=' num2str(alpha_constrained)]))
    dir_info = dir;
    dir_info = dir_info(3:end);
    iter_counter = string(length(dir_info)+1);
    % change the title
    
    set(objectiveStruct.ax.Title,'String',strcat('Nodal positions for iter ', num2str(iter_counter),...
        [' c=' num2str(1000*c_constrained, '%0.3f') 'kPa'], [', alpha=' num2str(alpha_constrained, '%0.3f') ],...
        [', rf=' num2str(rf_constrained, '%0.3f')],         [', rt=' num2str(rt_constrained, '%0.3f')]));
    
    export_fig(['iter' iter_counter{:}], '-png', '-nocrop', '-trans')

    % write to diary, the set optimization values and the optimized residual
    cd(fullfile(febioModel_folder))
    diary opt.log
    
    % set optimization values
    fprintf(['Step ', num2str(iter_counter) '\n']);
    fprintf('SETTING MATERIAL PARAMETERS...\n');
    fprintf(['Proposed (norm.): [c alpha rf rt] = ',sprintf(repmat('%6.16e ',[1,numel(params_transIso_normalized)]),params_transIso_normalized) '\n']);
    fprintf(['Proposed        : [c alpha rf rt] = ',sprintf(repmat('%6.16e ',[1,numel(param_transIso)]),param_transIso) '\n']);
    fprintf(['Set (constr.)   : [c alpha rf rt] = ',sprintf(repmat('%6.16e ',[1,numel(param_transIso_constrained)]),param_transIso_constrained) '\n']);

    % objective values
    fprintf(['objective sum = ' num2str(sum(Fopt),'%6.16e') '\n']);
    fprintf(['objective sum normalized = ' num2str(sum(Fopt)/numel(Fopt),'%6.16e') '\n']);
    diary off

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Parameter Reformulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [param_values_transIso,param_names_transIso] = ortho_2_transIso(param_values_ortho,~,reform)
    
    % transverse isotropy (Guccione) from fung orthotropy
    bf  = param_values_ortho(1)/param_values_ortho(10);
    bt  = param_values_ortho(2)/param_values_ortho(10);
    bft = param_values_ortho(4)/param_values_ortho(10)*2;
    c   = param_values_ortho(10);
    
    % reformulation, so can optimize for alpha and anisotropy ratios separately
    alpha   = bf + bt + bft;
    rf      = bf/alpha;
    rt      = bt/alpha;
    % rft     = 1-rf-rt;

    if reform == 1
        % reformulated parameters
        param_names_transIso    = {'c','alpha','rf','rt'};
        param_values_transIso   = [c, alpha, rf, rt];
    else
        % non reformulated params
        param_names_transIso    = {'c','bf','bt','bft'};
        param_values_transIso   = [c, bf, bt, bft];
    end

end

function [param_values_ortho,param_names_ortho] = transIso_2_ortho(param_values_transIso,~,reform)
    if reform == 1
        c       = (param_values_transIso(1));
        alpha   = (param_values_transIso(2));
        rf      = (param_values_transIso(3));
        rt      = (param_values_transIso(4));
        rft     = 1-rf-rt;
        
        bf      = rf*alpha;
        bt      = rt*alpha;
        bft     = rft*alpha;
    else
        c       = (param_values_transIso(1));
        bf      = (param_values_transIso(2));
        bt      = (param_values_transIso(3));
        bft     = (param_values_transIso(4));    
        alpha   = bf + bt + bft;
    end

    % transverse isotropy (Guccione) to fung orthotropy
    param_values_ortho(1)   = bf*c;
    param_values_ortho(2)   = bt*c;
    param_values_ortho(3)   = bt*c;
    param_values_ortho(4)   = bft*c/2;
    param_values_ortho(5)   = bt*c/2;
    param_values_ortho(6)   = bft*c/2;
    param_values_ortho(7)   = 0;
    param_values_ortho(8)   = 0;
    param_values_ortho(9)   = 0;
    param_values_ortho(10)  = c;
    param_values_ortho(11)  = 100*(alpha*c);   

    param_names_ortho       = {'E1','E2','E3','G12','G23','G31','v12','v23','v31','c','k'};

end

