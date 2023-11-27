wirelessnetworkSupportPackageCheck
rng('default'); % Reset the random number generator
simParameters = []; % Clear simParameters variable
simParameters.NumFramesSim = 100; % Simulation time in terms of number of 10 ms frames
simParameters.SchedulingType = 1; % Set the value to 0 (slot-based scheduling) or 1 (symbol-based scheduling)
simParameters.NumUEs = 4;
% Assign position to the UEs assuming that the gNB is at (0, 0, 0). N-by-3
% matrix where 'N' is the number of UEs. Each row has (x, y, z) position of a
% UE (in meters)
simParameters.UEPosition = [100 0 0;
                            150 0 0;
                            300 0 0;
                            400 0 0];
% Validate the UE positions
validateattributes(simParameters.UEPosition,{'numeric'},{'nonempty','real','nrows',simParameters.NumUEs,'ncols',3, ...
    'finite'},'simParameters.UEPosition','UEPosition');
simParameters.NumRBs = 25;
simParameters.SCS = 30; % kHz  30 for urban 15 for rural
simParameters.DLBandwidth = 400e6; % Hz 400 for urban 100 for rural
simParameters.ULBandwidth = 400e6; % Hz 400 for urban 100 for rural
simParameters.DLCarrierFreq = 6e9; % Hz 6 for urban 1 for rural
simParameters.ULCarrierFreq = 6e9; % Hz 6 for urban 1 for rural
simParameters.DLULPeriodicity = 5; % Duration of the DL-UL pattern in ms
simParameters.NumDLSlots = 2; % Number of consecutive full DL slots at the beginning of each DL-UL pattern
simParameters.NumDLSyms = 8; % Number of consecutive DL symbols in the beginning of the slot following the last full DL slot
simParameters.NumULSyms = 4; % Number of consecutive UL symbols in the end of the slot preceding the first full UL slot
simParameters.NumULSlots = 2; % Number of consecutive full UL slots at the end of each DL-UL pattern
simParameters.SchedulerStrategy = 'PF'; % Supported scheduling strategies: 'PF', 'RR' and 'BestCQI'
simParameters.TTIGranularity = 4;
simParameters.RBAllocationLimitUL = 15; % For PUSCH
simParameters.RBAllocationLimitDL = 15; % For PDSCH
simParameters.BSRPeriodicity = 1; % Buffer status report transmission periodicity (in ms)
simParameters.PUSCHPrepTime = 200; % In microseconds
simParameters.ChannelUpdatePeriodicity = 0.2; % In sec
simParameters.CQIDelta = 1;
% Mapping between distance from gNB (first column in meters) and maximum achievable UL CQI value (second column)
simParameters.CQIvsDistance = [ 
    200  15;
    300  12;
    500  10;    
    1000  8;
    1200  7];
simParameters.DMRSTypeAPosition = 2; % Type-A DM-RS position as 2 or 3
% PUSCH DM-RS configuration
simParameters.PUSCHDMRSAdditionalPosTypeB = 0;
simParameters.PUSCHDMRSAdditionalPosTypeA = 0;
simParameters.PUSCHDMRSConfigurationType = 1;
% PDSCH DM-RS configuration
simParameters.PDSCHDMRSAdditionalPosTypeB = 0;
simParameters.PDSCHDMRSAdditionalPosTypeA = 0;
simParameters.PDSCHDMRSConfigurationType = 1;
% Set the periodic DL and UL application traffic pattern for UEs
dlAppDataRate = 16e4*ones(simParameters.NumUEs,1); % DL application data rate in kilo bits per second (kbps)
ulAppDataRate = 16e4*ones(simParameters.NumUEs,1); % UL application data rate in kbps
% Validate the DL application data rate
validateattributes(dlAppDataRate, {'numeric'},{'nonempty','vector','numel',simParameters.NumUEs,'finite','>',0}, ...
    'dlAppDataRate','dlAppDataRate');
% Validate the UL application data rate
validateattributes(ulAppDataRate,{'numeric'},{'nonempty','vector','numel',simParameters.NumUEs,'finite','>',0}, ...
    'ulAppDataRate', 'ulAppDataRate');
simParameters.CQIVisualization = false;
simParameters.RBVisualization = false;
enableTraces = true;
simParameters.NumMetricsSteps = 20;
parametersLogFile = 'simParameters'; % For logging the simulation parameters
simulationLogFile = 'simulationLogs'; % For logging the simulation traces
simulationMetricsFile = 'simulationMetrics'; % For logging the simulation metrics
simParameters.DuplexMode = 1; % FDD (Value as 0) or TDD (Value as 1)
simParameters.NCellID = 1; % Physical cell ID
simParameters.Position = [0 0 0]; % Position of gNB in (x,y,z) coordinates
numSlotsSim = (simParameters.NumFramesSim * 10 * simParameters.SCS)/15;
if simParameters.SchedulingType % Symbol-based scheduling
    simParameters.PUSCHMappingType = 'B';
    simParameters.PDSCHMappingType = 'B';
else % Slot-based scheduling
    simParameters.PUSCHMappingType = 'A';
    simParameters.PDSCHMappingType = 'A';
end
simParameters.MetricsStepSize = ceil(numSlotsSim / simParameters.NumMetricsSteps);
numLogicalChannels = 1;
simParameters.LCHConfig.LCID = 4;
simParameters.RLCConfig.EntityType = 2;
lchInfo = repmat(struct('RNTI',[],'LCID',[],'EntityDir',[]), [simParameters.NumUEs 1]);
for idx = 1:simParameters.NumUEs
    lchInfo(idx).RNTI = idx;
    lchInfo(idx).LCID = simParameters.LCHConfig.LCID;
    lchInfo(idx).EntityDir = simParameters.RLCConfig.EntityType;
end
rlcChannelConfigStruct.LCGID = 1; % Mapping between logical channel and logical channel group ID
rlcChannelConfigStruct.Priority = 1; % Priority of each logical channel
rlcChannelConfigStruct.PBR = 8; % Prioritized bitrate (PBR), in kilobytes per second, of each logical channel
rlcChannelConfigStruct.BSD = 10; % Bucket size duration (BSD), in ms, of each logical channel
rlcChannelConfigStruct.EntityType = simParameters.RLCConfig.EntityType;
rlcChannelConfigStruct.LogicalChannelID = simParameters.LCHConfig.LCID;
simParameters.maxRLCSDULength = 9000;
maxUECQIs = zeros(simParameters.NumUEs, 1); % To store the maximum achievable CQI value for UEs
for ueIdx = 1:simParameters.NumUEs
    % Based on the distance of the UE from gNB, find matching row in CQIvsDistance mapping
    matchingRowIdx = find(simParameters.CQIvsDistance(:, 1) > simParameters.UEPosition(ueIdx,1));
    if isempty(matchingRowIdx)
        maxUECQIs(ueIdx) = simParameters.CQIvsDistance(end, 2);
    else
        maxUECQIs(ueIdx) = simParameters.CQIvsDistance(matchingRowIdx(1), 2);
    end
end
simParameters.InitialChannelQualityUL = zeros(simParameters.NumUEs, simParameters.NumRBs); % To store current UL CQI values on the RBs for different UEs
simParameters.InitialChannelQualityDL = zeros(simParameters.NumUEs, simParameters.NumRBs); % To store current DL CQI values on the RBs for different UEs
for ueIdx = 1:simParameters.NumUEs
    % Assign random CQI values for the RBs, limited by the maximum achievable CQI value
    simParameters.InitialChannelQualityUL(ueIdx, :) = randi([1 maxUECQIs(ueIdx)], 1, simParameters.NumRBs);
    % Initially, DL and UL CQI values are assumed to be equal
    simParameters.InitialChannelQualityDL(ueIdx, :) = simParameters.InitialChannelQualityUL(ueIdx, :);
end
gNB = hNRGNB(simParameters); 
switch(simParameters.SchedulerStrategy)
    case 'RR' % Round-robin scheduler
        scheduler = hNRSchedulerRoundRobin(simParameters);
    case 'PF' % Proportional fair scheduler
        scheduler = hNRSchedulerProportionalFair(simParameters);
    case 'BestCQI' % Best CQI scheduler
        scheduler = hNRSchedulerBestCQI(simParameters);
end
addScheduler(gNB, scheduler); % Add scheduler to gNB

gNB.PhyEntity = hNRGNBPassThroughPhy(simParameters); % Add passthrough PHY
configurePhy(gNB, simParameters);
setPhyInterface(gNB); % Set the interface to PHY layer
UEs = cell(simParameters.NumUEs, 1);
for ueIdx = 1:simParameters.NumUEs
    simParameters.Position = simParameters.UEPosition(ueIdx, :); % Position of the UE
    UEs{ueIdx} = hNRUE(simParameters, ueIdx);
    UEs{ueIdx}.PhyEntity = hNRUEPassThroughPhy(simParameters, ueIdx); % Add passthrough PHY
    configurePhy(UEs{ueIdx}, simParameters);
    setPhyInterface(UEs{ueIdx}); % Set the interface to PHY layer

    % Initialize the UL CQI values at gNB scheduler
    channelQualityInfoUL = struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityUL(ueIdx, :));
    updateChannelQualityUL(gNB.MACEntity.Scheduler, channelQualityInfoUL);
    
    % Initialize the DL CQI values at gNB scheduler
    channelQualityInfoDL = struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityDL(ueIdx, :));
    updateChannelQualityDL(gNB.MACEntity.Scheduler, channelQualityInfoDL);
 
    % Initialize the DL CQI values at UE for packet error probability estimation
    updateChannelQualityDL(UEs{ueIdx}.MACEntity, channelQualityInfoDL);

    % Setup logical channel at gNB for the UE
    configureLogicalChannel(gNB, ueIdx, rlcChannelConfigStruct);
    % Setup logical channel at UE
    configureLogicalChannel(UEs{ueIdx}, ueIdx, rlcChannelConfigStruct);

    % Create an object for On-Off network traffic pattern and add it to the
    % specified UE. This object generates the uplink (UL) data traffic on the UE
    ulApp = networkTrafficOnOff('GeneratePacket', true, ...
        'OnTime', simParameters.NumFramesSim/100, 'OffTime', 0, 'DataRate', ulAppDataRate(ueIdx));
    UEs{ueIdx}.addApplication(ueIdx, simParameters.LCHConfig.LCID, ulApp);

    % Create an object for On-Off network traffic pattern for the specified
    % UE and add it to the gNB. This object generates the downlink (DL) data
    % traffic on the gNB for the UE
    dlApp = networkTrafficOnOff('GeneratePacket', true, ...
        'OnTime', simParameters.NumFramesSim/100, 'OffTime', 0, 'DataRate', dlAppDataRate(ueIdx));
    gNB.addApplication(ueIdx, simParameters.LCHConfig.LCID, dlApp);
end
% Initialize wireless network simulator
nrNodes = [{gNB}; UEs];
networkSimulator = hWirelessNetworkSimulator(nrNodes);
if enableTraces

    % RLC metrics are logged for every 1 slot duration
    simRLCLogger = hNRRLCLogger(simParameters, lchInfo, networkSimulator, gNB, UEs);
    
    simSchedulingLogger = hNRSchedulingLogger(simParameters, networkSimulator, gNB, UEs);
    
    % Create an object for CQI and RB grid visualization
    if simParameters.CQIVisualization || simParameters.RBVisualization
        gridVisualizer = hNRGridVisualizer(simParameters, 'MACLogger', simSchedulingLogger);
    end
end
metricsVisualizer = hNRMetricsVisualizer(simParameters, 'EnableSchedulerMetricsPlots', true, ...
    'EnableRLCMetricsPlots', true, 'LCHInfo', lchInfo, 'NetworkSimulator', networkSimulator, 'GNB', gNB, 'UEs', UEs);
% Calculate the simulation duration (in seconds) from 'NumFramesSim'
simulationTime = simParameters.NumFramesSim * 1e-2;
% Run the simulation
run(networkSimulator, simulationTime);
metrics = getMetrics(metricsVisualizer);
save(simulationMetricsFile, 'metrics'); % Save simulation metrics in a MAT-file
if enableTraces
    % Read the logs and write them to MAT-files
    % Get the logs
    simulationLogs = cell(1,1);
    logInfo = struct('TimeStepLogs',[], 'SchedulingAssignmentLogs',[] ,'RLCLogs', []);
    [logInfo.TimeStepLogs] = getSchedulingLogs(simSchedulingLogger);
    logInfo.SchedulingAssignmentLogs = getGrantLogs(simSchedulingLogger); % Scheduling assignments log
    logInfo.RLCLogs = getRLCLogs(simRLCLogger); % RLC statistics logs
    simulationLogs{1} = logInfo;

    save(simulationLogFile, 'simulationLogs'); % Save simulation logs in a MAT-file
    save(parametersLogFile, 'simParameters'); % Save simulation parameters in a MAT-file
end