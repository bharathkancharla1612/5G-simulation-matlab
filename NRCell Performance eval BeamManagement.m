wirelessnetworkSupportPackageCheck
rng('default');                  % Reset the random number generator
simParameters = [];              % Clear the simParameters variable
simParameters.NumFramesSim = 10; % Simulation time in terms of number of 10 ms frames
enableBeamforming = true;
if enableBeamforming
    % Number of SSBs
    numSSBs = 4;
    % Total azimuthal range for all SSBs corresponding to a 120 degrees sectorized cell
    azSweepRange = [-60, 60];
    validateattributes(azSweepRange,{'numeric'},{'nonempty','real', 'vector', ...
        '>=',-60,'<=',60,'finite'}, ...
        'azSweepRange','azimuthal range');
    % Sweep range for non-overlapping SSB directions
    ssbSweepRange = diff(azSweepRange)/numSSBs;
    % Azimuthal angles for SSB beam sweeping in the azimuthal range
    ssbTxAngles = [-60 -20 20 60];
    % Elevation range for all SSBs sweeping in the azimuthal angles
    elSweepRange = [-45, 0];
end
simParameters.NumUEs = 3;
% Position of gNB in (x,y,z) coordinates
simParameters.GNBPosition = [0, 0, 0];
% Specify UE position in spherical coordinates (r,azimuth,elevation)
% relative to gNB. Here 'r' is the distance, 'azimuth' is the angle in the
% horizontal plane and elevation is the angle in the vertical plane. When
% beamforming is enabled, all the UE positions must lie with in the azimuth
% and elevation range limits of SSBs as defined above. It is a matrix of
% size N-by-3 where N is the number of UEs. Row 'i' of the matrix contains
% the coordinates of UE with RNTI 'i'
ueRelPosition = [500, -5, -5;
    500, -60, -30;
    500, -10, -20];
if enableBeamforming
    validateattributes(ueRelPosition(:,2),{'numeric'},{'nonempty','real','vector', ...
        '>=',azSweepRange(1),'<=',azSweepRange(2),'finite'}, ...
        'ueRelPosition(:,2)','azimuth angles');
    validateattributes(ueRelPosition(:,3),{'numeric'},{'nonempty','real','vector', ...
        '>=',elSweepRange(1),'<=',elSweepRange(2),'finite'}, ...
        'ueRelPosition(:,3)','elevation angles');
end
% Convert Spherical to Cartesian coordinates considering gNB position as origin
[xPos, yPos, zPos] = sph2cart(deg2rad(ueRelPosition(:, 2)),deg2rad(ueRelPosition(:, 3)), ...
    ueRelPosition(:, 1));
% Global coordinates of UEs
simParameters.UEPosition = [xPos, yPos, zPos] + simParameters.GNBPosition;
% Validate the UE positions
validateattributes(simParameters.UEPosition, {'numeric'},{'nonempty','real', ...
    'nrows',simParameters.NumUEs,'ncols',3,'finite'}, ...
    'simParameters.UEPosition','UEPosition');
simParameters.NumRBs = 30;             % Moderate number of resource blocks
simParameters.SCS = 30;                % Increase the subcarrier spacing for better spectral efficiency (kHz)
simParameters.DLBandwidth = 10e6;      % Increase the downlink bandwidth for higher capacity (Hz)
simParameters.ULBandwidth = 10e6;      % Increase the uplink bandwidth for higher capacity (Hz)
simParameters.DLCarrierFreq = 700e6;   % Use a lower downlink carrier frequency suitable for rural coverage (Hz)
simParameters.ULCarrierFreq = 600e6;   % Use a lower uplink carrier frequency suitable for rural coverage (Hz)
gNBTxArraySize = [2 4 2];
ueRxArraySize = repmat([1 1 2],simParameters.NumUEs,1);
lambda = physconst('LightSpeed')/simParameters.DLCarrierFreq;
simParameters.TxAntPanel = phased.NRRectangularPanelArray('ElementSet', ...
    repmat({phased.NRAntennaElement},1,gNBTxArraySize(3)), ...
    'Size',[gNBTxArraySize(1:2) 1 1],'Spacing', ...
    [0.5*lambda 0.5*lambda 1 1]);
csirs = cell(1,simParameters.NumUEs);
csirsOffset = [5 6 10 11];
for csirsIdx = 1:simParameters.NumUEs
    csirs{csirsIdx} = nrCSIRSConfig('NID',1,'NumRB',simParameters.NumRBs, ...
        'RowNumber',11,'SubcarrierLocations',[1 3 5 7], ...
        'SymbolLocations',0,'CSIRSPeriod',[20 csirsOffset(csirsIdx)]);
end
simParameters.CSIRSConfig = csirs;
if enableBeamforming
    simParameters.NumCSIRSBeams = 4;
    csirsConfigRSRP = cell(1,numSSBs);
    for ssbIdx = 1:numSSBs
        csirsConfig = nrCSIRSConfig;
        csirsConfig.NID = 1;
        csirsConfig.CSIRSType = repmat({'nzp'},1,simParameters.NumCSIRSBeams);
        csirsConfig.CSIRSPeriod = [40 0];
        csirsConfig.Density = repmat({'one'},1,simParameters.NumCSIRSBeams);
        csirsConfig.RowNumber = repmat(2,1,simParameters.NumCSIRSBeams);
        csirsConfig.SymbolLocations = {1,5,6,7};
        csirsConfig.SubcarrierLocations = repmat({ssbIdx-1},1, ...
            simParameters.NumCSIRSBeams);
        csirsConfig.NumRB = simParameters.NumRBs;
        csirsConfigRSRP{ssbIdx} = csirsConfig;
    end
    simParameters.CSIRSConfigRSRP = csirsConfigRSRP;
end
simParameters.DownlinkSINR90pc = [-0.4600 4.5400 9.5400 14.0500 16.540 19.0400 ...
    20.5400 23.0400 25.0400 27.4300 29.9300 30.4300 ...
    32.4300 35.4300 38.4300];
simParameters.GNBTxPower = 34; % In dBm
channelModelDL = cell(1,simParameters.NumUEs);
waveformInfo = nrOFDMInfo(simParameters.NumRBs,simParameters.SCS);
dlchannel = nrCDLChannel('DelayProfile','CDL-D','DelaySpread',300e-9);
for ueIdx = 1:simParameters.NumUEs
    cdl = hMakeCustomCDL(dlchannel);
    % Configure the channel seed based on the UE number
    % Independent fading for each UE
    cdl.Seed = cdl.Seed + (ueIdx - 1);
    % Compute the LOS angle from gNB to UE
    [~,depAngle] = rangeangle(simParameters.UEPosition(ueIdx, :)', ...
        simParameters.GNBPosition');
    % Configure the azimuth and zenith angle offsets for this UE
    cdl.AnglesAoD(:) = cdl.AnglesAoD(:) + depAngle(1);
    % Convert elevation angle to zenith angle
    cdl.AnglesZoD(:) = cdl.AnglesZoD(:) - cdl.AnglesZoD(1) + (90 - depAngle(2));
    % Compute the range angle from UE to gNB
    [~,arrAngle] = rangeangle(simParameters.GNBPosition',...
        simParameters.UEPosition(ueIdx, :)');
    % Configure the azimuth and zenith arrival angle offsets for this UE
    cdl.AnglesAoA(:) = cdl.AnglesAoA(:) - cdl.AnglesAoA(1) + arrAngle(1);
    % Convert elevation angle to zenith angle
    cdl.AnglesZoA(:) = cdl.AnglesZoA(:) - cdl.AnglesZoA(1) + (90 - arrAngle(2));
    cdl.CarrierFrequency = simParameters.DLCarrierFreq;
    cdl.TransmitAntennaArray = simParameters.TxAntPanel;
    cdl.ReceiveAntennaArray.Size = [ueRxArraySize(ueIdx, :) 1 1];
    cdl.SampleRate = waveformInfo.SampleRate;
    channelModelDL{ueIdx} = cdl;
end
% Number of transmit antennas
numTxAntennas = prod(gNBTxArraySize);
if enableBeamforming
    enableVerticalBeamSweep = true;
    simParameters.BeamWeightTable = zeros(numTxAntennas, ...
        simParameters.NumCSIRSBeams*numSSBs);
    txArrayStv = phased.SteeringVector('SensorArray',simParameters.TxAntPanel, ...
        'PropagationSpeed',physconst('LightSpeed'), ...
        'IncludeElementResponse',true);
    azBW = beamwidth(simParameters.TxAntPanel,simParameters.DLCarrierFreq, ...
        'Cut','Azimuth');
    elBW = beamwidth(simParameters.TxAntPanel,simParameters.DLCarrierFreq, ...
        'Cut','Elevation');
    for beamIdx = 1:numSSBs
        % Get the azimuthal sweep range for every SSB based on the SSB transmit
        % beam direction and its beamwidth in azimuth plane
        azSweepRangeSSB = ssbTxAngles(beamIdx) + [-ssbSweepRange/2 ssbSweepRange/2];
        % Get the azimuth and elevation angle pairs for all CSI-RS transmit beams
        csirsBeamAng = hGetBeamSweepAngles(simParameters.NumCSIRSBeams,azSweepRangeSSB, ...
            elSweepRange,azBW,elBW,enableVerticalBeamSweep);
        for csirsBeamIdx = 1:simParameters.NumCSIRSBeams
            simParameters.BeamWeightTable(:, simParameters.NumCSIRSBeams*(beamIdx-1) ...
                + csirsBeamIdx) = txArrayStv(simParameters.DLCarrierFreq, ...
                csirsBeamAng(:,csirsBeamIdx))/sqrt(numTxAntennas);
        end
    end
    
    % Initialize the beam indices for all UEs with highest L1-RSRP (for visualization purpose)
    beamIndices = -1*ones(1,simParameters.NumUEs);
end
% Supported scheduling strategies: 'PF', 'RR', and 'BestCQI'
simParameters.SchedulerStrategy = 'PF';
simParameters.RBAllocationLimitDL = 15; % For PDSC
simParameters.CQIVisualization = true;
simParameters.RBVisualization = true;
enableTraces = true;
simParameters.NumMetricsSteps = 10;
parametersLogFile = 'simParameters';         % For logging the simulation parameters
simulationLogFile = 'simulationLogs';        % For logging the simulation traces
simulationMetricsFile = 'simulationMetrics'; % For logging the simulation metrics
% DL application data rate in kilo bits per second (kbps)
dlAppDataRate = 40e3*ones(1,simParameters.NumUEs);
% Validate the DL application data rate
validateattributes(dlAppDataRate,{'numeric'},{'nonempty','vector','numel', ...
    simParameters.NumUEs,'finite','>',0},'dlAppDataRate', ...
    'dlAppDataRate');
simParameters.DuplexMode = 0;                    % FDD (Value as 0) or TDD (Value as 1)
simParameters.SchedulingType = 0;                % Slot-based scheduling
simParameters.NCellID = 1;                       % Physical cell ID
simParameters.GNBTxAnts = numTxAntennas;
simParameters.UERxAnts = prod(ueRxArraySize,2);
% [N1 N2] as per 3GPP TS 38.214 Table 5.2.2.2.1-2
csiReportConfig.PanelDimensions = [8 1];
csiReportConfig.CQIMode = 'Wideband';              % 'Wideband' or 'Subband'
% Set codebook mode as 1 or 2. It is applicable only when the number of
% transmission layers is 1 or 2 and number of CSI-RS ports is greater than 2
csiReportConfig.CodebookMode = 1;
simParameters.CSIReportConfig = {csiReportConfig};
if enableBeamforming
    simParameters.SSBIndex = computeSSBToUE(simParameters,ssbTxAngles,ssbSweepRange,ueRelPosition);
end
numSlotsSim = (simParameters.NumFramesSim * 10 * simParameters.SCS)/15;
simParameters.MetricsStepSize = ceil(numSlotsSim/simParameters.NumMetricsSteps);
numLogicalChannels = 1;
simParameters.LCHConfig.LCID = 4;
simParameters.RLCConfig.EntityType = 2;
% Mapping between logical channel and logical channel group ID
rlcChannelConfigStruct.LCGID = 1;
% Priority of each logical channel
rlcChannelConfigStruct.Priority = 1;
% Prioritized bitrate (PBR), in kilobytes per second, of each logical channel
rlcChannelConfigStruct.PBR = 8;
% Bucket size duration (BSD), in ms, of each logical channel
rlcChannelConfigStruct.BSD = 10;
rlcChannelConfigStruct.EntityType = simParameters.RLCConfig.EntityType;
rlcChannelConfigStruct.LogicalChannelID = simParameters.LCHConfig.LCID;
simParameters.Position = simParameters.GNBPosition;
% Create gNB node
gNB = hNRGNB(simParameters);
% Create scheduler
switch(simParameters.SchedulerStrategy)
    % Round robin scheduler
    case 'RR'
        scheduler = hNRSchedulerRoundRobin(simParameters);
    % Proportional fair scheduler
    case 'PF'
        scheduler = hNRSchedulerProportionalFair(simParameters);
    % Best CQI scheduler
    case 'BestCQI'
        scheduler = hNRSchedulerBestCQI(simParameters);
end

% Add scheduler to gNB
addScheduler(gNB,scheduler);
% Create the PHY instance
gNB.PhyEntity = hNRGNBPhy(simParameters);
% Configure the PHY
configurePhy(gNB,simParameters);
% Set the interface to PHY
setPhyInterface(gNB);

% Create the set of UE nodes
UEs = cell(simParameters.NumUEs,1);
ueParam = simParameters;
for ueIdx=1:simParameters.NumUEs
    % Position of the UE
    ueParam.Position = simParameters.UEPosition(ueIdx,:);
    ueParam.UERxAnts = simParameters.UERxAnts(ueIdx);
    if enableBeamforming
        ueParam.SSBIdx = simParameters.SSBIndex(ueIdx);
    end
    % Assuming same CSI report configuration for all UEs
    ueParam.CSIReportConfig = simParameters.CSIReportConfig{1};
    ueParam.ChannelModel = channelModelDL{ueIdx};
    UEs{ueIdx} = hNRUE(ueParam,ueIdx);
    % Create the PHY instance
    UEs{ueIdx}.PhyEntity = hNRUEPhy(ueParam,ueIdx);
    % Configure the PHY
    configurePhy(UEs{ueIdx}, ueParam);
    % Set up the interface to PHY
    setPhyInterface(UEs{ueIdx});

    % Set up logical channel at gNB for the UE
    configureLogicalChannel(gNB,ueIdx,rlcChannelConfigStruct);

    % Set up logical channel at UE
    configureLogicalChannel(UEs{ueIdx},ueIdx,rlcChannelConfigStruct);

    % Set up application traffic
    % Create an object for on-off network traffic pattern for the specified
    % UE and add it to the gNB. This object generates the downlink data
    % traffic on the gNB for the UE
    dlApp = networkTrafficOnOff('GeneratePacket',true, ...
        'OnTime',simParameters.NumFramesSim*10e-3,'OffTime',0, ...
        'DataRate',dlAppDataRate(ueIdx));
    addApplication(gNB,ueIdx,simParameters.LCHConfig.LCID,dlApp);
end
% Initialize wireless network simulator
nrNodes = [{gNB}; UEs];
networkSimulator = hWirelessNetworkSimulator(nrNodes);
nodes = struct('UEs', {UEs}, 'GNB', gNB);
linkDir = 0; % Indicates DL
if enableTraces
    % Create an object for MAC traces logging
    simSchedulingLogger = hNRSchedulingLogger(simParameters, networkSimulator, gNB, UEs, linkDir);

    % Create an object for PHY traces logging
    simPhyLogger = hNRPhyLogger(simParameters, networkSimulator, gNB, UEs);

    % Create an object for CQI and RB grid visualization
    if simParameters.CQIVisualization || simParameters.RBVisualization
        gridVisualizer = hNRGridVisualizer(simParameters, 'MACLogger', simSchedulingLogger, 'VisualizationFlag', linkDir);
    end
end
metricsVisualizer = hNRMetricsVisualizer(simParameters, 'EnableSchedulerMetricsPlots', true, 'EnablePhyMetricsPlots', true, 'VisualizationFlag', linkDir, ...
    'NetworkSimulator', networkSimulator, 'GNB', gNB, 'UEs', UEs);
% Calculate the simulation duration (in seconds) from 'NumFramesSim'
simulationTime = simParameters.NumFramesSim * 1e-2;
% Run the simulation
run(networkSimulator, simulationTime);
if enableBeamforming
    hPlotBeamformingPattern(simParameters, gNB, beamIndices);
end
metrics = getMetrics(metricsVisualizer);
save(simulationMetricsFile,'metrics');
displayPerformanceIndicators(metricsVisualizer);
if enableTraces
    simulationLogs = cell(1,1);
    if simParameters.DuplexMode == 0 % FDD
        logInfo = struct('DLTimeStepLogs',[],'SchedulingAssignmentLogs',[],'BLERLogs',[],'AvgBLERLogs',[]);
        [logInfo.DLTimeStepLogs,~] = getSchedulingLogs(simSchedulingLogger);
    else % TDD
        logInfo = struct('TimeStepLogs',[],'SchedulingAssignmentLogs',[], ...
            'BLERLogs',[],'AvgBLERLogs',[]);
        logInfo.TimeStepLogs = getSchedulingLogs(simSchedulingLogger);
    end
    % BLER logs
    [logInfo.BLERLogs,logInfo.AvgBLERLogs] = getBLERLogs(simPhyLogger);
    % Scheduling assignments log
    logInfo.SchedulingAssignmentLogs = getGrantLogs(simSchedulingLogger);
    simulationLogs{1} = logInfo;
    % Save simulation parameters in a MAT-file
    save(parametersLogFile,'simParameters');
    % Save simulation logs in a MAT-file
    save(simulationLogFile,'simulationLogs');
end

function ssbIndex = computeSSBToUE(simParameters,ssbTxAngles,ssbSweepRange,ueRelPosition)
%   SSBINDEX = computeSSBtoUE(SIMPARAMETERS,SSBTXANGLES,SSBSWEEPRANGE,UERELPOSITION)
%   computes an SSB to UE based on SSBTXANGLES and SSBSWEEPRANGE. If a UE
%   falls within the sweep range of an SSB then the corresponding SSBINDEX
%   is assigned. If a UE does not fall under any of the SSBs then the 
%   function approximates UE to its nearest SSB.
ssbIndex = zeros(1,simParameters.NumUEs);
for i=1:simParameters.NumUEs
    if (ueRelPosition(i,2) > (ssbTxAngles(1) - ssbSweepRange/2)) && ...
            (ueRelPosition(i,2) <= (ssbTxAngles(1) + ssbSweepRange/2))
        ssbIndex(i) = 1;
    elseif (ueRelPosition(i,2) > (ssbTxAngles(2) - ssbSweepRange/2)) && ...
            (ueRelPosition(i,2) <= (ssbTxAngles(2) + ssbSweepRange/2))
        ssbIndex(i) = 2;
    elseif (ueRelPosition(i,2) > (ssbTxAngles(3) - ssbSweepRange/2)) && ...
            (ueRelPosition(i,2) <= (ssbTxAngles(3) + ssbSweepRange/2))
        ssbIndex(i) = 3;
    elseif (ueRelPosition(i,2) > (ssbTxAngles(4) - ssbSweepRange/2)) && ...
            (ueRelPosition(i,2) <= (ssbTxAngles(4) + ssbSweepRange/2))
        ssbIndex(i) = 4;
    else % Approximate UE to the nearest SSB
        % Initialize an array to store the difference between the UE azimuth and
        % SSB azimuth beam sweep angles
        angleDiff = zeros(1,length(ssbTxAngles));
        for angleIdx = 1: length(ssbTxAngles)
            angleDiff(angleIdx) = abs(ssbTxAngles(angleIdx) - ueRelPosition(i,2));
        end
        % Minimum azimuth angle difference
        [~,idx] = min(angleDiff);
        ssbIndex(i) = idx;
    end
end
end

