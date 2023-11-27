%% 

wirelessnetworkSupportPackageCheck
rng("default") % Reset the random number generator
numFrameSimulation = 50; % Simulation time in terms of number of 10 ms frames
networkSimulator = wirelessNetworkSimulator.init;
gNB = nrGNB(DuplexMode="FDD",CarrierFrequency=6e9,ChannelBandwidth=400e6,SubcarrierSpacing=120e3,ReceiveGain=11);

% For Rural
% SCS - 30
% Carrier frequency - 1GHz
% Bandwidth - 100MHz

% For Urban
% SCS - 120
% Carrier frequency - 6GHz
% Bandwidth - 400MHz

configureScheduler(gNB,Scheduler="RoundRobin",ResourceAllocationType=0);
uePositions = [100 0 0; 250 0 0; 700 0 0; 750 0 0];
ueNames = "UE-" + (1:size(uePositions,1));
UEs = nrUE(Name=ueNames,Position=uePositions,ReceiveGain=11);
load("NRFDDAppConfig.mat")
% Validate the host device type for the applications configured
validateattributes(AppConfig.HostDevice,{'numeric'},{"nonempty","integer",">=",0,"<=",1},"AppConfig.HostDevice", ...
    "HostDevice");
load("NRFDDRLCChannelConfig.mat")
rlcBearerConfig = cell(1, length(UEs));
for rlcBearerInfoIdx = 1:size(RLCChannelConfig, 1)
    rlcBearerConfigStruct = table2struct(RLCChannelConfig(rlcBearerInfoIdx, 2:end));
    ueIdx = RLCChannelConfig.RNTI(rlcBearerInfoIdx);

    % Create an RLC bearer configuration object with the specified logical
    % channel ID
    rlcBearerObj = nrRLCBearerConfig(LogicalChannelID=rlcBearerConfigStruct.LogicalChannelID);
    % Set the other properties of the configuration object with specified values
    rlcBearerObj.LogicalChannelGroup = rlcBearerConfigStruct.LogicalChannelGroup;
    rlcBearerObj.SNFieldLength = rlcBearerConfigStruct.SNFieldLength;
    rlcBearerObj.BufferSize = rlcBearerConfigStruct.BufferSize;
    rlcBearerObj.ReassemblyTimer=rlcBearerConfigStruct.ReassemblyTimer;
    rlcBearerObj.Priority=rlcBearerConfigStruct.Priority;
    rlcBearerObj.PrioritizedBitRate=rlcBearerConfigStruct.PrioritizedBitRate;
    rlcBearerObj.BucketSizeDuration=rlcBearerConfigStruct.BucketSizeDuration;
    rlcBearerObj.RLCEntityType=rlcBearerConfigStruct.RLCEntityType;

    rlcBearerConfig{ueIdx} = [rlcBearerConfig{ueIdx} rlcBearerObj];
end
for ueIdx = 1:length(UEs)
    connectUE(gNB,UEs(ueIdx),BSRPeriodicity=5,RLCBearerConfig=rlcBearerConfig{ueIdx})
end
for appIdx = 1:size(AppConfig,1)
    
    % Create an object for On-Off network traffic pattern
    app = networkTrafficOnOff(PacketSize=AppConfig.PacketSize(appIdx),GeneratePacket=true, ...
            OnTime=numFrameSimulation/100,OffTime=0,DataRate=AppConfig.DataRate(appIdx));

    if AppConfig.HostDevice(appIdx) == 0
        % Add traffic pattern that generates traffic on downlink
        addTrafficSource(gNB,app,DestinationNode=UEs(AppConfig.RNTI(appIdx)),LogicalChannelID=AppConfig.LogicalChannelID(appIdx))
    else
        % Add traffic pattern that generates traffic on uplink
        addTrafficSource(UEs(AppConfig.RNTI(appIdx)),app,LogicalChannelID=AppConfig.LogicalChannelID(appIdx))
    end
end
addNodes(networkSimulator,gNB)
addNodes(networkSimulator,UEs)
enableTraces = true;
cqiVisualization = true;
rbVisualization = true;
if enableTraces
    % Create an object for RLC traces logging
    simRLCLogger = helperNRRLCLogger(numFrameSimulation,gNB,UEs);
    % Create an object for scheduler traces logging
    simSchedulingLogger = helperNRSchedulingLogger(numFrameSimulation,gNB,UEs);
    % Create an object for CQI and RB grid visualization
    gridVisualizer = helperNRGridVisualizer(numFrameSimulation,gNB,UEs,CQIGridVisualization=cqiVisualization, ...
        ResourceGridVisualization=rbVisualization,SchedulingLogger=simSchedulingLogger);
end
numMetricsSteps = 20;
metricsVisualizer = helperNRMetricsVisualizer(gNB,UEs,NumMetricsSteps=numMetricsSteps, ...
    PlotSchedulerMetrics=true,PlotRLCMetrics=true);
simulationLogFile = "simulationLogs"; % For logging the simulation traces
% Calculate the simulation duration (in seconds)
simulationTime = numFrameSimulation*1e-2;
% Run the simulation
run(networkSimulator,simulationTime)
gNBStats = statistics(gNB);
ueStats = statistics(UEs);
displayPerformanceIndicators(metricsVisualizer)
if enableTraces
    simulationLogs = cell(1,1);
    if gNB.DuplexMode == "FDD"
        logInfo = struct("DLTimeStepLogs",[],"ULTimeStepLogs",[],"SchedulingAssignmentLogs",[],"RLCLogs",[]);
        [logInfo.DLTimeStepLogs,logInfo.ULTimeStepLogs] = getSchedulingLogs(simSchedulingLogger);
    else % TDD
        logInfo = struct("TimeStepLogs",[],"SchedulingAssignmentLogs",[],"RLCLogs",[]);
        logInfo.TimeStepLogs = getSchedulingLogs(simSchedulingLogger);
    end
    % Get the scheduling assignments log
    logInfo.SchedulingAssignmentLogs = getGrantLogs(simSchedulingLogger);
    % Get the RLC logs
    logInfo.RLCLogs = getRLCLogs(simRLCLogger);
    % Save simulation logs in a MAT-file
    simulationLogs{1} = logInfo;
    save(simulationLogFile,"simulationLogs")
end