wirelessnetworkSupportPackageCheck
rng("default") % Reset the random number generator
numFrameSimulation = 10; % Simulation time in terms of number of 10 ms frames
networkSimulator = wirelessNetworkSimulator.init;
gNB = nrGNB(CarrierFrequency=3.5e9, ...
    ChannelBandwidth=100e6, ...
    SubcarrierSpacing=60e3,...
    NumTransmitAntennas=32, ...
    NumReceiveAntennas=16, ...
    ReceiveGain=14);
uePositions = [7 0 0; 10 0 0; 12 0 0; 20 0 0; 40 0 0; 50 0 0; 52 0 0];
ueNames = "UE-" + (1:size(uePositions,1));
UEs = nrUE(Name=ueNames,Position=uePositions,NumTransmitAntennas=4,NumReceiveAntennas=2,ReceiveGain=11);
rlcBearerConfig = nrRLCBearerConfig(SNFieldLength=6,BucketSizeDuration=10);
connectUE(gNB,UEs,RLCBearerConfig=rlcBearerConfig)
appDataRate = 40e3; % Application data rate in kilo bits per second (kbps)
for ueIdx = 1:length(UEs)
    % Install DL application traffic on gNB for the UE
    dlApp = networkTrafficOnOff(GeneratePacket=true,OnTime=numFrameSimulation*10e-3,...
        OffTime=0,DataRate=appDataRate);
    addTrafficSource(gNB,dlApp,DestinationNode=UEs(ueIdx))

    % Install UL application traffic on UE for the gNB
    ulApp = networkTrafficOnOff(GeneratePacket=true,OnTime=numFrameSimulation*10e-3,...
        OffTime=0,DataRate=appDataRate);
    addTrafficSource(UEs(ueIdx),ulApp)
end
addNodes(networkSimulator,gNB)
addNodes(networkSimulator,UEs)
channelConfig = struct("DelayProfile","CDL-C","DelaySpread",300e-9);
channels = createCDLChannels(channelConfig,gNB,UEs);
customChannelModel = hNRCustomChannelModel(channels);
addChannelModel(networkSimulator,@customChannelModel.applyChannelModel)
enableTraces = true;
if enableTraces
    % Create an object for scheduler traces logging
    simSchedulingLogger = helperNRSchedulingLogger(numFrameSimulation,gNB,UEs);
    % Create an object for PHY traces logging
    simPhyLogger = helperNRPhyLogger(numFrameSimulation,gNB,UEs);
end
numMetricsSteps = 10;
metricsVisualizer = helperNRMetricsVisualizer(gNB,UEs,NumMetricsSteps=numMetricsSteps,...
    PlotSchedulerMetrics=true,PlotPhyMetrics=true);
simulationLogFile = "simulationLogs"; % For logging the simulation traces
% Calculate the simulation duration (in seconds)
simulationTime = numFrameSimulation * 1e-2;
% Run the simulation
run(networkSimulator,simulationTime);
gNBStats = statistics(gNB);
ueStats = statistics(UEs);
displayPerformanceIndicators(metricsVisualizer)
