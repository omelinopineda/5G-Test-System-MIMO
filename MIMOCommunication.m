% This promgration is part of a study in the 5G network article, supported by
% Research Tecchnology privated Company.
% //Mutual Coupling on MIMO Communication//*

%% System Parameters QPSK modulated Alamouti OSTBC is simulated over(2x2) quasi-static frequency-flat Rayleigh.
fc = 60e9;         % Center frequency
Nt = 2;             % Number of Tx antennas
Nr = 2;             % Number of Rx antennas
blkLen = 2;         % Alamouti code block length
snr = 0:10;         % SNR range
maxNumErrs = 70e2;   % Maximum number of errors
maxNumBits = 66e4;   % Maximum number of bits

qpskMod = comm.QPSKModulator;
qpskDemod = comm.QPSKDemodulator; 
alamoutiEnc = comm.OSTBCEncoder( ...
    'NumTransmitAntennas', Nt);
alamoutiDec = comm.OSTBCCombiner( ...
    'NumTransmitAntennas', Nt, ...
    'NumReceiveAntennas',  Nr);
awgnChanNC = comm.AWGNChannel( ... % For no coupling case
    'NoiseMethod', 'Signal to noise ratio (SNR)',...
    'SignalPower', 1);
berCalcNC = comm.ErrorRate;       % For no coupling case

% Clone objects for mutual coupling case 
awgnChanMC = clone(awgnChanNC); 
berCalcMC  = clone(berCalcNC);

% Antenna Arrays and Coupling Matrices
txSpacing = 0.5;
rxSpacing = 0.1;
lambda = physconst('lightspeed')/fc;
antElement = dipole( ...
    'Length', lambda/2, ...
    'Width',  lambda/100);
txArray = linearArray( ...
    'Element',        antElement,...
    'NumElements',    Nt,...
    'ElementSpacing', txSpacing*lambda);
rxArray = linearArray( ...
    'Element',        antElement,...
    'NumElements',    Nr,...
    'ElementSpacing', rxSpacing*lambda);
txMCMtx = helperCalculateCouplingMatrix(txArray, fc, [1 Nt]);
rxMCMtx = helperCalculateCouplingMatrix(rxArray, fc, [1 Nr]);

% Spatial Correlation Matrices
txCorrMtx = eye(2);
rxCorrMtx = [1 0.9; 0.9 1];
combCorrMtx = kron(txCorrMtx, rxCorrMtx);
txMCCorrMtx = txMCMtx * txCorrMtx * txMCMtx';
rxMCCorrMtx = rxMCMtx * rxCorrMtx * rxMCMtx';
txSqrtCorrMtx = txMCMtx * sqrtm(txCorrMtx);
rxSqrtCorrMtx = rxMCMtx * sqrtm(rxCorrMtx);
combMCCorrMtx = kron(txSqrtCorrMtx, rxSqrtCorrMtx);
combMCCorrMtx = combMCCorrMtx * combMCCorrMtx';

% MIMO Channel Modeling
mimoChanNC = comm.MIMOChannel( ...  % For no coupling case 
    'MaximumDopplerShift',             0, ...
    'SpatialCorrelationSpecification', 'Combined', ...
    'SpatialCorrelationMatrix',        combCorrMtx,...
    'PathGainsOutputPort',             true);

% Clone objects for mutual coupling case 
mimoChanMC = clone(mimoChanNC);
mimoChanMC.SpatialCorrelationMatrix = combMCCorrMtx;

% Simulations
% Set up a figure to visualize BER results
h1 = figure; grid on; hold on;
ax = gca;
ax.YScale = 'log';
xlim([snr(1), snr(end)]); ylim([1e-3 1]);
xlabel('SNR (dB)'); ylabel('BER'); 
h1.NumberTitle = 'off';
h1.Name = 'Orthogonal Space-Time Block Coding';
h1.Renderer = 'zbuffer';
title('Alamouti-coded 2x2 System - High Coupling, High Correlation');

s = rng(108);  % For repeatability
[berNC, berMC] = deal(zeros(3,length(snr)));

% Loop over SNR values
for idx = 1:length(snr)
    awgnChanNC.SNR = snr(idx); 
    awgnChanMC.SNR = snr(idx); 
    reset(berCalcNC); 
    reset(berCalcMC);    
    
    while min(berNC(2,idx),berMC(2,idx)) <= maxNumErrs && (berNC(3,idx) <= maxNumBits)    
        % Generate random data
        txData = randi([0 3], blkLen, 1);
        
        % Perform QPSK modulation and Alamouti encoding
        txSig = alamoutiEnc(qpskMod(txData)); 
        
        % Pass through MIMO channel
        reset(mimoChanNC); reset(mimoChanMC);
        [chanOutNC, estChanNC] = mimoChanNC(txSig);
        [chanOutMC, estChanMC] = mimoChanMC(txSig);
        
        % Add AWGN
        rxSigNC = awgnChanNC(chanOutNC);
        rxSigMC = awgnChanMC(chanOutMC);
        
        % Perform Alamouti decoding with known channel state information
        decSigNC = alamoutiDec(rxSigNC, squeeze(estChanNC));
        decSigMC = alamoutiDec(rxSigMC, squeeze(estChanMC));
                        
        % Perform QPSK demodulation 
        rxDataNC = qpskDemod(decSigNC);
        rxDataMC = qpskDemod(decSigMC);
        
        % Update BER
        berNC(:, idx) = berCalcNC(txData, rxDataNC);
        berMC(:, idx) = berCalcMC(txData, rxDataMC);
    end 

    % Plot results
    semilogy(snr(1:idx), berNC(1,1:idx), 'r*');
    semilogy(snr(1:idx), berMC(1,1:idx), 'bo');
    legend({'Channel Without Coupling', 'Channel With Coupling'});
    drawnow;
end

% Perform curve fitting
fitBERNC = berfit(snr, berNC(1,:));
fitBERMC = berfit(snr, berMC(1,:));
semilogy(snr, fitBERNC, 'r', snr, fitBERMC, 'b');
legend({'Channel Without Coupling', 'Channel With Coupling'});
rng(s); % Restore RNG


