For analyisis Matlab R2018a or R2021a was used.
Installing Matlab 2021a: download the installer for your platform and MATLAB release from MathWorks Downloads, installation takes 30-45 minutes.
System requirements: 
https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2021a-windows.pdf
https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2021a-macintosh.pdf

*************************************************************************************************************************************
What does the sData file contain?

Information about the session:
sData.sessionInfo

Information about the animal, surgery details:
sData.mouseInfo

Raw behavioral data:
sData.daqdata: recorded by LabView
The sampling rate is sData.daqdata.meta.fs = 3000, I downsampled it to the imaging sampling rate, which is sData.imdata.meta.fps = 30.9406. 
sData.daqdata.frameSignal and sData.daqdata.frameIndex shows in which sample the imaging scan started. Using this we could align behavioral and imaging data samples.

Optogenetic-stimulation protocol: 
sData.stimProtocols: 
0: control trials = no opto-stim, only masking light (20 Hz sinusoid red LED flickering directed to the eyes)
1: opto-stim protocol (the intensity: voltage Range Percentage has to be converted to mW/mm2, sData does not contain this information.)
2: in some sessions more than one intensity of opto-stimulation was tested, labelled by differnt trial types (2,3,...). In most of the sessions intensities were pooled together for analysis
	If more intensities were tested, details can be found in sData.behavior.optoMoreProts struct, where sData.behavior.optoMoreProts.OptoStimProtTrials shows which trial belongs to which type (1: control, 2: intensity #1, 3: intensity #2, ...)
	In the sData.behavior.opto struct, all opto-stim trials are pooled together, even if there were more than one intensities tested.
	Related imaging data for more intensities can be found in sData.imdata.binned.optoMoreProts struct.
 
Downsampled behavioral and opto-stimulation data:
sData.behavior: contain information about the animal position on the wheel, animal speed, if the animal licks, when the reward was given, when the optogenetic stimulation was applied.
sData.behavior.binning: I spatially binned the behavioral data (speed, licks/cm), rows: trials, columns: bins. I used sData.behavior.meta.nBins=80 spatial bins, bin size was sData.behavior.meta.binSize=1.9635 cm. 
sData.behavior.opto: optogenetic stimulation parameters, trials are sotred based on the opto-stimulation (off, on, after). Binned lick and speed data also sorted here. 
sData.behavior.opto.OptoOnSignalDS: value is one if the stimulation is at its peak (for Chrimson it is sinusoid waveform, for Arch it is a square pulse). 
sData.behavior.performance and sData.behavior.performance2lick: I calculated different behavioral descriptors from the lick and speed profile to check if there is change in the behavior during optogenetic stimulation.

Imaging data: 
Calcium signals were recorded with SciScan and ROIs were detected and signal was extracted using NANSEN (Hennestad et al., submitted).
sData.imdata.roiSignals
sData.imdata.roiSignals(2).roif: raw Ca signals for each RIOs (rows) and each frame (columns)
sData.imdata.roiSignals(2).npilf: raw neurophil signals around each RIOs (rows) and each frame (columns)
sData.imdata.roiSignals(2).dff: dF/F signals for each RIOs (rows) and each frame (columns) calculated by NANSEN (neurophil is subtracted)
sData.imdata.roiSignals(2).deconv: deconvolved dF/F signals for each RIOs (rows) and each frame (columns) calculated by NANSEN (using CalmAn)
(For all analysis the sData.imdata.roiSignals(2).dff was used.) 

sData.imdata also contain information about recording parameters (sData.imdata.meta), signal extraction options (sData.imdata.signalExtractionOptions),
ROI statistics (sData.imdata.roiStat), starting position of cues on the wheel (sData.imdata.cues).
binned Ca data (sData.imdata.binned) and sorted based on optogenetic stimulation.

Place cell detection based on Mao et al., 2017:
sData.imdata.MaoPC_Opto_dff struct contain parameters for the detection and identity (ROI ID) of cells detected as place cells on different stimulation conditions (off, on, after).
Detected one peak place cells' ROI ID can be found here (position tuned cells having multiple peaks were excluded from analysis): sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3
Detected landmark/cue cells' ROI ID can be found here (excluded from analysis): sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF 

************************************************************************************************************************************
TEST SESSIONS
m8059-20200623-00 is a session recorded in RSC, VIP cells express ChrimsonR
m8060-20200622-00 is a session recorded in RSC, VIP cells express ArchT

************************************************************************************************************************************
BASIC ANALYSIS (all codes were written by Nora Lenkey)

Codes used to get sData.behavior and sData.imaging from raw data:
use '..._sData_raw' file (open in Matlab)
expected output should be the same as '-sData_analyzed' file
Expected run time for a 'normal' desktop computer: each scripts takes a few minutes (1-10)

Before running the codes set a save path in the sData file, here: sData.sessionInfo.savePath
The default is 'C:\MATLAB\SAVE\sessionID'

run these scripts (scripts has to be in the Matlab Path):

%% calculating behavioral data from raw LabView data
%Set parameters:
nBins = 80; % number of bins
IsOpto = 1; % was it an optical stimulated session? 0:no, 1:yes
OptoStimLimitMs = 15000; % in ms, if the opto-stimulation longer, labelled this trial as 'failed stimulation'. (Trial is discarded from analysis.)
DiscardCmBeginningCm = 10; % discard data in the beginnig of the track.
OptoSensitivity = 5; % sensitivity to detect opto-stimulation. Setting to '5' for full trial sinus and square stimulation works fine. 
sData = CalcBehav2(sData,nBins,IsOpto,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm); % basic analysis of behavioral data from raw LabView data   
LapsTested = 20;
sData = BehavFirstSecondHalf(sData,LapsTested); % compare behavior in the first and last 20 laps
close all

%% first processing of calcium data from raw data
%Set parameters:
VelMin = 0.1; % set a velocity limit (for pyramidal cells I used 0.1 cm/s). Below this value the animal considered 'sitting still', and data were discarded.
IsDeconv = 1; % 0: no deconvolved data, 1: use NANSEN deconvolution
gol = 5; % type of the task. I used for all experiments a goal oriented task (No #5), cues were at fixed positions.
sData = calcCaData2(sData,VelMin,IsDeconv,gol); % basic analysis of imaging data from raw data
close all

%% sorting and analysing Calcium data based on the trial type ( if it was optically stimulated or not: opto-off, opto-on, after-opto)
FigVisible = 'off'; % set the figure visibility 'off' if your computer tends to slow down during the script (plotting too many figures)
sData = CaDataOptoSortingdFFsData2(sData,FigVisible); 

%% place cell detection separately in opto-off and opto-on conditions (code is based on Mao et al., 2017)
datatype = 0; % datatype = 0 use dff data, datatype = 1 use deconvolved data (I used always dF/F for the analysis of this paper)
FigVisible = 'off';
sData = placecellMaosDataOptoOn2(sData,datatype,FigVisible); % place cell detection
sData = LandmarkCellDetection3(sData); % landmark cell detection; one-peak place cell detection (among place cells)

%% compare the same place cells in different types of protocol within one session (opto-on, opto-off)
datatype = 0; % datatype = 0 use dff , datatype = 1 use deconv
GaussFilter = 5;
sData = placecellComparison4_1stim(sData,datatype,GaussFilter); 

************************************************************************************************************************************
FURTHER ANALYSIS 
(Don't forget to set a folder for figure saving in the beginning of scripts. Default is = 'C:\MATLAB\SAVE\...')

sData = effectOnPCCtrOnePF3(sData); % checks average postition tuning curve changes for each place ROI during optical stimulation 
sData = EffectPCinCtr3(sData); % calculations if effect was significant on cell level, regression line calculation for each place cell (opto-on vs opto-off)
sData = gainModulationInOutRatioPCinCtr(sData); % calcuting effect of VIP modulation inside the place field (In) and outside the place field (Out)
