
BCcutoff = 10; %%%%%Read threshold for barcode
UMIcutoff=10;
UMIlim = 1000;
VFreqCutoff = 0.6;
JFreqCutoff=0.6;
[bName,fPath,FilterIndex] = uigetfile('*.txt', 'TCR Assemble batch file')
bFile = [fPath bName];
fid = fopen(bFile,'r');
sampleList = textscan(fid,'%s%s','Delimiter','\t');
fclose(fid);
for samp=2:length(sampleList{1})
tic;
disp('Parsing Sequencing data')
clearvars -except sampleList samp BCcutoff UMIcutoff BClim fPath UMIlim VFreqCutoff JFreqCutoff
BAMfile0 = sampleList{1}{samp};
BCFile=sampleList{2}{samp};
if exist(BCFile,'file')
    bcLim=1;
fid=fopen(BCFile,'r');
BCList=textscan(fid,'%s');
BCList=BCList{1};
BCList=randsample(BCList,length(BCList));

%%%%Make 1 hamming variants
rChar = repmat(permute(reshape([BCList{:}],length(BCList{1}),length(BCList)),[2 1]),1,1,length(BCList{1})*4);
bases={'A','T','C','G'};
for x=1:size(rChar,3)
    p=1+floor((x-1)/4);
    b=bases{1+rem(x-1,4)};
    rChar(:,p,x)=b;
end
BCVar=num2cell(rChar,2);
BCVar=permute(BCVar,[3 1 2]);
BCVar=BCVar(:);
else
    bcLim=0;
end


[bPath bName bExt] = fileparts(BAMfile0);

disp(BAMfile0)
disp('Loading sequencing data...')
startdir = bPath; % directory from which Enumerator is started
startdrive = regexp(startdir, '([A-Z]:).*', 'tokens'); % drive from which it is started
if isempty(startdrive)==0
    [~, netuse] = system('net use'); % detailed map of network drives
    drivemap = regexp(netuse, '\n[OK]* *([A-Z]:) *([^ ]*)', 'tokens'); % actual map
    for i = 1:length(drivemap)
        drivemap{i}{2} = regexprep(drivemap{i}{2}, '\\', '\\\\'); % needed later on
        if strcmp(drivemap{i}{1},startdrive{1}{1})
            startdir = regexprep(startdir, startdrive{1}{1}, drivemap{i}{2});
        end
    end
    if isempty(strfind(startdir,'.local'))==0
        startdir = ['\\PHILEO' startdir(strfind(startdir,'.local')+6:end)];
    end
    BAMfile = [startdir '\' bName bExt];
end

startdir = fPath; % directory from which Enumerator is started
startdrive = regexp(startdir, '([A-Z]:).*', 'tokens'); % drive from which it is started
if isempty(startdrive)==0
    [~, netuse] = system('net use'); % detailed map of network drives
    drivemap = regexp(netuse, '\n[OK]* *([A-Z]:) *([^ ]*)', 'tokens'); % actual map
    for i = 1:length(drivemap)
        drivemap{i}{2} = regexprep(drivemap{i}{2}, '\\', '\\\\'); % needed later on
        if strcmp(drivemap{i}{1},startdrive{1}{1})
            startdir = regexprep(startdir, startdrive{1}{1}, drivemap{i}{2});
        end
    end
    if isempty(strfind(startdir,'.local'))==0
        startdir = ['\\PHILEO' startdir(strfind(startdir,'.local')+6:end)];
    end
    fPath = startdir;
end


startDir = pwd; % drive from which ClusterPrep is started; ProcessImages must be in the same network folder
startDrive = regexp(startDir, '([A-Z]:).*', 'tokens');
if ~isempty(startDrive)
    for i = 1:length(drivemap)
        if strcmp(drivemap{i}{1},startDrive{1}{1})
            startDirNet = regexprep(startDir, startDrive{1}{1}, drivemap{i}{2});
        end
    end
end
pd = startDirNet;

addpath([pd '\bin']);

                    


if exist([BAMfile '.bai'],'file')==0 %%%index bam if hasn't been done
    if exist([BAMfile(1:end-4) 'Sort.bam'],'file')
        BAMfile0 = [BAMfile0(1:end-4) 'Sort.bam'];
        BAMfile = [BAMfile(1:end-4) 'Sort.bam'];
        if exist([BAMfile '.bai'],'file')==0
            disp('Indexing bam file')
            message = ['C:\cygwin64\bin\bash -login -c "samtools index ""' BAMfile0 '"""'];
            system(message);
        end
    else   
    disp('Sorting and indexing bam file')
    sortBAM = [BAMfile0(1:end-4) 'Sort.bam'];
    message = ['C:\cygwin64\bin\bash -login -c "samtools sort -@ 12 -o ""' sortBAM '"" ""' BAMfile0 '"""'];
    system(message);
    %delete(BAMfile0);
    BAMfile0 = sortBAM;
    message = ['C:\cygwin64\bin\bash -login -c "samtools index ""' BAMfile0 '"""'];
    system(message);
    BAMfile = [BAMfile(1:end-4) 'Sort.bam'];
    end
end
[foo fName ext] = fileparts(BAMfile);

regions = {'TRBV','TRBJ 5UTR','TRBJ','TRBC','TRAV','TRAJ 5UTR','TRAJ','TRAC'};
m=BioMap(BAMfile);
message = ['C:\cygwin64\bin\bash -login -c "samtools view -H ""' BAMfile0 '"" | grep -o ''MOUSE\|mouse\|mm10\|mm9''"'];
[foo moM] = system(message);
message = ['C:\cygwin64\bin\bash -login -c "samtools view -H ""' BAMfile0 '"" | grep -o ''HUMAN\|human\|hg19\|hg38''"'];
 [foo huM] = system(message);
 message = ['C:\cygwin64\bin\bash -login -c "samtools view -H ""' BAMfile0 '"" | grep -o ''macfas''"'];
 [foo macfasM] = system(message);
   mo=0;
 hu=0;
 macfas=0;
 d=m.SequenceDictionary;
 d2=strfind(d,'TRBC');
 d3=cellfun(@isempty,d2)==0;
 if sum(d3)>0        %%%%Aligned directly to TCR segments
     if isempty(moM) ==0
         genome=1;
         ref=[pd '\bin\mouseTCR_Ctrunc_TRJ_UTR_1allele.fa'];
         disp('Mouse TCR Analysis');
         mo=1;
     elseif isempty(macfasM) ==0
         genome=3;
         ref=[pd '\bin\CynoTCR.fa'];
         disp('MacFas TCR Analysis');
         macfas=1;
     else
         
         genome=2;
         disp('Human TCR Analysis');
         ref=[pd '\bin\humanTCR_Ctrunc_TRJ_UTR_1allele.fa'];
         hu=1;
     end
     ind = find(filterByFlag(m,'unmappedQuery',0));
     ind2 = [];
 end

 if isempty(moM) ==0 && mo==0
     q9=strfind(moM,'mm9');
     
     if isempty(q9)==0
        mo=1;
     end
     q10=strfind(moM,'mm10');
     
     if  isempty(q10)==0
        mo=2;
     end
     if mo==0
         mo=input('Mouse mm9 (1), or mm10 (2) genome alignment?:');
     end
 elseif isempty(huM) ==0 && hu==0
     q9=strfind(huM,'hg19');
     
     if isempty(q9)==0
        hu=1;
     end
     q10=strfind(huM,'hg38');
     
     if isempty(q10)==0
        hu=2;
     end
    if hu==0
         hu=input('Human hg19 (1), or hg38 (2) genome alignment?:');
     end
 end   
 if sum([mo;hu;macfas])==0
     spec=input('Human hg19 (1), human hg38 (2), mouse mm9 (3) or mouse mm10 (4) dataset?:');
     if spec == 1
         hu=1;
     elseif spec ==2
         hu=2;
     elseif spec ==3
         mo = 1;
     elseif spec ==4
         mo=2;
     end
 end

dirInd = filterByFlag(m,'StrandQuery',0); %%%%Both TCRa and TCRb are on forward strand for human and mouse
secInd = filterByFlag(m,'readIsSecond',1);
altInd = filterByFlag(m,'alnNotPrimary',0);
nReads = m.NSeqs/2;
ind=[ind;ind2];

L=false(size(dirInd));
L(ind)=1;

indF = find(all([dirInd secInd altInd L],2));
forReads=length(indF);
revReads = sum(all([dirInd==0 L],2));

disp('Reads mapped to correct strand of TCR locus:')
disp(forReads);
disp('Reads mapped to opposite strand of TCR locus:')
disp(revReads);

toc

disp('Extracting barcodes')
TCRHeader=m.Header(indF);
indBC=find(secInd==0);
try
indAll = ismember(m.Header(indBC),TCRHeader);
catch
    h=m.Header(indBC);
    h2=cellfun(@isempty,h);
    h(h2)=repmat({'A'},sum(h2),1);
    indAll = ismember(h,TCRHeader);
end

indBC=indBC(indAll);

if bcLim==1
BC = m.Sequence(indBC);
% clear TCRHeader flags ind ind2 indAll indF dirInd L
BCLen = length(BC{1});
BCchar = permute(reshape([BC{:}],BCLen,length(BC)),[2 1]);
beadBCChar = BCchar(:,1:12);

beadBC = mat2cell(beadBCChar,repmat(1,size(beadBCChar,1),1),12);
[BCHit BCI] =ismember(beadBC,BCVar);
BCICol=1+floor((BCI-1)/(4*length(beadBC{1})));

count = histc(BCICol,1:length(BCList));
b = BCList(count>BCcutoff);

BInd=1:length(BCList);
BInd=BInd(count>BCcutoff);
nBC = length(b);

else
BC = m.Sequence(indBC);
% clear TCRHeader flags ind ind2 indAll indF dirInd L
BCLen = length(BC{1});
BCchar = permute(reshape([BC{:}],BCLen,length(BC)),[2 1]);
beadBCChar = BCchar(:,1:12);

beadBC = mat2cell(beadBCChar,repmat(1,size(beadBCChar,1),1),12);


[beadList BCICol] = umiCollapseBC(beadBC);
count = histc(BCICol,1:max(BCICol));
b = beadList(count>BCcutoff);
[b ord]=sort(b);
BInd=1:length(beadList);
BInd=BInd(count>BCcutoff);
BInd=BInd(ord);
nBC = length(b);
end


mkdir([fPath 'BCFastQ']);
mkdir([fPath 'Alignments']);
mkdir([fPath 'Assemble']);
mkdir([fPath 'Clones/']);


mess=['Number of BC to process: ' num2str(nBC)];
disp(mess)
[~, coreInfo] = system('node listcores');
coreNames = regexp(coreInfo, '\n([^-][^ ]*) - [0-9]* *([^-][^ ]*)', 'tokens');
nCores = 0;
for iCore = 1:length(coreNames)
    if ~strcmp(coreNames{iCore}{2},'Offline')
        nCores = nCores + 1;
    end
end
nCores = ceil(length(b)/ceil(length(b)/nCores));

toc
h = waitbar(0,'Distributing barcode reads..');
for x = 1:nCores
    waitbar(x/nCores)
    indexCoreFile = [fPath fName 'CoreIndex' num2str(x) '.txt'];
    fid = fopen(indexCoreFile,'w');
    if x<nCores
        BCLoc = 1+(x-1)*ceil(nBC/nCores):(x)*ceil(nBC/nCores);
    else
        BCLoc = 1+(x-1)*ceil(nBC/nCores):nBC;
    end
    for y = 1:length(BCLoc)
        bamIndex = indBC(BCICol==BInd(BCLoc(y)));
        out = [repmat(BCLoc(y),1,length(bamIndex));bamIndex'];
        fprintf(fid,'%u\t%u\n',out);
        
        
    end
    fclose(fid);
end
close(h);


toc
disp('Distributed TCR assembly...')
%%%%%Distribute TCR assembly across cluster
jobname = 'TCRAssembly';
exeName = 'matlab.exe';
passFile = [fPath fName ext];
info.passfile = passFile;
info.refFile = ref;
info.BAMfile = BAMfile;
info.g = genome;
info.nBC=nBC;
info.nCores = nCores;
info.m = m;
info.pd=pd;
info.UMIcutoff=UMIcutoff;
info.UMIlim=UMIlim;
info.b=b;
info.VFreqCutoff=VFreqCutoff;
info.JFreqCutoff=JFreqCutoff;
infoFile = [fPath fName 'info.mat'];
save(infoFile,'info');

[~, job] = system(['job new /jobname:',jobname,' /nodegroup:ComputeNodes']); % create a new job from default template
job_id = regexprep(job,'.* ([0-9]*)\n','$1'); % get the unique job id


exeFile = [pd '\bin\TCRmapv12_4.exe'];
%  str = ['job add ', job_id, ' /parametric:1:' num2str(nCores) ':1 ', 'start \"\" /b /LOW /WAIT "' exeFile '" "' passFile '" "' FASTQfile '" ' num2str(genome) ' ' num2str(nBC) ' ' num2str(nCores) ' *'];
 str = ['job add ', job_id, ' /parametric:1:' num2str(nCores) ':1 ', 'start \"\" /b /LOW /WAIT "' exeFile '" \"' infoFile '\" *'];

system(str); % add the task to the job;

system(strcat('job submit /id:', job_id)); % finally, submit the job

% now, wait for the cluster to finish it
finished = false;
while ~finished
    [~, jobview] = system(['job view ', job_id]); % get detailed info about the job
    progress = str2double(regexprep(jobview, '.*Progress[ ]*: ([0-9]*)%.*', '$1')); % percentage completed
    if progress == 100
        finished = true;
    end
    pause(1);
end
fclose('all');
toc
disp('Finalizing results...')
%%%%%Collate results

cloneDir = [fPath 'Clones/'];
summaryFile = [fPath fName 'MappingSummary.txt'];
fidM = fopen(summaryFile,'wt');
fprintf(fidM,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','UMI_Name', 'BC', ...
    'UMI','BCNum','UMINum','nReads','TopVregion','TopJregion','TopJ_5UTR','CDR3','CDR3nuc','TopVReadCount','TopVFreq','TopJReadCount', ... 
    'TopJFreq','TopJ_5UTRCount','CDR3ReadCount','CDR3Freq',regions{1},regions{2},regions{3},regions{4},regions{5},regions{6},regions{7},regions{8});
fclose(fidM);
TCRFinal = [];


alignDir = [fPath 'Alignments\'];
for x=1:nCores
    alignResults = [alignDir fName 'AlignmentResults' num2str(x) '.txt'];
    depthResults = [alignDir fName 'DepthResults' num2str(x) '.mat'];
    system(['type "' alignResults '" >> "' summaryFile '"']);
    load(depthResults)
    if x==1
        depACDR3F=depACDR3;
        depAConF=depACon;
        depATruncF=depATrunc;
        depBCDR3F=depBCDR3;
        depBConF=depBCon;
        depBTruncF=depBTrunc;
        hTRBC1ACDR3F=hTRBC1ACDR3;
        hTRBC2ACDR3F=hTRBC2ACDR3;
        hTRBC1AConF=hTRBC1ACon;
        hTRBC2AConF=hTRBC2ACon;
        hTRBC1ATrunF=hTRBC1ATrun;
        hTRBC2ATrunF=hTRBC2ATrun;
        concensusSeq=C;
    else
        depACDR3F=[depACDR3F depACDR3];
        depAConF=[depAConF depACon];
        depATruncF=[depATruncF depATrunc];
        depBCDR3F=[depBCDR3F depBCDR3];
        depBConF=[depBConF depBCon];
        depBTruncF=[depBTruncF depBTrunc];
        hTRBC1ACDR3F=[hTRBC1ACDR3F hTRBC1ACDR3];
        hTRBC2ACDR3F=[hTRBC2ACDR3F hTRBC2ACDR3];
        hTRBC1AConF=[hTRBC1AConF hTRBC1ACon];
        hTRBC2AConF=[hTRBC2AConF hTRBC2ACon];
        hTRBC1ATrunF=[hTRBC1ATrunF hTRBC1ATrun];
        hTRBC2ATrunF=[hTRBC2ATrunF hTRBC2ATrun];
        concensusSeq=[concensusSeq;C];
    end

        
        
end
depACDR3F=depACDR3F(:,sum(depACDR3F,1)>0);
depAConF=depAConF(:,sum(depAConF,1)>0);
depATruncF=depATruncF(:,sum(depATruncF,1)>0);
depBCDR3F=depBCDR3F(:,sum(depBCDR3F,1)>0);
depBConF=depBConF(:,sum(depBConF,1)>0);
depBTruncF=depBTruncF(:,sum(depBTruncF,1)>0);

hTRBC1ACDR3F=hTRBC1ACDR3F(:,sum(hTRBC1ACDR3F,1)>0);
hTRBC2ACDR3F=hTRBC2ACDR3F(:,sum(hTRBC2ACDR3F,1)>0);
hTRBC2AConF=hTRBC2AConF(:,sum(hTRBC2AConF,1)>0);
hTRBC1AConF=hTRBC1AConF(:,sum(hTRBC1AConF,1)>0);
hTRBC2ATrunF=hTRBC2ATrunF(:,sum(hTRBC2ATrunF,1)>0);
hTRBC1ATrunF=hTRBC1ATrunF(:,sum(hTRBC1ATrunF,1)>0);

fid=fopen(summaryFile,'r');
dHeaders=textscan(fid,repmat('%s',1,26),1);
d=textscan(fid,'%s%s%s%u%u%u%s%s%s%s%s%u%.2f%u%.2f%u%u%.2f%u%u%u%u%u%u%u%u','Delimiter','\t');
fclose(fid);
TCRFinal=d;
[uniqBC L J] = unique(TCRFinal{2});
TCRCount = [];
CDR3 = [];
BCCount = zeros(length(uniqBC),1);
for x=1:length(uniqBC)
    BCHit = TCRFinal{10}(J==x);
    BCHit=BCHit(cellfun(@isempty,BCHit)==0);
    [uniq a bb] = unique(BCHit);
    cLoc = histc(bb,(1:length(uniq)))';
    TCRCount = [TCRCount;cLoc'];
    CDR3= [CDR3;uniq];
    BCCount(x) = length(uniq);
end
[uCDR3 foo uCI]=unique(CDR3);
CDRc=histcounts(uCI,1:length(uCDR3)+1);
CDR3L=cellfun(@length,uCDR3);

BCSumFile = [fPath fName 'TCRperBC.txt'];
fid=fopen(BCSumFile,'wt');
fprintf(fid,'%s\t%s\n','Barcode','nTCR');
for x=1:length(BCCount)
 fprintf(fid,'%s\t%u\n',uniqBC{x},BCCount(x)); 
end
fclose(fid);

TCRcFile = [fPath fName 'TCRCounts.txt'];
fid=fopen(TCRcFile,'wt');
fprintf(fid,'%s\t%s\n','TCR','nCells');
for x=1:length(CDRc)
 fprintf(fid,'%s\t%u\n',uCDR3{x},CDRc(x)); 
end
fclose(fid);

%%%%%Summary figures
if exist([fPath 'SummaryPlots'],'dir')==0
    mkdir([fPath 'SummaryPlots']);
end
summaryFolder = [fPath 'SummaryPlots\' fName '\'];
if exist(summaryFolder,'dir')==0
    mkdir(summaryFolder)
end


%%%%%Cells with same CDR3
count = histc(CDRc,1:max(CDRc));
f=figure;
bar(count,.5,'FaceColor','r','EdgeColor','k','LineWidth',1.5);
ylabel('Count')
xlabel('# Cells/CDR3')
saveas(f,[summaryFolder fName '_CellsperCDR3.png']);

%%%%%TCRperBC
count = histc(BCCount,0:max(BCCount));
f=figure;
bar(0:max(BCCount),count,.5,'FaceColor','r','EdgeColor','k','LineWidth',1.5);
ylabel('Count')
xlabel('# TCR/Barcode')
saveas(f,[summaryFolder fName '_TCRPerBC.png']);

%%%%CDR3 Length
count = histc(CDR3L,min(CDR3L):max(CDR3L));
f=figure;
bar(min(CDR3L):max(CDR3L),count,.5,'FaceColor','r','EdgeColor','k','LineWidth',1.5);
ylabel('Counts')
xlabel('CDR3 length')
saveas(f,[summaryFolder fName '_CDR3Length.png']);


%%%%ReadsPerUMI

f=figure;
bar(TCRFinal{6},'FaceColor','r','Horizontal','on');
set(gca,'xscale','log');
ylabel('UMI')
xlabel('reads')
saveas(f,[summaryFolder fName '_ReadsperUMI.png']);

%%%%UMIPerTCR
count = histc(TCRCount,1:max(TCRCount));
f=figure;
bar(1:max(TCRCount),count,.5,'FaceColor','r','EdgeColor','k','LineWidth',1.5);
ylabel('CDR3')
xlabel('UMI')
saveas(f,[summaryFolder fName '_UMIperTCR.png']);

%%%%ReadsMapping
trunReads = TCRFinal{6};
trunReads(trunReads>1000)=1000; %%%%Only analyzed 1000 reads of any UMI
tReads = sum(trunReads);
fracMapping = [sum(TCRFinal{18})/tReads sum(TCRFinal{19})/tReads sum(TCRFinal{20})/tReads sum(TCRFinal{21})/tReads sum(TCRFinal{22})/tReads sum(TCRFinal{23})/tReads sum(TCRFinal{24})/tReads sum(TCRFinal{25})/tReads ]*100;
f=figure;
bar(fracMapping,'FaceColor','r','EdgeColor','k','LineWidth',1.5)   
set(gca,'xticklabel',regions)
ylabel('% of Reads')
xlabel('TCR Region') 
saveas(f,[summaryFolder fName '_TCRMappingFrequencies.png']);     
         
%%%%TCRMapping 
f=figure;
bar([forReads/nReads*100;revReads/nReads*100],'FaceColor','r','EdgeColor','k','LineWidth',1.5)   
set(gca,'xticklabel',{'TCR PS/Total','TCR NS/Total'})
ylabel('% of Total or Mapped Reads')
xlabel('TCR Locus Positive or Negative strand mapping') 
saveas(f,[summaryFolder fName '_TCRLocusMapping.png']);  

% %%%ConstantRegionMapping
% f=figure;
% bar(mean(depACDR3F(end-4830:end,:),2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_ConsMappingACDR3UMI.png']); 
% f=figure;
% bar(mean(depBCDR3F,2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_ConsMappingBCDR3UMI.png']);
% 
% f=figure;
% bar(mean(depAConF(end-4830:end,:),2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_ConsMappingANoCDR3UMI.png']);
% f=figure;
% bar(mean(depBConF,2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_ConsMappingBNoCDR3UMI.png']);
% f=figure;
% bar(mean(depATruncF(end-4830:end,:),2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_ConsMappingATruncUMI.png']);
% f=figure;
% bar(mean(depBTruncF,2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_ConsMappingBTruncUMI.png']);
% 
% %%%%VariableRegionMapping
% 
% f=figure;
% 
% bar(mean(depAConF(1:end-4830,:),2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingANoCDR3UMI.png']);
% 
% f=figure;
% 
% bar(mean(depATruncF(1:end-4830,:),2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingATrunUMI.png']);
% 
% f=figure;
% 
% bar(mean(depACDR3F(1:end-4830,:),2),'FaceColor','r','EdgeColor','r')   
% 
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingACDR3UMI.png']);
% 
% 
% 
% f=figure;
% bar(sum(hTRBC1AConF,2),'FaceColor','r','EdgeColor','r')   
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingB1NoCDR3UMI.png']);
% 
% f=figure;
% bar(sum(hTRBC1ATrunF,2),'FaceColor','r','EdgeColor','r')   
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingB1TrunUMI.png']);
% 
% f=figure;
% bar(sum(hTRBC1ACDR3F,2),'FaceColor','r','EdgeColor','r')   
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingB1CDR3UMI.png']);
% 
% f=figure;
% bar(sum(hTRBC2AConF,2),'FaceColor','r','EdgeColor','r')   
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingB2NoCDR3UMI.png']);
% 
% f=figure;
% bar(sum(hTRBC2ATrunF,2),'FaceColor','r','EdgeColor','r')   
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingB2TrunUMI.png']);
% 
% f=figure;
% bar(sum(hTRBC2ACDR3F,2),'FaceColor','r','EdgeColor','r')   
% ylabel('Frac. of Reads')
% xlabel('Constant Region') 
% saveas(f,[summaryFolder fName '_UpStreamMappingB2CDR3UMI.png']);

end
        
            
     

                



