

function SeqWellTCRAnalysisAggregatev1(infoFile)
disp('Finalizing results...')
%%%%%Collate results
load(infoFile);

[fPath,fName,~]=fileparts(info.passfile);
fPath=[fPath '/'];
regions = {'TRBV','TRBJ 5UTR','TRBJ','TRBC','TRAV','TRAJ 5UTR','TRAJ','TRAC'};
summaryFile = [fPath fName 'MappingSummary.txt'];
fidM = fopen(summaryFile,'wt');
fprintf(fidM,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','UMI_Name', 'BC', ...
    'UMI','BCNum','UMINum','nReads','TopVregion','TopJregion','TopJ_5UTR','CDR3','CDR3nuc','TopVReadCount','TopVFreq','TopJReadCount', ... 
    'TopJFreq','TopJ_5UTRCount','CDR3ReadCount','CDR3Freq',regions{1},regions{2},regions{3},regions{4},regions{5},regions{6},regions{7},regions{8});
fclose(fidM);



alignDir = [fPath 'Alignments/'];
for x=1:info.nCores
    alignResults = [alignDir fName 'AlignmentResults' num2str(x) '.txt'];
    
    system(['cat "' alignResults '" >> "' summaryFile '"']);
        
end

fid=fopen(summaryFile,'r');

d=textscan(fid,'%s%s%s%u%u%u%s%s%s%s%s%u%.2f%u%.2f%u%u%.2f%u%u%u%u%u%u%u%u','Delimiter','\t','headerLines', 1);
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
summaryFolder = [fPath 'SummaryPlots/' fName '/'];
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
bar([info.forReads/info.nReads*100;info.revReads/info.nReads*100],'FaceColor','r','EdgeColor','k','LineWidth',1.5)   
set(gca,'xticklabel',{'TCR PS/Total','TCR NS/Total'})
ylabel('% of Total or Mapped Reads')
xlabel('TCR Locus Positive or Negative strand mapping') 
disp('pre-Done')
saveas(f,[summaryFolder fName '_TCRLocusMapping.png']);  
disp('Done')
quit force


end
        
            
     

                



