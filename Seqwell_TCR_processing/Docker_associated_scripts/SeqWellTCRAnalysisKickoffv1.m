function SeqWellTCRAnalysisKickoffv1(BAMfile,species,BCFile,bFile,settingsFile, nCores)

fid = fopen(settingsFile,'r');
settings = textscan(fid,'%s%s','Delimiter','\t');
fclose(fid);
for L=1:length(settings{1})
    if strcmp(settings{1}{L},'BCcutoff')
        BCcutoff = str2double(settings{2}{L}); %%%%%Read threshold for barcode
    elseif strcmp(settings{1}{L},'UMIcutoff')
        UMIcutoff = str2double(settings{2}{L});
    elseif strcmp(settings{1}{L},'UMIlim')
        UMIlim = str2double(settings{2}{L});
    elseif strcmp(settings{1}{L},'VFreqCutoff')
        VFreqCutoff = str2double(settings{2}{L});
    elseif strcmp(settings{1}{L},'JFreqCutoff')
        JFreqCutoff = str2double(settings{2}{L});
    elseif strcmp(settings{1}{L},'TCRBin')
        JFreqCutoff = str2double(settings{2}{L});    
    end
end

[fPath,bName,ext]=fileparts(bFile);
fPath=[fPath '/'];



tic;
disp('Parsing Sequencing data')


if exist(BCFile,'file') %%%%Already defined list of barcodes to focus anaysis on
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


disp(BAMfile)
disp('Loading sequencing data...')



if exist([BAMfile '.bai'],'file')==0 %%%index bam if hasn't been done
    if exist([BAMfile(1:end-4) 'Sort.bam'],'file')
        BAMfile = [BAMfile(1:end-4) 'Sort.bam'];
        if exist([BAMfile '.bai'],'file')==0
            disp('Indexing bam file')
            message = ['samtools index ' BAMfile ];
            system(message);
        end
    else
        disp('Sorting and indexing bam file')
        sortBAM = [BAMfile(1:end-4) 'Sort.bam'];
        message = ['samtools sort -@ 12 -o ' sortBAM ' ' BAMfile ];
        system(message);
        %delete(BAMfile0);
        BAMfile = sortBAM;
        message = ['samtools index ' BAMfile ];
        system(message);
    end
end
[foo fName ext] = fileparts(BAMfile);


m=BioMap(BAMfile);

ref = settings{2}{all([cellfun(@isempty,strfind(lower(settings{1}),species))==0 cellfun(@isempty,strfind(lower(settings{1}),'ref'))==0],2)};
CD3baseFile = settings{2}{all([cellfun(@isempty,strfind(lower(settings{1}),species))==0 cellfun(@isempty,strfind(lower(settings{1}),'cdr3base'))==0],2)};

UTRCref='';
if strcmpi(species,'human')
    UTRCref=settings{2}{strcmp('humanUTRcRef',settings{1})};
end
fid = fopen(CD3baseFile,'r');
CD3base = textscan(fid,'%s%u','Delimiter','\t');
fclose(fid);

ind = find(filterByFlag(m,'unmappedQuery',0));
ind2 = [];


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


mkdir([fPath 'BCFASTQ']);
mkdir([fPath 'Alignments']);
mkdir([fPath 'Assemble']);
mkdir([fPath 'Clones/']);


mess=['Number of BC to process: ' num2str(nBC)];
disp(mess)

toc

for x = 1:nCores
    
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



toc
disp('Distributed TCR assembly...')
%%%%%Distribute TCR assembly across cluster

passFile = [fPath fName ext];
info.passfile = passFile;
info.refFile = ref;
info.UTRC = UTRCref;
info.BAMfile = BAMfile;
info.species=species;
info.nBC=nBC;
info.nCores = nCores;
info.m = m;
info.UMIcutoff=UMIcutoff;
info.UMIlim=UMIlim;
info.b=b;
info.VFreqCutoff=VFreqCutoff;
info.JFreqCutoff=JFreqCutoff;
info.CDR3base=CD3base;
info.forReads=forReads;
info.revReads=revReads;
info.nReads=nReads;
infoFile = [fPath fName 'info.mat'];
save(infoFile,'info');
out=['infoFile=' infoFile];
disp(out)
end
        



