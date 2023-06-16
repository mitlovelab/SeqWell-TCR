

       function TCRmapv12_3(infoFile,iSlice)
tic
load(infoFile);
warning('off','bioinfo:saminfo:InvalidTagField')
% addpath([info.pd '\bin']);
cutoff = info.UMIcutoff;        %%%%%Number of reads required to consider UMI reads
[fPath,fName,ext] = fileparts(info.passfile);
if strcmp(fPath(end),'\')==0
    fPath = [fPath '\'];
end
regions = {'TRBV','TRBJ','TRBC','TRAV','TRAJ','TRAC'};
g=info.g;
nBC = info.nBC;
nCores = info.nCores;
vCut=info.VFreqCutoff;
jCut=info.JFreqCutoff;
%%reference Sequences
refFile = info.refFile;
fid=fopen(refFile,'r');
ref=textscan(fid,'%s');
fclose(fid);
ref=ref{1};
UMIlim=info.UMIlim;
%%%%Parse reference sequences
gR = strfind(ref,'>');
gL = cellfun(@isempty, gR)==0;
gI = find(gL);  %%Each alignment
TCRsequences=cell(length(gI),1);
TCRRegions = cell(length(gI),1);
c=1;
cU=1;
for x=1:length(gI)-1
    i=ref(gI(x)+1:gI(x+1)-1);
    
    [s e]=regexp(ref{gI(x)},'TR[ABGD][VDJC][0-9]*[A-Z]*-*[0-9]*');
    [sU eU]=regexp(ref{gI(x)},'UTR');
    if isempty(sU)
        TCRsequences{c} = [i{:}];
        TCRRegions{c}=ref{gI(x)}(s:e);
        c=c+1;
    else
        UTRRegions{cU}=ref{gI(x)}(s:eU);
        cU=cU+1;
    end
end

x=length(gI);
i=ref(gI(x)+1:end);

[s e]=regexp(ref{gI(x)},'TR[ABGD][VDJC][0-9]*[A-Z]*-*[0-9]*');
[sU eU]=regexp(ref{gI(x)},'UTR');
if isempty(sU)
    TCRsequences{c} = [i{:}];
    TCRRegions{c}=ref{gI(x)}(s:e);
    c=c+1;
else
    UTRRegions{cU}=ref{gI(x)}(s:eU);
    cU=cU+1;
end


rem=cellfun(@isempty,TCRRegions)==0;
TCRRegions=TCRRegions(rem);
TCRsequences=TCRsequences(rem);
[TCRRegions mm foo]=unique(TCRRegions);
TCRsequences = TCRsequences(mm);

%%%Determine start and finish of CDR3 region in each V and J region
%%%sequence
CDR3base=zeros(size(TCRsequences));
for x=[1,4] %%%%find CDR3 start base number in V regions
    we = strfind(TCRRegions,regions{x});
    ve = find(cellfun(@isempty, we)==0);%%V sequences
    for z=1:length(ve)
        s=TCRsequences{ve(z)};
        aa=nt2aa(s);
        C=strfind(aa,'C');
        
        if isempty(C) || ((strcmp(TCRRegions{ve(z)},'TRBV17') || strcmp(TCRRegions{ve(z)},'TRBV26')) && g==2)  || (strcmp(TCRRegions{ve(z)},'TRAV33-1') || strcmp(TCRRegions{ve(z)},'TRAV37-1') && g==3)%%%Do not have C in CDR3
           if strcmp(TCRRegions{ve(z)},'TRAV33-1') && g==3 
               CDR3base(ve(z))=274;
           elseif strcmp(TCRRegions{ve(z)},'TRAV37-1') && g==3 
               CDR3base(ve(z))=253;
           else
            C=strfind(aa,'SS');
            CDR3base(ve(z))=C(end)*3-2;
           end
            
        else
            CDR3base(ve(z))=C(end)*3-2;
        end
    end
end

for x=[2,5] %%%%find CDR3 end base number in J regions
    we = strfind(TCRRegions,regions{x});
    ve = find(cellfun(@isempty, we)==0);%%V sequences
    for z=1:length(ve)
        s=TCRsequences{ve(z)};
        for zz=1:3
            aa=nt2aa(s(zz:end));
            FG=strfind(aa,'FG');
            if isempty(FG)==0
                CDR3base(ve(z))=FG(end)*3+zz-1;
                break
            end
        end
        if CDR3base(ve(z))==0
            for zz=1:3
                aa=nt2aa(s(zz:end));
                if isempty(strfind(aa,'*'))
                    KG=strfind(aa,'KG');
                    if isempty(KG)==0
                        CDR3base(ve(z))=KG(1)*3+zz-1;
                        break
                    end
                    WG=strfind(aa,'WG');
                    if isempty(WG)==0
                        CDR3base(ve(z))=WG(1)*3+zz-1;
                        break
                    end
                    FS=strfind(aa,'FS');
                    if isempty(FS)==0
                        CDR3base(ve(z))=FS(1)*3+zz-1;
                        break
                    end
                    LG=strfind(aa,'LG');
                    if isempty(LG)==0
                        CDR3base(ve(z))=LG(1)*3+zz-1;
                        break
                    end
                    WS=strfind(aa,'WS');
                    if isempty(WS)==0
                        CDR3base(ve(z))=WS(1)*3+zz-1;
                        break
                    end
                    TG=strfind(aa,'TG');
                    if isempty(TG)==0
                        CDR3base(ve(z))=TG(1)*3+zz-1;
                        break
                    end
                    FA=strfind(aa,'FA');
                    if isempty(FA)==0
                        CDR3base(ve(z))=FA(1)*3+zz-1;
                        break
                    end
                    WQ=strfind(aa,'WQ');
                    if isempty(WQ)==0
                        CDR3base(ve(z))=WQ(1)*3+zz-1;
                        break
                    end
                    CDR3base(ve(z))=30;
                    
                end
            end
        end
    end
end



m=BioMap(info.BAMfile);

g
if ischar(iSlice)
    iSlice = str2double(iSlice);
end

if iSlice<nCores
    BCInd = 1+(iSlice-1)*ceil(nBC/nCores):(iSlice)*ceil(nBC/nCores);
else
    BCInd = 1+(iSlice-1)*ceil(nBC/nCores):nBC;
end
indexCoreFile = [fPath fName 'CoreIndex' num2str(iSlice) '.txt'];
fid = fopen(indexCoreFile,'r');
i=textscan(fid,'%u%u');
BCNum = i{1};
bamInd = i{2};
fclose(fid);
delete(indexCoreFile);


alignDir = [fPath 'Alignments\'];
assemDir=[fPath 'Assemble\'];
cloneDir = [fPath 'Clones\'];
if g == 1
    TCRBLASTIndex = [info.pd '\bin\mouseTCR_Ctrunc_TRJ_UTR_1allele.fa'];
    chrPosB=142498525:142500432;
    chrA='chr14';
    depBTrunc=zeros(length(chrPosB),1);
    depBCDR3=zeros(length(chrPosB),1);
    depBCon=zeros(length(chrPosB),1);
    chrPosA=22090057:23021075;
    depATrunc=zeros(length(chrPosA),1);
    depACDR3=zeros(length(chrPosA),1);
    depACon=zeros(length(chrPosA),1);
    chrB='chr6';
elseif g == 2
    TCRBLASTIndex = [info.pd '\bin\humanTCR_Ctrunc_TRJ_UTR_1allele.fa'];
    TRBCupstBLASTIndex = [info.pd '\bin\upstream_hTRBC.fa'];
    chrA='chr14';
    chrB='chr7';
    chrPosB=142498525:142500432;
    depBTrunc=zeros(length(chrPosB),1);
    depBCDR3=zeros(length(chrPosB),1);
    depBCon=zeros(length(chrPosB),1);
    chrPosA=22090057:23021075;
    depATrunc=zeros(length(chrPosA),1);
    depACDR3=zeros(length(chrPosA),1);
    depACon=zeros(length(chrPosA),1);
elseif g==3
    TCRBLASTIndex = [info.pd '\bin\CynoTCR.fa'];
    TRBCupstBLASTIndex = [info.pd '\bin\upstream_hTRBC.fa'];
    chrA='chr14';
    chrB='chr7';
    chrPosB=142498525:142500432;
    depBTrunc=zeros(length(chrPosB),1);
    depBCDR3=zeros(length(chrPosB),1);
    depBCon=zeros(length(chrPosB),1);
    chrPosA=22090057:23021075;
    depATrunc=zeros(length(chrPosA),1);
    depACDR3=zeros(length(chrPosA),1);
    depACon=zeros(length(chrPosA),1);
end
alignResults = [alignDir fName 'AlignmentResults' num2str(iSlice) '.txt'];
fidM = fopen(alignResults,'wt');
hTRBC1ACDR3=zeros(52250,1);
hTRBC1ACon=zeros(52250,1);
hTRBC1ATrun=zeros(52250,1);
hTRBC2ACDR3=zeros(7891,1);
hTRBC2ACon=zeros(7891,1);
hTRBC2ATrun=zeros(7891,1);
BCList=info.b;
C=cell(length(BCInd),100);
a=0;
for x=BCInd
    disp(x)
    a=a+1;
    locBamInd = bamInd(BCNum==x);
    BC = m.Sequence(locBamInd);
    BCLen = length(BC{1});
    BCChar = permute(reshape([BC{:}],BCLen,length(BC)),[2 1]);
    BCSeq=BCList{x};
    UMI = BCChar(:,13:20);
    UMI = mat2cell(UMI,repmat(1,size(UMI,1),1),8);
    [UMIList U] = umiCollapse(UMI);
    count = histc(U,1:max(U));
    UMI = UMIList(count>cutoff);
    UMIInd = 1:length(UMIList);
    UMIInd = UMIInd(count>cutoff);
    
    if ~isempty(UMI)
        for y=1:length(UMI)
            disp(y)
            UMISeqInd = locBamInd(U==UMIInd(y))+1;
            UMIReads = length(UMISeqInd);
            if UMIReads>UMIlim
                UMISeqInd = UMISeqInd(sort(randsample(1:UMIReads,UMIlim)));
            end
            UMIBam = getSubset(m,UMISeqInd);
            %        q{y}=UMIBam.Quality;
            UMIBamFile = [fPath  'BCFASTQ\' fName 'Bam' num2str(x) '_' num2str(y) '.bam'];
            UMIFASTAFile = [fPath  'BCFASTQ\' fName 'Bam' num2str(x) '.fa'];
            
            UMIDepth = [fPath  'BCFASTQ\' fName 'Depth' num2str(x) '_' num2str(y) '.fa'];
            if exist(UMIBamFile,'file')
                delete(UMIBamFile);
            end
            write(UMIBam,UMIBamFile);
            
            message = ['C:\cygwin64\bin\bash -login -c "samtools fasta -2 ""/\' UMIFASTAFile '"" ""/\' UMIBamFile '"""'];
            system(message);
            
            
            message = ['C:\cygwin64\bin\bash -login -c "samtools index ""/\' UMIBamFile '"""'];
            system(message);
            message = ['C:\cygwin64\bin\bash -login -c "samtools depth ""/\' UMIBamFile '"" > ""/\' UMIDepth '"""'];
            system(message);
            fid = fopen(UMIDepth,'r');
            dd=textscan(fid,'%s%u%u');
            fclose (fid);
            %%%%Make concensus sequence
            s=fastaread(UMIFASTAFile);
            C{a,y} = seqconsensus(s,'gaps','all');
            %%%%BLAST against TCR regions
            %%%%alignment
            bOut = [fPath  'BCFASTQ\' fName 'BC' num2str(x) 'UMI' num2str(y) 'BLAST.txt'];
            message = ['blastn -task blastn -evalue .001 -query "' UMIFASTAFile '" -db "' TCRBLASTIndex '" -out "' bOut '"'];
            [foo foo2]=system(message);
            if g ==2 %%%%Blast TRBC upstream regions to hg38
                
                bOutUp = [fPath  'BCFASTQ\' fName 'BC' num2str(x) 'UMI' num2str(y) 'upstreamBLAST.txt'];
                message = ['blastn -task blastn -evalue .001 -query "' UMIFASTAFile '" -db "' TRBCupstBLASTIndex '" -out "' bOutUp '"'];
                [foo foo2]=system(message);
                
            end
            
            
            
            
            %%%%Parse BLAST file
            fid = fopen(bOut,'r');
            results = textscan(fid,'%s','Delimiter','\n');
            results = results{1};
            fclose(fid);
            %%%Parse BLAST results
            entryR = strfind(results,'Query=');
            entryL = cellfun(@isempty, entryR)==0;
            entryI = find(entryL); %%%Start of each read BLAST results
            
            %%%Parse BLAST Alignments
            entryR = strfind(results,'Query');
            entryA = all([cellfun(@isempty, entryR)==0 entryL==0],2);
            entryAI = find(entryA); %%%Start of each alignment
            
            
            gR = strfind(results,'>');
            gL = cellfun(@isempty, gR)==0;
            gI = find(gL);  %%Each alignment
            
            gUTR = strfind(results,'UTR');
            gUTRL=cellfun(@isempty, gUTR)==0;
            gUTRR= find(gUTRL);
            
            scoreR = strfind(results,'Score =');
                    scoreL = cellfun(@isempty, scoreR)==0;
                    scoreI = find(scoreL);  %%Each alignment
            
            regHit = zeros(6,1);
            regTopHit = cell(6,1);
            topHitC= zeros(6,1);
            topHitF=zeros(6,1);
            regUTRHit = zeros(6,1);
            regUTRTopHit = cell(6,1);
            topHitUTRC= zeros(6,1);
            for r=1:6
                w = strfind(results,regions{r});
                v = find(all([gL gUTRL==0 cellfun(@isempty, w)==0],2));%%Alignments of a given region
                vUTR = find(all([gL gUTRL cellfun(@isempty, w)==0],2));%%Alignments of a UTR of given region
                if length(v)>0
                    align=all(cat(3,repmat(entryI',length(v),1)-repmat(v,1,length(entryI))<0,repmat([entryI(2:end);length(results)]',length(v),1)-repmat(v,1,length(entryI))>0),3);
                    align2=align(:,sum(align,1)>0);
                    topHits=zeros(size(align2,2),1);
                    for z=1:size(align2,2)
                        if isempty(find(align2(:,z),1))==0
                            topHits(z)=find(align2(:,z),1);
                        end
                    end
                    topHitsInd=v(topHits);
                    aa=results(topHitsInd);
                    exactRegion=cell(length(aa),1);
                    for z=1:length(aa)
                        [s e]=regexp(aa{z},[regions{r} '[0-9]*[A-Z]*-*[0-9]*']);
                        exactRegion{z}=aa{z}(s:e);
                    end
                    regHit(r)=length(exactRegion);
                    [u foo uI]=unique(exactRegion);
                    uC=histcounts(uI,1:length(u)+1);
                    h=find(uC==max(uC),1);
                    regTopHit{r}=u{h};
                    topHitC(r)=max(uC);
                    uC2=uC;
                    uC2(h)=0;
                    topHitF(r)=max(uC)/sum(max(uC)+max(uC2));
                end
                if length(vUTR)>0 %%%%Aggregate UTR alignments
                    align=all(cat(3,repmat(entryI',length(vUTR),1)-repmat(vUTR,1,length(entryI))<0,repmat([entryI(2:end);length(results)]',length(vUTR),1)-repmat(vUTR,1,length(entryI))>0),3);
                    align2=align(:,sum(align,1)>0);
                    topHits=zeros(size(align2,2),1);
                    for z=1:size(align2,2)
                        if isempty(find(align2(:,z),1))==0
                            topHits(z)=find(align2(:,z),1);
                        end
                    end
                    topHitsInd=vUTR(topHits);
                    aa=results(topHitsInd);
                    exactRegion=cell(length(aa),1);
                    for z=1:length(aa)
                        [s e]=regexp(aa{z},[regions{r} '[0-9]*[A-Z]*-*[0-9]*']);
                        exactRegion{z}=aa{z}(s:e);
                    end
                    regUTRHit(r)=length(exactRegion);
                    [u foo uI]=unique(exactRegion);
                    uC=histcounts(uI,1:length(u)+1);
                    h=find(uC==max(uC),1);
                    regUTRTopHit{r}=[u{h} '_5UTR'];
                    topHitUTRC(r)=max(uC);
                end
            end
            
            CDR3hit=0;
            vFinal=[];
            jFinal=[];
            jUTRFinal=[];
            vCFinal=0;
            jCFinal=0;
            jUTRCFinal=0;
            vFFinal=0;
            jFFinal=0;
            finalCDR3=[];
            finalCDR3nuc=[];
            finalCDR3c=0;
            finalCDR3f=0;
            chain=0;
            if sum(topHitC(2))>sum(topHitC(5)) && (topHitF(1)>vCut || strcmp(regTopHit{1},'TRBV7-2')) && topHitF(2)>jCut && sum(topHitC(1:2)>4)==2 %%%UMI is TRB and passes QC filters
                CDR3hit=1;
                chain=1;
                jUTRFinal=regUTRTopHit{2};
                jUTRCFinal=topHitUTRC(2);
                vFinal=regTopHit{1};
                jFinal=regTopHit{2};
                cFinal='TRBC';
                vCFinal=topHitC(1);
                jCFinal=topHitC(2);
                
                vFFinal=topHitF(1);
                jFFinal=topHitF(2);
                w = strfind(results,regTopHit{2});
                v = find(all([gL gUTRL==0 cellfun(@isempty, w)==0],2));%%Alignments of top TRBJ
                vT=results(v);
                exactRegion=cell(length(vT),1);
                for z=1:length(vT)
                    [s e]=regexp(vT{z},[regions{2} '[0-9]*[A-Z]*-*[0-9]*']);
                    exactRegion{z}=vT{z}(s:e);
                end
                vOK=strcmp(exactRegion,regTopHit{2});
                v=v(vOK);
                vJ=v;
                alignJ=sum(all(cat(3,repmat(entryI',length(v),1)-repmat(v,1,length(entryI))<0,repmat([entryI(2:end);length(results)]',length(v),1)-repmat(v,1,length(entryI))>0),3))>0;
                w = strfind(results,regTopHit{1});
                v = find(all([gL gUTRL==0 cellfun(@isempty, w)==0],2));%%Alignments of top TRBV
                vT=results(v);
                exactRegion=cell(length(vT),1);
                for z=1:length(vT)
                    [s e]=regexp(vT{z},[regions{1} '[0-9]*[A-Z]*-*[0-9]*']);
                    exactRegion{z}=vT{z}(s:e);
                end
                vOK=strcmp(exactRegion,regTopHit{1});
                v=v(vOK);
                vV=v;
                alignV=sum(all(cat(3,repmat(entryI',length(v),1)-repmat(v,1,length(entryI))<0,repmat([entryI(2:end);length(results)]',length(v),1)-repmat(v,1,length(entryI))>0),3))>0;
                CDRLines=entryI(sum([uint8(alignJ);uint8(alignV)])==2);%%%First line of results for reads with V and J region
                
                CDRreads=UMISeqInd(sum([uint8(alignJ);uint8(alignV)])==2); %%%Index of reads with V and J region in original bam file
            elseif sum(topHitC(2))<sum(topHitC(5)) && topHitF(4)>vCut && topHitF(5)>jCut && sum(topHitC(4:5)>4)==2 %%%UMI is TRA and passes QC filters
                CDR3hit=1;
                chain=2;
                jUTRFinal=regUTRTopHit{5};
                jUTRCFinal=topHitUTRC(5);
                vFinal=regTopHit{4};
                jFinal=regTopHit{5};
                cFinal='TRAC';
                vCFinal=topHitC(4);
                jCFinal=topHitC(5);
                vFFinal=topHitF(4);
                jFFinal=topHitF(5);
                w = strfind(results,regTopHit{5});
                v = find(all([gL gUTRL==0 cellfun(@isempty, w)==0],2));%%Alignments of top TRAJ
                vT=results(v);
                exactRegion=cell(length(vT),1);
                for z=1:length(vT)
                    [s e]=regexp(vT{z},[regions{5} '[0-9]*[A-Z]*-*[0-9]*']);
                    exactRegion{z}=vT{z}(s:e);
                end
                vOK=strcmp(exactRegion,regTopHit{5});
                v=v(vOK);
                vJ=v;
                alignJ=sum(all(cat(3,repmat(entryI',length(v),1)-repmat(v,1,length(entryI))<0,repmat([entryI(2:end);length(results)]',length(v),1)-repmat(v,1,length(entryI))>0),3))>0;
                w = strfind(results,regTopHit{4});
                v = find(all([gL gUTRL==0 cellfun(@isempty, w)==0],2));%%Alignments of top TRAV
                vT=results(v);
                exactRegion=cell(length(vT),1);
                for z=1:length(vT)
                    [s e]=regexp(vT{z},[regions{4} '[0-9]*[A-Z]*-*[0-9]*']);
                    exactRegion{z}=vT{z}(s:e);
                end
                vOK=strcmp(exactRegion,regTopHit{4});
                v=v(vOK);
                vV=v;
                alignV=sum(all(cat(3,repmat(entryI',length(v),1)-repmat(v,1,length(entryI))<0,repmat([entryI(2:end);length(results)]',length(v),1)-repmat(v,1,length(entryI))>0),3))>0;
                CDRLines=entryI(sum([uint8(alignJ);uint8(alignV)])==2);
                CDRreads=UMISeqInd(sum([uint8(alignJ);uint8(alignV)])==2);
            elseif sum(regHit(1:3))>sum(regHit(4:6)) || sum(regHit(2))>sum(regHit(5))
                if topHitF(1)>0.6
                    vFinal=regTopHit{1};
                end
                if topHitF(2)>0.6
                    jFinal=regTopHit{2};
                end
                jUTRFinal=regUTRTopHit{2};
                jUTRCFinal=topHitUTRC(2);
                vCFinal=topHitC(1);
                jCFinal=topHitC(2);
                vFFinal=topHitF(1);
                jFFinal=topHitF(2);
                chain=1;
            elseif sum(regHit(1:3))<sum(regHit(4:6)) || sum(regHit(2))<sum(regHit(5))
                if topHitF(4)>0.6
                    vFinal=regTopHit{4};
                end
                if topHitF(5)>0.6
                    jFinal=regTopHit{5};
                end
                jUTRFinal=regUTRTopHit{5};
                jUTRCFinal=topHitUTRC(5);
                vCFinal=topHitC(4);
                jCFinal=topHitC(5);
                vFFinal=topHitF(4);
                jFFinal=topHitF(5);
                chain=2;
            end
            if CDR3hit==1 && isempty(CDRreads)==0
                
                readsSeq=m.Sequence(CDRreads);
                
                vSeq=TCRsequences{strcmp(TCRRegions,vFinal)};
                vCDR3=CDR3base(strcmp(TCRRegions,vFinal));
                jSeq=TCRsequences{strcmp(TCRRegions,jFinal)};
                jCDR3=CDR3base(strcmp(TCRRegions,jFinal));
                jLen=length(jSeq);
                cSeq=TCRsequences{strcmp(TCRRegions,cFinal)};
                vJAlign=repmat(vJ,1,length(CDRLines))-repmat(CDRLines',length(vJ),1);
                vJAlign(vJAlign<=0)=10000; %%%min value in each column is best alignment of J region for each CDR3read
                vVAlign=repmat(vV,1,length(CDRLines))-repmat(CDRLines',length(vV),1);
                vVAlign(vVAlign<=0)=10000;%%%min value in each column is best alignment of V region for each CDR3read
                CDR3=cell(length(CDRLines),1);
                CDR3nuc=cell(length(CDRLines),1);
                for z=1:length(CDRLines)
                    read=readsSeq{z};
                    
                    vJA=vJ(vJAlign(:,z)==min(vJAlign(:,z)));%%%%Line number for best J alignment
                    
                    scoreLoc=scoreI-vJA;
                                if sum(scoreLoc>0)>1
                                    scoreLoc(scoreLoc<0)=10000;
                                    sC=scoreI(find(scoreLoc==min(scoreLoc))+1);
                                    endJAL=entryAI<sC;
                                    endJA=entryAI(find(endJAL,1,'last'));
                                else
                                    endJA = entryAI(end);
                                end
                    [JSQ,JEQ]=regexp(results{endJA},'[0-9]*'); %%Location of Base numbers in sequencing read for beginning and end of alignment
                    JreadE=str2double(results{endJA}(JSQ(2):JEQ(2))); %Base numbers in sequencing read for beginning and end of alignment
                    
                    [JSS,JES]=regexp(results{endJA+2},'[0-9]*'); %%Location of Base numbers in J sequence for beginning and end of alignment
                    JJE=str2double(results{endJA+2}(JSS(2):JES(2))); %Base numbers in J sequence for beginning and end of alignment
                    
                    vVA=vV(vVAlign(:,z)==min(vVAlign(:,z)));%%%%Line number for best V alignment
                    
                                scoreLoc=scoreI-vVA;
                                if sum(scoreLoc>0)>1
                                    scoreLoc(scoreLoc<0)=10000;
                                    sC=scoreI(find(scoreLoc==min(scoreLoc))+1);
                                    endVAL=entryAI<sC;
                                    endVA=entryAI(find(endVAL,1,'last'));
                                else
                                    endVA = entryAI(end);
                                end

                    [VSQ,VEQ]=regexp(results{endVA},'[0-9]*'); %%Location of Base numbers in sequencing read for beginning and end of alignment
                    VreadE=str2double(results{endVA}(VSQ(2):VEQ(2)));
                    [VSS,VES]=regexp(results{endVA+2},'[0-9]*'); %%Location of Base numbers in V sequence for beginning and end of alignment
                    VVE=str2double(results{endVA+2}(VSS(2):VES(2)));
                    
                    missingJ=JJE-jCDR3; %%%Check to make sure whole J region in alignment
                    if JreadE>VreadE
                        if missingJ>0
                            finalSeq=[vSeq(1:VVE) read(VreadE+1:JreadE)];
                            CDR3{z}=nt2aa(finalSeq(vCDR3:end-(JJE-jCDR3)),'ACGTOnly','false');
                            CDR3nuc{z}=finalSeq(vCDR3:end-(JJE-jCDR3));
                            
                        else
                            finalSeq=[vSeq(1:VVE) read(VreadE+1:end)];
                            CDR3{z}=[nt2aa(finalSeq(vCDR3:end),'ACGTOnly','false') '_'];
                            CDR3nuc{z}=finalSeq(vCDR3:end);
                        end
                    end
                    
                end
                CDR3nuc=CDR3nuc(cellfun(@isempty,CDR3)==0);
                CDR3=CDR3(cellfun(@isempty,CDR3)==0);
                
                if isempty(CDR3)==0
                    trunc=strfind(CDR3,'_');
                    nonTrunc=cellfun(@isempty,trunc);
                    if sum(nonTrunc)>0
                        CDR3=CDR3(nonTrunc);
                    end
                    [uCDR3 foo CI]=unique(CDR3);
                    cCDR3=histcounts(CI,1:length(uCDR3)+1);
                    CDI=find(cCDR3==max(cCDR3),1);
                    cc=uCDR3{CDI};
                    %if isempty(strfind(cc,'*'))
                        finalCDR3=uCDR3{CDI};
                        SelCDR3nuc=CDR3nuc{CI==CDI};
                        finalCDR3nuc=seqconsensus(SelCDR3nuc);
                        %Collect 1 hamming dist CDR3 for stats
                        L=cellfun(@length,uCDR3);
                        pHamInd=find(L==length(finalCDR3));
                        pHam=uCDR3(pHamInd);
                        hamL=false(length(pHam),1);
                        for zz=1:length(pHam)
                            hamDis=sum((int8(pHam{zz})-int8(finalCDR3))~=0);
                            if hamDis<=1
                                hamL(zz)=1;
                            end
                        end
                        finalCDR3c=sum(cCDR3(hamL));
                        finalCDR3f=finalCDR3c/length(CI);
                    %end
                end
            end
            if chain ==2
                chr=strcmp(chrA,dd{1});
                [alI al]=ismember(chrPosA,dd{2}(chr));
                al=al(alI);
                alF=double(dd{3}(al));
                alF=alF/sum(alF);
                
                if isempty(finalCDR3)==0
                    depACDR3=[depACDR3 zeros(size(depACDR3,1),1)];
                    depACDR3(alI,end)=alF;
                elseif jUTRCFinal>15
                    depATrunc=[depATrunc zeros(size(depATrunc,1),1)];
                    depATrunc(alI,end)=alF;
                else
                    depACon=[depACon zeros(size(depACon,1),1)];
                    depACon(alI,end)=alF;
                end
            elseif chain ==1
                if g==2
                    
                    locTRBC1A=zeros(52250,1);
                    locTRBC2A=zeros(7891,1);
                    while exist(bOutUp,'file')==0
                        pause(1)
                    end
                    %%%%Parse BLAST file
                    fid = fopen(bOutUp,'r');
                    results = textscan(fid,'%s','Delimiter','\n');
                    results = results{1};
                    fclose(fid);
                    %%%Parse BLAST results
                    entryR = strfind(results,'Query=');
                    entryL = cellfun(@isempty, entryR)==0;
                    entryI = find(entryL); %%%Start of each read BLAST results
                    
                    %%%Parse BLAST Alignments
                    entryR = strfind(results,'Sbjct');
                    entryA = all([cellfun(@isempty, entryR)==0 entryL==0],2);
                    entryAI = find(entryA); %%%Start of each alignment
                    
                    
                    gR = strfind(results,'>');
                    gL = cellfun(@isempty, gR)==0;
                    gI = find(gL);  %%Each alignment
                    
                    scoreR = strfind(results,'Score =');
                    scoreL = cellfun(@isempty, scoreR)==0;
                    scoreI = find(scoreL);  %%Each alignment
                    
                    for zz=1:length(entryI)-1
                        locI=gI(all([gI>entryI(zz) gI<entryI(zz+1)],2));
                        if isempty(locI)==0
                            for z3=1:length(locI)-1
                                q=strfind(results{locI(z3)},'TRBC1');
                                
                                
                                entries=entryAI(all([entryAI>locI(z3) entryAI<locI(z3+1)],2));
                                scoreLoc=scoreI-entries(1);
                                if sum(scoreLoc>0)>1
                                scoreLoc(scoreLoc<0)=10000;
                                sC=scoreI(find(scoreLoc==min(scoreLoc))+1);
                                entries=entries(entries<sC);
                                end
                                
                                [VSQ,VEQ]=regexp(results{entries(1)},'[0-9]*');
                                sH=str2double(results{entries(1)}(VSQ(1):VEQ(1)));
                                [VSQ,VEQ]=regexp(results{entries(end)},'[0-9]*');
                                sE=str2double(results{entries(end)}(VSQ(2):VEQ(2)));
                                if isempty(q)==0
                                    locTRBC1A(sH:sE)=locTRBC1A(sH:sE)+1;
                                else
                                    locTRBC2A(sH:sE)=locTRBC2A(sH:sE)+1;
                                end
                            end
                            
                            q=strfind(results{locI(end)},'TRBC1');
                            entries=entryAI(all([entryAI>locI(end) entryAI<entryI(zz+1)],2));
                            
                            [VSQ,VEQ]=regexp(results{entries(1)},'[0-9]*');
                            sH=str2double(results{entries(1)}(VSQ(1):VEQ(1)));
                            [VSQ,VEQ]=regexp(results{entries(end)},'[0-9]*');
                            sE=str2double(results{entries(end)}(VSQ(2):VEQ(2)));
                            if isempty(q)==0
                                    locTRBC1A(sH:sE)=locTRBC1A(sH:sE)+1;
                                else
                                    locTRBC2A(sH:sE)=locTRBC2A(sH:sE)+1;
                                end
                        end
                    end
                    
                    if isempty(finalCDR3)==0
                        hTRBC1ACDR3=[hTRBC1ACDR3 locTRBC1A];
                        hTRBC2ACDR3=[hTRBC2ACDR3 locTRBC2A];
                    elseif jUTRCFinal>15
                        hTRBC1ATrun=[hTRBC1ATrun locTRBC1A];
                        hTRBC2ATrun=[hTRBC2ATrun locTRBC2A];
                    else
                        hTRBC1ACon=[hTRBC1ACon locTRBC1A];
                        hTRBC2ACon=[hTRBC2ACon locTRBC2A];
                    end
                    
                    
                    
                end
                
                
                
                
                chr=strcmp(chrB,dd{1});
                [alI al]=ismember(chrPosB,dd{2}(chr));
                al=al(alI);
                alF=double(dd{3}(al));
                alF=alF/sum(alF);
                
                if isempty(finalCDR3)==0
                    depBCDR3=[depBCDR3 zeros(size(depBCDR3,1),1)];
                    depBCDR3(alI,end)=alF;
                elseif jUTRCFinal>15
                    depBTrunc=[depBTrunc zeros(size(depBTrunc,1),1)];
                    depBTrunc(alI,end)=alF;
                else
                    depBCon=[depBCon zeros(size(depBCon,1),1)];
                    depBCon(alI,end)=alF;
                end
                
            end
            
            fprintf(fidM,'%s\t%s\t%s\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%u\t%.2f\t%u\t%.2f\t%u\t%u\t%.2f\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n',['BC' num2str(x) 'UMI' num2str(y)],BCSeq,UMI{y},x,y,UMIReads, ...
                vFinal,jFinal,jUTRFinal,finalCDR3,finalCDR3nuc,vCFinal,vFFinal,jCFinal,jFFinal,jUTRCFinal,finalCDR3c,finalCDR3f,regHit(1),regUTRHit(2),regHit(2),regHit(3),regHit(4),regUTRHit(5),regHit(5),regHit(6));
            
            
            delete(UMIFASTAFile)
            delete(UMIBamFile)
            
            
            
            %             delete(bOut)
            
        end
    end
    toc
end
depthResults = [alignDir fName 'DepthResults' num2str(iSlice) '.mat'];
save(depthResults,'depACon','depATrunc','depACDR3','depBCon','depBTrunc','depBCDR3','hTRBC1ACDR3','hTRBC1ATrun','hTRBC1ACon','hTRBC2ACDR3','hTRBC2ATrun','hTRBC2ACon','C');
fclose(fidM);


%  end