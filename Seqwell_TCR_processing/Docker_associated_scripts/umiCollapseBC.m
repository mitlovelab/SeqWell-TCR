 function [finalBC readInd] =umiCollapseBC(u) 
            [rUMI, foo, umiInd] = unique(u);    %%%Unique UMI
            rUMIInd=(1:length(rUMI))';
            count = [];
            count(:,1) = histc(umiInd,1:max(umiInd));      %%%Number of each unique UMI
            [count,ord] = sort(count,'descend'); %%%ordering UMI by count
            rUMI=rUMI(ord);
            rUMIInd=rUMIInd(ord);
            [foo reOrd]=sort(rUMIInd);
            rChar = uint8(permute(reshape([rUMI{:}],length(rUMI{1}),length(rUMI)),[2 1])); %%%%Make UMI character array to allow comparison for UMI collapse by directional adjacency network analysis  (UMI-Tools Smith T, Sudbury,I)
            bcLen=size(rChar,2);
            finalBC={};
            finalInd=zeros(size(rChar,1),1,'uint32');
            c=1;
            h = waitbar(0,'Collapsing barcodes..');
            nS=length(rUMI);
            finished=0;
 while finished==0
     nTest = floor(10000000000/(size(rChar,1)*bcLen)); %%%%How many BC can be screened within memory contraints
     waitbar((nS-size(rChar,1))/nS)
     if nTest >= size(rChar,1)
         nTest=size(rChar,1);
         finished = 1;
        rCharTest = double(uint16(rChar));
        rUMITest=rUMI;
        indNum=uint32(repmat((1:nTest)',1,nTest));
     hDist=(pdist2(rCharTest,rCharTest,'hamming')*bcLen)>1;
     indNum(hDist)=nTest+1;
     indNum=min(indNum);
     indVal=unique(indNum);
     finalBC=[finalBC;rUMITest(indVal)];
     rCharTest=rCharTest(indVal,:);
     cE=c+size(rCharTest,1)-1;
     indNum=repmat(uint32((c:cE)'),1,size(rChar,1));
     hDist=(pdist2(rCharTest,double(uint16(rChar)),'hamming')*bcLen)>1;
     indNum(hDist)=cE+1;
     indNum=min(indNum);
     indNum(indNum>cE)=0;
     i=rUMIInd(indNum>0);
     finalInd(i)=indNum(indNum>0);
     else
     
     rCharTest = double(uint16(rChar(1:nTest,:)));
     rUMITest=rUMI(1:nTest);
     indNum=uint32(repmat((1:nTest)',1,nTest));
     hDist=(pdist2(rCharTest,rCharTest,'hamming')*bcLen)>1;
     indNum(hDist)=nTest+1;
     indNum=min(indNum);
     indVal=unique(indNum);
     finalBC=[finalBC;rUMITest(indVal)];
     rCharTest=rCharTest(indVal,:);
     cE=c+size(rCharTest,1)-1;
     indNum=repmat(uint32((c:cE)'),1,size(rChar,1));
     hDist=(pdist2(rCharTest,double(uint16(rChar)),'hamming')*bcLen)>1;
     indNum(hDist)=cE+1;
     indNum=min(indNum);
     indNum(indNum>cE)=0;
     i=rUMIInd(indNum>0);
     finalInd(i)=indNum(indNum>0);
     rChar=rChar(indNum==0,:);
     rUMI=rUMI(indNum==0);
     rUMIInd=rUMIInd(indNum==0);
     c=cE+1;
     end
 end
 
 readInd=finalInd(umiInd);
 
 close(h);