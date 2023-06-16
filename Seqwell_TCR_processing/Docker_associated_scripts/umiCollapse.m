%%%%%%Collapse of UMI sequences based on directional adjacency network analysis  (UMI-Tools Smith T, Sudbury,I) 
    function [tUMI umiInd] =umiCollapse(u) 
            [rUMI, foo, umiInd] = unique(u);    %%%Unique UMI
            rUMIInd=(1:length(rUMI))';
            count = [];
            count(:,1) = histc(umiInd,1:max(umiInd));      %%%Number of each unique UMI
            [count,ord] = sort(count,'descend'); %%%ordering UMI by count
            rUMI=rUMI(ord);
            rUMIInd=rUMIInd(ord);
            uTemp = cell(1);
            test=0;
            c=0;
            tUMI=cell(1);
            rChar = uint8(permute(reshape([rUMI{:}],length(rUMI{1}),length(rUMI)),[2 1])); %%%%Make UMI character array to allow comparison for UMI collapse by directional adjacency network analysis  (UMI-Tools Smith T, Sudbury,I)
            while test==0
                if c==0
                    tUMI=rUMI(1);
                else
                    tUMI=[tUMI rUMI(1)];
                end
                c=c+1;
                hamDistRep = sum(rChar~=repmat(rChar(1,:),size(rChar,1),1),2);
                %%%%Create Del variants
                rCharDel = zeros(length(rChar(1,:))-1,length(rChar(1,:)),'uint8');
                for zz=1:length(rChar(1,:))-2
                    if zz==1
                        rCharDel(zz,:) = [rChar(1,2:end) uint8('T')];
                    else
                        rCharDel(zz,:) = [rChar(1,1:zz-1) rChar(1,zz+1:end) uint8('T')];
                    end
                end
                rCharDel(end,:) = [rChar(1,1:end-1) uint8('T')];
                hamDistDel = min(sum(repmat(permute(rCharDel,[3 2 1]),size(rChar,1),1,1)~=repmat(rChar,1,1,size(rCharDel,1)),2),[],3);
                
                %%%%Create Insert variants
                rCharIn = zeros(length(rChar(1,:))-1,length(rChar(1,:)),'uint8');
                for zz=1:length(rChar(1,:))-2
                    if zz==1
                        rCharIn(zz,:) = [uint8('N') rChar(1,1:end-1)];
                    else
                        rCharIn(zz,:) = [rChar(1,1:zz-1) uint8('N') rChar(1,zz:end-1)];
                    end
                end
                rCharIn(end,:) = [rChar(1,1:end-1) uint8('N')];
                hamDistIn = min(sum(repmat(permute(rCharIn,[3 2 1]),size(rChar,1),1,1)~=repmat(rChar,1,1,size(rCharIn,1)),2),[],3);
                hamDistTot = any([hamDistRep==1 hamDistIn==1 hamDistDel==0],2);
                
                
                relNodes = all([2*count-1<count(1) hamDistTot],2);   %%%%Logical for UMIs 1 change away
                relNodes(1)=1;
                relNodesSeq=rChar(relNodes,:);           %%%%UMIs 1 change away
                relNodesCount = count(relNodes);
                indRelNodes = rUMIInd(relNodes);
                umiInd(ismember(umiInd,indRelNodes))=c*1000;
                if sum(relNodes==0)==0
                    test=1;
                elseif sum(relNodes)==1 %%%No matching UMI with hamming dist of 1
                    rChar = rChar(relNodes==0,:);
                    count = count(relNodes==0);
                    rUMI = rUMI(relNodes==0);
                    rUMIInd = rUMIInd(relNodes==0);
                else    %%%%test whether UMI with 2 changes related to UMIs with 1 change
                    a=1;
                    while a==1
                        a=0;
                        unRelNodesSeq = rChar(relNodes==0,:);
                        if isempty(unRelNodesSeq)
                            break
                        end
                        unRelNodesCount= count(relNodes==0);
                        unRelNodesInd = find(relNodes==0);
                        hamDist2=sum(repmat(relNodesSeq,1,1,size(unRelNodesSeq,1))~=repmat(permute(unRelNodesSeq,[3 2 1]),size(relNodesSeq,1),1,1),2);
                        countComp = 2*repmat(relNodesCount,1,1,size(unRelNodesSeq,1))-1>repmat(permute(unRelNodesCount,[3 2 1]),size(relNodesSeq,1),1,1);
                        relNodesEdge2 = squeeze(sum(all([countComp hamDist2==1],2)))>0;
                        if sum(relNodesEdge2)>0
                            a=1;
                            relNodes(unRelNodesInd(relNodesEdge2))=1;
                            relNodesSeq=rChar(relNodes,:);           %%%%UMIs 1 change away
                            relNodesCount = count(relNodes);
                            indRelNodes = rUMIInd(relNodes);
                            umiInd(ismember(umiInd,indRelNodes))=c*1000;
                        end
                    end
                    unRelNodesCount=[];
                    if sum(relNodes==0)==0
                        test=1;
                    else
                        rChar = rChar(relNodes==0,:);
                        count = count(relNodes==0);
                        rUMI = rUMI(relNodes==0);
                        rUMIInd = rUMIInd(relNodes==0);

                    end
                    
                end
            end
            umiInd = umiInd/1000;
        