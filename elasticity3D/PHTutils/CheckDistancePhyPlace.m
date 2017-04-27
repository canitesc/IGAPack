function [ ] = CheckDistancePhyPlace( PHUTelem, GIFTmesh,patchBoundaries)
numPts = 21; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);
wref = linspace(-1,1,numPts);

%figure
view(-37.50,30)

patchA= cell2mat(patchBoundaries(1));
patchB=cell2mat(patchBoundaries(2));
faceA=cell2mat(patchBoundaries(3));
faceB=cell2mat(patchBoundaries(4));
coordFaceA=[];
coordFaceB=[];
for patch=1:2

    if patch==1
        patchIndex=patchA;
        face=faceA;
    else
        patchIndex=patchB;
        face=faceB;
    end
    indexElem=1;
    
    
    pxmin = PHUTelem{patchIndex}(indexElem).vertex(1);
    pymin = PHUTelem{patchIndex}(indexElem).vertex(2);
    pzmin = PHUTelem{patchIndex}(indexElem).vertex(3);
    pxmax = PHUTelem{patchIndex}(indexElem).vertex(4);
    pymax = PHUTelem{patchIndex}(indexElem).vertex(5);
    pzmax = PHUTelem{patchIndex}(indexElem).vertex(6);
    
%     face_patchA=zeros(numPts,numPts,3);
%     face_patchB=zeros(numPts,numPts,3);
    
    for i=1:numPts
        for j=1:numPts
        
            switch face
                
                case 1 %front
                    u = uref(i);
                    v = -1;
                    w = wref(j);
                    
                case 2 %right
                    u = 1;
                    v = vref(i);
                    w = wref(j);
                    
                case 3 %back
                    u = uref(i);
                    v = 1;
                    w = wref(j);
                    
                case 4 %left
                    u = -1;
                    v = vref(i);
                    w = wref(j);
                    
                case 5 %down
                    u = uref(i);
                    v = vref(j);
                    w = -1;
                    
                case 6 %up
                    u = uref(i);
                    v = vref(j);
                    w = 1;
            end
            
            if patch==1
               % face_patchA(i,j,:)= paramMap3D(GIFTmesh{patchIndex}, u, v, w, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                facePatchA=paramMap3D(GIFTmesh{patchIndex}, u, v, w, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                plot3(facePatchA(1),facePatchA(2),facePatchA(3),'*c')
                coordFaceA=[coordFaceA;facePatchA];
                hold on
              % pause
            else
               % face_patchB(i,j,:)= paramMap3D(GIFTmesh{patchIndex}, u, v, w, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                facePatchB= paramMap3D(GIFTmesh{patchIndex}, u, v, w, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                plot3(facePatchB(1),facePatchB(2),facePatchB(3),'.m')
                coordFaceB=[coordFaceB;facePatchB];
                hold on 
              % pause
            end
            
            
        end
    end
    
%     coordFaceA
%     coordFaceB
end

end

