function  plotPHTMesh3DMP_faceColored( PHTelem, GIFTmesh,indexPatch,indexElem )
%plots selected patch
%surfaces marked with dots with different colors
%front -red
%right - blue
%back -black
%left - green
%down - magenta
%up - yellow

numPts = 21; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);
wref = linspace(-1,1,numPts);

figure
view(-37.50,30)

for patchIndex = indexPatch
    
    for i=1:length(PHTelem{patchIndex})
       
        %if isempty(PHTelem{patchIndex}(indexElem).children)
        
        %store the corners of the element in parameter space for easy
        %access
        pxmin = PHTelem{patchIndex}(indexElem).vertex(1);
        pymin = PHTelem{patchIndex}(indexElem).vertex(2);
        pzmin = PHTelem{patchIndex}(indexElem).vertex(3);
        pxmax = PHTelem{patchIndex}(indexElem).vertex(4);
        pymax = PHTelem{patchIndex}(indexElem).vertex(5);
        pzmax = PHTelem{patchIndex}(indexElem).vertex(6);
        
        %determine the edges of the element in the physical space
        %             coord_up_left = zeros(numPts,3);
        %             coord_down_left = zeros(numPts,3);
        %             coord_up_right = zeros(numPts,3);
        %             coord_down_right = zeros(numPts,3);
        %             coord_up_front = zeros(numPts,3);
        %             coord_down_front = zeros(numPts,3);
        %             coord_up_back = zeros(numPts,3);
        %             coord_down_back = zeros(numPts,3);
        %             coord_left_front = zeros(numPts,3);
        %             coord_right_front = zeros(numPts,3);
        %             coord_left_back = zeros(numPts,3);
        %             coord_right_back = zeros(numPts,3);
        face_up=zeros(numPts,numPts,3);
        face_down=zeros(numPts,numPts,3);
        face_left=zeros(numPts,numPts,3);
        face_right=zeros(numPts,numPts,3);
        face_front=zeros(numPts,numPts,3);
        face_back=zeros(numPts,numPts,3);
        
        for k=1:numPts
            for j=1:numPts
                face_front(k,j,:)= paramMap3D(GIFTmesh{patchIndex}, uref(k), -1, wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                face_right(k,j,:)= paramMap3D(GIFTmesh{patchIndex}, 1, vref(k), wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                face_back(k,j,:)= paramMap3D(GIFTmesh{patchIndex}, uref(k), 1, wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                face_left(k,j,:)= paramMap3D(GIFTmesh{patchIndex}, -1, vref(k), wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                face_down(k,j,:)= paramMap3D(GIFTmesh{patchIndex}, uref(k), vref(j), -1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                face_up(k,j,:)= paramMap3D(GIFTmesh{patchIndex}, uref(k), vref(j), 1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                
            end
        end
        
        
        plot3(face_front(:,:,1),face_front(:,:,2),face_front(:,:,3),'.c')
        hold on
        plot3(face_right(:,:,1),face_right(:,:,2),face_right(:,:,3),'.y')
        hold on
        plot3(face_back(:,:,1),face_back(:,:,2),face_back(:,:,3),'.b')
        hold on
        plot3(face_left(:,:,1),face_left(:,:,2),face_left(:,:,3),'.g')
        hold on
        plot3(face_down(:,:,1),face_down(:,:,2),face_down(:,:,3),'.m')
        hold on
        plot3(face_up(:,:,1),face_up(:,:,2),face_up(:,:,3),'.r')
        hold on
        
        
        %  end
    end
end

axis tight
hold on
drawnow
