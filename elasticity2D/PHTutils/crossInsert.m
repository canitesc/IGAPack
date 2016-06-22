function [ PHTelem, dimBasis ] = crossInsert( PHTelem, elm_index, dimBasis, p, q )
%inserts cross in PHTelem at elm_index
%dimBasis: dimension of the PHT space


%define corner indices
[ sw, south, se, west, ~, east, nw, north, ne ] = getCornerIndices( p,q );

for e=elm_index    
    lastElem = length(PHTelem);

    %collect information about new basis vertices and T junctions    
    [ newBasisVert,  ~, RTjunct, knotUl, knotUr, knotVd, knotVu, warningFlag ] = checkNeighbors( PHTelem, e );
    if warningFlag
        disp(['Skipping refinement in element ', num2str(e)])
        continue
    end
     
    %update parent element                 
    PHTelem(e).children = lastElem+1:lastElem+4;                    
    
    xmin = PHTelem(e).vertex(1);
    xmax = PHTelem(e).vertex(3);
    ymin = PHTelem(e).vertex(2);
    ymax = PHTelem(e).vertex(4);
    
    %add children elements and their easy neighbors
    %lower-left child
    PHTelem(lastElem+1).parent = e;
    PHTelem(lastElem+1).children = [];
    PHTelem(lastElem+1).vertex = [xmin, ymin, (xmin+xmax)/2, (ymin+ymax)/2];
    PHTelem(lastElem+1).level = PHTelem(e).level+1;
    PHTelem(lastElem+1).neighbor_up = lastElem+3;
    PHTelem(lastElem+1).neighbor_right = lastElem+2;    
    
    %lower-right child
    PHTelem(lastElem+2).parent = e;
    PHTelem(lastElem+2).children = [];
    PHTelem(lastElem+2).vertex = [(xmin+xmax)/2, ymin, xmax, (ymin+ymax)/2];
    PHTelem(lastElem+2).level = PHTelem(e).level+1;
    PHTelem(lastElem+2).neighbor_up = lastElem+4;
    PHTelem(lastElem+2).neighbor_left = lastElem+1;
    
    %upper-left child
    PHTelem(lastElem+3).parent = e;
    PHTelem(lastElem+3).children = [];
    PHTelem(lastElem+3).vertex = [xmin, (ymin+ymax)/2, (xmin+xmax)/2, ymax];
    PHTelem(lastElem+3).level = PHTelem(e).level+1;
    PHTelem(lastElem+3).neighbor_down = lastElem+1;
    PHTelem(lastElem+3).neighbor_right = lastElem+4;
    
    %upper-right child
    PHTelem(lastElem+4).parent = e;
    PHTelem(lastElem+4).children = [];
    PHTelem(lastElem+4).vertex = [(xmin+xmax)/2, (ymin+ymax)/2, xmax, ymax];
    PHTelem(lastElem+4).level = PHTelem(e).level+1;
    PHTelem(lastElem+4).neighbor_down = lastElem+2;
    PHTelem(lastElem+4).neighbor_left = lastElem+3;
    
    %add the neighbors outside the refined element, taking into account
    %T-junctions
    
    if (length(PHTelem(e).neighbor_down)==2)
        PHTelem(lastElem+1).neighbor_down = PHTelem(e).neighbor_down(1);
        PHTelem(lastElem+2).neighbor_down = PHTelem(e).neighbor_down(2);
    else
        PHTelem(lastElem+1).neighbor_down = PHTelem(e).neighbor_down;
        PHTelem(lastElem+2).neighbor_down = PHTelem(e).neighbor_down;
    end
    
    if (length(PHTelem(e).neighbor_right)==2)
        PHTelem(lastElem+2).neighbor_right = PHTelem(e).neighbor_right(1);
        PHTelem(lastElem+4).neighbor_right = PHTelem(e).neighbor_right(2);
    else
        PHTelem(lastElem+2).neighbor_right = PHTelem(e).neighbor_right;
        PHTelem(lastElem+4).neighbor_right = PHTelem(e).neighbor_right;
    end
        
    if (length(PHTelem(e).neighbor_up)==2)
        PHTelem(lastElem+3).neighbor_up = PHTelem(e).neighbor_up(1);
        PHTelem(lastElem+4).neighbor_up = PHTelem(e).neighbor_up(2);
    else
        PHTelem(lastElem+3).neighbor_up = PHTelem(e).neighbor_up;
        PHTelem(lastElem+4).neighbor_up = PHTelem(e).neighbor_up;
    end
    
    if (length(PHTelem(e).neighbor_left)==2)
        PHTelem(lastElem+1).neighbor_left = PHTelem(e).neighbor_left(1);
        PHTelem(lastElem+3).neighbor_left = PHTelem(e).neighbor_left(2);
    else
        PHTelem(lastElem+1).neighbor_left = PHTelem(e).neighbor_left;
        PHTelem(lastElem+3).neighbor_left = PHTelem(e).neighbor_left;
    end
    
    %handle the removed T-junctions
    dimBasisTemp = dimBasis;
    for ijunct = RTjunct
                
        %update the neighbors of the neighbor with children elements so we
        %get correct knotUln, knotUrn, knotVdn, knotVun
        switch ijunct
            case 1                
                down_neighbor = PHTelem(e).neighbor_down(1);
                parent_neighbor = PHTelem(down_neighbor).parent;
                PHTelem(parent_neighbor).neighbor_up = [lastElem+1, lastElem+2];   
                                
            case 2               
                right_neighbor = PHTelem(e).neighbor_right(1);
                parent_neighbor = PHTelem(right_neighbor).parent;
                PHTelem(parent_neighbor).neighbor_left = [lastElem+2, lastElem+4];
            case 3
                up_neighbor = PHTelem(e).neighbor_up(1);
                parent_neighbor = PHTelem(up_neighbor).parent;
                PHTelem(parent_neighbor).neighbor_down = [lastElem+3, lastElem+4];
            case 4
                left_neighbor = PHTelem(e).neighbor_left(1);
                parent_neighbor = PHTelem(left_neighbor).parent;
                PHTelem(parent_neighbor).neighbor_right = [lastElem+1, lastElem+3];
        end
                  
        [ newBasisVertn,  ~, RTjunctn, knotUln, knotUrn, knotVdn, knotVun ] = checkNeighbors( PHTelem, parent_neighbor );
        %update the Bezier extraction operators and element_nod indices        
        C_n = PHTelem(parent_neighbor).C;
        elmIn_n = PHTelem(parent_neighbor).nodes;
        knotU1n = PHTelem(parent_neighbor).vertex(1);
        knotU2n = PHTelem(parent_neighbor).vertex(3);
        knotV1n = PHTelem(parent_neighbor).vertex(2);
        knotV2n = PHTelem(parent_neighbor).vertex(4);
        
        [ Ce1t, Ce2t, Ce3t, Ce4t] = deCasteljau2dai( C_n, knotU1n, knotU2n, knotV1n, knotV2n, p, q, [RTjunctn, newBasisVertn, 5], knotUln, knotUrn, knotVdn, knotVun, elmIn_n, dimBasisTemp );
        
         switch ijunct
            case 1
                [newBasisSet, dimBasisTemp] = newBasisIndices( 3, elmIn_n, p, q, dimBasisTemp );
                down_neighbors = PHTelem(e).neighbor_down;
                PHTelem(down_neighbors(1)).C = Ce3t;
                PHTelem(down_neighbors(2)).C = Ce4t;
                PHTelem(down_neighbors(1)).nodes(ne) = newBasisSet{3}(ne);
                PHTelem(down_neighbors(1)).nodes(north) = newBasisSet{3}(north);
                PHTelem(down_neighbors(2)).nodes(nw) = newBasisSet{4}(nw);
                PHTelem(down_neighbors(2)).nodes(north) = newBasisSet{4}(north);
                
            case 2
                [newBasisSet, dimBasisTemp] = newBasisIndices( 4, elmIn_n, p, q, dimBasisTemp );
                right_neighbors = PHTelem(e).neighbor_right;
                PHTelem(right_neighbors(1)).C = Ce1t;
                PHTelem(right_neighbors(2)).C = Ce3t;
                PHTelem(right_neighbors(1)).nodes(nw) = newBasisSet{1}(nw);
                PHTelem(right_neighbors(1)).nodes(west) = newBasisSet{1}(west);
                PHTelem(right_neighbors(2)).nodes(sw) = newBasisSet{3}(sw);
                PHTelem(right_neighbors(2)).nodes(west) = newBasisSet{3}(west);
             case 3
                 [newBasisSet, dimBasisTemp] = newBasisIndices( 1, elmIn_n, p, q, dimBasisTemp );
                 up_neighbors = PHTelem(e).neighbor_up;
                 PHTelem(up_neighbors(1)).C = Ce1t;
                 PHTelem(up_neighbors(2)).C = Ce2t;
                 PHTelem(up_neighbors(1)).nodes(se) = newBasisSet{1}(se);
                 PHTelem(up_neighbors(1)).nodes(south) = newBasisSet{1}(south);
                 PHTelem(up_neighbors(2)).nodes(sw) = newBasisSet{2}(sw);
                 PHTelem(up_neighbors(2)).nodes(south) = newBasisSet{2}(south);
             case 4
                 [newBasisSet, dimBasisTemp] = newBasisIndices( 2, elmIn_n, p, q, dimBasisTemp );
                 left_neighbors = PHTelem(e).neighbor_left;
                 PHTelem(left_neighbors(1)).C = Ce2t;
                 PHTelem(left_neighbors(2)).C = Ce4t;
                 PHTelem(left_neighbors(1)).nodes(ne) = newBasisSet{2}(ne);
                 PHTelem(left_neighbors(1)).nodes(east) = newBasisSet{2}(east);
                 PHTelem(left_neighbors(2)).nodes(se) = newBasisSet{4}(se);
                 PHTelem(left_neighbors(2)).nodes(east) = newBasisSet{4}(east);
         end   
    end
         
    %calculate the new Bezier extraction operators and element_nod indices
    %of the children elements
    C_temp = PHTelem(e).C;
    elmIn = PHTelem(e).nodes;
        
    [ Ce1, Ce2, Ce3, Ce4, in1, in2, in3, in4, dimBasis] = deCasteljau2dai( C_temp, xmin, xmax, ymin, ymax, p, q, [RTjunct, newBasisVert, 5], knotUl, knotUr, knotVd, knotVu, elmIn, dimBasis);        
    
    PHTelem(lastElem+1).C = Ce1;
    PHTelem(lastElem+2).C = Ce2;
    PHTelem(lastElem+3).C = Ce3;
    PHTelem(lastElem+4).C = Ce4;
    PHTelem(lastElem+1).nodes = in1;
    PHTelem(lastElem+2).nodes = in2;
    PHTelem(lastElem+3).nodes = in3;
    PHTelem(lastElem+4).nodes = in4;
                   
    %update the neighbors of the neighbors with self
    for ichild = 1:4
        if (length(PHTelem(lastElem+ichild).neighbor_down)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_down).neighbor_up = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_down).neighbor_up, lastElem+ichild]),e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_right)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left, lastElem+ichild]),e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_up)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_up).neighbor_down = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_up).neighbor_down, lastElem+ichild]),e);
        end
        
        if (length(PHTelem(lastElem+ichild).neighbor_left)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right, lastElem+ichild]),e);
        end
    end                    
end