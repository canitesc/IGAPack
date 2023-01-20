function [ PHTelem, dimBasis ] = crossInsert1D( PHTelem, elm_index, dimBasis, p)
%inserts cross in PHTelem at elm_index
%dimBasis: dimension of the PHT space

for e=elm_index    
    lastElem = length(PHTelem);   
     
    %update parent element                 
    PHTelem(e).children = lastElem+1:lastElem+2;                    
    
    xmin = PHTelem(e).vertex(1);
    xmax = PHTelem(e).vertex(2);

    
    %add children elements
    %left child
    PHTelem(lastElem+1).parent = e;
    PHTelem(lastElem+1).children = [];
    PHTelem(lastElem+1).vertex = [xmin, (xmin+xmax)/2];
    PHTelem(lastElem+1).level = PHTelem(e).level+1;    
    PHTelem(lastElem+1).neighbor_right = lastElem+2;    
    
    %right child
    PHTelem(lastElem+2).parent = e;
    PHTelem(lastElem+2).children = [];
    PHTelem(lastElem+2).vertex = [(xmin+xmax)/2, xmax];
    PHTelem(lastElem+2).level = PHTelem(e).level+1;    
    PHTelem(lastElem+2).neighbor_left = lastElem+1;
        
    %add the neighbors outside the refined element           
    if ~isempty(PHTelem(e).neighbor_right)
        PHTelem(lastElem+2).neighbor_right = PHTelem(e).neighbor_right;        
    end
            
    if ~isempty(PHTelem(e).neighbor_left)
        PHTelem(lastElem+1).neighbor_left = PHTelem(e).neighbor_left;            
    end
    
         
    %calculate the new Bezier extraction operators and element_nod indices
    %of the children elements
    C_temp = PHTelem(e).C;
    elmIn = PHTelem(e).nodes;   
        
    [ Ce1, Ce2, in1, in2, dimBasis] = deCasteljau1dai( C_temp, xmin, xmax, p, elmIn, dimBasis);            
    PHTelem(lastElem+1).C = Ce1;
    PHTelem(lastElem+2).C = Ce2;
  
    PHTelem(lastElem+1).nodes = in1;
    PHTelem(lastElem+2).nodes = in2;
                   
    %update the neighbors of the neighbors with self
    for ichild = 1:2
               
        if (length(PHTelem(lastElem+ichild).neighbor_right)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_right).neighbor_left, lastElem+ichild]),e);
        end                
        
        if (length(PHTelem(lastElem+ichild).neighbor_left)==1)
            PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right = setdiff(unique([PHTelem(PHTelem(lastElem+ichild).neighbor_left).neighbor_right, lastElem+ichild]),e);
        end
    end                    
end