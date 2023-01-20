function [W,Q] = quadrature( quadorder, qt, sdim )

% The function quadrature returns a n x 1 column vector W of quadrature
% weights and a n x dim matrix of quadrature points, where n is the
% number of quadrature points.  The function is called as follows:
%
%   [W,Q]=quadrature( nint, type, dim )
%
% nint is the quadrature order, type is the type of quadrature
% (i.e. gaussian, triangular, etc.. ) and dim is the number of spacial
% dimentions of the problem.  The default for type is GAUSS and the
% default for dim is unity.
%
% wrQ=quadrature(nint,'TRIANGULAR',2);itten by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering 
% Northwestern University


  if ( nargin < 3 )   % set default arguments
    if ( strcmp(qt,'GAUSS') == 1 )
      dim = 1;
    else
      dim = 2;
    end
  end

  if ( nargin < 2 )
    type = 'GAUSS';
  end

  if ( strcmp(qt,'GAUSS') == 1 ) 

    if ( quadorder > 16 )  % check for valid quadrature order
      disp('Order of quadrature too high for Gaussian Quadrature'); 
      quadorder =16;
    end
    
    quadpoint=zeros(quadorder^sdim ,sdim);
    quadweight=zeros(quadorder^sdim,1);
  
    r1pt=zeros(quadorder,1); r1wt=zeros(quadorder,1);

    switch ( quadorder ) 
      case 1
        r1pt(1) = 0.000000000000000;
        r1wt(1) = 2.000000000000000;

      case 2
        r1pt(1) = -0.577350269189626;
        r1pt(2) = 0.577350269189626;

        r1wt(1) = 1.000000000000000; 
        r1wt(2) = 1.000000000000000;         

      case 3
        r1pt(1) = 0.774596669241483;
        r1pt(2) =-0.774596669241483;
        r1pt(3) = 0.000000000000000;

        r1wt(1) = 0.555555555555556;
        r1wt(2) = 0.555555555555556; 
        r1wt(3) = 0.888888888888889;   

      case 4
        r1pt(1) = 0.861136311594053;
        r1pt(2) =-0.861136311594053;
        r1pt(3) = 0.339981043584856;
        r1pt(4) =-0.339981043584856;

        r1wt(1) = 0.347854845137454;
        r1wt(2) = 0.347854845137454; 
        r1wt(3) = 0.652145154862546;
        r1wt(4) = 0.652145154862546;  

      case 5
        r1pt(1) = 0.906179845938664;
        r1pt(2) =-0.906179845938664;
        r1pt(3) = 0.538469310105683;
        r1pt(4) =-0.538469310105683;
        r1pt(5) = 0.000000000000000;

        r1wt(1) = 0.236926885056189;
        r1wt(2) = 0.236926885056189;
        r1wt(3) = 0.478628670499366;
        r1wt(4) = 0.478628670499366;  
        r1wt(5) = 0.568888888888889;  

      case 6
        r1pt(1) = 0.932469514203152;
        r1pt(2) =-0.932469514203152;
        r1pt(3) = 0.661209386466265;
        r1pt(4) =-0.661209386466265;
        r1pt(5) = 0.238619186003152;
        r1pt(6) =-0.238619186003152;

        r1wt(1) = 0.171324492379170;
        r1wt(2) = 0.171324492379170;
        r1wt(3) = 0.360761573048139;
        r1wt(4) = 0.360761573048139;   
        r1wt(5) = 0.467913934572691; 
        r1wt(6) = 0.467913934572691;
	
      case 7
        r1pt(1) =  0.949107912342759;
        r1pt(2) = -0.949107912342759;
        r1pt(3) =  0.741531185599394;
        r1pt(4) = -0.741531185599394;
        r1pt(5) =  0.405845151377397;
        r1pt(6) = -0.405845151377397;
        r1pt(7) =  0.000000000000000;

        r1wt(1) = 0.129484966168870;
        r1wt(2) = 0.129484966168870;
        r1wt(3) = 0.279705391489277;
        r1wt(4) = 0.279705391489277;
        r1wt(5) = 0.381830050505119;
        r1wt(6) = 0.381830050505119;
        r1wt(7) = 0.417959183673469;

      case 8
        r1pt(1) =  0.960289856497536;
        r1pt(2) = -0.960289856497536;
        r1pt(3) =  0.796666477413627;
        r1pt(4) = -0.796666477413627;
        r1pt(5) =  0.525532409916329;
        r1pt(6) = -0.525532409916329;
        r1pt(7) =  0.183434642495650;
        r1pt(8) = -0.183434642495650;

        r1wt(1) = 0.101228536290376;
        r1wt(2) = 0.101228536290376;
        r1wt(3) = 0.222381034453374;
        r1wt(4) = 0.222381034453374;
        r1wt(5) = 0.313706645877887;
        r1wt(6) = 0.313706645877887;
        r1wt(7) = 0.362683783378362;
        r1wt(8) = 0.362683783378362;
        
        case 9
            r1pt  = [0.9681602395076260898355762, ...
            -0.9681602395076260898355762,...
             0.8360311073266357942994298, ...
            -0.8360311073266357942994298,...
             0.6133714327005903973087020, ...
            -0.6133714327005903973087020, ...
             0.3242534234038089290385380, ...
            -0.3242534234038089290385380, ...
                0]';
           
          
            r1wt = [0.0812743883615744119718922, ...
            0.0812743883615744119718922, ...
            0.1806481606948574040584720, ...
            0.1806481606948574040584720, ...
            0.2606106964029354623187429, ...
            0.2606106964029354623187429, ...
            0.3123470770400028400686304, ...
            0.3123470770400028400686304, ...
            0.3302393550012597631645251]';
                        
            
            
        case 10
            r1pt = [0.9739065285171717200779640, ...
            -0.9739065285171717200779640, ...
             0.8650633666889845107320967, ...
            -0.8650633666889845107320967, ...
            0.6794095682990244062343274, ...
            -0.6794095682990244062343274, ...
            0.4333953941292471907992659, ...
             -0.4333953941292471907992659, ...
            0.1488743389816312108848260, ...
              -0.1488743389816312108848260]';
             
             
            r1wt = [0.0666713443086881375935688, ...
             0.0666713443086881375935688, ...
             0.1494513491505805931457763, ...
             0.1494513491505805931457763, ...
             0.2190863625159820439955349, ...
             0.2190863625159820439955349, ...
             0.2692667193099963550912269, ...
             0.2692667193099963550912269, ...             
             0.2955242247147528701738930, ...
             0.2955242247147528701738930]';
         
        case 11
            
            r1pt = [0.9782286581460570,...
            -0.9782286581460570, ...
            0.8870625997680953, ...
            -0.8870625997680953, ...
            0.7301520055740494, ...
            -0.7301520055740494, ...
            0.5190961292068118 ,...
            -0.5190961292068118, ...
            0.2695431559523450, ...
            -0.2695431559523450, ...
            0.0000000000000000]';
        
            r1wt = [0.0556685671161737, ...
            0.0556685671161737, ...
            0.1255803694649046, ...
            0.1255803694649046, ...
            0.1862902109277343, ...
            0.1862902109277343, ...
            0.2331937645919905, ...
            0.2331937645919905, ...
            0.2628045445102467, ...
            0.2628045445102467, ...
            0.2729250867779006]';
        
        case 12
            r1pt = [0.9815606342467192, ...
                -0.9815606342467192, ...
                0.9041172563704749, ...
                -0.9041172563704749, ...
                0.7699026741943047, ...
                -0.7699026741943047, ...
                0.5873179542866175, ...
                -0.5873179542866175, ...
                0.3678314989981802, ...
                -0.3678314989981802, ...
                0.1252334085114689, ...
                -0.1252334085114689]';
            
            r1wt = [0.0471753363865118	, ...
                0.0471753363865118, ...
                0.1069393259953184, ...
                0.1069393259953184, ...
                0.1600783285433462, ...
                0.1600783285433462, ...
                0.2031674267230659, ...
                0.2031674267230659, ...
                0.2334925365383548, ...
                0.2334925365383548, ...
                0.2491470458134028, ...
                0.2491470458134028]';
                
        case 13
            
            r1pt = [0.9841830547185881,...
                -0.9841830547185881,...
                0.9175983992229779, ... 
                -0.9175983992229779, ...
                0.8015780907333099, ...
                -0.8015780907333099, ...
                0.6423493394403402, ...
                -0.6423493394403402, ...
                0.4484927510364469, ...
                -0.4484927510364469, ...
                0.2304583159551348, ...
                -0.2304583159551348, ...
                0.0000000000000000]';
            
            r1wt = [0.0404840047653159, ...
                0.0404840047653159, ...
                0.0921214998377285, ...
                0.0921214998377285, ...
                0.1388735102197872, ...
                0.1388735102197872, ...
                0.1781459807619457, ...
                0.1781459807619457, ...
                0.2078160475368885, ...
                0.2078160475368885, ...
                0.2262831802628972, ...
                0.2262831802628972, ...
                0.2325515532308739]';
            
            
        case 14
            
            r1pt = [0.9862838086968123,...
                -0.9862838086968123, ...
                0.9284348836635735, ...
                -0.9284348836635735, ...
                0.8272013150697650, ...
                -0.8272013150697650, ...
                0.6872929048116855, ...
                -0.6872929048116855, ...
                0.5152486363581541, ...
                -0.5152486363581541, ...
                0.3191123689278897, ...
                -0.3191123689278897, ...
                0.1080549487073437, ...
                -0.1080549487073437]';
            
            r1wt = [0.0351194603317519, ...
                0.0351194603317519, ...
                0.0801580871597602, ...
                0.0801580871597602, ...
                0.1215185706879032, ...
                0.1215185706879032, ...
                0.1572031671581935, ...
                0.1572031671581935, ...
                0.1855383974779378, ...
                0.1855383974779378, ...
                0.2051984637212956, ...
                0.2051984637212956, ...
                0.2152638534631578, ...
                0.2152638534631578]';
               
            
        case 15
            r1pt = [0.9879925180204854, ...
                -0.9879925180204854, ...
                0.9372733924007060, ...
                -0.9372733924007060, ...
                0.8482065834104272, ...
                -0.8482065834104272, ...
                0.7244177313601701, ...
                -0.7244177313601701, ...
                0.5709721726085388, ...
                -0.5709721726085388, ...
                0.3941513470775634, ...
                -0.3941513470775634, ...
                0.2011940939974345, ...
                -0.2011940939974345, ...
                0.0000000000000000]';
            
            r1wt = [0.0307532419961173, ...
                0.0307532419961173, ...
                0.0703660474881081, ...
                0.0703660474881081, ...
                0.1071592204671719, ...
                0.1071592204671719, ...
                0.1395706779261543, ...
                0.1395706779261543, ...
                0.1662692058169939, ...
                0.1662692058169939, ...
                0.1861610000155622, ...
                0.1861610000155622, ...
                0.1984314853271116, ...
                0.1984314853271116, ...
                0.2025782419255613]';
            
        case 16
            r1pt = [0.9894009349916499, ...
                -0.9894009349916499, ...
                0.9445750230732326, ...
                -0.9445750230732326, ...
                0.8656312023878318, ...
                -0.8656312023878318, ...
                0.7554044083550030, ...
                -0.7554044083550030, ...
                0.6178762444026438, ...
                -0.6178762444026438, ...
                0.4580167776572274, ...
                -0.4580167776572274, ...
                0.2816035507792589, ...
                -0.2816035507792589, ...
                0.0950125098376374, ...
                -0.0950125098376374]';
            
            r1wt = [0.0271524594117541, ...
                0.0271524594117541, ...
                0.0622535239386479, ...
                0.0622535239386479, ...
                0.0951585116824928, ...
                0.0951585116824928, ...
                0.1246289712555339, ...
                0.1246289712555339, ...
                0.1495959888165767, ...
                0.1495959888165767, ...
                0.1691565193950025, ...
                0.1691565193950025, ...
                0.1826034150449236, ...
                0.1826034150449236, ...
                0.1894506104550685, ...
                0.1894506104550685]';
            
      otherwise
        disp('Order of quadrature to high for Gaussian Quadrature'); 
	
    end  % end of quadorder switch

    n=1;
     
    if ( sdim == 1 ) 
      for i = 1:quadorder
        quadpoint(n,:) = [ r1pt(i) ];           
        quadweight(n) = r1wt(i); 
        n = n+1;
      end
    
    elseif ( sdim == 2 ) 
      for i = 1:quadorder
        for j = 1:quadorder
          quadpoint(n,:) = [ r1pt(i), r1pt(j)];           
          quadweight(n) = r1wt(i)*r1wt(j); 
          n = n+1;
        end
      end
  
    else % sdim == 3
      for i = 1:quadorder
        for j = 1:quadorder
          for k = 1:quadorder
            quadpoint(n,:) = [ r1pt(i), r1pt(j), r1pt(k) ];           
            quadweight(n) = r1wt(i)*r1wt(j)*r1wt(k); 
            n = n+1;
          end
	end
      end
      
    end
    
    Q=quadpoint;
    W=quadweight;
  % END OF GAUSSIAN QUADRATURE DEFINITION
  
  elseif ( strcmp(qt,'TRIANGULAR') == 1 ) 
    
    if ( sdim == 3 )  %%% TETRAHEDRA
      
      if ( quadorder ~= 1 &  quadorder ~= 2 &  quadorder ~= 3  ) 
        % check for valid quadrature order
        disp('Incorect quadrature order for triangular quadrature');
        quadorder = 1;
      end
      
      if  ( quadorder == 1 )
        quadpoint = [ 0.25 0.25 0.25 ];
        quadweight = 1;
        
      elseif ( quadorder == 2 ) 
        quadpoint = [ 0.58541020  0.13819660  0.13819660;
                      0.13819660  0.58541020  0.13819660;
                      0.13819660  0.13819660  0.58541020;
                      0.13819660  0.13819660  0.13819660];
        quadweight = [1; 1; 1; 1]/4;
        
      elseif ( quadorder == 3 ) 
        quadpoint = [ 0.25  0.25  0.25;
                      1/2   1/6   1/6;
                      1/6   1/2   1/6;
                      1/6   1/6   1/2;
                      1/6   1/6   1/6];
        quadweight = [-4/5 9/20 9/20 9/20 9/20]';
        
      end
      
      Q=quadpoint;
      W=quadweight/6;
         
    else  %%% TRIANGLES
      
      if ( quadorder > 7 ) % check for valid quadrature order
        disp('Quadrature order too high for triangular quadrature');
        quadorder = 1;
      end
      
      if ( quadorder == 1 )   % set quad points and quadweights
        quadpoint = [ 0.3333333333333, 0.3333333333333 ];
        quadweight = 1;
        
      elseif ( quadorder == 2 ) 
        quadpoint = zeros( 3, 2 );
        quadweight = zeros( 3, 1 );
        
        quadpoint(1,:) = [ 0.1666666666667, 0.1666666666667 ];
        quadpoint(2,:) = [ 0.6666666666667, 0.1666666666667 ];
        quadpoint(3,:) = [ 0.1666666666667, 0.6666666666667 ]; 
        
        quadweight(1) = 0.3333333333333; 
        quadweight(2) = 0.3333333333333; 
        quadweight(3) = 0.3333333333333;   
        
      elseif ( quadorder <= 5 ) 
        quadpoint = zeros( 7, 2 );
        quadweight = zeros( 7, 1 );
        
        quadpoint(1,:) = [ 0.1012865073235, 0.1012865073235 ];
        quadpoint(2,:) = [ 0.7974269853531, 0.1012865073235 ];
        quadpoint(3,:) = [ 0.1012865073235, 0.7974269853531 ]; 
        quadpoint(4,:) = [ 0.4701420641051, 0.0597158717898 ];
        quadpoint(5,:) = [ 0.4701420641051, 0.4701420641051 ];
        quadpoint(6,:) = [ 0.0597158717898, 0.4701420641051 ]; 
        quadpoint(7,:) = [ 0.3333333333333, 0.3333333333333 ];
        
        quadweight(1) = 0.1259391805448; 
        quadweight(2) = 0.1259391805448; 
        quadweight(3) = 0.1259391805448; 
        quadweight(4) = 0.1323941527885;
        quadweight(5) = 0.1323941527885;
        quadweight(6) = 0.1323941527885;
        quadweight(7) = 0.2250000000000;  
        
      else
        quadpoint = zeros( 13, 2 );
        quadweight = zeros( 13, 1 );
        
        quadpoint(1 ,:) = [ 0.0651301029022, 0.0651301029022 ];
        quadpoint(2 ,:) = [ 0.8697397941956, 0.0651301029022 ];
        quadpoint(3 ,:) = [ 0.0651301029022, 0.8697397941956 ];
        quadpoint(4 ,:) = [ 0.3128654960049, 0.0486903154253 ];
        quadpoint(5 ,:) = [ 0.6384441885698, 0.3128654960049 ];
        quadpoint(6 ,:) = [ 0.0486903154253, 0.6384441885698 ];
        quadpoint(7 ,:) = [ 0.6384441885698, 0.0486903154253 ];
        quadpoint(8 ,:) = [ 0.3128654960049, 0.6384441885698 ];
        quadpoint(9 ,:) = [ 0.0486903154253, 0.3128654960049 ];
        quadpoint(10,:) = [ 0.2603459660790, 0.2603459660790 ];
        quadpoint(11,:) = [ 0.4793080678419, 0.2603459660790 ];
        quadpoint(12,:) = [ 0.2603459660790, 0.4793080678419 ];
        quadpoint(13,:) = [ 0.3333333333333, 0.3333333333333 ];
        
        quadweight(1 ) = 0.0533472356088;
        quadweight(2 ) = 0.0533472356088; 
        quadweight(3 ) = 0.0533472356088;
        quadweight(4 ) = 0.0771137608903;
        quadweight(5 ) = 0.0771137608903;
        quadweight(6 ) = 0.0771137608903;
        quadweight(7 ) = 0.0771137608903;
        quadweight(8 ) = 0.0771137608903;
        quadweight(9 ) = 0.0771137608903;
        quadweight(10) = 0.1756152576332; 
        quadweight(11) = 0.1756152576332; 
        quadweight(12) = 0.1756152576332;
        quadweight(13) =-0.1495700444677; 
        
      end
      
      Q=quadpoint;
      W=quadweight/2;   % ATTENTION ATTENTION WHY DIVIDE TO 2?????
    end
    
  end  % end of TRIANGULAR initialization
  
% END OF FUNCTION
