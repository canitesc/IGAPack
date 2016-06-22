dRdxi=[];
numGaussX = 4;
numGaussY = 4;
noGpEle = numGaussX*numGaussY;

%calculate the number of actual elements (i.e., non-refined, without children)
numElem = 0;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
        numElem = numElem+1;
    end
    end
end

noElems =  numElem;
x     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
u     = zeros(3,noGpEle,noElems);  % displacements of Gauss points
sigma = zeros(3,noGpEle,noElems);  % stresses      of Gauss points
errorU = zeros(2,noGpEle,noElems); % error in displacement field
errorStress = zeros(3,noGpEle,noElems); %error in stress field

id = 1; 
xx=zeros(noGpEle*noElems,2);

invC = inv(Cmat);

errDisp = zeros(2,1);
[gwx, gpx]=quadrature(numGaussX, 'GAUSS', 1);
[gwy, gpy]=quadrature(numGaussY, 'GAUSS', 1);

gpx=gpx';
gpy=gpy';

%define the 2D Bernstein polynomials
[B_u, dB_u] = bernstein_basis(gpx,p);
[B_v, dB_v] = bernstein_basis(gpy,q);
B_uv = zeros(numGaussX, numGaussY,(p+1)*(q+1));
dBdu = zeros(numGaussX, numGaussY,(p+1)*(q+1));
dBdv = zeros(numGaussX, numGaussY,(p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        B_uv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

e = 1;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
        xmin = PHTelem{patchIndex}(i).vertex(1);
        xmax = PHTelem{patchIndex}(i).vertex(3);
        ymin = PHTelem{patchIndex}(i).vertex(2);
        ymax = PHTelem{patchIndex}(i).vertex(4);
        
        
        %jacobian of the transformation from reference [-1,1]x[-1,1]
        %element to the local element in parameter
        scalefac = (xmax - xmin)*(ymax - ymin)/4;
        
                       
        J = zeros(numGaussY, numGaussX);        
        
        
        nument = size(PHTelem{patchIndex}(i).C,1); %number of basis functions with support on current knotspan                
       
        scrt = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
        scrt_x = 2*scrt-1;
        scrt_y = 2*scrt;
        dscrtx = reshape([2*scrt-1; 2*scrt],1,2*nument);      
        
        gp = 1;
                        
        for jj=1:numGaussY
            for ii=1:numGaussX
                
                [ coord, dxdxi]  = paramMap( GIFTmesh{patchIndex}, gpx(ii), gpy(jj), xmin, ymin, xmax, ymax);
                
                
                  [theta, r] = cart2pol(coord(1)-0.5,coord(2)-1);
%                     
%                     cracklength = 0.5;
%                     Kexact1 = KI;
%                     Kexact2 = KII;
%                     xCr = [0 1;cracklength 1];
%                     adv = xCr(2,:) - xCr(1,:);
%                     xTip = [cracklength 1];
%                     [disp_x, disp_y] = exact_Griffith3(coord,Emod,nu,stressState, Kexact1, Kexact2,xTip,adv);
%                     disp_ex = [disp_x, disp_y];
%                     
%                     stress_ex=exact_stresses2(r, theta, Kexact1, Kexact2);
                
               % disp_ex = holeu_d([coord(1), coord(2)], rad, Emod, nu, tx);
               % stress_ex = ghole(coord(1), coord(2), rad, tx);                      
                
                J(jj, ii) = det(dxdxi);
                cR = PHTelem{patchIndex}(i).C * squeeze(B_uv(ii,jj,:));
                
                dR(1,:) = PHTelem{patchIndex}(i).C * squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                dR(2,:) = PHTelem{patchIndex}(i).C * squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                
                % Solve for first derivatives in global coordinates
                dR = dxdxi\dR;
                
                B = zeros(2*nument,3);
                B(1:2:2*nument-1,1) = dR(1,:);
                B(2:2:2*nument,2) = dR(2,:);
                B(1:2:2*nument-1,3) = dR(2,:);
                B(2:2:2*nument,3) = dR(1,:);
                                
                %calculate displacement values                
                disp_x = cR'*sol0(scrt_x);
                disp_y = cR'*sol0(scrt_y);    
                
                %calculate the error in stress values
                stressvect = Cmat*B'*sol0(dscrtx);
                
                x(1,gp,e)     = coord(1); x(2,gp,e)     = coord(2);
                xx(id,:)        = coord;
                u(1,gp,e)     = disp_x; u(2,gp,e)     = disp_y;
                sigma(:,gp,e)   = stressvect;
%                 errorStress(:,gp,e) = stress_ex - stressvect;
%                 errDisp(1)= disp_ex(1)-disp_x; errDisp(2)= disp_ex(2)-disp_y;
%                 errorU(:,gp,e) = errDisp; 
                id = id + 1;    gp = gp +1 ;  
                
            end
        end
        
        e = e + 1;
                       
    end
end

end

%Energy_TSpline
msh_to_vtu (x, sigma, u, errorU,errorStress,[numGaussX numGaussY], vtuFile);