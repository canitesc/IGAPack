function index = findVertex(vertices, vertex)
%     Find if a vertex exist in the vertices list up to tol_eq tolerance
% 
%     Parameters
%     ----------
%     vertices : (2D array with 3 columns)
%         list of vertices
%     vertex : (1D array of length 3)
%         vertex to check 
%     tol_eq : float, optional
%         Tolerance for equality. The default is 1e-10.
% 
%     Returns
%     -------
%     index : (int) index of vertex if it exists in vertices, None otherwise
tolEq = 1e-10;
index = NaN;
for i = 1:size(vertices,1)
    testArr = zeros(1,3, 'logical');
    for j = 1:3
        testArr(j) = (abs(vertices(i,j)-vertex(j))<tolEq);
    end
    if all(testArr)
        index = i;
        break
    end
end

