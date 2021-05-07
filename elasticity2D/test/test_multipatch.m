function tests = test_multipatch
    tests = functiontests(localfunctions);
end

function test_corners(testCase)
    deg = [3,2];
    index_dl = getCornerNode2D(deg, 'down_left');
    index_dr = getCornerNode2D(deg, 'down_right');
    index_ul = getCornerNode2D(deg, 'up_left');
    index_ur = getCornerNode2D(deg, 'up_right');
    verifyEqual(testCase, 1, index_dl)
    verifyEqual(testCase, 4, index_dr)
    verifyEqual(testCase, 9, index_ul)
    verifyEqual(testCase, 12, index_ur)
end

function test_fv(testCase)
    vertices = [-0.5       ,  0.        ,  0.     
               -0.5       ,  0.04587369,  0.       
               -0.48280077,  0.13881616,  0.       
                0.        ,  0.5       ,  0.        
               -0.9395114 ,  0.11542752,  0.        
               -0.9405877 ,  0.34786882,  0.        ];
    indx1 = findVertex(vertices, [0, 0, 0]);
    indx2 = findVertex(vertices, [0, 0.5, 0]);
        
    verifyEqual(testCase, NaN, indx1)
    verifyEqual(testCase, 4, indx2)    
end
