function Sx = skew_symmetric_matrix(x)
    Sx = [ 0    -x(3)  x(2);
           x(3)  0    -x(1);
          -x(2)  x(1)  0];
end

