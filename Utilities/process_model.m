function xk_1 = process_model(xk, Ts)
    tauf = 200;
    taug = 300;
    ek = xk(1:3);
    q4k = xk(4);
    omegak = xk(5:7);
    bfk = xk(8:10);
    bgk = xk(14:16);
    
    xk_1 = xk;
    xk_1(1:3) = xk_1(1:3) + Ts/2*(-cross(omegak, ek) + q4k*omegak);
    xk_1(4) = xk_1(4) + Ts/2*(-omegak'*ek);
    xk_1(8:10) =  xk_1(8:10) + Ts*(-bfk/tauf);
    xk_1(14:16) =  xk_1(14:16) + Ts*(-bgk/taug);
    xk_1(1:4) = xk_1(1:4)/norm(xk_1(1:4));
    if xk_1(4) < 0
        xk_1(1:4) = -xk_1(1:4);
    end
end