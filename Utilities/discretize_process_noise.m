function Qd = discretize_process_noise(A, Q, Ts)
    num_states = size(A, 1);
    tempMat = expm(Ts*[-A Q; zeros(num_states) A']);
    Ad = tempMat(num_states+1:num_states*2, num_states+1:num_states*2)';
    Qd = Ad*tempMat(1:num_states, num_states+1:num_states*2);
end

