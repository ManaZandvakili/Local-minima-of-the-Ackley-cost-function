function a = policy_action(obs, NN)
%POLICY_ACTION  The trained network: normalise the observation, two tanh layers, linear output.
%   a = policy_action(obs, NN), 6 numbers in [-1, 1] (the deterministic action of the PPO
%   policy). NN from data/controllers/<residual|plain>_<case>.mat (struct net).
%   Computed in single precision, like PyTorch.
%#codegen

o = (obs - NN.obs_mean) ./ sqrt(NN.obs_var + NN.epsilon);
o = single(min(max(o, -NN.clip_obs), NN.clip_obs));
h = tanh(single(NN.W1) * o + single(NN.b1));
h = tanh(single(NN.W2) * h + single(NN.b2));
a = single(NN.W3) * h + single(NN.b3);
a = min(max(double(a), -1), 1);
end
