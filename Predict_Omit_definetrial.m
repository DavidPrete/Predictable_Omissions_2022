function trl = Predict_Omit_definetrial(cfg)

trl    = zeros(length(cfg.sample),3);
sample = cfg.sample;

% determine the number of samples before and after the trigger
pretrig  = round(cfg.prestim  * cfg.fsample);
posttrig = round(cfg.poststim * cfg.fsample);

for jj = 1:(length(sample))
    trl_begin = sample(jj) + pretrig;       
    trl_end   = sample(jj) + posttrig;       
    offset    = pretrig;
    trl(jj,:) = [trl_begin trl_end offset];
    
end