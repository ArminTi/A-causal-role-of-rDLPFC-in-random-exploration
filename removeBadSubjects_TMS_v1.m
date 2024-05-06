function [sub, sub_bad] = removeBadSubjects_TMS_v1(sub)

% [sub, sub_bad] = removeBadSubjects_E1_v1(sub)

for sn = 1:length(sub)
    dum = sub(sn).co(:,:);
    sub(sn).fC = nanmean(dum(:));
end
thresh = 0.55;
ind_good = [sub.fC]>=thresh;
sub_bad = sub(~ind_good);
sub = sub(ind_good);