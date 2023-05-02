function candidateName = get_name(mName)
%SAVE_RESULTS Gets the results file name based on the owner script and
%availability

[path, name, ~] = fileparts(mName);
if ~exist('results', 'dir')
    mkdir('results')
end
candidateName = fullfile(path, 'results', [name '_0.mat']);
i = 1;
while exist(candidateName, "file")
    candidateName = fullfile(path, 'results', [name sprintf('_%d.mat', i)]);
    i = i + 1;
end

end

