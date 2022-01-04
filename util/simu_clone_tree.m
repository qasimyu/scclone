function [tree, muta_assignments, clone_states] = simu_clone_tree(num_clone, num_muta)
% simulate a subclonal tree
% inputs: num_clone, number of subclones
%         num_muta, number of mutations
% outputs: tree, subclonal tree
%          muta_assignments, mutation assignments to edges of the tree
%          clone_states, mutation profiles of subclones

tree = zeros(1, num_clone);
if num_clone > 1
    tree(2) = 1;
    for n = 3:num_clone
        k = randperm(n-2, 1);
        tree(n) = k+1;
    end
end

muta_assignments = randsrc(1,num_muta,[2:num_clone; ones(1,num_clone-1)/(num_clone-1)]);
clone_states = get_clone_states(tree, muta_assignments);

end

function clone_states = get_clone_states(tree, muta_assignments)

num_clone = length(tree);
num_muta = length(muta_assignments);
clone_states = zeros(num_clone, num_muta);

visited = zeros(1,num_clone);
visited(1) = 1; %root

for i = 2:num_clone
    if visited(i) == 1
        continue;
    end
    p = tree(i); %parent node
    if p == 0
        continue;
    end
    nodes = i;
    while p ~= 0 % not root
        nodes = [nodes p];
        if visited(p) == 1 % the state of node p has been evaluated
            break;
        end
        visited(p) = 1;
        p = tree(p);
    end
    k = length(nodes);
    for n = k-1:-1:1
        indx = nodes(n);
        indx_p = nodes(n+1);
        clone_states(indx,:) = clone_states(indx_p,:);
        for j = 1:num_muta
            if muta_assignments(j) == indx
                clone_states(indx, j) = 1;
            end
        end
    end
end

end