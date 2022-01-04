function [E, cell_assignments] = simu_mutation_data(num_cell, clone_states)
% simulate mutation data given a subclonal tree
% inputs: num_cell, number of cells
%         clone_states, mutation profiles of subclones
% outputs: E, mutation profiles of single cells
%          cell_assignments, cell assignments to subclones

num_clone = size(clone_states, 1);
if num_cell < num_clone
    disp('number of cells should be larger than number of subclones');
    E = [];
    cell_assignments = [];
else
    cell_assignments = 1:num_clone;
    for i = num_clone+1:num_cell
        k = randsrc(1,1,[cell_assignments; ones(1,length(cell_assignments))/length(cell_assignments)]);
        cell_assignments = [cell_assignments k];
    end
    E = clone_states(cell_assignments,:);
end

end