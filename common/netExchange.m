function [idx, coef] = netExchange(rxns, expr)
% Reactions and coefficients of a flux written as one reaction id or as a signed sum,
% e.g. '+Liver_EX_glc_D(e)_[bc] -Liver_EX_glc_D[bpL]_[bp]'. A signed sum is used for an
% organ that exchanges a metabolite with more than one blood compartment, with each term
% signed so that the sum is the net release from the organ.
%
% INPUT:  rxns  reaction ids of the model
%         expr  one reaction id, or terms '+id' / '-id' separated by spaces
% OUTPUT: idx, coef  reaction indices and coefficients (empty if a reaction is missing)
%
% Anna Sheehy 2026.

expr = strtrim(char(expr));
if any(expr(1) == '+-')
    tok = regexp(expr, '([+-])(\S+)', 'tokens');
    ids = cellfun(@(t) t{2}, tok, 'UniformOutput', false);
    coef = cellfun(@(t) 1 - 2 * strcmp(t{1}, '-'), tok);
else
    ids = {expr};  coef = 1;
end
[ok, idx] = ismember(ids, rxns);
if ~all(ok), idx = [];  coef = []; end
idx = idx(:);  coef = coef(:);
end
