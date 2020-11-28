function [px_rs,py_rs,ai_rs,li_rs,Li_rs] = reshape_state(z)
global n;
sz = [n,n,9*n,9*n,9*9*n];
en = 0;

st = en+1; en = sum(sz(1));
px_rs = z(st:en,1);

st = en+1; en = sum(sz(1:2));
py_rs = z(st:en,1);

st = en+1; en = sum(sz(1:3));
ai_rs = reshape(z(st:en,1),9,n);

st = en+1; en = sum(sz(1:4));
li_rs = reshape(z(st:en,1),9,n);

st = en+1; en = sum(sz(1:5));
Li_rs = reshape(z(st:en,1),9,9,n);

% all(x0==px_rs)
% all(y0==py_rs)
% all(ai==ai_rs)
% all(li==li_rs)
% all(Li==Li_rs)
end

