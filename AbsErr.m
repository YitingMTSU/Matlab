function [error] = AbsErr(A,B)
ERR = abs(A - B);
error = max(max(ERR));
