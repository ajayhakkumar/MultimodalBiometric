function [alpha,C,fs]=splitsequence(bseq)

global n1 n2;

seq1 = sprintf('%d',bseq(1:n1));
n22 = n1+n2;

seq2 = sprintf('%d',bseq(n1+1:n22));
seq3 = bseq(n22+1:end);

alphas =[0.1 1];
dval=bin2dec(seq1);
alpha = (alphas(2)-alphas(1))*( dval/(2^n1-1) ) + alphas(1);

cmx =[0.1 10];
dval= bin2dec(seq2)+1;
C = (cmx(2)-cmx(1))*( dval/(2^n2-1) ) + cmx(1);

fs = logical(seq3);

