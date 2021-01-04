function [ path ] = getexproot( )
%GETEXPROOT gets the root of the file

path = strrep(mfilename('fullpath'), mfilename(), '');

end

