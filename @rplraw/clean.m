function [obj, varargout] = clean(obj,varargin)
%@rplraw/clean Remove line noise in rplraw object.
%   OBJ = clean(OBJ) replaces the data in the object with line noise removed

Args = struct('LineNoiseFreq',50);
Args.flags = {};
[Args,varargin2] = getOptArgs(varargin,Args);

obj.data.analogData = nptRemoveLineNoise(obj.data.analogData,Args.LineNoiseFreq,obj.data.analogInfo.SampleRate);
