function [xlimits,ylimits] = queryXYLimits(axHandle)
xlimits = get(axHandle,'XLim');
ylimits = get(axHandle,'YLim');