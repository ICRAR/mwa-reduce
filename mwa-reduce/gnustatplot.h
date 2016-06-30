#ifndef GNUSTATPLOT_H
#define GNUSTATPLOT_H

#include "uvector.h"
#include "gnuplot.h"
#include <vector>

class GNUStatPlot
{
public:
	GNUStatPlot() :
		_drawMinLine(false), _drawMaxLine(false),
		_drawUpperStddev(true), _drawLowerStddev(true)
		{ }
	bool HasXValues() const { return !_xValues.empty(); }
	void SetXValues(const ao::uvector<double>& xValues)
	{
		_xValues = xValues;
		_values.assign(xValues.size(), ao::uvector<double>());
	}
	void AddYSet(const ao::uvector<double>& yValues)
	{
		for(size_t xi=0; xi!=yValues.size(); ++xi)
		{
			_values[xi].push_back(yValues[xi]);
		}
	}
	void SetDrawExtremes(bool drawMinLine, bool drawMaxLine)
	{
		_drawMinLine = drawMinLine;
		_drawMaxLine = drawMaxLine;
	}
	void SetDrawStdDev(bool drawUpperStddev, bool drawLowerStddev)
	{
		_drawUpperStddev = drawUpperStddev;
		_drawLowerStddev = drawLowerStddev;
	}
	void Plot(GNUPlot& plot)
	{
		const std::string filenamePrefix = plot.FilenamePrefix();
		GNUPlot::Line
			*meanLine = plot.AddLine(filenamePrefix+"-mean.txt", "Mean"),
			*q1Line = plot.AddLine(filenamePrefix+"-q1.txt", "Quartile1"),
			*medLine = plot.AddLine(filenamePrefix+"-med.txt", "Median"),
			*q3Line = plot.AddLine(filenamePrefix+"-q3.txt", "Quartile3"),
			*minLine = 0, *maxLine = 0,
			*stddev1Line = 0, *stddev2Line = 0;
		if(_drawUpperStddev)
			stddev1Line = plot.AddLine(filenamePrefix+"stddev1.txt", "+stddev");
		if(_drawLowerStddev)
			stddev2Line = plot.AddLine(filenamePrefix+"stddev2.txt", "-stddev");
		if(_drawMinLine)
			minLine = plot.AddLine(filenamePrefix+"-min.txt", "Minimum");
		if(_drawMaxLine)
			maxLine = plot.AddLine(filenamePrefix+"-max.txt", "Maximum");
		for(size_t xi=0; xi!=_values.size(); ++xi)
		{
			ao::uvector<double> val;
			for(size_t pi=0; pi!=_values[xi].size(); ++pi)
			{
				if(std::isfinite(_values[xi][pi]))
					val.push_back(_values[xi][pi]);
			}
			
			if(!val.empty())
			{
				std::sort(val.begin(), val.end());
				double minVal = val.front(), maxVal = val.back();
				double medVal;
				if(val.size()%2 == 0)
					medVal = (val[val.size()/2-1] + val[val.size()/2])*0.5;
				else
					medVal = val[val.size()/2];
				double q1val = val[val.size()/4], q3val = val[val.size()*3/4];
				double sum = 0.0;
				for(size_t i=0; i!=val.size(); ++i)
					sum += val[i];
				double mean = sum/val.size();
				double stddev = 0.0;
				for(size_t i=0; i!=val.size(); ++i)
					stddev += (val[i] - mean) * (val[i] - mean);
				stddev = sqrt(stddev / val.size());
				
				meanLine->AddPoint(_xValues[xi], mean);
				medLine->AddPoint(_xValues[xi], medVal);
				q1Line->AddPoint(_xValues[xi], q1val);
				q3Line->AddPoint(_xValues[xi], q3val);
				if(_drawUpperStddev)
					stddev1Line->AddPoint(_xValues[xi], mean+stddev);
				if(_drawLowerStddev)
					stddev2Line->AddPoint(_xValues[xi], mean-stddev);
				if(_drawMinLine)
					minLine->AddPoint(_xValues[xi], minVal);
				if(_drawMaxLine)
					maxLine->AddPoint(_xValues[xi], maxVal);
			}
		}
	}
	
private:
	bool _drawMinLine, _drawMaxLine, _drawUpperStddev, _drawLowerStddev;

	ao::uvector<double> _xValues;
	std::vector<ao::uvector<double>> _values;
};

#endif
