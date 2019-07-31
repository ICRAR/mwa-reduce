#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

class GNUPlot
{
public:
	class Line
	{
	public:
		friend class GNUPlot;
		void AddPoint(double x, double y)
		{
			_file.precision(16);
			_file << x << '\t' << y << '\n';
		}
		void AddPoint(double x, double y, double yErr)
		{
			_file.precision(16);
			_file << x << '\t' << y << '\t' << yErr << '\n';
		}
	private:
		explicit Line(const std::string& filename) : _file(filename)
		{ }
		
		std::ofstream _file;
	};
	
	GNUPlot(const std::string& filenamePrefix, const std::string& xLabel, const std::string& yLabel, bool xLog=false, bool yLog=false) : _file(filenamePrefix + ".plt"), _filenamePrefix(filenamePrefix)
	{
		_file <<
			"set terminal pdfcairo enhanced color font 'Helvetica,16'\n";
		if(xLog && yLog)
			_file << "set logscale xy\n";
		else if(xLog)
			_file << "set logscale x\n";
		else if(yLog)
			_file << "set logscale y\n";
		_file <<
			"#set yrange [0.1:]\n"
			"set output \"" << filenamePrefix << ".ps\"\n"
			"#set key bottom left\n"
			"set xlabel \"" << xLabel << "\"\n"
			"set ylabel \"" << yLabel << "\"\n";
	}
	
	void SetTitle(const std::string& title)
	{
		_file << "set title \"" << title << "\"\n";
	}
	
	void SetXRange(double min, double max)
	{
		_file << "set xrange [" << min << ":" << max << "]\n";
	}
	void SetXRangeMinimum(double min)
	{
		_file << "set xrange [" << min << ":]\n";
	}
	void SetXRangeMaximum(double max)
	{
		_file << "set xrange [:" << max << "]\n";
	}
	
	void SetYRange(double min, double max)
	{
		_file << "set yrange [" << min << ":" << max << "]\n";
	}
	void SetYRangeMinimum(double min)
	{
		_file << "set yrange [" << min << ":]\n";
	}
	void SetYRangeMaximum(double max)
	{
		_file << "set yrange [:" << max << "]\n";
	}
	
	void SetKeyBottomLeft()
	{
		_file << "set key bottom left\n";
	}
	
	~GNUPlot()
	{
		_file << '\n';
		for(std::vector<Line*>::iterator i=_lines.begin(); i!=_lines.end(); ++i)
			delete *i;
	}
	
	Line* AddLine(const std::string& filename, const std::string& caption, double width=2.0)
	{
		if(_lines.empty())
			_file << "plot \\\n";
		else
			_file << ",\\\n";
		_file << "\"" << filename << "\" using 1:2 with lines title '" << caption << "' lw " << width;
		_lines.push_back(new Line(filename));
		return _lines.back();
	}
	
	Line* AddErrorSet(const std::string& filename, const std::string& caption, double sigmaLevel=1.0, double width=1.0)
	{
		if(_lines.empty())
			_file << "plot \\\n";
		else
			_file << ",\\\n";
		_file << "\"" << filename << "\" using 1:2:";
		if(sigmaLevel==1.0)
			_file << "3";
		else
			_file << "(column(3)*" << sigmaLevel << ")";
		_file << " with errorbars title '' lc rgb \"#000000\" lw " << width;
		_file << ",\\\n\"" << filename << "\" using 1:2 with lines title '" << caption << "' lw " << width;
		_lines.push_back(new Line(filename));
		return _lines.back();
	}
	
	void AddLineFromExistingFile(const std::string& filename, const std::string& caption, double width=2.0)
	{
		if(_lines.empty())
			_file << "plot \\\n";
		else
			_file << ",\\\n";
		_file << "\"" << filename << "\" using 1:2 with lines title '" << caption << "' lw " << width;
		_lines.push_back(0);
	}
	
	void AddFunction(const std::string& funcStr, const std::string& caption = "")
	{
		if(_lines.empty())
			_file << "plot \\\n";
		else
			_file << ",\\\n";
		_file << funcStr << " with lines lw 1.0 title '" << caption << "'";
		_lines.push_back(0);
	}
	
	Line* AddPointSet(const std::string& filename, const std::string& caption)
	{
		if(_lines.empty())
			_file << "plot \\\n";
		else
			_file << ",\\\n";
		_file << "\"" << filename << "\" using 1:2 with points title '" << caption << "' lw 2.0";
		_lines.push_back(new Line(filename));
		return _lines.back();
	}
	
	void DrawLine(double x1, double y1, double x2, double y2)
	{
		if(!_lines.empty())
			throw std::runtime_error("Wrong call to DrawLine");
		_file << "set arrow from " << x1 << ',' << y1 << " to " << x2 << ',' << y2 << " nohead\n";
	}
	
	void VerticalLine(double x, size_t lt=1)
	{
		_file << "set arrow from " << x << ",graph(0,0) to " << x << ",graph(1,1) nohead lt " << lt << "\n";
	}
	
	const std::string& FilenamePrefix() const { return _filenamePrefix; }
private:
	std::vector<Line*> _lines;
	std::ofstream _file;
	std::string _filenamePrefix;
};

#endif
