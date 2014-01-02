#ifndef TRIANGLE_INTERPOLATOR_H
#define TRIANGLE_INTERPOLATOR_H

#include <cmath>
#include <cstring>

#include "delaunay.h"

class TriangleInterpolator
{
	typedef std::size_t size_t;
public:
	void Interpolate(double* grid, size_t gridWidth, size_t gridHeight,
		double x1, double y1, double v1,
		double x2, double y2, double v2,
		double x3, double y3, double v3)
	{
		double d1 = PointToLineDistance(x1, y1, x2, y2, x3, y3);
		double d2 = PointToLineDistance(x2, y2, x3, y3, x1, y1);
		double d3 = PointToLineDistance(x3, y3, x1, y1, x2, y2);
		double minY = floor(std::min(std::min(y1, y2), y3));
		double maxY = ceil(std::max(std::max(y1, y2), y3));
		if(minY >= gridHeight || maxY < 0.0) 
			return;
		size_t
			minYi = size_t(std::max(minY, 0.0)),
			maxYi = size_t(std::min(maxY, double(gridHeight)));
		for(size_t yi=minYi; yi!=maxYi; ++yi)
		{
			double y = yi;
			double leftX1, leftX2, leftX3;
			bool before1 = y <= y1, before2 = y <= y2, before3 = y <= y3;
			bool after1 = y >= y1, after2 = y >= y2, after3 = y >= y3;
			bool active1 = (before2 && after3) || (after2 && before3);
			bool active2 = (before3 && after1) || (after3 && before1);
			bool active3 = (before1 && after2) || (after1 && before2);
			int actCt = 0;
			if(active1) {
				leftX1 = GetXForYCrossing(y, x2, y2, x3, y3);
				++actCt;
			} else leftX1 = 1e100;
			if(active2) {
				leftX2 = GetXForYCrossing(y, x3, y3, x1, y1);
				++actCt;
			} else leftX2 = 1e100;
			if(active3) {
				leftX3 = GetXForYCrossing(y, x1, y1, x2, y2);
				++actCt;
			} else leftX3 = 1e100;
			if(leftX1 > leftX2) std::swap(leftX1, leftX2);
			if(leftX2 > leftX3) std::swap(leftX2, leftX3);
			if(leftX1 > leftX2) std::swap(leftX1, leftX2);
			if(actCt == 3)
				leftX2 = leftX3;
			else if(actCt <= 1) continue;
			size_t leftX = size_t(std::max(leftX1, 0.0));
			size_t rightX = size_t(ceil(std::max(std::min(leftX2, double(gridWidth)), 0.0)));
			//std::cout << yi << ':' << leftX << '-' << rightX << ' ' << std::flush;
			for(size_t xi=leftX; xi<rightX; ++xi)
			{
				double x = xi;
				double w1 = PointToLineDistance(x, y, x2, y2, x3, y3);
				double w2 = PointToLineDistance(x, y, x3, y3, x1, y1);
				double w3 = PointToLineDistance(x, y, x1, y1, x2, y2);
				double sum = w1 / d1 + w2 / d2 + w3 / d3;
				//std::cout << sum << ' ';
				double val = v1 * w1 / d1 + v2 * w2 / d2 + v3 * w3 / d3;
				grid[xi + yi*gridWidth] = val / sum;
			}
		}
	}
	
private:
	static double PointToLineDistance(
		double px, double py,
		double lx1, double ly1,
		double lx2, double ly2)
	{
		double dx = lx2 - lx1, negDy = ly1 - ly2;
		return fabs(negDy*px + dx*py + (lx1*ly2 - lx2*ly1)) / sqrt(dx*dx + negDy*negDy);
	}
	
	static double GetXForYCrossing(
		double y,
		double lx1, double ly1,
		double lx2, double ly2)
	{
		return (y - ly1) * (lx2 - lx1) / (ly2 - ly1) + lx1;
	}
};

#endif
