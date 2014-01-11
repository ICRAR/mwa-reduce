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
		double x3, double y3, double v3) const
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
	
	void InterpolateConvexHullVertex(double* grid, size_t gridWidth, size_t gridHeight,
		double x1, double y1,
		double x2, double y2, double v2,
		double x3, double y3) const
	{
		if((x1 == x2 && y1 == y2) || (x2 == x3 && y2 == y3))
			return;
		double
			yd1 = x1-x2, xd1 = y2-y1,
			yd3 = x2-x3, xd3 = y3-y2;
		double yStart, yEnd;
		if(xd1 < 0.0)
			yStart = GetXForYCrossing(0.0, y2, x2, y2+yd1, x2+xd1);
		else if(xd1 > 0.0)
			yStart = GetXForYCrossing(gridWidth, y2, x2, y2+yd1, x2+xd1);
		else
			yStart = y2;
		if(xd3 < 0.0)
			yEnd = GetXForYCrossing(0.0, y2, x2, y2+yd3, x2+xd3);
		else if(xd3 > 0.0)
			yEnd = GetXForYCrossing(gridWidth, y2, x2, y2+yd3, x2+xd3);
		else
			yEnd = y2;
		if(yEnd < yStart) std::swap(yStart, yEnd);
		if(yStart < 0.0 || !std::isfinite(yStart)) yStart = 0.0;
		if(yEnd < 0.0) yEnd = 0.0;
		if(yEnd > gridHeight || !std::isfinite(yEnd)) yEnd = gridHeight;
		if(yEnd < y2) yEnd = y2;
		if(yStart > y2) yStart = y2;
		size_t yiStart = size_t(yStart), yiEnd = size_t(ceil(yEnd));
		for(size_t yi=yiStart; yi<yiEnd; ++yi)
		{
			double y = yi;
			double xc1, xc3;
			xc1 = GetXForYCrossing(y, x2, y2, x2+xd1, y2+yd1);
			xc3 = GetXForYCrossing(y, x2, y2, x2+xd3, y2+yd3);
			if(xc3 < xc1) std::swap(xc1, xc3);
			if(yd1 >= 0.0 && yd3 < 0.0) {
				xc3 = xc1;
				xc1 = 0.0;
			}
			if(yd1 < 0.0 && yd3 >= 0.0) {
				xc1 = xc3;
				xc3 = gridWidth;
			}
			if(xc1 < 0.0) xc1 = 0.0;
			if(xc3 < 0.0) xc3 = 0.0;
			if(xc3 > gridWidth) xc3 = gridWidth;
			size_t xStart = size_t(xc1), xEnd = size_t(ceil(xc3));
			for(size_t xi=xStart; xi<xEnd; ++xi)
			{
				grid[xi + yi*gridWidth] = v2;
			}
		}
	}
	
	void InterpolateConvexHullEdge(double* grid, size_t gridWidth, size_t gridHeight,
		double x1, double y1, double v1,
		double x2, double y2, double v2) const
	{
		double dx = y2-y1, dy = x1-x2;
		double yStart, yEnd;
		if(dx < 0.0) {
			yStart = GetXForYCrossing(0.0, y1, x1, y1+dy, x1+dx);
			yEnd = GetXForYCrossing(0.0, y2, x2, y2+dy, x2+dx);
		}
		else {
			yStart = GetXForYCrossing(gridWidth, y1, x1, y1+dy, x1+dx);
			yEnd = GetXForYCrossing(gridWidth, y2, x2, y2+dy, x2+dx);
		}
		if(yStart > yEnd) std::swap(yStart, yEnd);
		double minY = std::min(y1, y2), maxY = std::max(y1, y2);
		if(yStart > minY) yStart = minY;
		if(yEnd < maxY) yEnd = maxY;
		if(yStart < 0.0 || !std::isfinite(yStart)) yStart = 0.0;
		if(yEnd < 0.0) yEnd = 0.0;
		if(yEnd > gridHeight || !std::isfinite(yEnd)) yEnd = gridHeight;
		size_t yiStart = size_t(yStart), yiEnd = size_t(ceil(yEnd));
		for(size_t yi=yiStart; yi<yiEnd; ++yi)
		{
			double y = yi;
			double cx1 = GetXForYCrossing(y, x1, y1, x1+dx, y1+dy);
			double cx2 = GetXForYCrossing(y, x2, y2, x2+dx, y2+dy);
			if(cx1 > cx2) std::swap(cx1, cx2);
			if(cx1 < 0.0) cx1 = 0.0;
			if(cx2 < 0.0) cx2 = 0.0;
			else if(cx2 > gridWidth) cx2 = gridWidth;
			double cm = GetXForYCrossing(y, x1, y1, x2, y2);
			if(cm < 0.0) cm = 0.0;
			else if(cm > gridWidth) cm = gridWidth;
			if(dx >= 0.0 && cm > cx1)
				cx1 = cm;
			else if(dx <= 0.0 && cm < cx2)
				cx2 = cm;
			size_t xStart = size_t(cx1), xEnd = size_t(ceil(cx2));
			for(size_t xi=xStart; xi<xEnd; ++xi)
			{
				double x = xi;
				double d1 = PointToLineDistance(x, y, x1, y1, x1+dx, y1+dy);
				double d2 = PointToLineDistance(x, y, x2, y2, x2+dx, y2+dy);
				grid[yi*gridWidth + xi] = d2/(d1+d2)*v1 + d1/(d1+d2)*v2;
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
