#ifndef AREA_SET_H
#define AREA_SET_H

#include <vector>
#include <string>

#include "units/imagecoordinates.h"
#include "units/radeccoord.h"

class SkyAreaElement
{
	public:
		enum TypeEnum { Circle, Line, Box };
		
		SkyAreaElement() : 
			_type(Circle), _circleRadius(-1.0),
			_circleCentreRA(0.0), _circleCentreDec(0.0),
			_boxUpperRA(0.0), _boxLowerRA(0.0),
			_boxUpperDec(0.0), _boxLowerDec(0.0)
		{
		}
		
		void SetToCircle(double radius, double raCentre, double decCentre)
		{
			_type = Circle;
			_circleRadius = radius;
			_circleCentreRA = raCentre;
			_circleCentreDec = decCentre;
		}
		
		void SetToBox(double upperRA, double lowerRA, double upperDec, double lowerDec)
		{
			_type = Box;
			if(upperRA < lowerRA) std::swap(upperRA, lowerRA);
			if(upperDec < lowerDec) std::swap(upperDec, lowerDec);
			_boxUpperRA = upperRA;
			_boxLowerRA = lowerRA;
			_boxUpperDec = upperDec;
			_boxLowerDec = lowerDec;
		}
		
		template<typename NumType>
		bool IsIn(NumType ra, NumType dec) const
		{
			switch(_type)
			{
				case Circle:
					return ImageCoordinates::AngularDistance<NumType>(ra, dec, _circleCentreRA, _circleCentreDec) <= _circleRadius;
				case Box:
					return ra >= _boxLowerRA && ra <= _boxUpperRA && dec >= _boxLowerDec && dec <= _boxUpperDec;
				default:
					return false;
			}
		}
		
		void Save(std::ostream &stream) const
		{
			stream << "  element {\n    shape ";
			switch(_type)
			{
				case Circle:
					stream << "circle deg " << (_circleRadius * 180.0 / M_PI) << ' ' <<
						RaDecCoord::RAToString(_circleCentreRA) << ' ' <<
						RaDecCoord::DecToString(_circleCentreDec);
						break;
				case Box:
					stream << "box " <<
						RaDecCoord::RAToString(_boxUpperRA) << ' ' <<
						RaDecCoord::DecToString(_boxUpperDec) << ' ' <<
						RaDecCoord::RAToString(_boxLowerRA) << ' ' <<
						RaDecCoord::DecToString(_boxLowerDec);
						break;
				default:
					stream << "?";
			}
			stream << "\n  }\n";
		}

	private:
		TypeEnum _type;
		double _circleRadius, _circleCentreRA, _circleCentreDec;
		double _boxUpperRA, _boxLowerRA, _boxUpperDec, _boxLowerDec;
};

class SkyArea
{
	public:
		SkyArea() : _weight(0.0), _allowCleaning(true)
		{
		}
		SkyArea(const SkyArea &source) :
			_elements(source._elements),
			_name(source._name),
			_weight(source._weight),
			_allowCleaning(source._allowCleaning)
		{
		}
		void operator=(const SkyArea &source)
		{
			_elements = source._elements;
			_name = source._name;
			_weight = source._weight;
			_allowCleaning = source._allowCleaning;
		}
		template<typename NumType>
		bool IsIn(NumType ra, NumType dec) const
		{
			for(std::vector<SkyAreaElement>::const_iterator i=_elements.begin();
					i!=_elements.end(); ++i)
			{
				if(i->IsIn(ra, dec))
					return true;
			}
			return false;
		}
		void AddElement(const SkyAreaElement &element)
		{
			_elements.push_back(element);
		}
		void SetName(const std::string &name) { _name = name; }
		void SetWeight(double weight) { _weight = weight; }
		void SetAllowCleaning(bool allowCleaning) { _allowCleaning = allowCleaning; }
		
		bool AllowCleaning() const { return _allowCleaning; }
		
		void Save(std::ostream &stream) const
		{
			stream << "area {\n  name \"" << _name << "\"\n  allow-cleaning " << (_allowCleaning ? "True" : "False") <<
			'\n';
			if(_weight != 0.0)
				stream << "weight " << _weight << '\n';
			for(std::vector<SkyAreaElement>::const_iterator i=_elements.begin(); i!=_elements.end(); ++i)
			{
				i->Save(stream);
			}
			stream << "}\n";
		}
	private:
		std::vector<SkyAreaElement> _elements;
		std::string _name;
		double _weight;
		bool _allowCleaning;
};

class AreaSet
{
	public:
		void AddArea(const SkyArea &area)
		{
			_areas.push_back(area);
		}
		
		template<typename NumType>
		void FindAreas(std::vector<const SkyArea*> &foundAreas, NumType ra, NumType dec) const
		{
			for(std::vector<SkyArea>::const_iterator i=_areas.begin(); i!=_areas.end(); ++i)
				if(i->IsIn(ra, dec))
					foundAreas.push_back(&*i);
		}
		
		template<typename NumType>
		bool IsInArea(NumType ra, NumType dec) const
		{
			for(std::vector<SkyArea>::const_iterator i=_areas.begin(); i!=_areas.end(); ++i)
				if(i->IsIn(ra, dec)) return true;
			return false;
		}
		
		void FindAreasInImage(std::vector<const SkyArea*> &foundAreas, size_t x, size_t y) const
		{
			long double l, m, ra, dec;
			ImageCoordinates::XYToLM<long double>(x, y, _pixelSizeX, _pixelSizeY, _imageWidth, _imageHeight, l, m);
			ImageCoordinates::LMToRaDec<long double>(l, m, _imageRA, _imageDec, ra, dec);
			FindAreas(foundAreas, ra, dec);
		}
		
		template<typename NumType>
		bool AllowCleaning(NumType ra, NumType dec) const
		{
			std::vector<SkyArea*> areas;
			FindAreas<NumType>(areas, ra, dec);
			// TODO sort on weight / size
			if(areas.empty())
				return false;
			else
				return (*areas.rbegin())->AllowCleaning();
		}
		
		bool AllowCleaningInImage(size_t x, size_t y) const
		{
			std::vector<const SkyArea*> areas;
			FindAreasInImage(areas, x, y);
			// TODO sort on weight / size
			if(areas.empty())
				return false;
			else
				return (*areas.rbegin())->AllowCleaning();
		}
		
		void SetImageProperties(double pixelSizeX, double pixelSizeY, double ra, double dec, size_t width, size_t height)
		{
			_pixelSizeX = pixelSizeX;
			_pixelSizeY = pixelSizeY;
			_imageRA = ra;
			_imageDec = dec;
			_imageWidth = width;
			_imageHeight = height;
		}
		
		void Save(std::ostream &stream) const
		{
			stream << "skyareas fileformat 1.0\n";
			for(std::vector<SkyArea>::const_iterator i=_areas.begin(); i!=_areas.end(); ++i)
				i->Save(stream);
		}
	private:
		std::vector<SkyArea> _areas;
		double _pixelSizeX, _pixelSizeY, _imageRA, _imageDec;
		size_t _imageWidth, _imageHeight;
};

#endif
