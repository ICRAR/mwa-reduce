#ifndef CLEANABLE_IMAGE_SET_H
#define CLEANABLE_IMAGE_SET_H

#include "../cachedimageset.h"
#include "../imagebufferallocator.h"

#include <vector>

namespace clean_algorithms {
		
	class SingleImageSet {
	public:
		struct Value {
			double value;
			Value() { }
			Value(double _value) : value(_value) { }
			double GetValue(size_t i) { 
				return value;
			}
			static Value Zero() { return Value(0.0); }
		};
		
		SingleImageSet(size_t size, SingleImageSet& prototype) :
			image(prototype._allocator->Allocate(size)),
			_allocator(prototype._allocator)
		{
		}
		
		SingleImageSet(size_t size, ImageBufferAllocator<double>& allocator) :
			image(allocator.Allocate(size)),
			_allocator(&allocator)
		{
		}
		
		~SingleImageSet()
		{
			_allocator->Free(image);
		}
		
		void Load(CachedImageSet& set, PolarizationEnum pol, size_t freqIndex)
		{
			set.Load(image, pol, freqIndex, false);
		}
		
		void Store(CachedImageSet& set, PolarizationEnum pol, size_t freqIndex) const
		{
			set.Store(image, pol, freqIndex, false);
		}
		
		Value Get(size_t pixelIndex) const
		{
			return Value(image[pixelIndex]);
		}
		
		double JoinedValue(size_t pixelIndex) const
		{
			return image[pixelIndex];
		}
		
		double JoinedValueNormalized(size_t pixelIndex) const
		{
			return image[pixelIndex];
		}
		
		double AbsJoinedValue(size_t pixelIndex) const
		{
			return fabs(image[pixelIndex]);
		}
		
		bool IsComponentNegative(size_t pixelIndex) const
		{
			return image[pixelIndex]<0.0;
		}
		
		void AddComponent(const SingleImageSet& source, size_t pixelIndex, double factor)
		{
			image[pixelIndex] += source.image[pixelIndex] * factor;
		}
		
		size_t ImageCount() const { return 1; }
		
		static size_t StaticImageCount() { return 1; }
		
		double* GetImage(size_t imageIndex)
		{
			return image;
		}
		static size_t PSFIndex(size_t imageIndex)
		{
			return 0;
		}
		ImageBufferAllocator<double>* Allocator() { return _allocator; }
		
		double* Data() { return image; }
	private:
		double *image;
		
		ImageBufferAllocator<double> *_allocator;
	};
	
	class PolarizedImageSet {
	public:
		struct Value {
			double xx, xyr, xyi, yy;
			double GetValue(size_t i) { 
				switch(i) {
					default:
					case 0: return xx;
					case 1: return xyr;
					case 2: return xyi;
					case 3: return yy;
				}
			}
			static Value Zero() {
				Value zero;
				zero.xx = 0.0; zero.xyr = 0.0;
				zero.xyi = 0.0; zero.yy = 0.0;
				return zero;
			}
		};
		
		PolarizedImageSet(size_t size, PolarizedImageSet& prototype) :
			xx(prototype._allocator->Allocate(size)),
			xyr(prototype._allocator->Allocate(size)),
			xyi(prototype._allocator->Allocate(size)),
			yy(prototype._allocator->Allocate(size)),
			_allocator(prototype._allocator)
		{
		}
		
		PolarizedImageSet(size_t size, ImageBufferAllocator<double>& allocator) :
			xx(allocator.Allocate(size)),
			xyr(allocator.Allocate(size)),
			xyi(allocator.Allocate(size)),
			yy(allocator.Allocate(size)),
			_allocator(&allocator)
		{
		}
		
		~PolarizedImageSet()
		{
			_allocator->Free(xx);
			_allocator->Free(xyr);
			_allocator->Free(xyi);
			_allocator->Free(yy);
		}
		
		void LoadLinear(CachedImageSet& set, size_t freqIndex)
		{
			set.Load(xx, PolarizationEnum::XX, freqIndex, false);
			set.Load(xyr, PolarizationEnum::XY, freqIndex, false);
			set.Load(xyi, PolarizationEnum::XY, freqIndex, true);
			set.Load(yy, PolarizationEnum::YY, freqIndex, false);
		}
		
		void LoadStokes(CachedImageSet& set, size_t freqIndex)
		{
			set.Load(xx, PolarizationEnum::StokesI, freqIndex, false);
			set.Load(xyr, PolarizationEnum::StokesQ, freqIndex, false);
			set.Load(xyi, PolarizationEnum::StokesU, freqIndex, false);
			set.Load(yy, PolarizationEnum::StokesV, freqIndex, false);
		}
		
		void StoreLinear(CachedImageSet& set, size_t freqIndex) const
		{
			set.Store(xx, PolarizationEnum::XX, freqIndex, false);
			set.Store(xyr, PolarizationEnum::XY, freqIndex, false);
			set.Store(xyi, PolarizationEnum::XY, freqIndex, true);
			set.Store(yy, PolarizationEnum::YY, freqIndex, false);
		}
		
		void StoreStokes(CachedImageSet& set, size_t freqIndex) const
		{
			set.Store(xx, PolarizationEnum::StokesI, freqIndex, false);
			set.Store(xyr, PolarizationEnum::StokesQ, freqIndex, false);
			set.Store(xyi, PolarizationEnum::StokesU, freqIndex, false);
			set.Store(yy, PolarizationEnum::StokesV, freqIndex, false);
		}
		
		Value Get(size_t index) const
		{
			Value v;
			v.xx = xx[index];
			v.xyr = xyr[index];
			v.xyi = xyi[index];
			v.yy = yy[index];
			return v;
		}
		
		double JoinedValue(size_t index) const
		{
			return SquaredSum(index);
		}
		
		double JoinedValueNormalized(size_t index) const
		{
			return sqrt(SquaredSum(index));
		}
		
		double AbsJoinedValue(size_t index) const
		{
			return SquaredSum(index);
		}
		
		double SquaredSum(size_t index) const
		{
			return
				xx[index]*xx[index] +
				xyr[index]*xyr[index] + xyi[index]*xyi[index] +
				yy[index]*yy[index];
		}
		
		bool IsComponentNegative(size_t index) const
		{
			return xx[index]<0.0 || yy[index]<0.0;
		}
		
		void AddComponent(const PolarizedImageSet& source, size_t index, double factor)
		{
			xx[index] += source.xx[index] * factor;
			xyr[index] += source.xyr[index] * factor;
			xyi[index] += source.xyi[index] * factor;
			yy[index] += source.yy[index] * factor;
		}
		
		size_t ImageCount() const { return 4; }
		
		static size_t StaticImageCount() { return 4; }
		
		double* GetImage(size_t imageIndex)
		{
			double* vals[4] = { xx, xyr, xyi, yy };
			return vals[imageIndex];
		}
		static size_t PSFIndex(size_t imageIndex)
		{
			return 0;
		}
		ImageBufferAllocator<double>* Allocator() { return _allocator; }
	private:
		double *xx, *xyr, *xyi, *yy;
		
		ImageBufferAllocator<double> *_allocator;
	};
	
	class MultiImageSet {
	public:
		struct Value {
			std::vector<PolarizedImageSet::Value> values;
			double GetValue(size_t i)
			{
				return values[i/4].GetValue(i%4);
			}
			static Value Zero() { return Value(); }
		};
		
		MultiImageSet(size_t imageSize, MultiImageSet& prototype)
		{
			for(size_t i=0; i!=prototype._sets.size(); ++i)
			{
				_sets.push_back(new PolarizedImageSet(imageSize, *prototype.Allocator()));
			}
		}
		
		MultiImageSet(size_t imageSize, size_t count, ImageBufferAllocator<double>& allocator)
		{
			for(size_t i=0; i!=count; ++i)
			{
				_sets.push_back(new PolarizedImageSet(imageSize, allocator));
			}
		}
		
		~MultiImageSet()
		{
			for(std::vector<PolarizedImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				delete *i;
			}
		}
		
		void LoadLinear(CachedImageSet& set, size_t i)
		{
			_sets[i]->LoadLinear(set, i);
		}
		
		void StoreLinear(CachedImageSet& set, size_t i) const
		{
			_sets[i]->StoreLinear(set, i);
		}
		
		void LoadStokes(CachedImageSet& set, size_t i)
		{
			_sets[i]->LoadStokes(set, i);
		}
		
		void StoreStokes(CachedImageSet& set, size_t i) const
		{
			_sets[i]->StoreStokes(set, i);
		}
		
		double JoinedValue(size_t index) const
		{
			double val = 0.0;
			for(std::vector<PolarizedImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				val += (*i)->JoinedValueNormalized(index);
			}
			return val;
		}
		
		double JoinedValueNormalized(size_t index) const
		{
			return JoinedValue(index) / _sets.size();
		}
		
		double AbsJoinedValue(size_t index) const
		{
			return JoinedValue(index);
		}
		
		bool IsComponentNegative(size_t index) const
		{
			for(std::vector<PolarizedImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				if((*i)->IsComponentNegative(index)) return true;
			}
			return false;
		}
		
		void AddComponent(const MultiImageSet& source, size_t index, double factor)
		{
			for(size_t i=0; i!=_sets.size(); ++i)
			{
				_sets[i]->AddComponent(*source._sets[i], index, factor);
			}
		}
		
		Value Get(const size_t index)
		{
			Value v;
			v.values.resize(_sets.size());
			for(size_t i=0; i!=_sets.size(); ++i)
				v.values[i] = _sets[i]->Get(index);
			return v;
		}
		
		size_t ImageCount() const { return PolarizedImageSet::StaticImageCount() * _sets.size(); }
		
		double* GetImage(size_t imageIndex)
		{
			return _sets[imageIndex/4]->GetImage(imageIndex%4);
		}
		
		double* GetImage(size_t polIndex, size_t freqIndex)
		{
			return _sets[freqIndex]->GetImage(polIndex);
		}
		
		static size_t PSFIndex(size_t imageIndex)
		{
			return imageIndex/4;
		}
		ImageBufferAllocator<double>* Allocator()
		{ 
			return _sets.front()->Allocator();
		}
	private:
		std::vector<PolarizedImageSet*> _sets;
	};	
}

#endif
