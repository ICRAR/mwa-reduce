/***************************************************************************
 *   Copyright (C) 2007 by Andre Offringa   *
 *   offringa@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "pngfile.h"

#include <stdexcept>

PngFile::PngFile(const std::string &filename, unsigned width, unsigned height) : _filename(filename), _width(width), _height(height), _pixelSize(4)
{
}

PngFile::~PngFile()
{
}

void PngFile::BeginWrite()
{
	_fp = fopen(_filename.c_str(), "wb");
	if(!_fp)
		throw std::runtime_error("Can not open file");
	
	_png_ptr =
		png_create_write_struct(PNG_LIBPNG_VER_STRING, (png_voidp) NULL, NULL, NULL);
	
	if(!_png_ptr) {
		fclose(_fp);
		throw std::runtime_error("Can not create png write structure");
	}
	
	_info_ptr = png_create_info_struct(_png_ptr);
	if(!_info_ptr) {
		png_destroy_write_struct(&_png_ptr, (png_infopp) NULL);
		fclose(_fp);
		throw std::runtime_error("Can not write info structure to file");
	}
	
	if (setjmp(png_jmpbuf(_png_ptr)))
	{
		png_destroy_write_struct(&_png_ptr, &_info_ptr);
		fclose(_fp);
		throw std::runtime_error("Unknown error occured during writing of png file");
	}
	
	png_init_io(_png_ptr, _fp);
	
	png_set_IHDR(_png_ptr, _info_ptr, _width, _height, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	
	_row_pointers = (png_bytep*) png_malloc(_png_ptr, _height*sizeof(png_bytep));
	
	for(unsigned i=0; i<_height;i++)
		_row_pointers[i]=(png_bytep) png_malloc(_png_ptr, _width*_pixelSize);
}

void PngFile::Close() {
	png_set_rows(_png_ptr, _info_ptr, _row_pointers);
	png_write_png(_png_ptr, _info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
	png_write_end(_png_ptr, _info_ptr);
	
	for(unsigned i=0; i<_height;i++)
		png_free(_png_ptr, _row_pointers[i]);
	png_free(_png_ptr, _row_pointers);
	
	png_destroy_write_struct(&_png_ptr, &_info_ptr);
	fclose(_fp);
}

void PngFile::Clear(int colorR, int colorG, int colorB, int colorA) throw()
{
	for(unsigned y=0;y<_height;y++) {
		int xa = 0;
		for(unsigned x=0;x<_width;x++) {
			_row_pointers[y][xa++] = colorR;
			_row_pointers[y][xa++] = colorG;
			_row_pointers[y][xa++] = colorB;
			_row_pointers[y][xa++] = colorA;
		}
	}
}

void PngFile::PlotDatapoint(unsigned x, unsigned y, int colorR, int colorG, int colorB, int colorA)
{
	int width = 3;
	for(int xi=(signed) x-width/2;xi<=(signed) x+(width-1)/2;xi++) {
		if(xi >= 0 && xi < (signed) _width) {
			_row_pointers[y][xi*_pixelSize] = colorR;
			_row_pointers[y][xi*_pixelSize+1] = colorG;
			_row_pointers[y][xi*_pixelSize+2] = colorB;
			_row_pointers[y][xi*_pixelSize+3] = colorA;
		}
	}
	for(int yi=(signed) y-width/2;yi<=(signed) y+(width-1)/2;yi++) {
		if(yi >= 0 && yi < (signed) _height) {
			_row_pointers[yi][x*_pixelSize] = colorR;
			_row_pointers[yi][x*_pixelSize+1] = colorG;
			_row_pointers[yi][x*_pixelSize+2] = colorB;
			_row_pointers[yi][x*_pixelSize+3] = colorA;
		}
	}
}
