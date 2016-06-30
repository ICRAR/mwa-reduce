#ifndef RADECCOORD_H
#define RADECCOORD_H

#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cmath>

class RaDecCoord
{
	private:
		static bool isRASeparator(char c) { return c==':' || c==' '; }
		static bool isDecSeparator(char c) { return c=='.' || c==' '; }
	public:
		static long double ParseRA(const std::string &str)
		{
			char *cstr;
			bool sign = false;
			for(size_t i=0; i!=str.size(); ++i)
			{
				if(str[i] == '-')
				{
					sign = true;
					break;
				} else if(str[i] != ' ')
				{
					sign = false;
					break;
				}
			}
			long double secs=0.0, mins=0.0,
				hrs = strtol(str.c_str(), &cstr, 10);
			// Parse format '00h00m00.0s'
			if(*cstr == 'h')
			{
				++cstr;
				mins = strtol(cstr, &cstr, 10);
				if(*cstr == 'm')
				{
					++cstr;
					secs = strtold(cstr, &cstr);
					if(*cstr == 's')
						++cstr;
					else throw std::runtime_error("Missing 's'");
				} else throw std::runtime_error("Missing 'm'");
			}
			// Parse format '00:00:00.0'
			else if(isRASeparator(*cstr))
			{
				++cstr;
				mins = strtol(cstr, &cstr, 10);
				if(isRASeparator(*cstr))
				{
					++cstr;
					secs = strtold(cstr, &cstr);
				} else throw std::runtime_error("Missing ':' after minutes");
			}
			else throw std::runtime_error("Missing 'h' or ':'");
			if(*cstr != 0)
				throw std::runtime_error("Could not parse RA (string contains more tokens than expected)");
			if(sign)
				return (hrs/24.0 - mins/(24.0*60.0) - secs/(24.0*60.0*60.0))*2.0*M_PI;
			else
				return (hrs/24.0 + mins/(24.0*60.0) + secs/(24.0*60.0*60.0))*2.0*M_PI;
		}
		
		static long double ParseDec(const std::string &str)
		{
			char *cstr;
			bool sign = false;
			for(size_t i=0; i!=str.size(); ++i)
			{
				if(str[i] == '-')
				{
					sign = true;
					break;
				} else if(str[i] != ' ')
				{
					sign = false;
					break;
				}
			}
			long double secs=0.0, mins=0.0,
				degs = strtol(str.c_str(), &cstr, 10);
			// Parse format '00d00m00.0s'
			if(*cstr == 'd')
			{
				++cstr;
				mins = strtol(cstr, &cstr, 10);
				if(*cstr == 'm')
				{
					++cstr;
					secs = strtold(cstr, &cstr);
					if(*cstr == 's')
						++cstr;
					else throw std::runtime_error("Missing 's'");
				} else throw std::runtime_error("Missing 'm'");
			}
			// Parse format '00.00.00.0'
			else if(isDecSeparator(*cstr))
			{
				++cstr;
				mins = strtol(cstr, &cstr, 10);
				if(isDecSeparator(*cstr))
				{
					++cstr;
					secs = strtold(cstr, &cstr);
				} else throw std::runtime_error("Missing '.' after minutes");
			}
			else throw std::runtime_error("Missing 'd' or '.' after degrees");
			if(*cstr != 0)
				throw std::runtime_error("Could not parse Dec (string contains more tokens than expected)");
			else if(sign)
				return (degs/360.0 - mins/(360.0*60.0) - secs/(360.0*60.0*60.0))*2.0*M_PI;
			else
				return (degs/360.0 + mins/(360.0*60.0) + secs/(360.0*60.0*60.0))*2.0*M_PI;
		}
		
		static std::string RAToString(long double ra)
		{
			const long double partsPerHour = 60.0L*60.0L*1000.0L;
			long double hrs = fmodl(ra * (12.0L / M_PIl), 24.0L);
			hrs = roundl(hrs * partsPerHour)/partsPerHour;
			std::stringstream s;
			if(hrs < 0) {
				hrs = -hrs;
				s << '-';
			}
			hrs = (roundl(hrs*partsPerHour)+0.5) / partsPerHour;
			int hrsInt = int(floorl(hrs)), minInt = int(floorl(fmodl(hrs,1.0L)*60.0L)),
				secInt = int(floorl(fmodl(hrs*60.0L,1.0L)*(60.0L))),
				subSecInt = int(floorl(fmodl(hrs*(60.0L*60.0L),1.0L)*(1000.0L)));
			s << (char) ((hrsInt/10)+'0') << (char) ((hrsInt%10)+'0') << 'h'
				<< (char) ((minInt/10)+'0') << (char) ((minInt%10)+'0') << 'm'
				<< (char) ((secInt/10)+'0') << (char) ((secInt%10)+'0');
			s << '.' << (char) (subSecInt/100+'0');
			if(subSecInt%100!=0)
			{
				s << (char) ((subSecInt%100)/10+'0');
				if(subSecInt%10!=0)
					s << (char) ((subSecInt%10)+'0');
			}
			s << 's';
			return s.str();
		}
		
		static std::string RAToString(long double ra, char delimiter)
		{
			const long double partsPerHour = 60.0L*60.0L*1000.0L;
			long double hrs = fmodl(ra * (12.0L / M_PIl), 24.0L);
			hrs = roundl(hrs * partsPerHour)/partsPerHour;
			std::stringstream s;
			if(hrs < 0) {
				hrs = -hrs;
				s << '-';
			}
			hrs = (roundl(hrs*partsPerHour)+0.5) / partsPerHour;
			int hrsInt = int(floorl(hrs)), minInt = int(floorl(fmodl(hrs,1.0L)*60.0L)),
				secInt = int(floorl(fmodl(hrs*60.0L,1.0L)*(60.0L))),
				subSecInt = int(floorl(fmodl(hrs*(60.0L*60.0L),1.0L)*(1000.0L)));
			s << (char) ((hrsInt/10)+'0') << (char) ((hrsInt%10)+'0') << delimiter
				<< (char) ((minInt/10)+'0') << (char) ((minInt%10)+'0') << delimiter
				<< (char) ((secInt/10)+'0') << (char) ((secInt%10)+'0');
			s << '.' << (char) (subSecInt/100+'0');
			if(subSecInt%100!=0)
			{
				s << (char) ((subSecInt%100)/10+'0');
				if(subSecInt%10!=0)
					s << (char) ((subSecInt%10)+'0');
			}
			return s.str();
		}
		
		static void RAToHMS(long double ra, int& hrs, int& min, double& sec)
		{
			const long double partsPerHour = 60.0L*60.0L*100.0L;
			long double hrsf = fmodl(ra * (12.0L / M_PIl), 24.0L);
			hrsf = roundl(hrsf * partsPerHour)/partsPerHour;
			bool negate = hrsf < 0;
			if(negate) {
				hrsf = -hrsf;
			}
			hrsf = (roundl(hrsf*partsPerHour)+0.5) / partsPerHour;
			if(negate)
				hrs = -int(floorl(hrsf));
			else
				hrs = int(floorl(hrsf));
			min = int(floorl(fmodl(hrsf,1.0L)*60.0L));
			sec = floorl(100.0L*fmodl(hrsf*60.0L,1.0L)*60.0L)/100.0L;
		}
		
		static std::string DecToString(long double dec)
		{
			const long double partsPerDeg = 60.0L*60.0L*1000.0L;
			long double deg = dec * (180.0L / M_PIl);
			deg = round(deg * partsPerDeg)/partsPerDeg;
			std::stringstream s;
			if(deg < 0) {
				deg = -deg;
				s << '-';
			}
			deg = (round(deg*partsPerDeg)+0.5L) / partsPerDeg;
			int
				degInt = int(floorl(deg)),
				minInt = int(floorl(fmodl(deg,1.0L)*60.0L)),
				secInt = int(floorl(fmodl(deg,1.0L/60.0L)*(60.0L*60.0L))),
				subSecInt = int(floorl(fmodl(deg,1.0L/60.0L/60.0L)*(60.0L*60.0L*1000.0L)));
			s << (char) ((degInt/10)+'0') << (char) ((degInt%10)+'0') << 'd'
				<< (char) ((minInt/10)+'0') << (char) ((minInt%10)+'0') << 'm'
				<< (char) ((secInt/10)+'0') << (char) ((secInt%10)+'0');
			if(subSecInt!=0)
			{
				s << '.' << (char) (subSecInt/100+'0');
				if(subSecInt%100!=0)
				{
					s << (char) ((subSecInt/10)%10+'0');
					if(subSecInt%10!=0)
						s << (char) (subSecInt%10+'0');
				}
			}
			s << 's';
			return s.str();
		}
		
		static std::string DecToString(long double dec, char delimiter)
		{
			const long double partsPerDeg = 60.0L*60.0L*1000.0L;
			long double deg = dec * (180.0L / M_PIl);
			deg = roundl(deg * partsPerDeg)/partsPerDeg;
			std::stringstream s;
			if(deg < 0) {
				deg = -deg;
				s << '-';
			}
			deg = (roundl(deg*partsPerDeg)+0.5L) / partsPerDeg;
			int
				degInt = int(floorl(deg)),
				minInt = int(floorl(fmodl(deg,1.0L)*60.0L)),
				secInt = int(floorl(fmodl(deg,1.0L/60.0L)*(60.0L*60.0L))),
				subSecInt = int(floorl(fmodl(deg,1.0L/60.0L/60.0L)*(60.0L*60.0L*1000.0L)));
			s << (char) ((degInt/10)+'0') << (char) ((degInt%10)+'0') << delimiter
				<< (char) ((minInt/10)+'0') << (char) ((minInt%10)+'0') << delimiter
				<< (char) ((secInt/10)+'0') << (char) ((secInt%10)+'0');
			if(subSecInt!=0)
			{
				s << '.' << (char) (subSecInt/100+'0');
				if(subSecInt%100!=0)
				{
					s << (char) ((subSecInt/10)%10+'0');
					if(subSecInt%10!=0)
						s << (char) (subSecInt%10+'0');
				}
			}
			return s.str();
		}
		
		static void DecToDMS(long double dec, int& deg, int& min, double& sec)
		{
			const long double partsPerDeg = 60.0L*60.0L*100.0L;
			long double degf = dec * (180.0 / M_PI);
			degf = round(degf * partsPerDeg)/partsPerDeg;
			bool negate = degf < 0;
			if(negate) {
				degf = -degf;
			}
			degf = (round(degf*partsPerDeg)+0.5) / partsPerDeg;
			
			if(negate)
				deg = -int(floor(degf));
			else
				deg = int(floor(degf));
			min = int(floor(fmod(degf,1.0)*60.0)),
			sec = floor(100.*fmod(degf,1.0/60.0)*(60.0*60.0))/100.0;
		}
	private:
		RaDecCoord() { }
};

#endif
