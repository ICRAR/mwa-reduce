#include <iostream>

using namespace std;

int count = 0;

bool isnum(char c) { return (c>='0' && c<='9') || c=='.'; }

void outputComplex(string real, string imag)
{
	if(count == 0)
		cout << "// " << (count+50) << " - " << (count+69) << " MHz.\n";
	else
	{
		std::cout << ", ";
		if((count%2)==0)
		{
			cout << '\n';
			if((count%20)==0) {
				if(count+69 < 500)
					cout << "// " << (count+50) << " - " << (count+69) << " MHz.\n";
				else
					cout << "// " << (count+50) << " - 500 MHz.\n";
			}
		}
	}
	std::cout << "ctype(" << real << ',' << imag << ")";
	++count;
}

void parse(const string line)
{
	enum State { Begin, InReal, AfterReal, InImag } state = Begin;
	size_t pos = 0;
	std::string real, imag;
	for(size_t i=0; i!=line.size(); ++i)
	{
		char c = line[i];
		switch(state)
		{
			case Begin:
				if(isnum(c))
				{
					state = InReal;
					real += c;
				}
				break;
			case InReal:
				if(isnum(c))
					real += c;
				else
					state = AfterReal;
				break;
			case AfterReal:
				if(isnum(c))
				{
					state = InImag;
					imag += c;
				}
				break;
			case InImag:
				if(isnum(c))
					imag += c;
				else {
					state = Begin;
					outputComplex(real, imag);
					real = string(); imag = string();
				}
		}
	}
	if(state == InImag)
		outputComplex(real, imag);
}

int main(int argc, char *argv[])
{
	while(cin.good())
	{
		string line;
		getline(cin, line);
		parse(line);
	}
	std::cout << '\n';
}