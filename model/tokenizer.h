#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <fstream>
#include <stdexcept>

class Tokenizer
{
	public:
		std::string getString()
		{
			std::string token;
			getToken(token);
			if(token.size()<2 || token[0]!='\"' || token[token.size()-1]!='\"')
				throw std::runtime_error("Expecting string");
			return token.substr(1, token.size()-2);
		}
		double getTokenAsDouble()
		{
			std::string token;
			getToken(token);
			return atof(token.c_str());
		}
		
		bool getTokenAsBool()
		{
			std::string token;
			getToken(token);
			if(token == "1" || token == "true" || token == "True")
				return true;
			else if(token == "0" || token == "false" || token == "False")
				return false;
			else throw std::runtime_error("Expecting boolean");
		}
		
		bool getToken(std::string &token)
		{
			if(!_stream->good()) return false;
			
			std::stringstream s;
			bool finished = false;
			enum { StateStart, StateInToken, StateInString, StateEscaped, StateInLineComment } state = StateStart;
			do {
				char c;
				_stream->get(c);
				if(!_stream->fail())
				{
					switch(state)
					{
					case StateStart:
						if(!(c == ' ' || c == '\t' || c == '\n'))
						{
							if(c == '/')
								state = StateEscaped;
							else
							{
								s << c;
								if(c == '\"')
									state = StateInString;
								else
									state = StateInToken;
							}
						}
						break;
					case StateInString:
						s << c;
						if(c == '\"')
							finished = true;
						break;
					case StateInToken:
						if(c == ' ' || c == '\t' || c == '\n' || c == ';')
							finished = true;
						else
							s << c;
						break;
					case StateEscaped:
						if(c == '/')
							state = StateInLineComment;
						else
							throw std::runtime_error("Incorrect /");
						break;
					case StateInLineComment:
						if(c == '\n')
							state = StateStart;
						break;
					}
				}
				else {
					token = s.str();
					return !token.empty();
				}
			} while(!finished);
			token = s.str();
			return true;
		}
		
		void SetStream(std::ifstream &stream)
		{
			_stream = &stream;
		}
		std::ifstream &Stream() { return *_stream; }
		const std::ifstream &Stream() const { return *_stream; }
		
		std::ifstream *_stream;
};

#endif
