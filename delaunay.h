#ifndef DELAUNAY_H
#define DELAUNAY_H

//#define CHECK_DELAUNAY

#include <algorithm>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <stack>

class Delaunay
{
private:
	class TriangleNode;
	class Vertex
	{
	public:
		Vertex() : inConvexHull(false), outsideConvexHull(false)
		{ }
		
		double x, y;
		bool inConvexHull, outsideConvexHull;
		std::vector<TriangleNode*> triangles;
		void *userData;
		
		void Print()
		{
			std::cout << this << '(';
			if(inConvexHull) std::cout << "convex)";
			else if(outsideConvexHull) std::cout << "outside)";
			else std::cout << "inside)";
		}
		
		bool operator<(const Vertex& other) const
		{
			if(x < other.x) return true;
			else if(x == other.x) return y < other.y;
			else return false;
		}
		static bool is_smaller(const Vertex* a, const Vertex* b)
		{
			return *a < *b;
		}
		
		bool IsLeftTurn(const std::vector<Vertex*>& vertices) const
		{
			if(vertices.size() < 2) return false;
			
			std::vector<Vertex*>::const_reverse_iterator i=vertices.rbegin();
			const Vertex* mid = *i;
			++i;
			const Vertex* left = *i;
			return IsLeftTurn(left, mid, this);
		}
		static bool IsLeftTurn(const Vertex* a, const Vertex* b, const Vertex* c)
		{
			return (b->x - a->x)*(c->y - a->y) - (b->y - a->y)*(c->x - a->x) > 0.0;
		}
		bool IsRightTurn(const std::vector<Vertex*>& vertices) const
		{
			if(vertices.size() < 2) return false;
			
			std::vector<Vertex*>::const_reverse_iterator i=vertices.rbegin();
			const Vertex* mid = *i;
			++i;
			const Vertex* left = *i;
			return IsRightTurn(left, mid, this);
		}
		static bool IsRightTurn(const Vertex* a, const Vertex* b, const Vertex* c)
		{
			return (b->x - a->x)*(c->y - a->y) - (b->y - a->y)*(c->x - a->x) < 0.0;
		}
		static double TurnDirection(const Vertex* a, const Vertex* b, double x, double y)
		{
			return (b->x - a->x)*(y - a->y) - (b->y - a->y)*(x - a->x);
		}
		void Disconnect(TriangleNode* t)
		{
			for(std::vector<TriangleNode*>::iterator i=triangles.begin(); i!=triangles.end(); ++i)
			{
				if(*i == t) {
					triangles.erase(i);
					break;
				}
			}
		}
		double Distance(const Vertex& other) const
		{
			double dx = other.x - x, dy = other.y - y;
			return sqrt(dx*dx + dy*dy);
		}
	};
	
	class Edge
	{
	public:
		class Vertex *vertices[2];
		
		Edge(Vertex* v1, Vertex* v2)
		{
			vertices[0] = v1; vertices[1] = v2;
		}
		bool HasSameVertex(const Edge& other, unsigned v1, unsigned v2)
		{
			for(v1=0; v1!=2; ++v1)
			{
				for(v2=0; v2!=2; ++v2)
					if(vertices[v1] == other.vertices[v2])
						return true;
			}
			//v1=0; v2=0;
			return false;
		}
		bool Crosses(const Vertex* a, const Vertex* b)
		{
			if(a == vertices[0] || a == vertices[1] ||
				 b == vertices[0] || b == vertices[1]) return false;
			double
				x0 = vertices[0]->x,
				y0 = vertices[0]->y,
				x1 = vertices[1]->x,
				y1 = vertices[1]->y,
				s1_x = x1 - x0,
				s1_y = y1 - y0,
				s2_x = b->x - a->x,
				s2_y = b->y - a->y;
			double
				s = (-s1_y * (x0 - a->x) + s1_x * (y0 - a->y)) / (-s2_x * s1_y + s1_x * s2_y),
				t = ( s2_x * (y0 - a->y) - s2_y * (x0 - a->x)) / (-s2_x * s1_y + s1_x * s2_y);
			return s > 0.0 && s < 1.0 && t > 0.0 && t < 1.0;
		}
	};
	
	class TriangleNode
	{
	public:
		TriangleNode()
		{
			for(size_t i=0; i!=3; ++i)
			{
				vertices[i] = 0;
				subNodes[i] = 0;
			}
			parent = 0;
		}
		~TriangleNode()
		{
			for(size_t i=0; i!=3; ++i)
			{
				delete subNodes[i];
				if(vertices[i] != 0)
					vertices[i]->Disconnect(this);
			}
			if(parent != 0)
			{
				for(size_t i=0; i!=3; ++i)
				{
					if(parent->subNodes[i] == this)
						parent->subNodes[i] = 0;
				}
			}
		}
		void Print()
		{
			vertices[0]->Print();
			for(size_t i=1; i!=3; ++i) {
				std::cout << ' ';
				vertices[i]->Print();
			}
			if(vertices[0]==vertices[1] || vertices[1]==vertices[2] || vertices[2]==vertices[0])
				std::cout << " degenerate";
			if(HasSubNodes())
				std::cout << " parent node";
			else
				std::cout << " leaf node";
			if(Vertex::IsRightTurn(vertices[0], vertices[1], vertices[2]))
				std::cout << " CW";
			else
				std::cout << " ACW";
			std::cout << '\n';
		}
		bool IsIn(double x, double y) const
		{
			for(size_t i=0; i!=3; ++i)
			{
				double d = Vertex::TurnDirection(vertices[i], vertices[(i+1)%3], x, y);
				if(d < 0.0) { // right turn
					return false;
				}
			}
			return true;
		}
		bool IsOutsideConvexHull() const
		{
			size_t outsideCount = 0;
			for(size_t i=0; i!=3; ++i)
			{
				if(vertices[i]->outsideConvexHull)
					++outsideCount;
			}
			return outsideCount>=2;
		}
		bool HasConvexEdge() const
		{
			size_t convexCount = 0;
			for(size_t i=0; i!=3; ++i)
			{
				if(vertices[i]->inConvexHull)
					++convexCount;
			}
			return convexCount>=2;
		}
		bool HasOutsideEdge() const
		{
			size_t outsideCount = 0;
			for(size_t i=0; i!=3; ++i)
			{
				if(vertices[i]->outsideConvexHull)
					++outsideCount;
			}
			return outsideCount>=2;
		}
		bool HasConvexVertex() const
		{
			return
				vertices[0]->inConvexHull ||
				vertices[1]->inConvexHull ||
				vertices[2]->inConvexHull;
		}
		bool HasOutsideVertex() const
		{
			return
				vertices[0]->outsideConvexHull ||
				vertices[1]->outsideConvexHull ||
				vertices[2]->outsideConvexHull;
		}
		bool HasInsideVertex() const
		{
			for(size_t i=0; i!=3; ++i)
				if(!vertices[i]->outsideConvexHull && !vertices[i]->inConvexHull)
					return true;
			return false;
		}
		bool HasVertex(Vertex* v) const
		{
			return vertices[0] == v || vertices[1] == v || vertices[2] == v;
		}
		bool HasSubNodes() const
		{
			return subNodes[0] != 0 || subNodes[1] != 0 || subNodes[2] != 0;
		}
		
		TriangleNode* GetNode(double x, double y)
		{
			for(size_t i=0; i!=3; ++i)
			{
				if(subNodes[i]!=0 && subNodes[i]->IsIn(x, y))
					return subNodes[i]->GetNode(x, y);
			}
			if(HasSubNodes())
				std::cout << "Didn't find a leaf node in GetNode()!\n";
			return this;
		}
		bool IsClockwise() const
		{
			return Vertex::IsRightTurn(vertices[0], vertices[1], vertices[2]);
		}
		void TurnAround()
		{
			std::swap(vertices[0], vertices[2]);
		}
		Vertex* vertices[3];
		TriangleNode* subNodes[3];
		TriangleNode* parent;
	};
	
public:
	struct Triangle
	{
		double x[3], y[3];
		void* userData[3];
	};
	struct ConvexVertex
	{
		double x, y;
		void* userData;
	};
	
	Delaunay() : _boundTop(0), _topNode(0)
	{
	}
	
	~Delaunay() {
		Clear();
	}
	
	void Clear()
	{
		delete _topNode;
		for(std::vector<Vertex*>::iterator i=_vertices.begin(); i!=_vertices.end(); ++i)
			delete *i;
		_vertices.clear();
		for(std::vector<Edge*>::iterator i=_convexHull.begin(); i!=_convexHull.end(); ++i)
			delete *i;
		_convexHull.clear();
		_boundTop = 0;
		_triangles.clear();
	}
	
	void AddVertex(double ra, double dec, void* userData=0)
	{
		Vertex* vertex = new Vertex();
		vertex->x = ra;
		vertex->y = dec;
		vertex->userData = userData;
		_vertices.push_back(vertex);
	}
	
	size_t TriangleCount() const { return _triangles.size(); }
	
	Triangle GetTriangle(size_t triangleIndex) const
	{
		TriangleNode* node = _triangles[triangleIndex];
		Triangle t;
		for(size_t i=0; i!=3; ++i)
		{
			const Vertex& v = *node->vertices[i];
			t.x[i] = v.x;
			t.y[i] = v.y;
			t.userData[i] = v.userData;
		}
		return t;
	}
	
	size_t ConvexVerticesCount() const { return _convexHull.size(); }
	
	ConvexVertex GetConvexVertex(size_t vertexIndex) const
	{
		ConvexVertex cv;
		Vertex* v = _convexHull[vertexIndex]->vertices[0];
		cv.x = v->x;
		cv.y = v->y;
		cv.userData = v->userData;
		return cv;
	}
	
	void Triangulate()
	{
		if(_vertices.size() < 2) return;
		
		for(size_t i=0; i!=_vertices.size(); ++i)
		{
			_vertices[i]->x = fmod(_vertices[i]->x, M_PI*2.0);
			if(_vertices[i]->x < 0.0)
				_vertices[i]->x += M_PI*2.0;
		}
		
		std::sort(_vertices.begin(), _vertices.end(), &Vertex::is_smaller);
		
		// Might phase wrap: find the largest distance and assume that's where
		// it is wrapped if it is > PI.
		double deltaX = 0.0;
		size_t phaseWrappingVertex = 0;
		for(size_t i=0; i!=_vertices.size()-1; ++i)
		{
			double thisDeltaX = _vertices[i+1]->x - _vertices[i]->x;
			if(thisDeltaX > deltaX)
			{
				deltaX = thisDeltaX;
				phaseWrappingVertex = i+1;
			}
		}
		if(deltaX > M_PI*0.5)
		{
#ifdef CHECK_DELAUNAY
			std::cout << "Detected phase wrapping: dist=" << deltaX << ", vertex=" << phaseWrappingVertex << ", coord=" << kvisCoord(*_vertices[phaseWrappingVertex]) << '\n';
#endif
			for(size_t i=phaseWrappingVertex; i!=_vertices.size(); ++i)
				_vertices[i]->x -= 2.0*M_PI;
			std::sort(_vertices.begin(), _vertices.end(), &Vertex::is_smaller);
		}
#ifdef CHECK_DELAUNAY
		else {
			std::cout << "No phase wrapping: max dist=" << deltaX << '\n';
		}
#endif
		
		// Start by building the convex hull
		std::vector<Vertex*> convexVertices;
		convexVertices.push_back(_vertices[0]);
		convexVertices.push_back(_vertices[1]);
		for(size_t i=2; i!=_vertices.size(); ++i)
		{
			while(_vertices[i]->IsLeftTurn(convexVertices))
				convexVertices.pop_back();
			convexVertices.push_back(_vertices[i]);
		}
		for(std::vector<Vertex*>::const_reverse_iterator i=_vertices.rbegin()+1; i!=_vertices.rend(); ++i)
		{
			while((*i)->IsLeftTurn(convexVertices))
				convexVertices.pop_back();
			convexVertices.push_back(*i);
		}
		// Initial point will have been added twice
		convexVertices.pop_back();
		
		convexVertices.front()->inConvexHull = true;
		for(std::vector<Vertex*>::const_iterator i=convexVertices.begin()+1; i!=convexVertices.end(); ++i)
		{
			(*i)->inConvexHull = true;
			
			addConvex(*(i-1), *i);
		}
		addConvex(convexVertices.back(), convexVertices.front());
		
		// Make a bounding triangle
		Vertex *topVertex = _vertices[0], *bottomVertex = _vertices[0];
		size_t topVertexIndex = 0;
		for(size_t i=1; i!=_vertices.size(); ++i)
		{
			Vertex* v= _vertices[i];
			if(v->y > topVertex->y)
			{
				topVertex = v;
				topVertexIndex = i;
			}
			if(v->y < bottomVertex->y)
				bottomVertex = v;
		}
		Vertex *rightOfTop = 0;
		for(size_t i=topVertexIndex+1; i<_vertices.size(); ++i)
		{
			if((rightOfTop == 0 || rightOfTop->y < _vertices[i]->y) && _vertices[i]->inConvexHull)
				rightOfTop = _vertices[i];
		}
		if(rightOfTop == 0) rightOfTop = topVertex;
		double
			margin = 1.0,
			eTopX = topVertex->x,
			eTopY = topVertex->y,
			eTopXR = rightOfTop->x,
			eTopYR = rightOfTop->y,
			eBottomY = bottomVertex->y - margin,
			eRightX = eTopX + (eBottomY-eTopY)*(eTopXR-eTopX)/(eTopYR-eTopY) + margin;
		
#ifdef CHECK_DELAUNAY
		Vertex* nearestBottom = 0;
#endif
		double smallestAngle = 1e100;
		for(size_t i=0; i!=_vertices.size(); ++i)
		{
			if(_vertices[i]->inConvexHull)
			{
				double coef = (_vertices[i]->x - eRightX) / (_vertices[i]->y - eBottomY);
				if(smallestAngle > coef)
				{
#ifdef CHECK_DELAUNAY
					nearestBottom = _vertices[i];
#endif
					smallestAngle = coef;
				}
			}
		}
#ifdef CHECK_DELAUNAY
		std::cout << "Left and right: " <<
			nearestBottom->x << ',' <<
			nearestBottom->y << "  " <<
			rightOfTop->x << ',' <<
			rightOfTop->y << '\n';
#endif
		_boundTop = topVertex;
		_boundLeft.x = eRightX + (eTopY - eBottomY) * smallestAngle - margin;
		_boundLeft.y = eTopY;
		_boundLeft.outsideConvexHull = true;
		_boundRight.x = eRightX; _boundRight.y = eBottomY;
		_boundRight.outsideConvexHull = true;
#ifdef CHECK_DELAUNAY
		std::cout << "Bounding triangle: "
			<< _boundLeft.x << ',' << _boundLeft.y << "  "
			<< _boundRight.x << ',' << _boundRight.y << "  "
			<< _boundTop->x << ',' << _boundTop->y << '\n';
#endif
		_topNode = makeTriangle(0, &_boundLeft, &_boundRight, _boundTop);
#ifdef CHECK_DELAUNAY
		if(_topNode->IsClockwise())
			std::cout << "Top node is CW!\n";
#endif
		for(std::vector<Vertex*>::iterator v=_vertices.begin(); v!=_vertices.end(); ++v)
		{
			if(*v != _boundTop)
			{
				TriangleNode* node = _topNode->GetNode((*v)->x, (*v)->y);
				for(size_t i=0; i!=3; ++i)
				{
					node->subNodes[i] =
						makeTriangle(node, node->vertices[i], node->vertices[(i+1)%3], *v);
#ifdef CHECK_DELAUNAY
					if(node->subNodes[i]->IsClockwise())
						std::cout << "Made CW node!\n";
#endif
				}
			}
		}
		
#ifdef CHECK_DELAUNAY
		SaveTriangulationAsKvis("notconvex.ann");
#endif
		makeLeafTriangleList();
		
		std::vector<Edge> edges;
		getEdges(edges);
	
		std::stack<Edge*> nonDelaunay;
		for(std::vector<Edge>::iterator e=edges.begin(); e!=edges.end(); ++e)
		{
			if(!isLocallyDelaunay(*e))
				nonDelaunay.push(new Edge(*e));
		}
#ifdef CHECK_DELAUNAY
		std::cout << "Non-delaunay edges: " << nonDelaunay.size() << '\n';
		size_t nFlips = 0;
#endif
		while(!nonDelaunay.empty())
		{
			Edge* edge = nonDelaunay.top();
			nonDelaunay.pop();
			if(!isLocallyDelaunay(*edge))
			{
				// Replace by the edge connecting the respective
				// third vertices of the two incident triangles
				Vertex *vertices[4];
				bool isFlipped = flip(*edge, vertices);
				
				// Push other four edges of the two triangles into the stack if unmarked
				if(isFlipped)
				{
#ifdef CHECK_DELAUNAY
					++nFlips;
#endif
					for(size_t i=0; i!=4; ++i)
					{
						nonDelaunay.push(new Edge(vertices[i], vertices[(i+1)%4]));
					}
				}
			}
			delete edge;
		}
#ifdef CHECK_DELAUNAY
		std::cout << "Number of flips: " << nFlips << '\n';		
		SaveTriangulationAsKvis("flipped.ann");
#endif
		
		removeBoundingVertices();
		makeLeafTriangleList();
	}
	
	void SaveConvexHullAsKvis(const std::string& filename)
	{
		std::ofstream file(filename.c_str());
		std::ofstream file2((filename+".data").c_str());
		file << "COLOR RED\n";
		for(std::vector<Edge*>::const_iterator e=_convexHull.begin(); e!=_convexHull.end(); ++e)
		{
			file << "LINE " << kvisCoord(*(*e)->vertices[0]) << ' '
				<< kvisCoord(*(*e)->vertices[1]) << '\n';
			file2 << (*e)->vertices[0]->x << '\t' << (*e)->vertices[0]->y << '\n';
			file2 << (*e)->vertices[1]->x << '\t' << (*e)->vertices[1]->y << "\n\n";
		}
	}
	
	void SaveTriangulationAsKvis(const std::string& filename)
	{
		std::ofstream file(filename.c_str());
		std::ofstream file2((filename+".data").c_str());
		file << "COLOR BLUE\n";
		std::vector<TriangleNode*> nodeStack;
		nodeStack.push_back(_topNode);
		while(!nodeStack.empty())
		{
			TriangleNode* node = nodeStack.back();
			nodeStack.pop_back();
			if(node->HasSubNodes())
			{
				for(size_t i=0; i!=3; ++i)
				{
					if(node->subNodes[i] != 0)
						nodeStack.push_back(node->subNodes[i]);
				}
			}
			else {
				for(size_t i=0; i!=3; ++i) {
					file << "LINE " << kvisCoord(*node->vertices[i]) << ' '
						<< kvisCoord(*node->vertices[(i+1)%3]) << '\n';
					file2 << node->vertices[i]->x << '\t'
						<< node->vertices[i]->y << "\n";
				}
				file2 << node->vertices[0]->x << '\t'
					<< node->vertices[0]->y << "\n";
				file2 << '\n';
			}
		}
	}
	private:
	std::string kvisCoord(const Vertex& vertex)
	{
		double ra = vertex.x*180.0/M_PI;
		if(ra < 0.0) ra += 360.0;
		std::ostringstream str;
		str << fmod(ra, 360.0) << ' ' << (vertex.y*180.0/M_PI);
		return str.str();
	}
	
	static TriangleNode* makeTriangle(TriangleNode* parent, Vertex* a, Vertex* b, Vertex* c)
	{
		TriangleNode* node = new TriangleNode();
		node->vertices[0] = a;
		node->vertices[1] = b;
		node->vertices[2] = c;
		a->triangles.push_back(node);
		b->triangles.push_back(node);
		c->triangles.push_back(node);
		node->parent = parent;
		return node;
	}
	
	bool isBound(const Vertex* v) const
	{
		return v == &_boundLeft || v == &_boundRight || v == _boundTop;
	}
	
	bool isLocallyDelaunay(const Edge& edge)
	{
		Vertex *a = edge.vertices[0], *b = edge.vertices[1];
		// Find the involved triangles
		size_t index = 0;
		TriangleNode* triangles[2];
		for(size_t i=0; i!=a->triangles.size(); ++i)
		{
			TriangleNode* t = a->triangles[i];
			if(t->HasVertex(b) && !t->HasSubNodes()) {
				triangles[index] = t;
				++index;
				if(index == 2) break;
			}
		}
		if(index != 2) return true;
		
		// Find third and fourth point
		Vertex *d=0;
		for(size_t i=0; i!=3; ++i)
		{
//			if(triangles[0]->vertices[i] != a &&
//				triangles[0]->vertices[i] != b)
//				c = triangles[0]->vertices[i];
			if(triangles[1]->vertices[i] != a &&
				triangles[1]->vertices[i] != b)
				d = triangles[1]->vertices[i];
		}
		if((isBound(a) && isBound(b)) || d==0)
			return true;
		/*if(a->outsideConvexHull || b->outsideConvexHull ||
			c->outsideConvexHull || d->outsideConvexHull)
		{
			if(a == &_boundLeft || b == &_boundLeft)
				return false;
			if(c == &_boundLeft || d == &_boundLeft)
				return true;
			else if(a == &_boundRight || b == &_boundRight)
				return false;
			else if(c == &_boundRight || d == &_boundRight)
				return true;
		}*/
		
		return isLocallyDelaunay(
			triangles[0]->vertices[2], triangles[0]->vertices[1], triangles[0]->vertices[0],
			d);
	}
	
	bool flip(const Edge& edge, Vertex** vertices)
	{
		// Find the involved triangles
		size_t index = 0;
		TriangleNode* triangles[2];
		Vertex *a = edge.vertices[0], *b = edge.vertices[1];
		for(size_t i=0; i!=a->triangles.size(); ++i)
		{
			TriangleNode* t = a->triangles[i];
			if(t->HasVertex(b) && !t->HasSubNodes()) {
				triangles[index] = t;
				++index;
				if(index == 2) break;
			}
		}
		if(index != 2)
			return false;
		// Find third and fourth point
		Vertex *c=0, *d=0;
		for(size_t i=0; i!=3; ++i)
		{
			if(triangles[0]->vertices[i] != a &&
				triangles[0]->vertices[i] != b)
				c = triangles[0]->vertices[i];
			if(triangles[1]->vertices[i] != a &&
				triangles[1]->vertices[i] != b)
				d = triangles[1]->vertices[i];
		}
#ifdef CHECK_DELAUNAY
		if(c == d) {
			std::cout << "Triangulation contains doubles!\n";
			return false;
		}
		if(c==0 || d==0) {
			std::cout << "Could not find c or d in flip()\n";
			return false;
		}
#endif
		for(size_t i=0; i!=3; ++i)
		{
			if(triangles[0]->vertices[i] == a)
			{
				triangles[0]->vertices[i]->Disconnect(triangles[0]);
				triangles[0]->vertices[i] = d;
				d->triangles.push_back(triangles[0]);
			}
			if(triangles[1]->vertices[i] == b)
			{
				triangles[1]->vertices[i]->Disconnect(triangles[1]);
				triangles[1]->vertices[i] = c;
				c->triangles.push_back(triangles[1]);
			}
		}
#ifdef CHECK_DELAUNAY
		if(triangles[0]->IsClockwise() || triangles[1]->IsClockwise())
			std::cout << "Created clockwise triangle!\n";
#endif
		
		vertices[0] = edge.vertices[0];
		vertices[1] = d;
		vertices[2] = edge.vertices[1];
		vertices[3] = c;
		return true;
	}
	
	static bool isLocallyDelaunay(const Vertex* a, const Vertex* b, const Vertex* c, const Vertex* d)
	{
		return circumDeterminant(a, b, c, d) >= 0.0;
	}

	static double circumDeterminant(const Vertex* a, const Vertex* b, const Vertex* c, const Vertex* d)
	{
		double
			ma = a->x - d->x,
			mb = a->y - d->y,
			mc = a->x*a->x - d->x*d->x + a->y*a->y - d->y*d->y,
			md = b->x - d->x,
			me = b->y - d->y,
			mf = b->x*b->x - d->x*d->x + b->y*b->y - d->y*d->y,
			mg = c->x - d->x,
			mh = c->y - d->y,
			mi = c->x*c->x - d->x*d->x + c->y*c->y - d->y*d->y;
		return
			(ma * me * mi + mb * mf * mg + mc * md * mh) -
			(mc * me * mg + mb * md * mi + ma * mf * mh);
	}
	
	void addConvex(Vertex* a, Vertex* b)
	{
		Edge *edge = new Edge(a, b);
		_convexHull.push_back(edge);
	}
	
	void makeLeafTriangleList()
	{
		_triangles.clear();
		std::vector<TriangleNode*> nodeStack;
		nodeStack.push_back(_topNode);
		while(!nodeStack.empty())
		{
			TriangleNode* node = nodeStack.back();
			nodeStack.pop_back();
			if(node->HasSubNodes())
			{
				for(size_t i=0; i!=3; ++i)
				{
					if(node->subNodes[i] != 0)
						nodeStack.push_back(node->subNodes[i]);
				}
			}
			else
				_triangles.push_back(node);
		}
	}
	
	void getEdges(std::vector<Edge>& edges)
	{
		for(std::vector<TriangleNode*>::const_iterator i=_triangles.begin(); i!=_triangles.end(); ++i)
		{
			TriangleNode* t = *i;
			edges.push_back(Edge(t->vertices[0], t->vertices[1]));
			edges.push_back(Edge(t->vertices[1], t->vertices[2]));
			edges.push_back(Edge(t->vertices[2], t->vertices[0]));
		}
	}
	
	void removeBoundingVertices()
	{
		Vertex* bounds[2] = { &_boundLeft, &_boundRight };
		bool didErase = false;
		do {
			for(size_t vIndex=0; vIndex!=2; ++vIndex)
			{
				Vertex* v = bounds[vIndex];
				std::stack<size_t> erasables;
				for(size_t i=0; i!=v->triangles.size(); ++i)
				{
					TriangleNode* t = v->triangles[i];
					if(!t->HasSubNodes()) {
#ifdef CHECK_DELAUNAY
						if(t->HasInsideVertex()) {
							std::cout << "Removing inside vertex: ";
							t->Print();
						}
#endif
						erasables.push(i);
					}
				}
				didErase = !erasables.empty();
				while(!erasables.empty())
				{
					delete v->triangles[erasables.top()];
					erasables.pop();
				}
			}
		} while(didErase);
	}
	
	// Remove copy constructor
	Delaunay(const Delaunay&) = delete;
	
	// Remove assignment operator
	void operator=(const Delaunay&) = delete;
	
	std::vector<Vertex*> _vertices;
	std::vector<Edge*> _convexHull;
	Vertex _boundLeft, *_boundTop, _boundRight;
	TriangleNode *_topNode;
	std::vector<TriangleNode*> _triangles;
};

#endif
