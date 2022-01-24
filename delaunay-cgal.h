#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

class Delaunay {
 public:
  struct ConvexVertex {
    double x, y;
    void* userData;
  };
  struct Triangle {
    double x[3], y[3];
    void* userData[3];
  };

  void Clear();
  void AddVertex(double ra, double dec, void* userData = 0) {
    if (!std::isfinite(ra) || !std::isfinite(dec))
      throw std::runtime_error("Non-finite value in vector coordinates");
    Vertex v;
    v.x = ra;
    v.y = dec;
    v.userData = userData;
    _vertices.push_back(v);
  }
  void Triangulate();
  void SaveConvexHullAsKvis(const std::string& filename);
  void SaveTriangulationAsKvis(const std::string& filename);

  size_t TriangleCount();
  Triangle GetTriangle(size_t triangleIndex) const;
  size_t ConvexVerticesCount();
  ConvexVertex GetConvexVertex(size_t vertexIndex) const;

 private:
  struct Vertex {
    double x, y;
    void* userData;
  };
  std::vector<Vertex> _vertices;
};

#endif
