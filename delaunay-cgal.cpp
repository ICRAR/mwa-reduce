#include "delaunay.h"

#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
//#include <CGAL/Vertex_circulator.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<void*, K> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> DelaunayAlg;
typedef K::Point_2 Point;

void Delaunay::Triangulate() {
  DelaunayAlg dt;
  for (std::vector<Vertex>::const_iterator i = _vertices.begin();
       i != _vertices.end(); ++i) {
    Point p(i->x, i->y);
    DelaunayAlg::Vertex_handle v = dt.insert(p);
    v->info() = i->userData;
  }

  std::cout << "Vertices=" << dt.number_of_vertices() << std::endl;

  DelaunayAlg::Vertex_circulator incident_faces(dt.infinite_vertex());
  DelaunayAlg::Vertex_circulator iter(incident_faces);
  do {
    std::cout << *iter << '\n';
    ++iter;
  } while (incident_faces != iter);
  std::cout << "End.\n";
}

void Delaunay::Clear() {}

size_t Delaunay::ConvexVerticesCount() {}

Delaunay::ConvexVertex Delaunay::GetConvexVertex(size_t vertexIndex) const {}

Delaunay::Triangle Delaunay::GetTriangle(size_t triangleIndex) const {}

void Delaunay::SaveConvexHullAsKvis(const std::__cxx11::string& filename) {}

void Delaunay::SaveTriangulationAsKvis(const std::__cxx11::string& filename) {}

size_t Delaunay::TriangleCount() {}
