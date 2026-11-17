/*
 * Poly2Tri Copyright (c) 2009-2022, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "sweep_context.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include "advancing_front.h"

namespace p2t {

namespace {

// ---------------------------------------------------------------------------
// Deterministic vertex ordering and equality
//
// The original implementation used std::sort with the equivalence-only
// comparator p2t::cmp (equal-coordinate points are equivalent, so their
// relative order is implementation-defined) and an std::unordered_set for
// de-duplication (whose traversal/insertion layout depends on the standard
// library implementation and allocator addresses). Both make the output
// depend on the platform rather than only on the input.
//
// Explicit rules used here:
//   1. ordering: quantized coordinates (y, then x) ascending;
//   2. ties:     insertion order (guaranteed by std::stable_sort);
//   3. equality: quantized coordinates (see SameVertex for the exact-mode
//                compatibility rule).
//
// Quantization precision is configurable without changing any public
// signature through the P2T_QUANTIZATION_PRECISION environment variable:
//   unset/0 -> exact mode (default; output is vertex-identical to before);
//   1..17   -> coordinates are rounded to that many decimal places before
//              comparison, so vertices coinciding within the quantization
//              grid are treated as the same vertex.
// Illegal values make the first triangulation call throw std::runtime_error.
// ---------------------------------------------------------------------------

const int kMaxQuantizationPrecision = 17;

// Reads and validates the precision on every call so that different CDT
// instances in one process can use different settings. Any illegal value
// throws std::runtime_error, which propagates out of the triangulation.
int QuantizationPrecision()
{
  const char* raw = std::getenv("P2T_QUANTIZATION_PRECISION");
  if (raw == nullptr) {
    return 0;
  }
  const std::string value(raw);
  if (value.empty()) {
    throw std::runtime_error(
        "P2T_QUANTIZATION_PRECISION: empty value, expected an integer in "
        "[0, 17]");
  }
  std::size_t i = 0;
  if (value[i] == '+' || value[i] == '-') {
    ++i;
  }
  if (i == value.size()) {
    throw std::runtime_error(
        "P2T_QUANTIZATION_PRECISION: not an integer: '" + value + "'");
  }
  for (std::size_t k = i; k < value.size(); ++k) {
    if (value[k] < '0' || value[k] > '9') {
      throw std::runtime_error(
          "P2T_QUANTIZATION_PRECISION: not an integer: '" + value + "'");
    }
  }
  int parsed;
  try {
    parsed = std::stoi(value);
  } catch (const std::out_of_range&) {
    throw std::runtime_error(
        "P2T_QUANTIZATION_PRECISION: value out of range [0, 17]: '" + value
        + "'");
  }
  if (parsed < 0 || parsed > kMaxQuantizationPrecision) {
    throw std::runtime_error(
        "P2T_QUANTIZATION_PRECISION: value out of range [0, 17]: '" + value
        + "'");
  }
  return parsed;
}

// Exact-mode key: an unsigned 64-bit value whose integer ordering is
// identical to the floating-point ordering. IEEE-754 doubles are ordered by
// their raw bits interpreted as sign-magnitude; flipping the sign bit of
// positives and all bits of negatives yields an unsigned monotone key.
// -0.0 is normalized to +0.0 so both share one key, matching the original
// p2t::cmp equivalence class (it compared doubles with <, where -0.0 ==
// +0.0); the insertion-order tie break then resolves their relative order.
unsigned long long OrderPreservingKey(double v)
{
  unsigned long long bits;
  std::memcpy(&bits, &v, sizeof(bits));
  if (bits == 0x8000000000000000ULL) { // negative zero
    bits = 0ULL;
  }
  const unsigned long long mask = 0ULL - (bits >> 63);
  return bits ^ (mask | (1ULL << 63));
}

// Unsigned, monotone quantization key. In exact mode it is the raw ordered
// bit key. In grid mode signed llround values are shifted into unsigned
// order by a fixed bias (only meaningful for |coordinate * 10^p| < 2^63,
// which holds for every realistic geometric input).
unsigned long long QuantizedKey(double v)
{
  const int precision = QuantizationPrecision();
  if (precision == 0) {
    return OrderPreservingKey(v);
  }
  const long long scaled = std::llround(v * std::pow(10.0, precision));
  return static_cast<unsigned long long>(scaled) + (1ULL << 63);
}

// Strict weak ordering on quantized coordinates only. Equivalent elements
// keep their relative order under std::stable_sort, which supplies the
// insertion-order tie break.
struct DeterministicPointLess {
  bool operator()(const Point* a, const Point* b) const
  {
    const unsigned long long ay = QuantizedKey(a->y);
    const unsigned long long by = QuantizedKey(b->y);
    if (ay != by) {
      return ay < by;
    }
    const unsigned long long ax = QuantizedKey(a->x);
    const unsigned long long bx = QuantizedKey(b->x);
    return ax < bx;
  }
};

// Explicit vertex-equality rule. In exact mode (precision 0) the historical
// pointer-identity semantics are preserved, so default output stays
// vertex-identical; with a configured quantization grid, distinct objects
// whose coordinates quantize to the same cell count as one vertex.
bool SameVertex(const Point* a, const Point* b)
{
  if (QuantizationPrecision() == 0) {
    return a == b;
  }
  return QuantizedKey(a->x) == QuantizedKey(b->x)
         && QuantizedKey(a->y) == QuantizedKey(b->y);
}

bool AlreadyPresent(const std::vector<Point*>& points, const Point* candidate)
{
  for (const Point* p : points) {
    if (SameVertex(p, candidate)) {
      return true;
    }
  }
  return false;
}

} // namespace

SweepContext::SweepContext(std::vector<Point*> polyline) : points_(std::move(polyline)),
  front_(nullptr),
  head_(nullptr),
  tail_(nullptr),
  af_head_(nullptr),
  af_middle_(nullptr),
  af_tail_(nullptr)
{
  InitEdges(points_);
}

void SweepContext::AddHole(const std::vector<Point*>& polyline)
{
  InitEdges(polyline);
  for (auto i : polyline) {
    // Deterministic de-duplication: linear scan in insertion order instead
    // of an address-ordered hash set.
    if (!AlreadyPresent(points_, i)) {
      points_.push_back(i);
    }
  }
}

void SweepContext::AddPoint(Point* point) {
  if (AlreadyPresent(points_, point)) {
    throw std::runtime_error("Point already exists");
  }
  points_.push_back(point);
}

std::vector<Triangle*> &SweepContext::GetTriangles()
{
  return triangles_;
}

std::list<Triangle*> &SweepContext::GetMap()
{
  return map_;
}

void SweepContext::InitTriangulation()
{
  double xmax(points_[0]->x), xmin(points_[0]->x);
  double ymax(points_[0]->y), ymin(points_[0]->y);

  // Calculate bounds.
  for (auto& point : points_) {
    Point& p = *point;
    if (p.x > xmax)
      xmax = p.x;
    if (p.x < xmin)
      xmin = p.x;
    if (p.y > ymax)
      ymax = p.y;
    if (p.y < ymin)
      ymin = p.y;
  }

  double dx = kAlpha * (xmax - xmin);
  double dy = kAlpha * (ymax - ymin);
  head_ = new Point(xmin - dx, ymin - dy);
  tail_ = new Point(xmax + dx, ymin - dy);

  // Deterministic sweep order: quantized coordinates first, insertion order
  // for equivalent points. std::stable_sort makes the tie break explicit and
  // independent of any standard-library/allocator implementation details.
  std::stable_sort(points_.begin(), points_.end(), DeterministicPointLess());

}

void SweepContext::InitEdges(const std::vector<Point*>& polyline)
{
  size_t num_points = polyline.size();
  for (size_t i = 0; i < num_points; i++) {
    size_t j = i < num_points - 1 ? i + 1 : 0;
    edge_list.push_back(new Edge(*polyline[i], *polyline[j]));
  }
}

Point* SweepContext::GetPoint(size_t index)
{
  return points_[index];
}

void SweepContext::AddToMap(Triangle* triangle)
{
  map_.push_back(triangle);
}

Node* SweepContext::LocateNode(const Point& point)
{
  // TODO implement search tree
  return front_->LocateNode(point.x);
}

void SweepContext::CreateAdvancingFront()
{

  // Initial triangle
  Triangle* triangle = new Triangle(*points_[0], *head_, *tail_);

  map_.push_back(triangle);

  af_head_ = new Node(*triangle->GetPoint(1), *triangle);
  af_middle_ = new Node(*triangle->GetPoint(0), *triangle);
  af_tail_ = new Node(*triangle->GetPoint(2));
  front_ = new AdvancingFront(*af_head_, *af_tail_);

  // TODO: More intuitive if head is middles next and not previous?
  //       so swap head and tail
  af_head_->next = af_middle_;
  af_middle_->next = af_tail_;
  af_middle_->prev = af_head_;
  af_tail_->prev = af_middle_;
}

void SweepContext::RemoveNode(Node* node)
{
  delete node;
}

void SweepContext::MapTriangleToNodes(Triangle& t)
{
  for (int i = 0; i < 3; i++) {
    if (!t.GetNeighbor(i)) {
      Node* n = front_->LocatePoint(t.PointCW(*t.GetPoint(i)));
      if (n)
        n->triangle = &t;
    }
  }
}

void SweepContext::RemoveFromMap(Triangle* triangle)
{
  map_.remove(triangle);
}

void SweepContext::MeshClean(Triangle& triangle)
{
  std::vector<Triangle *> triangles;
  triangles.push_back(&triangle);

  while(!triangles.empty()){
	Triangle *t = triangles.back();
	triangles.pop_back();

    if (t != nullptr && !t->IsInterior()) {
      t->IsInterior(true);
      triangles_.push_back(t);
      for (int i = 0; i < 3; i++) {
        if (!t->constrained_edge[i])
          triangles.push_back(t->GetNeighbor(i));
      }
    }
  }
}

SweepContext::~SweepContext()
{

    // Clean up memory

    delete head_;
    delete tail_;
    delete front_;
    delete af_head_;
    delete af_middle_;
    delete af_tail_;

    for (auto ptr : map_) {
      delete ptr;
    }

    for (auto& i : edge_list) {
      delete i;
    }
}

} // namespace p2t
