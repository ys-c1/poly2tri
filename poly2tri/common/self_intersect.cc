/*
 * Poly2Tri Copyright (c) 2009-2018, Poly2Tri Contributors
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

#include "self_intersect.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace p2t {
namespace {

// ---------------------------------------------------------------------------
// Shared geometry helpers (used by both the sweep and the brute-force oracle,
// so that both compute intersection points with the exact same arithmetic).
// ---------------------------------------------------------------------------

struct Seg {
  // Original endpoints as taken from the polyline.
  double ax, ay, bx, by;
  // Ordered endpoints: (x1,y1) <= (x2,y2) lexicographically by (x,y).
  double x1, y1, x2, y2;
  double slope;    // dy/dx; +inf for vertical segments.
  bool vertical;   // x1 == x2 exactly.
  std::size_t edge; // Original edge index in the outline.

  // y at sweep position x. Snaps to the exact endpoint ordinates when x
  // coincides with an endpoint abscissa (events always do), which keeps the
  // status ordering consistent at event points.
  double y_at(double x) const
  {
    if (x == x1) return y1;
    if (x == x2) return y2;
    return y1 + (x - x1) * slope;
  }
};

inline double Cross(double ox, double oy, double px, double py, double qx, double qy)
{
  return (px - ox) * (qy - oy) - (py - oy) * (qx - ox);
}

// Intersection point of the two supporting lines (non-parallel segments).
inline void LineIntersect(const Seg& s, const Seg& t, double& ix, double& iy)
{
  const double rx = s.bx - s.ax, ry = s.by - s.ay;
  const double sx = t.bx - t.ax, sy = t.by - t.ay;
  const double den = rx * sy - ry * sx;
  const double tt = ((t.ax - s.ax) * sy - (t.ay - s.ay) * sx) / den;
  ix = s.ax + tt * rx;
  iy = s.ay + tt * ry;
}

inline bool OnSegment(const Seg& s, double px, double py)
{
  return px >= std::min(s.ax, s.bx) && px <= std::max(s.ax, s.bx) &&
         py >= std::min(s.ay, s.by) && py <= std::max(s.ay, s.by);
}

inline bool Collinear(const Seg& s, const Seg& t)
{
  return Cross(s.ax, s.ay, s.bx, s.by, t.ax, t.ay) == 0.0 &&
         Cross(s.ax, s.ay, s.bx, s.by, t.bx, t.by) == 0.0;
}

// Result of a segment-pair intersection test.
enum PairResult { kNone, kPoint, kOverlap };

// Full segment-segment intersection test with degenerate handling.
// On kPoint/kOverlap, (ix,iy) is the reported point (for kOverlap: the
// lexicographically smallest point of the shared interval).
PairResult IntersectPair(const Seg& s, const Seg& t, double& ix, double& iy)
{
  const double o1 = Cross(s.ax, s.ay, s.bx, s.by, t.ax, t.ay);
  const double o2 = Cross(s.ax, s.ay, s.bx, s.by, t.bx, t.by);
  const double o3 = Cross(t.ax, t.ay, t.bx, t.by, s.ax, s.ay);
  const double o4 = Cross(t.ax, t.ay, t.bx, t.by, s.bx, s.by);

  if (o1 == 0.0 && o2 == 0.0) {
    // Collinear: measure the 1D overlap along x (along y for vertical
    // segments), so the reported point is the lexicographically smallest
    // point of the shared interval.
    const bool use_x = (s.ax != s.bx);
    const double sa = use_x ? s.ax : s.ay, sb = use_x ? s.bx : s.by;
    const double ta = use_x ? t.ax : t.ay, tb = use_x ? t.bx : t.by;
    const double lo = std::max(std::min(sa, sb), std::min(ta, tb));
    const double hi = std::min(std::max(sa, sb), std::max(ta, tb));
    if (lo > hi) return kNone;
    // The interval start is the "leading" endpoint: the rightmost/uppermost
    // of the two left/lower endpoints. Report its exact coordinates.
    const Seg& lead = (std::min(sa, sb) >= std::min(ta, tb)) ? s : t;
    const bool a_is_min =
        use_x ? (lead.ax <= lead.bx) : (lead.ay <= lead.by);
    ix = a_is_min ? lead.ax : lead.bx;
    iy = a_is_min ? lead.ay : lead.by;
    return (lo < hi) ? kOverlap : kPoint;
  }

  if ((o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0) && o1 != 0.0 && o2 != 0.0 &&
      o3 != 0.0 && o4 != 0.0) {
    // Proper crossing.
    LineIntersect(s, t, ix, iy);
    return kPoint;
  }

  // Touching: an endpoint of one segment lying on the other (non-collinear,
  // so there is at most one such point).
  double px = 0.0, py = 0.0;
  bool found = false;
  const struct { const Seg* seg; double x, y; } ends[4] = {
      {&t, t.ax, t.ay}, {&t, t.bx, t.by}, {&s, s.ax, s.ay}, {&s, s.bx, s.by}};
  const Seg* hosts[4] = {&s, &s, &t, &t};
  for (int i = 0; i < 4; ++i) {
    const double cx = ends[i].x, cy = ends[i].y;
    const Seg& host = *hosts[i];
    if (Cross(host.ax, host.ay, host.bx, host.by, cx, cy) == 0.0 &&
        OnSegment(host, cx, cy)) {
      if (!found) { px = cx; py = cy; found = true; }
      // A shared endpoint is found twice; same point, ignore the duplicate.
    }
  }
  if (found) { ix = px; iy = py; return kPoint; }
  return kNone;
}

// Proper-crossing point only (used to schedule sweep events). Returns false
// if the segments do not properly cross.
bool ProperCross(const Seg& s, const Seg& t, double& ix, double& iy)
{
  const double o1 = Cross(s.ax, s.ay, s.bx, s.by, t.ax, t.ay);
  const double o2 = Cross(s.ax, s.ay, s.bx, s.by, t.bx, t.by);
  const double o3 = Cross(t.ax, t.ay, t.bx, t.by, s.ax, s.ay);
  const double o4 = Cross(t.ax, t.ay, t.bx, t.by, s.bx, s.by);
  if (o1 == 0.0 || o2 == 0.0 || o3 == 0.0 || o4 == 0.0) return false;
  if ((o1 > 0) == (o2 > 0) || (o3 > 0) == (o4 > 0)) return false;
  LineIntersect(s, t, ix, iy);
  return true;
}

// Build the edge list, dropping zero-length edges.
std::vector<Seg> BuildSegments(const std::vector<Point>& polyline, bool closed,
                               std::size_t& num_edges)
{
  const std::size_t n = polyline.size();
  num_edges = closed ? n : (n > 0 ? n - 1 : 0);
  std::vector<Seg> segs;
  for (std::size_t i = 0; i < num_edges; ++i) {
    const Point& a = polyline[i];
    const Point& b = polyline[(i + 1) % n];
    if (a.x == b.x && a.y == b.y) continue; // zero-length edge: ignored
    Seg s;
    s.ax = a.x; s.ay = a.y; s.bx = b.x; s.by = b.y;
    if ((a.x < b.x) || (a.x == b.x && a.y < b.y)) {
      s.x1 = a.x; s.y1 = a.y; s.x2 = b.x; s.y2 = b.y;
    } else {
      s.x1 = b.x; s.y1 = b.y; s.x2 = a.x; s.y2 = a.y;
    }
    s.vertical = (s.x1 == s.x2);
    s.slope = s.vertical ? std::numeric_limits<double>::infinity()
                         : (s.y2 - s.y1) / (s.x2 - s.x1);
    s.edge = i;
    segs.push_back(s);
  }
  return segs;
}

inline bool EdgesAdjacent(std::size_t i, std::size_t j, std::size_t num_edges,
                          bool closed)
{
  const std::size_t d = i > j ? i - j : j - i;
  if (d == 1) return true;
  return closed && d == num_edges - 1;
}

inline std::uint64_t PairKey(std::size_t i, std::size_t j)
{
  if (i > j) std::swap(i, j);
  return (static_cast<std::uint64_t>(i) << 32) | static_cast<std::uint64_t>(j);
}

void SortOutput(std::vector<SelfIntersection>& out)
{
  std::sort(out.begin(), out.end(),
            [](const SelfIntersection& a, const SelfIntersection& b) {
              if (a.x != b.x) return a.x < b.x;
              if (a.y != b.y) return a.y < b.y;
              if (a.edge_a != b.edge_a) return a.edge_a < b.edge_a;
              return a.edge_b < b.edge_b;
            });
}

// ---------------------------------------------------------------------------
// Brute-force reference (O(n^2)).
// ---------------------------------------------------------------------------

bool BruteForce(const std::vector<Point>& polyline, bool closed,
                std::size_t max_report, std::vector<SelfIntersection>& out,
                bool* truncated)
{
  std::size_t num_edges = 0;
  const std::vector<Seg> segs = BuildSegments(polyline, closed, num_edges);
  std::vector<SelfIntersection> all;
  for (std::size_t i = 0; i < segs.size(); ++i) {
    for (std::size_t j = i + 1; j < segs.size(); ++j) {
      if (EdgesAdjacent(segs[i].edge, segs[j].edge, num_edges, closed)) continue;
      double ix, iy;
      const PairResult r = IntersectPair(segs[i], segs[j], ix, iy);
      if (r == kNone) continue;
      SelfIntersection rec;
      rec.x = ix; rec.y = iy;
      rec.edge_a = std::min(segs[i].edge, segs[j].edge);
      rec.edge_b = std::max(segs[i].edge, segs[j].edge);
      rec.overlap = (r == kOverlap);
      all.push_back(rec);
    }
  }
  SortOutput(all);
  bool tr = false;
  if (all.size() > max_report) {
    all.resize(max_report);
    tr = true;
  }
  out = all;
  if (truncated) *truncated = tr;
  return !out.empty();
}

// ---------------------------------------------------------------------------
// Bentley-Ottmann sweep line.
// ---------------------------------------------------------------------------

struct Event {
  double x, y;
  int type;   // 0 = END, 1 = CROSS, 2 = START, 3 = DUMMY (vertical marker)
  int s1, s2; // segment slots; CROSS uses both, DUMMY none.
  bool operator<(const Event& o) const
  {
    if (x != o.x) return x < o.x;
    if (y != o.y) return y < o.y;
    if (type != o.type) return type < o.type;
    if (s1 != o.s1) return s1 < o.s1;
    return s2 < o.s2;
  }
};

class Sweep {
public:
  Sweep(const std::vector<Point>& polyline, bool closed, std::size_t max_report)
      : closed_(closed), max_report_(max_report)
  {
    segs_ = BuildSegments(polyline, closed, num_edges_);
    // Extra slot used as an ordering probe for range queries.
    Seg probe;
    probe.ax = probe.ay = probe.bx = probe.by = 0.0;
    probe.x1 = probe.y1 = probe.x2 = probe.y2 = 0.0;
    probe.slope = 0.0;
    probe.vertical = false;
    probe.edge = 0;
    segs_.push_back(probe);
    probe_slot_ = static_cast<int>(segs_.size() - 1);

    cmp_.segs = &segs_;
    cmp_.cur_x = &cur_x_;
    cmp_.cur_y = &cur_y_;
    cmp_.at_point = &at_point_;
    at_point_.assign(segs_.size(), 0);
    status_ = std::set<int, SegCmp>(cmp_);
    pos_.assign(segs_.size(), status_.end());
    in_status_.assign(segs_.size(), 0);

    for (int i = 0; i < static_cast<int>(segs_.size()) - 1; ++i) {
      const Seg& s = segs_[i];
      if (s.vertical) {
        verticals_[s.x1].push_back(i);
        // Dummy events only make the sweep visit this x coordinate.
        pq_.insert(Event{s.x1, s.y1, 3, i, -1});
        pq_.insert(Event{s.x2, s.y2, 3, i, -1});
      } else {
        pq_.insert(Event{s.x1, s.y1, 2, i, -1});
        pq_.insert(Event{s.x2, s.y2, 0, i, -1});
      }
    }
  }

  bool Run(std::vector<SelfIntersection>& out, bool* truncated)
  {
    double last_x = std::numeric_limits<double>::quiet_NaN();
    while (!pq_.empty() && !truncated_) {
      const double x = pq_.begin()->x;
      if (last_x != x) {
        last_x = x;
        cur_x_ = x;
        ProcessVerticalsAt(x);
        if (truncated_) break;
      }
      ProcessPoint();
    }
    SortOutput(out_);
    out = out_;
    if (truncated) *truncated = truncated_;
    return !out.empty();
  }

private:
  struct SegCmp {
    const std::vector<Seg>* segs;
    const double* cur_x;
    const double* cur_y;
    // Slots known to pass exactly through the current event point compare
    // with the exact event ordinate, so that their mutual order is decided
    // by the slope tie-break (the correct order just right of the event).
    const std::vector<char>* at_point;
    double YOf(int s) const
    {
      return (*at_point)[s] ? *cur_y : (*segs)[s].y_at(*cur_x);
    }
    bool operator()(int a, int b) const
    {
      if (a == b) return false;
      const double ya = YOf(a);
      const double yb = YOf(b);
      if (ya != yb) return ya < yb;
      const double sa = (*segs)[a].slope;
      const double sb = (*segs)[b].slope;
      if (sa != sb) return sa < sb;
      return a < b;
    }
  };

  bool Adjacent(std::size_t i, std::size_t j) const
  {
    return EdgesAdjacent(i, j, num_edges_, closed_);
  }

  // Returns false when the report cap was hit and the sweep must abort.
  bool Report(std::size_t ea, std::size_t eb, double x, double y, bool overlap)
  {
    if (Adjacent(ea, eb)) return true;
    if (!reported_.insert(PairKey(ea, eb)).second) return true;
    if (out_.size() >= max_report_) {
      truncated_ = true;
      reported_.erase(PairKey(ea, eb));
      return false;
    }
    SelfIntersection rec;
    rec.x = x; rec.y = y;
    rec.edge_a = std::min(ea, eb);
    rec.edge_b = std::max(ea, eb);
    rec.overlap = overlap;
    out_.push_back(rec);
    return true;
  }

  // Report a collinear pair (overlap or single-point touch) using exact
  // endpoint coordinates, matching the brute-force oracle. Returns false
  // when the report cap was hit.
  bool ReportCollinear(int a, int b)
  {
    const Seg& A = segs_[a];
    const Seg& B = segs_[b];
    const double lo = std::max(A.x1, B.x1);
    const double hi = std::min(A.x2, B.x2);
    if (lo > hi) return true;
    const Seg& lead = (A.x1 >= B.x1) ? A : B;
    return Report(A.edge, B.edge, lead.x1, lead.y1, /*overlap=*/lo < hi);
  }

  // Collinear overlap / touch report plus crossing-event scheduling for a
  // pair of non-vertical segments that just became neighbours in the status.
  void CheckPair(int a, int b)
  {
    const Seg& A = segs_[a];
    const Seg& B = segs_[b];
    if (Collinear(A, B)) {
      ReportCollinear(a, b);
      return;
    }
    if (Adjacent(A.edge, B.edge)) return;
    double ix, iy;
    if (!ProperCross(A, B, ix, iy)) return;
    // Only crossings strictly ahead of the current event point.
    if (ix < cur_x_ || (ix == cur_x_ && iy <= cur_y_)) return;
    const std::uint64_t key = PairKey(A.edge, B.edge);
    if (!scheduled_.insert(key).second) return;
    pq_.insert(Event{ix, iy, 1, a, b});
  }

  void InsertWithChecks(int s)
  {
    const Seg& S = segs_[s];
    const auto res = status_.insert(s);
    pos_[s] = res.first;
    in_status_[s] = 1;
    const auto it = res.first;
    // Walk left over the collinear run, reporting overlaps.
    auto itl = it;
    while (itl != status_.begin()) {
      auto prv = std::prev(itl);
      if (!Collinear(segs_[*prv], S)) {
        CheckPair(*prv, s);
        break;
      }
      CheckPair(*prv, s); // reports the overlap
      itl = prv;
    }
    // Walk right over the collinear run, reporting overlaps.
    auto itr = std::next(it);
    while (itr != status_.end()) {
      if (!Collinear(segs_[*itr], S)) {
        CheckPair(s, *itr);
        break;
      }
      CheckPair(s, *itr);
      ++itr;
    }
  }

  // Vertical segments at x against the active (non-vertical) status, plus
  // vertical-vertical pairs at this x.
  void ProcessVerticalsAt(double x)
  {
    auto it = verticals_.find(x);
    if (it == verticals_.end()) return;
    const std::vector<int>& verts = it->second;

    // Vertical vs. active status: range query on y_at(x).
    for (int vi : verts) {
      const Seg& V = segs_[vi];
      Seg& probe = segs_[probe_slot_];
      probe.x1 = x; probe.y1 = V.y1; probe.x2 = x + 1.0; probe.y2 = V.y1;
      probe.slope = -std::numeric_limits<double>::infinity();
      auto lo = status_.lower_bound(probe_slot_);
      for (auto sit = lo; sit != status_.end(); ++sit) {
        const Seg& S = segs_[*sit];
        const double y = S.y_at(x);
        if (y > V.y2) break;
        if (y < V.y1) continue;
        if (!Report(V.edge, S.edge, x, y, /*overlap=*/false)) return;
      }
    }

    // Vertical vs. vertical at the same x (all mutually collinear).
    std::vector<int> sorted = verts;
    std::sort(sorted.begin(), sorted.end(),
              [&](int a, int b) { return segs_[a].y1 < segs_[b].y1; });
    for (std::size_t i = 0; i < sorted.size(); ++i) {
      const Seg& A = segs_[sorted[i]];
      for (std::size_t j = i + 1; j < sorted.size(); ++j) {
        const Seg& B = segs_[sorted[j]];
        if (B.y1 > A.y2) break;
        const double lo = B.y1;
        const double hi = std::min(A.y2, B.y2);
        if (!Report(A.edge, B.edge, x, lo, /*overlap=*/lo < hi)) return;
      }
    }
  }

  // Process all events at one point (pq_.begin() determines the point).
  void ProcessPoint()
  {
    const double px = pq_.begin()->x;
    const double py = pq_.begin()->y;
    cur_x_ = px;
    cur_y_ = py;

    std::vector<int> starts, ends, cross;
    while (!pq_.empty() && pq_.begin()->x == px && pq_.begin()->y == py) {
      const Event e = *pq_.begin();
      pq_.erase(pq_.begin());
      if (e.type == 2) starts.push_back(e.s1);
      else if (e.type == 0) ends.push_back(e.s1);
      else if (e.type == 1) { cross.push_back(e.s1); cross.push_back(e.s2); }
      // DUMMY (vertical marker): nothing to do here.
    }

    // Absorb status segments passing through (px,py) in their interior:
    // they form a contiguous run of entries with y_at(px) == py. The run is
    // located with a small tolerance: the event ordinate and the segments'
    // y_at are computed with different floating-point expressions and can
    // differ by a few ulps when they are mathematically equal.
    std::vector<int> interior;
    {
      const double eps = 1e-12 * (1.0 + std::fabs(py));
      Seg& probe = segs_[probe_slot_];
      probe.x1 = px; probe.x2 = px + 1.0;
      probe.y1 = py - eps; probe.y2 = probe.y1;
      probe.slope = -std::numeric_limits<double>::infinity();
      auto lo = status_.lower_bound(probe_slot_);
      probe.y1 = py + eps; probe.y2 = probe.y1;
      probe.slope = std::numeric_limits<double>::infinity();
      auto hi = status_.lower_bound(probe_slot_);
      for (auto it = lo; it != hi; ++it) {
        const int s = *it;
        if (std::find(ends.begin(), ends.end(), s) == ends.end())
          interior.push_back(s);
      }
    }

    // All segments involved at this point.
    std::vector<int> block = starts;
    block.insert(block.end(), ends.begin(), ends.end());
    block.insert(block.end(), cross.begin(), cross.end());
    block.insert(block.end(), interior.begin(), interior.end());
    std::sort(block.begin(), block.end());
    block.erase(std::unique(block.begin(), block.end()), block.end());

    // Report every pair meeting at this point. Collinear pairs are reported
    // as overlaps (or single-point touches) at the shared interval start,
    // not as a crossing at this event point.
    for (std::size_t i = 0; i < block.size() && !truncated_; ++i) {
      for (std::size_t j = i + 1; j < block.size(); ++j) {
        if (Collinear(segs_[block[i]], segs_[block[j]])) {
          if (!ReportCollinear(block[i], block[j])) return;
        } else if (!Report(segs_[block[i]].edge, segs_[block[j]].edge, px, py,
                           /*overlap=*/false)) {
          return;
        }
      }
    }
    // Vertical segments spanning this point meet every block member.
    auto vit = verticals_.find(px);
    if (vit != verticals_.end()) {
      for (int vi : vit->second) {
        const Seg& V = segs_[vi];
        if (py < V.y1 || py > V.y2) continue;
        for (int s : block) {
          if (!Report(V.edge, segs_[s].edge, px, py, /*overlap=*/false)) return;
        }
      }
    }

    // Remove ending and crossing segments by stored iterator (no comparator
    // involved, so the pre-event physical order does not matter), then mark
    // the whole block and re-insert the continuing segments: while marked,
    // they compare with the exact event ordinate and their mutual order is
    // the slope tie-break, i.e. the correct order just right of the event.
    std::vector<int> remove = ends;
    remove.insert(remove.end(), cross.begin(), cross.end());
    remove.insert(remove.end(), interior.begin(), interior.end());
    std::sort(remove.begin(), remove.end());
    remove.erase(std::unique(remove.begin(), remove.end()), remove.end());
    for (int s : remove) {
      if (in_status_[s]) {
        status_.erase(pos_[s]);
        in_status_[s] = 0;
      }
    }
    for (int s : block) at_point_[s] = 1;

    std::vector<int> reinsert = cross;
    reinsert.insert(reinsert.end(), interior.begin(), interior.end());
    reinsert.insert(reinsert.end(), starts.begin(), starts.end());
    std::sort(reinsert.begin(), reinsert.end());
    reinsert.erase(std::unique(reinsert.begin(), reinsert.end()),
                   reinsert.end());
    for (int s : reinsert) {
      if (segs_[s].vertical) continue;
      InsertWithChecks(s);
    }

    if (reinsert.empty()) {
      // The segments above and below the removed block became neighbours.
      Seg& probe = segs_[probe_slot_];
      probe.x1 = px; probe.y1 = py; probe.x2 = px + 1.0; probe.y2 = py;
      probe.slope = std::numeric_limits<double>::infinity();
      auto above = status_.lower_bound(probe_slot_);
      if (above != status_.end() && above != status_.begin()) {
        auto below = std::prev(above);
        CheckPair(*below, *above);
      }
    }

    for (int s : block) at_point_[s] = 0;
  }

  bool closed_;
  std::size_t max_report_;
  std::size_t num_edges_ = 0;
  std::vector<Seg> segs_;
  int probe_slot_ = -1;
  double cur_x_ = 0.0, cur_y_ = 0.0;
  std::vector<char> at_point_;
  SegCmp cmp_;
  std::set<int, SegCmp> status_;
  // Stable iterators into the status, for comparator-free removal.
  std::vector<std::set<int, SegCmp>::iterator> pos_;
  std::vector<char> in_status_;
  std::set<Event> pq_;
  std::map<double, std::vector<int>> verticals_;
  std::set<std::uint64_t> reported_;
  std::set<std::uint64_t> scheduled_;
  std::vector<SelfIntersection> out_;
  bool truncated_ = false;
};

} // namespace

bool FindSelfIntersections(const std::vector<Point>& polyline, bool closed,
                           std::size_t max_report,
                           std::vector<SelfIntersection>& out, bool* truncated)
{
  if (max_report == 0) {
    // Still need the boolean answer: run with capacity for one record.
    std::vector<SelfIntersection> tmp;
    bool tr = false;
    const bool found = FindSelfIntersections(polyline, closed, 1, tmp, &tr);
    out.clear();
    if (truncated) *truncated = tr;
    return found;
  }
  Sweep sweep(polyline, closed, max_report);
  return sweep.Run(out, truncated);
}

bool FindSelfIntersectionsBruteForce(const std::vector<Point>& polyline,
                                     bool closed, std::size_t max_report,
                                     std::vector<SelfIntersection>& out,
                                     bool* truncated)
{
  if (max_report == 0) {
    std::vector<SelfIntersection> tmp;
    bool tr = false;
    const bool found =
        FindSelfIntersectionsBruteForce(polyline, closed, 1, tmp, &tr);
    out.clear();
    if (truncated) *truncated = tr;
    return found;
  }
  return BruteForce(polyline, closed, max_report, out, truncated);
}

} // namespace p2t
