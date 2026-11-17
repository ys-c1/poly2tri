// Serialization round-trip test for poly2tri triangulation results.
//
// Standalone test/script: not wired into CMake/meson, links against the
// in-repo poly2tri sources only (STL + poly2tri, no third-party deps).
//
// Build (from repo root):
//   g++ -std=c++14 -O2 -DP2T_STATIC_EXPORTS -I.
//       serialization/serialize_roundtrip.cpp
//       poly2tri/common/shapes.cc poly2tri/sweep/*.cc
//       -o artifacts/serialize_roundtrip
//
// Binary format (all integers little-endian):
//   Header (20 bytes):
//     char[4]  magic = "P2TS"
//     u16      format version (configurable at write time)
//     u16      flags (bit0 = payload is RLE-compressed; configurable)
//     u32      stored payload size (after compression)
//     u32      uncompressed payload size
//     u32      CRC32 (IEEE) of the *uncompressed* payload
//   Payload:
//     u32 n_vertices;  n * (f64 x, f64 y)
//     u32 n_triangles; n * (u32 a, u32 b, u32 c)   // indices into vertices
//     u32 n_edges;     n * (u32 a, u32 b)          // constrained edges
//     u32 n_polygons;  per polygon:
//       u32 outer_len, outer_len * u32 vertex indices
//       u32 n_holes;   per hole: u32 len, len * u32 vertex indices

#include <poly2tri/poly2tri.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace p2ts {

// ---------------------------------------------------------------- CRC32 (IEEE)
inline uint32_t Crc32(const uint8_t* data, size_t n) {
  static uint32_t table[256];
  static bool init = false;
  if (!init) {
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t c = i;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      table[i] = c;
    }
    init = true;
  }
  uint32_t c = 0xFFFFFFFFu;
  for (size_t i = 0; i < n; ++i) c = table[(c ^ data[i]) & 0xFF] ^ (c >> 8);
  return c ^ 0xFFFFFFFFu;
}

// ------------------------------------------------- PackBits-style byte RLE
inline std::vector<uint8_t> RleEncode(const std::vector<uint8_t>& in) {
  std::vector<uint8_t> out;
  size_t i = 0;
  while (i < in.size()) {
    size_t run = 1;
    while (i + run < in.size() && in[i + run] == in[i] && run < 128) ++run;
    if (run >= 3) {
      out.push_back(static_cast<uint8_t>(257 - run));  // run token
      out.push_back(in[i]);
      i += run;
    } else {
      size_t lit_start = i;
      i += run;
      while (i < in.size() && i - lit_start < 128) {
        size_t r = 1;
        while (i + r < in.size() && in[i + r] == in[i] && r < 128) ++r;
        if (r >= 3) break;
        i += r;
      }
      out.push_back(static_cast<uint8_t>(i - lit_start - 1));  // literal token
      out.insert(out.end(), in.begin() + lit_start, in.begin() + i);
    }
  }
  return out;
}

inline std::vector<uint8_t> RleDecode(const uint8_t* in, size_t n, size_t expected) {
  std::vector<uint8_t> out;
  out.reserve(expected);
  size_t i = 0;
  while (i < n) {
    uint8_t t = in[i++];
    if (t <= 127) {  // literal: t+1 bytes follow
      if (i + t + 1 > n) throw std::runtime_error("corrupt RLE stream (literal overrun)");
      out.insert(out.end(), in + i, in + i + t + 1);
      i += t + 1;
    } else if (t >= 129) {  // run: 257-t copies of next byte
      if (i >= n) throw std::runtime_error("corrupt RLE stream (run overrun)");
      out.insert(out.end(), 257 - t, in[i++]);
    }  // 128 = noop
  }
  if (out.size() != expected) throw std::runtime_error("corrupt RLE stream (size mismatch)");
  return out;
}

// ---------------------------------------------------------------- data model
struct Polygon {
  std::vector<uint32_t> outer;               // vertex indices, CCW
  std::vector<std::vector<uint32_t>> holes;  // vertex indices
};

struct Scene {
  std::vector<std::array<double, 2>> vertices;
  std::vector<Polygon> polygons;
  std::vector<std::array<uint32_t, 3>> triangles;  // triangulation result
  std::vector<std::array<uint32_t, 2>> edges;      // constrained edges
};

constexpr uint16_t kCurrentVersion = 1;

// ------------------------------------------------------------------ writer
class Writer {
 public:
  void U16(uint16_t v) { Put(v); }
  void U32(uint32_t v) { Put(v); }
  void F64(double v) {
    uint64_t bits;
    static_assert(sizeof(bits) == sizeof(v), "");
    std::memcpy(&bits, &v, 8);
    Put(bits);
  }
  std::vector<uint8_t> buf;

 private:
  template <typename T>
  void Put(T v) {
    uint8_t b[sizeof(T)];
    std::memcpy(b, &v, sizeof(T));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    std::reverse(b, b + sizeof(T));
#endif
    buf.insert(buf.end(), b, b + sizeof(T));
  }
};

inline std::vector<uint8_t> Serialize(const Scene& s, uint16_t version, bool compress) {
  Writer w;
  w.U32(static_cast<uint32_t>(s.vertices.size()));
  for (const auto& v : s.vertices) { w.F64(v[0]); w.F64(v[1]); }
  w.U32(static_cast<uint32_t>(s.triangles.size()));
  for (const auto& t : s.triangles) { w.U32(t[0]); w.U32(t[1]); w.U32(t[2]); }
  w.U32(static_cast<uint32_t>(s.edges.size()));
  for (const auto& e : s.edges) { w.U32(e[0]); w.U32(e[1]); }
  w.U32(static_cast<uint32_t>(s.polygons.size()));
  for (const auto& p : s.polygons) {
    w.U32(static_cast<uint32_t>(p.outer.size()));
    for (uint32_t i : p.outer) w.U32(i);
    w.U32(static_cast<uint32_t>(p.holes.size()));
    for (const auto& h : p.holes) {
      w.U32(static_cast<uint32_t>(h.size()));
      for (uint32_t i : h) w.U32(i);
    }
  }

  const uint32_t crc = Crc32(w.buf.data(), w.buf.size());
  const uint32_t raw_size = static_cast<uint32_t>(w.buf.size());
  std::vector<uint8_t> stored = compress ? RleEncode(w.buf) : std::move(w.buf);

  Writer h;
  h.buf = {'P', '2', 'T', 'S'};
  h.U16(version);
  h.U16(compress ? 1 : 0);
  h.U32(static_cast<uint32_t>(stored.size()));
  h.U32(raw_size);
  h.U32(crc);
  h.buf.insert(h.buf.end(), stored.begin(), stored.end());
  return h.buf;
}

// ------------------------------------------------------------------ reader
class Reader {
 public:
  Reader(const uint8_t* p, size_t n) : p_(p), n_(n) {}
  uint16_t U16() { return Get<uint16_t>(); }
  uint32_t U32() { return Get<uint32_t>(); }
  double F64() {
    uint64_t bits = Get<uint64_t>();
    double v;
    std::memcpy(&v, &bits, 8);
    return v;
  }
  bool AtEnd() const { return n_ == 0; }

 private:
  template <typename T>
  T Get() {
    if (n_ < sizeof(T)) throw std::runtime_error("unexpected end of data");
    uint8_t b[sizeof(T)];
    std::memcpy(b, p_, sizeof(T));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    std::reverse(b, b + sizeof(T));
#endif
    T v;
    std::memcpy(&v, b, sizeof(T));
    p_ += sizeof(T);
    n_ -= sizeof(T);
    return v;
  }
  const uint8_t* p_;
  size_t n_;
};

inline Scene Deserialize(const std::vector<uint8_t>& bytes, uint16_t expected_version) {
  if (bytes.size() < 20) throw std::runtime_error("buffer too small for header");
  char magic[4];
  std::memcpy(magic, bytes.data(), 4);
  if (std::memcmp(magic, "P2TS", 4) != 0) throw std::runtime_error("bad magic (not a P2TS stream)");
  Reader r(bytes.data(), bytes.size());
  r.U32();  // skip magic
  const uint16_t version = r.U16();
  const uint16_t flags = r.U16();
  const uint32_t stored_size = r.U32();
  const uint32_t raw_size = r.U32();
  const uint32_t crc = r.U32();

  if (version != expected_version)
    throw std::runtime_error("version mismatch: stream is v" + std::to_string(version) +
                             ", reader supports v" + std::to_string(expected_version));
  if (bytes.size() - 20 != stored_size)
    throw std::runtime_error("payload size mismatch (truncated or corrupt stream)");

  const uint8_t* payload = bytes.data() + 20;
  std::vector<uint8_t> raw;
  if (flags & 1) {
    raw = RleDecode(payload, stored_size, raw_size);
  } else {
    raw.assign(payload, payload + stored_size);
    if (raw.size() != raw_size) throw std::runtime_error("payload size mismatch");
  }
  if (Crc32(raw.data(), raw.size()) != crc)
    throw std::runtime_error("CRC32 mismatch: stream is corrupt");

  Scene s;
  Reader d(raw.data(), raw.size());
  const uint32_t nv = d.U32();
  s.vertices.reserve(nv);
  for (uint32_t i = 0; i < nv; ++i) s.vertices.push_back({d.F64(), d.F64()});
  const uint32_t nt = d.U32();
  s.triangles.reserve(nt);
  for (uint32_t i = 0; i < nt; ++i) s.triangles.push_back({d.U32(), d.U32(), d.U32()});
  const uint32_t ne = d.U32();
  s.edges.reserve(ne);
  for (uint32_t i = 0; i < ne; ++i) s.edges.push_back({d.U32(), d.U32()});
  const uint32_t np = d.U32();
  s.polygons.reserve(np);
  for (uint32_t i = 0; i < np; ++i) {
    Polygon p;
    const uint32_t ol = d.U32();
    for (uint32_t j = 0; j < ol; ++j) p.outer.push_back(d.U32());
    const uint32_t nh = d.U32();
    for (uint32_t j = 0; j < nh; ++j) {
      const uint32_t hl = d.U32();
      std::vector<uint32_t> hole(hl);
      for (auto& idx : hole) idx = d.U32();
      p.holes.push_back(std::move(hole));
    }
    s.polygons.push_back(std::move(p));
  }
  if (!d.AtEnd()) throw std::runtime_error("trailing bytes after payload");
  return s;
}

// ---------------------------------------------------- triangulation helpers
// Triangulate every polygon of a scene; returns triangles as index triples in
// the exact order poly2tri produces them.
inline std::vector<std::array<uint32_t, 3>> Triangulate(const Scene& s) {
  std::vector<std::array<uint32_t, 3>> out;
  for (const auto& poly : s.polygons) {
    std::vector<p2t::Point*> points(s.vertices.size(), nullptr);
    std::unordered_map<const p2t::Point*, uint32_t> index_of;
    auto make_point = [&](uint32_t idx) -> p2t::Point* {
      if (!points[idx]) {
        points[idx] = new p2t::Point(s.vertices[idx][0], s.vertices[idx][1]);
        index_of[points[idx]] = idx;
      }
      return points[idx];
    };
    std::vector<p2t::Point*> outer;
    for (uint32_t idx : poly.outer) outer.push_back(make_point(idx));
    p2t::CDT cdt(outer);
    for (const auto& hole : poly.holes) {
      std::vector<p2t::Point*> h;
      for (uint32_t idx : hole) h.push_back(make_point(idx));
      cdt.AddHole(h);
    }
    cdt.Triangulate();
    for (p2t::Triangle* t : cdt.GetTriangles())
      out.push_back({index_of[t->GetPoint(0)], index_of[t->GetPoint(1)], index_of[t->GetPoint(2)]});
    for (p2t::Point* p : points) delete p;
  }
  return out;
}

inline std::vector<std::array<uint32_t, 2>> ConstrainedEdges(const Scene& s) {
  std::vector<std::array<uint32_t, 2>> edges;
  auto add_loop = [&](const std::vector<uint32_t>& loop) {
    for (size_t i = 0; i < loop.size(); ++i)
      edges.push_back({loop[i], loop[(i + 1) % loop.size()]});
  };
  for (const auto& p : s.polygons) {
    add_loop(p.outer);
    for (const auto& h : p.holes) add_loop(h);
  }
  return edges;
}

}  // namespace p2ts

// ------------------------------------------------------------ test fixtures
using p2ts::Polygon;
using p2ts::Scene;

static Scene MakeSingleContour() {
  Scene s;
  s.vertices = {{0, 0}, {4, 0}, {5, 2}, {3, 1.5}, {4, 4}, {1, 3}, {0, 4}};
  Polygon p;
  p.outer = {0, 1, 2, 3, 4, 5, 6};
  s.polygons.push_back(p);
  return s;
}

static Scene MakeWithHole() {
  Scene s;
  s.vertices = {{0, 0}, {10, 0}, {10, 8}, {0, 8},        // outer
                {3, 3}, {7, 3},  {7, 5},  {3, 5}};       // hole
  Polygon p;
  p.outer = {0, 1, 2, 3};
  p.holes = {{4, 5, 6, 7}};
  s.polygons.push_back(p);
  return s;
}

static Scene MakeMultiContour() {
  Scene s;
  s.vertices = {// polygon 0: triangle
                {0, 0}, {3, 0}, {1.5, 2.5},
                // polygon 1: pentagon
                {6, 0}, {9, 0.5}, {9.5, 3}, {7, 4}, {5.5, 2},
                // polygon 2: quad with triangular hole
                {12, 0}, {16, 0}, {16, 4}, {12, 4}, {13.5, 1.5}, {14.5, 1.5}, {14, 2.5}};
  Polygon p0; p0.outer = {0, 1, 2};
  Polygon p1; p1.outer = {3, 4, 5, 6, 7};
  Polygon p2; p2.outer = {8, 9, 10, 11}; p2.holes = {{12, 13, 14}};
  s.polygons = {p0, p1, p2};
  return s;
}

static Scene MakeDegenerate() {
  // Collinear consecutive vertices along every edge (180-degree corners)
  // plus a very thin sliver hole.
  Scene s;
  s.vertices = {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0},          // bottom edge
                {4, 1}, {4, 2},                                   // right edge
                {3, 2}, {2, 2}, {1, 2}, {0, 2},                   // top edge
                {0, 1},                                           // left edge
                {1.5, 0.9}, {2.5, 0.9}, {2.5, 0.9 + 1e-4}, {1.5, 0.9 + 1e-4}};  // sliver hole
  Polygon p;
  p.outer = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  p.holes = {{12, 13, 14, 15}};
  s.polygons.push_back(p);
  return s;
}

// ------------------------------------------------------------------ driver
struct CaseResult {
  std::string name;
  size_t raw_bytes, rle_bytes;
  bool roundtrip_ok, retriangulate_ok;
  double write_us, read_us;
};

static CaseResult RunCase(const std::string& name, Scene scene) {
  using clock = std::chrono::steady_clock;
  CaseResult r;
  r.name = name;

  // 1. Original triangulation (ground truth, vertex order as produced).
  scene.triangles = p2ts::Triangulate(scene);
  scene.edges = p2ts::ConstrainedEdges(scene);
  const auto original = scene.triangles;

  // 2. Serialize (timed), both compression settings.
  auto t0 = clock::now();
  const auto bytes_raw = p2ts::Serialize(scene, p2ts::kCurrentVersion, /*compress=*/false);
  auto t1 = clock::now();
  r.write_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
  const auto bytes_rle = p2ts::Serialize(scene, p2ts::kCurrentVersion, /*compress=*/true);
  r.raw_bytes = bytes_raw.size();
  r.rle_bytes = bytes_rle.size();

  // 3. Deserialize (timed, compressed variant) and verify.
  t0 = clock::now();
  Scene back = p2ts::Deserialize(bytes_rle, p2ts::kCurrentVersion);
  t1 = clock::now();
  r.read_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

  // 4. Stored result must equal the original triangulation vertex-for-vertex.
  r.roundtrip_ok = (back.vertices == scene.vertices) && (back.edges == scene.edges) &&
                   (back.triangles == original);

  // 5. Re-triangulating the deserialized input must reproduce the stored result.
  const auto redone = p2ts::Triangulate(back);
  r.retriangulate_ok = (redone == back.triangles);
  return r;
}

int main() {
  std::vector<CaseResult> results;
  results.push_back(RunCase("single-contour", MakeSingleContour()));
  results.push_back(RunCase("with-hole", MakeWithHole()));
  results.push_back(RunCase("multi-contour", MakeMultiContour()));
  results.push_back(RunCase("degenerate", MakeDegenerate()));

  std::cout << "case            raw(B)  rle(B)  roundtrip  re-triangulate  write(us)  read(us)\n";
  bool all_ok = true;
  for (const auto& r : results) {
    char line[160];
    std::snprintf(line, sizeof(line), "%-15s %6zu  %6zu  %-9s  %-14s  %9.1f  %8.1f", r.name.c_str(),
                  r.raw_bytes, r.rle_bytes, r.roundtrip_ok ? "OK" : "FAIL",
                  r.retriangulate_ok ? "OK" : "FAIL", r.write_us, r.read_us);
    std::cout << line << '\n';
    all_ok &= r.roundtrip_ok && r.retriangulate_ok;
  }

  // Version mismatch must be reported as an error.
  const auto future = p2ts::Serialize(MakeSingleContour(), /*version=*/2, false);
  try {
    p2ts::Deserialize(future, /*expected_version=*/1);
    std::cout << "version-mismatch: FAIL (no error raised)\n";
    all_ok = false;
  } catch (const std::exception& e) {
    std::cout << "version-mismatch: OK, error = \"" << e.what() << "\"\n";
  }

  std::cout << (all_ok ? "ALL OK" : "FAILURES PRESENT") << '\n';
  return all_ok ? 0 : 1;
}
