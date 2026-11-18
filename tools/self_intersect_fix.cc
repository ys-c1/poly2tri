// 手绘轮廓自交检测与修复 + 三角化验证工具(独立脚本,不改动库代码与构建脚本)。
//
// 三段流程:
//   阶段1  扫描线剪枝的边对判交(sweep-line prune + 精确线段判交)
//   阶段2  求出全部交点(含共线重叠段端点)
//   阶段3  在交点处切分所有边,按位合并节点,重建平面图拓扑,
//          有向半边面遍历提取简单环,再交给 poly2tri CDT 三角化
//
// 可配参数(Options): 容差 tolerance;是否允许自动桥接重叠边 allow_auto_bridge。
// 重叠边不允许桥接时判定为无法修复,输出冲突边对。
//
// 编译运行见 tools/run.sh。

#include "poly2tri/poly2tri.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Vec {
  double x, y;
};
Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y}; }
Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.y + b.y}; }
Vec operator*(Vec a, double s) { return {a.x * s, a.y * s}; }
double Cross(Vec a, Vec b) { return a.x * b.y - a.y * b.x; }
double Dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y; }
double Len(Vec a) { return std::hypot(a.x, a.y); }

// ---------------- 可配参数 ----------------
struct Options {
  double tolerance = 1e-9;       // 顶点合并 / 交点吸附容差(坐标单位)
  bool allow_auto_bridge = true; // 重叠边自动桥接(去重);false 时重叠即报冲突
};

struct Intersection {
  Vec p;
  int edge_a, edge_b;
  bool overlap; // true: 共线重叠段的端点
};

struct FixResult {
  bool ok = true;
  std::vector<Intersection> intersections;
  std::vector<std::pair<int, int>> conflicts; // 冲突边对(原轮廓边下标)
  std::vector<std::vector<Vec>> loops;        // 重建出的简单环
};

struct Edge {
  Vec a, b;
};

using Clock = std::chrono::steady_clock;
struct Timings {
  double sweep_us = 0;   // 阶段1+2: 扫描判交 + 交点计算
  double rebuild_us = 0; // 阶段3: 切分 + 拓扑重建
  double tri_us = 0;     // CDT 三角化
};
double Micros(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::micro>(b - a).count();
}

double Shoelace(const std::vector<Vec>& p) {
  double s = 0;
  for (size_t i = 0; i < p.size(); ++i) {
    const Vec& a = p[i];
    const Vec& b = p[(i + 1) % p.size()];
    s += a.x * b.y - b.x * a.y;
  }
  return 0.5 * s;
}

// ---------------- 阶段1+2: 边对精确判交 ----------------
// adjacent: 两边共享轮廓端点(下标相邻),此时只判共线重叠,不判相交
void TestPair(int i, int j, const Edge& A, const Edge& B, bool adjacent, double tol,
              std::vector<std::vector<std::pair<double, Vec>>>& splits,
              std::vector<Intersection>& inters,
              std::set<std::pair<int, int>>& overlap_pairs) {
  const Vec r = A.b - A.a, s = B.b - B.a;
  const double lr = Len(r), ls = Len(s);
  if (lr < tol || ls < tol) return; // 零长边
  const Vec qp = B.a - A.a;
  const double denom = Cross(r, s);
  const double tol_dir = tol * std::max(lr, ls); // cross 量纲为 length^2
  if (std::fabs(denom) > tol_dir) {
    if (adjacent) return;
    const double t = Cross(qp, s) / denom; // A 上参数
    const double u = Cross(qp, r) / denom; // B 上参数
    const double tt = tol / lr, tu = tol / ls;
    if (t >= -tt && t <= 1.0 + tt && u >= -tu && u <= 1.0 + tu) {
      const double tc = std::min(1.0, std::max(0.0, t));
      const double uc = std::min(1.0, std::max(0.0, u));
      const Vec p = A.a + r * tc; // 同一个 p 对象推入两侧,保证节点按位一致
      inters.push_back({p, i, j, false});
      splits[i].push_back({tc, p});
      splits[j].push_back({uc, p});
    }
    return;
  }
  // 平行:检查是否共线(B 的端点到直线 A 的距离)
  if (std::fabs(Cross(qp, r)) > tol_dir) return;
  const double rr = Dot(r, r);
  const double t0 = Dot(B.a - A.a, r) / rr;
  const double t1 = Dot(B.b - A.a, r) / rr;
  const double lo = std::max(0.0, std::min(t0, t1));
  const double hi = std::min(1.0, std::max(t0, t1));
  if (hi - lo <= tol / lr) return; // 仅端点相接,不算重叠
  const double ss = Dot(s, s);
  overlap_pairs.insert({std::min(i, j), std::max(i, j)});
  for (double ta : {lo, hi}) {
    const Vec p = A.a + r * ta;
    splits[i].push_back({ta, p});
    splits[j].push_back({Dot(p - B.a, s) / ss, p});
    inters.push_back({p, i, j, true});
  }
}

// ---------------- 修复主流程 ----------------
FixResult FixContour(const std::vector<Vec>& contour, const Options& opt, Timings* tim) {
  FixResult res;
  const int n = static_cast<int>(contour.size());
  std::vector<Edge> edges(n);
  for (int i = 0; i < n; ++i) edges[i] = {contour[i], contour[(i + 1) % n]};

  const auto t0 = Clock::now();
  // 阶段1: 扫描线按 min-x 排序剪枝;阶段2: 对候选边对精确求交
  std::vector<std::vector<std::pair<double, Vec>>> splits(n);
  for (int i = 0; i < n; ++i) {
    splits[i].push_back({0.0, edges[i].a});
    splits[i].push_back({1.0, edges[i].b});
  }
  std::set<std::pair<int, int>> overlap_pairs;
  {
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return std::min(edges[i].a.x, edges[i].b.x) < std::min(edges[j].a.x, edges[j].b.x);
    });
    std::vector<int> active;
    for (int oi = 0; oi < n; ++oi) {
      const int i = order[oi];
      const double minx = std::min(edges[i].a.x, edges[i].b.x) - opt.tolerance;
      std::vector<int> kept;
      kept.reserve(active.size());
      for (int j : active)
        if (std::max(edges[j].a.x, edges[j].b.x) >= minx) kept.push_back(j);
      active.swap(kept);
      for (int j : active) {
        const bool adj = (j == (i + 1) % n) || (i == (j + 1) % n);
        TestPair(i, j, edges[i], edges[j], adj, opt.tolerance, splits, res.intersections,
                 overlap_pairs);
      }
      active.push_back(i);
    }
  }
  const auto t1 = Clock::now();
  tim->sweep_us = Micros(t0, t1);

  // 无交无重叠:原样透传,保证三角化输出逐顶点不变
  if (res.intersections.empty()) {
    res.loops.push_back(contour);
    return res;
  }

  // 阶段3: 在交点处切分,重建拓扑
  std::map<std::pair<double, double>, int> node_id; // 交点按位一致,精确合并
  std::vector<Vec> nodes;
  auto getNode = [&](Vec p) {
    const std::pair<double, double> key{p.x, p.y};
    const auto it = node_id.find(key);
    if (it != node_id.end()) return it->second;
    const int id = static_cast<int>(nodes.size());
    nodes.push_back(p);
    node_id[key] = id;
    return id;
  };
  std::map<std::pair<int, int>, std::vector<int>> seg_owners; // 无向子段 -> 原始边号
  for (int i = 0; i < n; ++i) {
    auto& sp = splits[i];
    std::sort(sp.begin(), sp.end(),
              [](const std::pair<double, Vec>& a, const std::pair<double, Vec>& b) {
                return a.first < b.first;
              });
    const double len = Len(edges[i].b - edges[i].a);
    const double tt = opt.tolerance / std::max(len, opt.tolerance);
    std::vector<std::pair<double, Vec>> uniq;
    for (const auto& s : sp) {
      if (!uniq.empty() && std::fabs(s.first - uniq.back().first) <= tt) continue;
      uniq.push_back(s);
    }
    for (size_t k = 0; k + 1 < uniq.size(); ++k) {
      const int u = getNode(uniq[k].second);
      const int v = getNode(uniq[k + 1].second);
      if (u == v) continue;
      seg_owners[{std::min(u, v), std::max(u, v)}].push_back(i);
    }
  }
  // 被多条原始边覆盖的子段 = 重叠冲突
  std::set<std::pair<int, int>> conflicts;
  for (const auto& kv : seg_owners) {
    const auto& o = kv.second;
    for (size_t a = 0; a < o.size(); ++a)
      for (size_t b = a + 1; b < o.size(); ++b)
        conflicts.insert({std::min(o[a], o[b]), std::max(o[a], o[b])});
  }
  res.conflicts.assign(conflicts.begin(), conflicts.end());
  if (!res.conflicts.empty() && !opt.allow_auto_bridge) {
    res.ok = false; // 无法自动修复:报告冲突边
    return res;
  }
  // 桥接:重复子段只保留一份,重建平面图邻接关系
  std::vector<std::vector<int>> adj(nodes.size());
  for (const auto& kv : seg_owners) {
    adj[kv.first.first].push_back(kv.first.second);
    adj[kv.first.second].push_back(kv.first.first);
  }
  // 有向半边面遍历:在每个节点沿"自base顺时针最小角"的边离开,左侧即所求面
  auto angleOf = [&](int from, int to) {
    const Vec d = nodes[to] - nodes[from];
    return std::atan2(d.y, d.x);
  };
  std::set<std::pair<int, int>> used;
  const double area_tol = opt.tolerance;
  for (const auto& kv : seg_owners) {
    for (int dir = 0; dir < 2; ++dir) {
      const int su = dir ? kv.first.second : kv.first.first;
      const int sv = dir ? kv.first.first : kv.first.second;
      if (used.count({su, sv})) continue;
      std::vector<int> poly;
      int cu = su, cv = sv;
      const int guard_max = 4 * static_cast<int>(seg_owners.size()) + 8;
      for (int guard = 0; guard < guard_max; ++guard) {
        used.insert({cu, cv});
        poly.push_back(cu);
        const double base = angleOf(cv, cu);
        int best = -1;
        double bestd = 0;
        for (int w : adj[cv]) {
          double d = base - angleOf(cv, w);
          while (d <= 0) d += 2 * kPi;
          while (d > 2 * kPi) d -= 2 * kPi;
          if (best < 0 || d < bestd) {
            best = w;
            bestd = d;
          }
        }
        const int nu = cv, nv = best;
        if (nu == su && nv == sv) break; // 闭合
        cu = nu;
        cv = nv;
      }
      if (poly.size() < 3) continue;
      std::vector<Vec> loop;
      loop.reserve(poly.size());
      for (int id : poly) loop.push_back(nodes[id]);
      if (Shoelace(loop) > area_tol) res.loops.push_back(loop); // 只保留正向有界面
    }
  }
  tim->rebuild_us = Micros(t1, Clock::now());
  return res;
}

// ---------------- 三角化(含嵌套环:内环作洞) ----------------
bool PointInLoop(const std::vector<Vec>& poly, Vec p) {
  bool inside = false;
  for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
    const Vec& a = poly[j];
    const Vec& b = poly[i];
    if ((a.y > p.y) != (b.y > p.y)) {
      const double xint = a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y);
      if (p.x < xint) inside = !inside;
    }
  }
  return inside;
}

struct TriOut {
  long triangles = 0;
  double area = 0;
  double expected = 0; // 按嵌套深度交错求和 = 三角化应得面积
  bool delaunay = true;
};

TriOut TriangulateLoops(const std::vector<std::vector<Vec>>& loops) {
  TriOut out;
  const int m = static_cast<int>(loops.size());
  std::vector<double> areas(m);
  std::vector<Vec> cent(m);
  for (int i = 0; i < m; ++i) {
    areas[i] = std::fabs(Shoelace(loops[i]));
    Vec c{0, 0};
    for (const Vec& v : loops[i]) c = c + v;
    cent[i] = c * (1.0 / loops[i].size());
  }
  std::vector<int> parent(m, -1); // 最内层包含者
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      if (i != j && areas[j] > areas[i] && PointInLoop(loops[j], cent[i]))
        if (parent[i] < 0 || areas[j] < areas[parent[i]]) parent[i] = j;
  std::vector<int> depth(m, 0);
  for (int i = 0; i < m; ++i) {
    int d = 0, p = parent[i];
    while (p >= 0) {
      ++d;
      p = parent[p];
    }
    depth[i] = d;
  }
  for (int i = 0; i < m; ++i) out.expected += (depth[i] % 2 == 0 ? 1.0 : -1.0) * areas[i];
  for (int i = 0; i < m; ++i) {
    if (depth[i] % 2 != 0) continue; // 奇数深度是洞,随父环一起三角化
    std::vector<p2t::Point*> poly;
    for (const Vec& v : loops[i]) poly.push_back(new p2t::Point(v.x, v.y));
    p2t::CDT cdt(poly);
    std::vector<std::vector<p2t::Point*>> holes;
    for (int j = 0; j < m; ++j)
      if (parent[j] == i) {
        holes.emplace_back();
        for (const Vec& v : loops[j]) holes.back().push_back(new p2t::Point(v.x, v.y));
        cdt.AddHole(holes.back());
      }
    cdt.Triangulate();
    const auto tris = cdt.GetTriangles();
    out.triangles += static_cast<long>(tris.size());
    for (p2t::Triangle* t : tris) {
      const p2t::Point *p0 = t->GetPoint(0), *p1 = t->GetPoint(1), *p2 = t->GetPoint(2);
      const double x1 = p1->x - p0->x, y1 = p1->y - p0->y;
      const double x2 = p2->x - p0->x, y2 = p2->y - p0->y;
      out.area += 0.5 * std::fabs(x1 * y2 - x2 * y1);
    }
    out.delaunay = out.delaunay && p2t::IsDelaunay(tris);
    for (p2t::Point* p : poly) delete p;
    for (auto& h : holes)
      for (p2t::Point* p : h) delete p;
  }
  return out;
}

// ---------------- 报告 ----------------
void PrintEdge(const std::vector<Vec>& c, int i) {
  const Vec a = c[i], b = c[(i + 1) % c.size()];
  std::printf("边%d(%.10g,%.10g)->(%.10g,%.10g)", i, a.x, a.y, b.x, b.y);
}

bool RunCase(const char* name, const std::vector<Vec>& contour, const Options& opt) {
  Timings tim;
  FixResult res = FixContour(contour, opt, &tim);
  std::printf("=== %s ===\n", name);
  std::printf("  配置: 容差=%.3g  自动桥接=%s\n", opt.tolerance,
              opt.allow_auto_bridge ? "开" : "关");
  if (res.intersections.empty()) {
    std::printf("  交点: 无\n");
  } else {
    for (const auto& ix : res.intersections)
      std::printf("  交点 (%.10g, %.10g)  边%d x 边%d%s\n", ix.p.x, ix.p.y, ix.edge_a,
                  ix.edge_b, ix.overlap ? "  [共线重叠段端点]" : "");
  }
  if (!res.ok) {
    std::printf("  无法自动修复,冲突边:\n");
    for (const auto& c : res.conflicts) {
      std::printf("    ");
      PrintEdge(contour, c.first);
      std::printf("  x  ");
      PrintEdge(contour, c.second);
      std::printf("\n");
    }
    return false;
  }
  if (!res.conflicts.empty()) {
    std::printf("  重叠冲突边(已自动桥接去重): ");
    for (const auto& c : res.conflicts) std::printf("边%dx边%d ", c.first, c.second);
    std::printf("\n");
  }
  const auto t2 = Clock::now();
  const TriOut tri = TriangulateLoops(res.loops);
  tim.tri_us = Micros(t2, Clock::now());
  const double dev = std::fabs(tri.expected - tri.area);
  const bool pass = dev < 1e-9 && tri.delaunay;
  std::printf("  重建环数: %zu  修复后三角面数: %ld\n", res.loops.size(), tri.triangles);
  std::printf("  面积: 原始鞋带=%.10g  重建环=%.10g  三角化=%.10g  偏差=%.3g  %s\n",
              Shoelace(contour), tri.expected, tri.area, dev,
              dev < 1e-9 ? "一致[OK]" : "不一致[FAIL]");
  std::printf("  Delaunay 校验: %s\n", tri.delaunay ? "通过" : "失败");
  std::printf("  耗时: 扫描判交 %.1fus  切分重建 %.1fus  三角化 %.1fus\n", tim.sweep_us,
              tim.rebuild_us, tim.tri_us);
  return pass;
}

// 逐顶点一致性:无交轮廓走修复管线 vs 直接 CDT,三角形顶点集合必须按位相同
bool IdentityCheck(const std::vector<Vec>& contour) {
  auto snapshot = [](const std::vector<p2t::Triangle*>& tris) {
    std::vector<std::pair<double, double>> v;
    for (p2t::Triangle* t : tris)
      for (int k = 0; k < 3; ++k) v.push_back({t->GetPoint(k)->x, t->GetPoint(k)->y});
    std::sort(v.begin(), v.end());
    return v;
  };
  std::vector<p2t::Point*> pts;
  for (const Vec& v : contour) pts.push_back(new p2t::Point(v.x, v.y));
  p2t::CDT cdt(pts);
  cdt.Triangulate();
  const auto direct = snapshot(cdt.GetTriangles());
  for (p2t::Point* p : pts) delete p;

  Timings tim;
  FixResult res = FixContour(contour, Options{}, &tim);
  bool same_loop = res.loops.size() == 1 && res.loops[0].size() == contour.size();
  if (same_loop)
    for (size_t k = 0; k < contour.size(); ++k)
      same_loop = same_loop && res.loops[0][k].x == contour[k].x &&
                  res.loops[0][k].y == contour[k].y;
  std::vector<p2t::Point*> pts2;
  for (const Vec& v : res.loops[0]) pts2.push_back(new p2t::Point(v.x, v.y));
  p2t::CDT cdt2(pts2);
  cdt2.Triangulate();
  const auto piped = snapshot(cdt2.GetTriangles());
  for (p2t::Point* p : pts2) delete p;

  const bool identical = same_loop && direct == piped;
  std::printf("=== 逐顶点一致性(无交轮廓: 修复管线 vs 直接 CDT) ===\n");
  std::printf("  环透传: %s  三角形顶点集合按位相同: %s\n", same_loop ? "是" : "否",
              direct == piped ? "是" : "否");
  return identical;
}

} // namespace

int main() {
  bool all = true;

  // 案例1: 无交(凸凹混合简单多边形)
  all &= RunCase("案例1: 无交轮廓",
                 {{0, 0}, {4, 0}, {5, 2}, {2, 4}, {0, 3}}, Options{});

  // 案例2: 单交(蝴蝶结,边0 x 边2 交于 (1,1))
  all &= RunCase("案例2: 单交轮廓(蝴蝶结)",
                 {{0, 0}, {2, 2}, {2, 0}, {0, 2}}, Options{});

  // 案例3: 多交(五角星,5 个交点)
  std::vector<Vec> star;
  for (int k = 0; k < 5; ++k) {
    const double a = (90.0 + 144.0 * k) * kPi / 180.0;
    star.push_back({std::cos(a), std::sin(a)});
  }
  all &= RunCase("案例3: 多交轮廓(五角星)", star, Options{});

  // 案例4: 重叠边(顶边被来回描了三遍),允许自动桥接
  const std::vector<Vec> overlap{{0, 0}, {4, 0}, {4, 2}, {1, 2}, {3, 2}, {0, 2}};
  all &= RunCase("案例4: 重叠边轮廓(自动桥接)", overlap, Options{});

  // 案例5: 同一重叠轮廓,禁止自动桥接 -> 必须报冲突边
  Options nobridge;
  nobridge.allow_auto_bridge = false;
  all &= !RunCase("案例5: 重叠边轮廓(禁止桥接,应报冲突)", overlap, nobridge);

  // 逐顶点一致性
  all &= IdentityCheck({{0, 0}, {4, 0}, {5, 2}, {2, 4}, {0, 3}});

  std::printf("%s\n", all ? "全部通过" : "存在失败项");
  return all ? 0 : 1;
}
