// Randomized check of liftfold::fold against a line-by-line port of the Python
// folds it replaces (_fold_to_intervals / _fold_to_blocks). Inputs where runs
// overlap in source are folded per chain instead, so those are checked for
// coverage and piece count rather than equality.
//   g++ -O2 -std=c++17 -I include tests/test_lift_fold.cpp -o /tmp/test_lift_fold
#include "pangenome_index/lift_fold.hpp"

#include <cstdio>
#include <random>
#include <tuple>

using liftfold::Kind;
using liftfold::Piece;
using liftfold::Run;
using Pts = std::vector<std::pair<int64_t, int64_t>>;

// _fold_to_intervals, one group.
static std::vector<Piece> py_intervals(Pts pts) {
    std::sort(pts.begin(), pts.end());
    std::vector<Piece> out;
    size_t i = 0, n = pts.size();
    while (i < n) {
        size_t j = i;
        int direction = 0;
        while (j + 1 < n) {
            int64_t dt = pts[j + 1].second - pts[j].second;
            if (dt == 0) { j++; continue; }
            int step = dt > 0 ? 1 : -1;
            if (direction == 0) direction = step;
            if (step != direction) break;
            j++;
        }
        int64_t lo = pts[i].second, hi = pts[i].second;
        for (size_t k = i; k <= j; k++) { lo = std::min(lo, pts[k].second); hi = std::max(hi, pts[k].second); }
        out.push_back({pts[i].first, pts[j].first + 1, lo, hi + 1, direction < 0 ? '-' : '+'});
        i = j + 1;
    }
    return out;
}

// _fold_to_blocks, one group.
static std::vector<Piece> py_blocks(Pts pts) {
    std::sort(pts.begin(), pts.end());
    std::vector<Piece> out;
    size_t i = 0, n = pts.size();
    while (i < n) {
        auto [s0, t0] = pts[i];
        size_t j = i;
        int direction = 0;
        while (j + 1 < n) {
            int64_t ds = pts[j + 1].first - pts[j].first;
            int64_t dt = pts[j + 1].second - pts[j].second;
            if (ds <= 0) break;
            int step;
            if (dt == ds) step = 1;
            else if (dt == -ds) step = -1;
            else break;
            if (direction == 0) direction = step;
            else if (step != direction) break;
            j++;
        }
        auto [s1, t1] = pts[j];
        int64_t lo = direction >= 0 ? t0 : t1, hi = direction >= 0 ? t1 : t0;
        out.push_back({s0, s1 + 1, lo, hi + 1, direction < 0 ? '-' : '+'});
        i = j + 1;
    }
    return out;
}

static bool same(const std::vector<Piece>& a, const std::vector<Piece>& b, bool with_source) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (std::tie(a[i].target_start, a[i].target_end, a[i].strand) !=
            std::tie(b[i].target_start, b[i].target_end, b[i].strand)) return false;
        if (with_source && std::tie(a[i].source_start, a[i].source_end) !=
                           std::tie(b[i].source_start, b[i].source_end)) return false;
    }
    return true;
}

// Every point lies in some piece (on its diagonal, for blocks) and every piece
// holds at least one point: the fold of overlapping runs loses and invents nothing.
static bool covers(const std::vector<Piece>& out, const Pts& pts, bool blocks) {
    std::vector<char> used(out.size(), 0);
    for (const auto& [s, t] : pts) {
        bool hit = false;
        for (size_t k = 0; k < out.size(); k++) {
            const Piece& p = out[k];
            if (s < p.source_start || s >= p.source_end || t < p.target_start || t >= p.target_end) continue;
            if (blocks) {
                int64_t want = p.strand == '-' ? p.target_end - 1 - (s - p.source_start)
                                               : p.target_start + (s - p.source_start);
                if (t != want) continue;
            }
            hit = true; used[k] = 1;
        }
        if (!hit) return false;
    }
    return std::all_of(used.begin(), used.end(), [](char u) { return u; });
}

// k colinear copies whose nodes interleave in source fold to exactly k pieces.
static size_t interleaved_failures(std::mt19937_64& rng) {
    auto R = [&](int64_t lo, int64_t hi) { return std::uniform_int_distribution<int64_t>(lo, hi)(rng); };
    size_t failures = 0;
    for (int it = 0; it < 20000; it++) {
        const int k = static_cast<int>(R(2, 6)), nodes = static_cast<int>(R(2, 200));
        std::vector<Run> runs;
        for (int c = 0; c < k; c++) {
            int d = R(0, 1) ? 1 : -1;
            int64_t s = 0, t = 100000 * (c + 1);
            for (int i = 0; i < nodes; i++) {
                int64_t n = R(1, 4);
                runs.push_back({s, t, n, d});
                s += n; t += d * n;
            }
        }
        std::shuffle(runs.begin(), runs.end(), rng);
        if (liftfold::fold(runs, Kind::INTERVALS).size() != static_cast<size_t>(k) ||
            liftfold::fold(runs, Kind::BLOCKS).size() != static_cast<size_t>(k)) {
            if (++failures <= 3) std::printf("FAIL interleaved case %d (k=%d)\n", it, k);
        }
    }
    return failures;
}

int main() {
    std::mt19937_64 rng(1);
    auto R = [&](int64_t lo, int64_t hi) { return std::uniform_int_distribution<int64_t>(lo, hi)(rng); };
    size_t cases = 0, overlapping = 0, failures = 0;
    for (int it = 0; it < 200000; it++) {
        // Node-like runs walking along the source, with SNP gaps, indels,
        // inversions, jumps, single-base nodes and occasional second copies.
        std::vector<Run> runs;
        int64_t s = R(0, 50), t = R(1000, 5000);
        int d = R(0, 4) == 0 ? -1 : 1;
        const int nruns = static_cast<int>(R(1, 12));
        for (int k = 0; k < nruns; k++) {
            Run r{s, t, R(1, 6), d};
            runs.push_back(r);
            int64_t gap = R(0, 3) == 0 ? R(1, 3) : 0;           // SNP: both skip
            int64_t indel = R(0, 5) == 0 ? R(-3, 3) : 0;        // offset change
            s += r.n + gap;
            t += d * (r.n + gap) + indel;
            if (R(0, 15) == 0) { d = -d; t += R(-20, 20); }    // inversion
            if (R(0, 30) == 0) t += R(-1000, 1000);             // jump
        }
        if (R(0, 10) == 0) {                                     // second copy on the path
            Run c = runs[R(0, static_cast<int64_t>(runs.size()) - 1)];
            c.t0 += R(-500, 500);
            if (R(0, 1)) c.d = -c.d;
            runs.push_back(c);
        }
        Pts pts;
        for (const Run& r : runs) for (int64_t i = 0; i < r.n; i++) pts.emplace_back(r.s0 + i, r.t0 + r.d * i);
        std::shuffle(runs.begin(), runs.end(), rng);

        auto fi = liftfold::fold(runs, Kind::INTERVALS);
        auto fb = liftfold::fold(runs, Kind::BLOCKS);
        cases++;
        std::vector<Run> sorted = runs;
        std::sort(sorted.begin(), sorted.end(), [](const Run& a, const Run& b) { return a.s0 < b.s0; });
        bool overlap = false;
        for (size_t i = 1; i < sorted.size(); i++) if (sorted[i].s0 < sorted[i-1].s0 + sorted[i-1].n) { overlap = true; break; }
        overlapping += overlap;
        bool ok = overlap ? covers(fi, pts, false) && covers(fb, pts, true) && fi.size() <= runs.size()
                          : same(fi, py_intervals(pts), true) && same(fb, py_blocks(pts), true);
        if (!ok && ++failures <= 5) {
            std::printf("FAIL case %d:", it);
            for (const Run& r : runs) std::printf(" (%ld,%ld,%ld,%d)", (long)r.s0, (long)r.t0, (long)r.n, r.d);
            std::printf("\n");
        }
    }
    std::printf("%zu cases (%zu with overlapping runs), %zu failures\n", cases, overlapping, failures);
    const size_t inter = interleaved_failures(rng);
    std::printf("20000 interleaved-copy cases, %zu failures\n", inter);
    failures += inter;
    return failures ? 1 : 0;
}
