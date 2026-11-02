#ifndef PANGENOME_INDEX_LIFT_FOLD_HPP
#define PANGENOME_INDEX_LIFT_FOLD_HPP

// Folds source->target correspondences into target intervals or chain blocks.
//
// This is the per-base fold the middleware used to run in Python
// (_fold_to_intervals / _fold_to_blocks), with the same output. It takes
// node-sized runs instead of single bases, so a run costs O(1) rather than one
// step per base; runs that overlap in source fall back to the per-base walk.

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace liftfold {

/// n correspondences: source s0 + i maps to target t0 + d * i (d is +1 or -1).
struct Run {
    int64_t s0 = 0, t0 = 0, n = 0;
    int d = 1;
};

/// One folded piece; both sides half-open.
struct Piece {
    int64_t source_start = 0, source_end = 0;
    int64_t target_start = 0, target_end = 0;
    char strand = '+';
};

/// INTERVALS: maximal runs of one target direction (_fold_to_intervals).
/// BLOCKS: maximal runs of constant source-target offset (_fold_to_blocks).
enum class Kind { INTERVALS, BLOCKS };

/// The per-base fold as a state machine over consecutive points.
class Folder {
public:
    Folder(Kind kind, std::vector<Piece>& out) : kind_(kind), out_(out) {}

    /// Add the next point in (source, target) order.
    void add(int64_t s, int64_t t) {
        if (!open_) { start(s, t); return; }
        const int64_t ds = s - s1_, dt = t - t1_;
        int step = 0;
        if (kind_ == Kind::INTERVALS) {
            if (dt == 0) { extend(s, t); return; }       // plateau
            step = dt > 0 ? 1 : -1;
        } else {
            if (ds <= 0) { restart(s, t); return; }
            if (dt == ds) step = 1;
            else if (dt == -ds) step = -1;
            else { restart(s, t); return; }               // indel
        }
        if (dir_ == 0) dir_ = step;
        if (step != dir_) { restart(s, t); return; }
        extend(s, t);
    }

    /// Add a run's points. By its third point the direction is the run's own,
    /// so every later step is accepted.
    void add(const Run& r) {
        for (int64_t i = 0; i < std::min<int64_t>(r.n, 3); i++) add(r.s0 + i, r.t0 + r.d * i);
        if (r.n >= 4) extend(r.s0 + r.n - 1, r.t0 + r.d * (r.n - 1));
    }

    void finish() {
        if (open_) emit();
        open_ = false;
    }

private:
    Kind kind_;
    std::vector<Piece>& out_;
    bool open_ = false;
    int dir_ = 0;
    int64_t s0_ = 0, t0_ = 0, s1_ = 0, t1_ = 0, lo_ = 0, hi_ = 0;

    void start(int64_t s, int64_t t) {
        open_ = true; dir_ = 0;
        s0_ = s1_ = s; t0_ = t1_ = lo_ = hi_ = t;
    }
    void extend(int64_t s, int64_t t) {
        s1_ = s; t1_ = t;
        lo_ = std::min(lo_, t); hi_ = std::max(hi_, t);
    }
    void restart(int64_t s, int64_t t) { emit(); start(s, t); }

    void emit() {
        Piece p;
        p.source_start = s0_;
        p.source_end = s1_ + 1;
        p.strand = dir_ < 0 ? '-' : '+';
        if (kind_ == Kind::INTERVALS) {
            p.target_start = lo_; p.target_end = hi_ + 1;
        } else {
            p.target_start = dir_ >= 0 ? t0_ : t1_;
            p.target_end = (dir_ >= 0 ? t1_ : t0_) + 1;
        }
        out_.push_back(p);
    }
};

/// Fold one alignment's runs (one target path), in source order.
inline std::vector<Piece> fold(std::vector<Run> runs, Kind kind) {
    std::vector<Piece> out;
    std::sort(runs.begin(), runs.end(), [](const Run& a, const Run& b) {
        return a.s0 != b.s0 ? a.s0 < b.s0 : a.t0 < b.t0;
    });
    bool overlap = false;
    for (size_t i = 1; i < runs.size() && !overlap; i++) {
        overlap = runs[i].s0 < runs[i - 1].s0 + runs[i - 1].n;
    }
    Folder f(kind, out);
    if (!overlap) {
        for (const Run& r : runs) f.add(r);
    } else {
        // Two copies on one path interleave in source order: walk every base,
        // sorted by (source, target), exactly as the per-base fold did.
        std::vector<std::pair<int64_t, int64_t>> pts;
        for (const Run& r : runs) {
            for (int64_t i = 0; i < r.n; i++) pts.emplace_back(r.s0 + i, r.t0 + r.d * i);
        }
        std::sort(pts.begin(), pts.end());
        for (const auto& p : pts) f.add(p.first, p.second);
    }
    f.finish();
    return out;
}

} // namespace liftfold

#endif // PANGENOME_INDEX_LIFT_FOLD_HPP
