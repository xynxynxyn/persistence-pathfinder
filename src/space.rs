use std::{
    collections::{HashMap, HashSet, VecDeque},
    f32::INFINITY,
    fmt::Display,
    hash::Hash,
    ops::IndexMut,
};

use anyhow::{Result, anyhow};
use itertools::Itertools;

use crate::grid::{Coord, Grid};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Ord, Eq, Hash)]
struct BiEdge {
    inner: (Coord, Coord),
}

impl BiEdge {
    // Use the constructor to create edges only so they are ordered internally and
    // ((a,b),(c,d)) is the same as ((c,d),(a,b))
    fn new(from: impl Into<Coord>, to: impl Into<Coord>) -> BiEdge {
        let from = from.into();
        let to = to.into();
        if from < to {
            BiEdge { inner: (from, to) }
        } else {
            BiEdge { inner: (to, from) }
        }
    }
}

impl Display for BiEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}<->{}", self.inner.0, self.inner.1)
    }
}

#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub struct CoordSig {
    c: Coord,
    s: EdgeSignature,
}

impl Display for CoordSig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/[{}]", self.c, self.s)
    }
}

// Represents the changes in H2-signature when traversing an edge
// When using cut-crossing approximation this is either going to be 0 or 1
// For a more general approach this should be f32 instead of i32
#[derive(Debug, Clone, PartialEq, PartialOrd, Ord, Eq, Hash)]
pub struct EdgeSignature {
    inner: Vec<i32>,
}

impl Display for EdgeSignature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{}]",
            self.inner.iter().map(|i| i.to_string()).join(",".into())
        )
    }
}

impl From<Vec<i32>> for EdgeSignature {
    fn from(inner: Vec<i32>) -> Self {
        EdgeSignature { inner }
    }
}

impl EdgeSignature {
    fn update(&self, other: &EdgeSignature) -> EdgeSignature {
        EdgeSignature {
            inner: self
                .inner
                .iter()
                .zip(other.inner.iter())
                .map(|(l, r)| (l + r) % 2)
                .collect_vec(),
        }
    }
}

#[derive(Debug)]
pub struct ThresholdedSpace {
    inner: Grid<bool>,
    obstacles: Vec<Coord>,
    edges: HashMap<BiEdge, EdgeSignature>,
    default_signature: EdgeSignature,
}

#[derive(Debug)]
pub struct ProbabilityMap {
    inner: Grid<f32>,
}

impl ProbabilityMap {
    pub fn empty(rows: usize, cols: usize) -> Self {
        Self {
            inner: Grid::new(rows, cols),
        }
    }

    pub fn update(&mut self, c: impl Into<Coord>, val: f32) -> Result<()> {
        let c = c.into();
        if !self.inner.in_bounds(c) {
            return Err(anyhow!("OOB"));
        }

        if 1f32 < val || val < 0f32 {
            return Err(anyhow!("invalid probability {}", val));
        }

        *self.inner.index_mut(c) = val;

        Ok(())
    }
}

impl ThresholdedSpace {
    fn sig(&self, edge: BiEdge) -> EdgeSignature {
        self.edges
            .get(&edge)
            .unwrap_or(&self.default_signature)
            .clone()
    }

    pub fn from_probability_map(prob_map: &ProbabilityMap, threshold: f32) -> Self {
        let b_map = prob_map.inner.map(|&p| p <= threshold);

        // Do a flood fill to find all possible obstacles
        let obstacles = b_map.filter(|b| !b).collect::<Vec<_>>();
        let mut obstacle_refs = vec![];
        let mut visited = HashSet::new();
        for o in obstacles {
            if visited.contains(&o) {
                continue;
            }

            // Found a new obstacle, save the initial point as a reference
            visited.insert(o);

            // do a flood fill and mark everything visited
            let mut queue = VecDeque::from([o]);
            // Take the smallest coord as the reference to be consistent
            let mut ref_coord = o;
            while let Some(cur) = queue.pop_front() {
                // Only consider neighbors that are also obstacles
                for n in b_map.filter_neighbors(cur, |b| !b) {
                    // skip if we already considered this cell
                    if visited.contains(&n) {
                        continue;
                    }
                    // Update the reference coordinate if it's better
                    if n < ref_coord {
                        ref_coord = n;
                    }
                    visited.insert(n);
                    queue.push_back(n);
                }
            }

            obstacle_refs.push(ref_coord);
        }

        // Compute the relevant edge signatures
        // If an edge has a 0-signature with no change it is not inserted
        let mut edges = HashMap::new();

        // Go over all coordinates that are not obstacles
        for from in b_map.filter(|&b| b) {
            for to in b_map.filter_neighbors(from, |&b| b) {
                let edge = BiEdge::new(from, to);
                // check if the edge crosses a boundary by doing the cross-cut approximation
                // for every obstacle ref
                let mut sig = vec![0; obstacle_refs.len()];
                let mut nonzero = false;
                for (i, o) in obstacle_refs.iter().enumerate() {
                    // Check if the coordinates are both to the left of the ref or in the same
                    // column
                    if o.diff(&to).1 > 0 || o.diff(&from).1 > 0 {
                        continue;
                    }
                    // Check if one is above and one below
                    if (o.row <= from.row) == (o.row <= to.row) {
                        continue;
                    }
                    // Now we know the edge crosses this boundary, add 1 to the sig
                    sig[i] = 1;
                    nonzero = true;
                }

                if nonzero {
                    let sig: EdgeSignature = sig.into();
                    edges.insert(edge, sig);
                }
            }
        }

        let default_signature = vec![0; obstacle_refs.len()].into();

        ThresholdedSpace {
            inner: b_map,
            obstacles: obstacle_refs,
            edges,
            default_signature,
        }
    }

    // Find an optimal path from start to goal with the given signature
    // The heuristic function is a simple euclidean distance from the current vertex to the goal
    // This is probably super inefficient, but later on you can reuse the paths computed at higher
    // layers for heuristics.
    pub fn a_star(
        &self,
        start: impl Into<Coord>,
        goal: impl Into<Coord>,
        sig: impl Into<EdgeSignature>,
    ) -> Option<Vec<Coord>> {
        let start = start.into();
        let goal = goal.into();
        let sig = sig.into();

        let h = |c: &CoordSig| {
            let (row_d, col_d) = goal.diff(&c.c);
            ((row_d * row_d + col_d * col_d) as f32).sqrt()
        };
        let mut came_from: HashMap<CoordSig, CoordSig> = HashMap::new();
        let mut g_score: HashMap<CoordSig, i32> = HashMap::new();
        let mut f_score: HashMap<CoordSig, f32> = HashMap::new();

        // The start always has a 0 signature
        let start = CoordSig {
            c: start,
            s: self.default_signature.clone(),
        };
        let goal = CoordSig { c: goal, s: sig };

        g_score.insert(start.clone(), 0);
        let start_f_score = g_score[&start] as f32 + h(&start);
        f_score.insert(start.clone(), start_f_score);

        let mut open: HashSet<CoordSig> = HashSet::new();
        open.insert(start);
        let mut current;

        while !open.is_empty() {
            current = open
                .iter()
                .min_by(|c1, c2| {
                    f_score
                        .get(c1)
                        .unwrap_or(&INFINITY)
                        .partial_cmp(f_score.get(c2).unwrap_or(&INFINITY))
                        .expect("floating point comparison failed")
                })
                .cloned()
                .expect("no min element");
            open.remove(&current);

            if current == goal {
                // Done, reconstruct the path
                let mut path = vec![goal.c];
                let mut cur = current;
                while let Some(parent) = came_from.get(&cur) {
                    //println!("cost @ {}: {}", cur, g_score[&cur]);
                    path.push(parent.c);
                    cur = parent.clone();
                }
                path.reverse();
                return Some(path);
            }

            let cur_loc = current.c;
            for neighbor in self.inner.filter_neighbors(cur_loc, |&b| b) {
                let edge_delta = BiEdge::new(cur_loc, neighbor);
                let new_sig = current.s.update(&self.sig(edge_delta));
                let neighbor = CoordSig {
                    c: neighbor,
                    s: new_sig,
                };

                // The cost to go to any adjacent node is always 1 in the grid
                let tentative_g_score = g_score[&current] + 1;
                if tentative_g_score < g_score.get(&neighbor).copied().unwrap_or(i32::max_value()) {
                    came_from.insert(neighbor.clone(), current.clone());
                    g_score.insert(neighbor.clone(), tentative_g_score);
                    f_score.insert(neighbor.clone(), tentative_g_score as f32 + h(&neighbor));

                    // Add the neighbor to the set of vertices to check
                    open.insert(neighbor);
                }
            }
        }

        None
    }

    pub fn print_with_path(&self, path: &[Coord]) {
        for (r, row) in self.inner.rows_iter().enumerate() {
            for (c, cell) in row.enumerate() {
                let coord = Coord::new(r, c);
                if let Some(&c) = path.first()
                    && coord == c
                {
                    print!("s");
                } else if let Some(&c) = path.last()
                    && coord == c
                {
                    print!("g");
                } else if path.contains(&coord) {
                    print!("*");
                } else if cell {
                    print!(".");
                } else {
                    print!("#");
                }
            }
            println!("");
        }
    }
}

impl Display for ThresholdedSpace {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.inner.rows_iter() {
            for cell in row {
                if cell {
                    write!(f, ".")?;
                } else {
                    write!(f, "#")?;
                }
            }
            write!(f, "\n")?;
        }

        Ok(())
    }
}

mod tests {
    #[test]
    fn one_obstacle() {
        use super::*;
        let mut pmap = ProbabilityMap::empty(10, 10);
        pmap.update((4, 4), 1f32).unwrap();
        pmap.update((5, 5), 1f32).unwrap();
        pmap.update((4, 5), 1f32).unwrap();
        pmap.update((5, 4), 1f32).unwrap();
        pmap.update((3, 4), 1f32).unwrap();

        let tspace = ThresholdedSpace::from_probability_map(&pmap, 0.5);
        println!("{:?}", tspace);
        assert!(tspace.obstacles.len() == 1);
        assert!(tspace.obstacles[0] == (3, 4).into());
    }

    #[test]
    fn two_obstacles() {
        use super::*;
        let mut pmap = ProbabilityMap::empty(5, 10);
        pmap.update((1, 4), 1.0).unwrap();
        pmap.update((1, 3), 1.0).unwrap();
        pmap.update((3, 7), 1.0).unwrap();
        pmap.update((3, 6), 1.0).unwrap();

        let tspace = ThresholdedSpace::from_probability_map(&pmap, 0.5);
        println!("{}", tspace);
        println!("{:?}", tspace.obstacles);
        assert!(tspace.obstacles.len() == 2);
        assert!(tspace.obstacles.contains(&(1, 3).into()));
        assert!(tspace.obstacles.contains(&(3, 6).into()));
    }
}
