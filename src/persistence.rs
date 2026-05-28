use bon::Builder;
use itertools::{Itertools, multizip};
use std::{
    collections::{HashMap, HashSet},
    f32::{INFINITY, NEG_INFINITY},
    fmt::Display,
    marker::PhantomData,
};

use crate::{
    grid::Coord,
    space::{EdgeSignature, ProbabilityMap, ThresholdedSpace},
};

type Interval = (f32, f32);
type Path = Vec<Coord>;

struct PathClass {
    paths: Vec<(f32, Path)>,
}

// Gives a single path per equivalenceclass and its valid epsilon interval
#[derive(Default)]
pub struct Barcode {
    path_classes: Vec<(PathClass, Interval)>,
    epsilons: Vec<f32>,
}

impl Barcode {
    fn new(epsilons: Vec<f32>) -> Barcode {
        Barcode {
            path_classes: vec![],
            epsilons,
        }
    }

    fn close_path_class(&mut self, paths: Vec<(f32, Path)>) {
        let interval = (
            paths.first().expect("empty path class").0,
            paths.last().expect("empty path class").0,
        );
        self.path_classes.push((PathClass { paths }, interval));
    }
}

impl Display for Barcode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, (paths, interval)) in self.path_classes.iter().enumerate() {
            writeln!(f, "path class {} [{:.1},{:.1}]", i, interval.0, interval.1)?;
            for (e, p) in &paths.paths {
                writeln!(
                    f,
                    "@{:.1}: {}",
                    e,
                    p.iter().map(Coord::to_string).join(",".into())
                )?;
            }
        }

        Ok(())
    }
}

#[derive(Builder)]
pub struct Solver {
    #[builder(default = 0.1)]
    delta: f32,
    pmap: ProbabilityMap,
    #[builder(into)]
    start: Coord,
    #[builder(into)]
    goal: Coord,
}

impl Solver {
    pub fn at_threshold(&self, threshold: f32) -> ThresholdedSpace {
        ThresholdedSpace::from_probability_map(&self.pmap, threshold)
    }

    pub fn barcode(&self) -> Barcode {
        // Determine all the epsilons
        let mut epsilons = (0..(1.0 / self.delta) as i32 + 1)
            .map(|i| i as f32 * self.delta)
            .collect_vec();
        // Append 1.0 if it's not already in
        if let Some(x) = epsilons.last()
            && *x != 1.0
        {
            epsilons.push(1.0)
        }

        let tspaces = epsilons.iter().map(|&e| self.at_threshold(e)).collect_vec();

        // Maintain a set of open path classes
        // An open class is a vec containing a tuple mapping between epsilon and a path at that
        // epsilon threshold
        // Whenever a path class is merged into another, close it and put it into
        // This is initialized to all the optimal paths at epsilon = 0.0
        let mut open_classes: Vec<Vec<(f32, Path)>> = vec![];
        let mut barcode = Barcode::default();

        for (e, tspace) in epsilons.iter().zip(&tspaces) {
            // Compute the optimal paths at this threshold as a mapping from signature to path
            let paths: HashMap<EdgeSignature, Path> =
                HashMap::from_iter(tspace.possible_signatures().filter_map(|sig| {
                    if let Some(p) = tspace.a_star(self.start, self.goal, sig.clone()) {
                        Some((sig, p))
                    } else {
                        None
                    }
                }));

            // For each open class
            // 1. Take the latest path
            // 2. project it into this tspace
            // 3. Record the computed signature
            // If two or more path classes are merged, compute hausdorff and only keep the closest
            // - The other path classes are closed and recorded in the barcode
            // - When this is done, remove the associated path from the paths map for the current
            //   threshold
            // - At the end record all remaining signature/paths in the map as a new open class

            // Map the signature of a projection to its index in the open_classes
            let mut projections: HashMap<EdgeSignature, Vec<usize>> =
                HashMap::from_iter(tspace.possible_signatures().map(|sig| (sig, vec![])));
            open_classes.iter().enumerate().for_each(|(i, pclass)| {
                let last_path = &pclass.last().expect("path class has no recorded paths").1;
                let sig = tspace
                    .project_path(last_path)
                    .expect("could not project path");
                projections.entry(sig).or_default().push(i);
            });

            let mut new_open_classes: Vec<Vec<(f32, Path)>> = vec![];
            // Go over each signature and do the respective action
            for (sig, mergers) in projections {
                // If no path exists for that signature it is not a possible path
                // For example an obstacle is next to a boundary making it impossible to go around
                if let Some(new_path) = paths.get(&sig).cloned() {
                    match mergers.len() {
                        // New pathclass, 0 paths from the prev tspace mapped to this one
                        0 => {
                            new_open_classes.push(vec![(*e, new_path)]);
                        }

                        // 1to1 mapping, simply extend the existing pathclass
                        1 => {
                            let mut class = open_classes[mergers[0]].clone();
                            class.push((*e, new_path));
                            new_open_classes.push(class);
                        }

                        // Merger, determine hausdorff distance to see which pathclass remains
                        // Close the other candidates
                        _ => {
                            let mut best_index = 0;
                            let mut best = INFINITY;
                            for &m in &mergers {
                                let candidate_path =
                                    open_classes[m].last().expect("empty path class").1.clone();
                                let h = hausdorff(&new_path, &candidate_path);
                                if h < best {
                                    best = h;
                                    best_index = m;
                                }
                            }
                            // Found the best path, extend the path class with it and close all others
                            let mut class = open_classes[best_index].clone();
                            class.push((*e, new_path));
                            new_open_classes.push(class);
                            for &m in mergers
                                .iter()
                                .enumerate()
                                .filter_map(|(i, e)| if i == best_index { None } else { Some(e) })
                            {
                                barcode.close_path_class(open_classes[m].clone());
                            }
                        }
                    }
                }
            }

            // replace the open_classes with the new vec
            open_classes = new_open_classes;
        }
        // Close all path classes that are still open
        for class in open_classes {
            barcode.close_path_class(class);
        }

        barcode
    }
}

// calculate the hausdorff distance via bruteforce between two sets of points
fn hausdorff(p1: &Path, p2: &Path) -> f32 {
    let hd1 = p1
        .iter()
        .map(|c1| {
            p2.iter()
                .map(|c2| c1.distance(c2))
                .fold(INFINITY, |acc, val| acc.min(val))
        })
        .fold(NEG_INFINITY, |acc, val| acc.max(val));
    let hd2 = p2
        .iter()
        .map(|c2| {
            p1.iter()
                .map(|c1| c2.distance(c1))
                .fold(INFINITY, |acc, val| acc.min(val))
        })
        .fold(NEG_INFINITY, |acc, val| acc.max(val));

    hd1.max(hd2)
}
