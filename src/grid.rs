use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use itertools::Itertools;

#[derive(Debug)]
pub struct Grid<T> {
    inner: Vec<T>,
    rows: usize,
    cols: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Coord {
    pub row: usize,
    pub col: usize,
}

impl Coord {
    pub fn new(row: usize, col: usize) -> Coord {
        Coord { row, col }
    }

    pub fn diff(&self, other: &Coord) -> (i64, i64) {
        (
            self.row as i64 - other.row as i64,
            self.col as i64 - other.col as i64,
        )
    }
}

impl Display for Coord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({},{})", self.row, self.col)
    }
}

impl From<(usize, usize)> for Coord {
    fn from((row, col): (usize, usize)) -> Self {
        Coord::new(row, col)
    }
}

impl<T: Default + Clone> Grid<T> {
    pub fn new(rows: usize, cols: usize) -> Grid<T> {
        Grid {
            inner: vec![T::default(); rows * cols],
            rows,
            cols,
        }
    }
}

impl<T: Clone> Grid<T> {
    pub fn rows_iter(&self) -> impl Iterator<Item = impl Iterator<Item = T>> {
        (0..self.rows).map(|i| {
            self.inner[i * self.cols..i * self.cols + self.cols]
                .iter()
                .cloned()
        })
    }
}

impl<T> Grid<T> {
    pub fn in_bounds(&self, c: impl Into<Coord>) -> bool {
        let c = c.into();
        c.row < self.rows && c.col < self.cols
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn map<U>(&self, f: impl Fn(&T) -> U) -> Grid<U> {
        Grid {
            inner: self.inner.iter().map(f).collect(),
            rows: self.rows,
            cols: self.cols,
        }
    }

    pub fn filter(&self, f: impl Fn(&T) -> bool) -> impl Iterator<Item = Coord> {
        self.inner.iter().enumerate().filter_map(move |(i, val)| {
            if f(&val) {
                Some(Coord::new(i / self.cols, i % self.cols))
            } else {
                None
            }
        })
    }

    pub fn neighbors(&self, coord: impl Into<Coord>) -> impl Iterator<Item = Coord> {
        let coord = coord.into();
        let (r, c) = (coord.row, coord.col);
        let check = [
            Coord::new(r.saturating_add(1), c),
            Coord::new(r.saturating_add(1), c.saturating_add(1)),
            Coord::new(r, c.saturating_add(1)),
            Coord::new(r.saturating_sub(1), c),
            Coord::new(r.saturating_sub(1), c.saturating_sub(1)),
            Coord::new(r, c.saturating_sub(1)),
            Coord::new(r.saturating_add(1), c.saturating_sub(1)),
            Coord::new(r.saturating_sub(1), c.saturating_add(1)),
        ];

        check.into_iter().dedup().filter(|&c| self.in_bounds(c))
    }

    pub fn filter_neighbors(
        &self,
        coord: impl Into<Coord>,
        f: impl Fn(&T) -> bool,
    ) -> impl Iterator<Item = Coord> {
        let coord = coord.into();
        self.neighbors(coord).filter(move |&c| f(&self[c]))
    }
}

impl<T> Index<Coord> for Grid<T> {
    type Output = T;
    fn index(&self, c: Coord) -> &Self::Output {
        &self.inner[c.row * self.cols + c.col]
    }
}

impl<T> IndexMut<Coord> for Grid<T> {
    fn index_mut(&mut self, c: Coord) -> &mut Self::Output {
        &mut self.inner[c.row * self.cols + c.col]
    }
}
