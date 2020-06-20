#![feature(const_generics)]
// #![feature(fn_traits)]
// #![feature(unboxed_closures)]
/// calculates erf using series expansion shown here:
/// http://people.math.sfu.ca/~cbm/aands/page_297.htm
// fn erf_series_0()
use std::f64;
use std::ops;

mod cheb {
    ///https://www.embeddedrelated.com/showarticle/152.php
    /// x is your normal x, u is the normalized coordinate
    fn u_from_x(x: f64, a: f64, b: f64) -> f64 {
        ((2.0 * x) - a - b) / (b - a)
    }

    fn x_from_u(u: f64, a: f64, b: f64) -> f64 {
        ((b - a) / 2.0) * u + ((a + b) / 2.0)
    }

    fn chebyshev_node(i: u64, n: u64) -> f64 {
        return (((2.0 * i as f64 - 1.0) * std::f64::consts::PI) / (2.0 * n as f64)).cos();
    }

    // fn chebyshev_nodes<const NODE_COUNT: usize>() -> [f64; NODE_COUNT]{
    //     let mut nodes = [0.0; NODE_COUNT];
    //     (1..=NODE_COUNT as u64).for_each(|i| nodes[i-1] = chebyshev_node(i, NODE_COUNT as u64));
    //     nodes
    // }

    fn chebyshev_nodes(node_count:u64) -> Vec<f64>{
        (1..=node_count).map(|i| chebyshev_node(i, node_count)).collect::<Vec<f64>>()
    }

    fn craft_approximator(x_range_min: f64, x_range_max: f64, c_values :Vec<f64>) -> impl FnOnce(f64) -> f64{
        let approximator = move |x: f64| -> f64{
            // let c_values = c_values;
            // let x_range_min= x_range_min;
            // let x_range_max= x_range_max;
            let u = u_from_x(x, x_range_min, x_range_max);
            let mut t_prev = 1.0;
            let mut t = u;
            let mut y = c_values[0];
            for i in 1..c_values.len() {
                y += t * c_values[i];
                let t_next = 2.0 * u * t - t_prev;
                t_prev = t;
                t = t_next;
            }
            y
        };
        approximator
    }

    pub fn gen_chebyshev_approximator(x_range_min: f64, x_range_max: f64, node_count:u64, fx: fn(f64) -> f64) -> impl FnOnce(f64) -> f64{
        let nodes = chebyshev_nodes(node_count);
        let x_values = nodes.iter().map(|u| x_from_u(*u, x_range_min, x_range_max)).collect::<Vec<f64>>();
        let y_values = x_values.iter().map(|x| (fx)(*x)).collect::<Vec<f64>>();
        println!("{:#?}",nodes);
        println!("{:#?}",x_values);
        println!("{:#?}",y_values);
        let mut c_values = vec![0.0; node_count as usize];
        let mut t_values_0 = vec![1.0; node_count as usize];
        let mut t_values_1 = nodes.clone();
        if node_count > 0 {
            c_values[0] = y_values.iter().sum::<f64>() / node_count as f64;
            if node_count > 1 {
                c_values[1] = 2.0 * t_values_1.iter().zip(&y_values).map(|(t, y)| (*t) * (*y)).sum::<f64>() / node_count as f64;
                for i in 2..node_count as usize {
                    //tn+1 = 2*u*tn - tn-1;
                    let new_t_values = if i & 1 == 0 {
                        for j in 0..node_count as usize {
                            let t_value_prev = t_values_0[j];
                            let t_value = t_values_1[j];
                            t_values_0[j] = (2.0 * nodes[j] * t_value) - t_value_prev;
                        }
                        &t_values_0
                    }else{
                        for j in 0..node_count as usize {
                            let t_value_prev = t_values_1[j];
                            let t_value = t_values_0[j];
                            t_values_1[j] = (2.0 * nodes[j] * t_value) - t_value_prev;
                        }
                        &t_values_1
                    };
                    c_values[i] = 2.0 * new_t_values.iter().zip(&y_values).map(|(t, y)| (*t) * (*y)).sum::<f64>() / node_count as f64;
                }
            }
        }
        craft_approximator(x_range_min, x_range_max, c_values)
    }


    #[derive(Debug)]
    pub struct Chebyshev {
        x_range_min: f64,
        x_range_max: f64,
        fx: fn(f64) -> f64,
        c_values: Vec<f64>

    }


    impl Chebyshev {
        pub fn new(x_range_min: f64, x_range_max: f64, node_count:u64, fx: fn(f64) -> f64) -> Chebyshev {


            let nodes = chebyshev_nodes(node_count);
            let x_values = nodes.iter().map(|u| x_from_u(*u, x_range_min, x_range_max)).collect::<Vec<f64>>();
            let y_values = x_values.iter().map(|x| (fx)(*x)).collect::<Vec<f64>>();
            let mut c_values = vec![0.0; node_count as usize];
            let mut t_values_0 = vec![1.0; node_count as usize];
            let mut t_values_1 = nodes.clone();
            if node_count > 0 {
                c_values[0] = y_values.iter().sum::<f64>() / node_count as f64;
                if node_count > 1 {
                    c_values[1] = 2.0 * t_values_1.iter().zip(&y_values).map(|(t, y)| (*t) * (*y)).sum::<f64>() / node_count as f64;
                    for i in 2..node_count as usize {
                        //tn+1 = 2*u*tn - tn-1;
                        let new_t_values = if i & 1 == 0 {
                            for j in 0..node_count as usize {
                                let t_value_prev = t_values_0[j];
                                let t_value = t_values_1[j];
                                t_values_0[j] = (2.0 * nodes[j] * t_value) - t_value_prev;
                            }
                            &t_values_0
                        }else{
                            for j in 0..node_count as usize {
                                let t_value_prev = t_values_1[j];
                                let t_value = t_values_0[j];
                                t_values_1[j] = (2.0 * nodes[j] * t_value) - t_value_prev;
                            }
                            &t_values_1
                        };
                        c_values[i] = 2.0 * new_t_values.iter().zip(&y_values).map(|(t, y)| (*t) * (*y)).sum::<f64>() / node_count as f64;
                    }
                }
            }
            let chebyshev = Chebyshev { x_range_min, x_range_max, fx, c_values};
            return chebyshev;
        }
        pub fn x_range_min(&self) -> f64 {
            self.x_range_min
        }
        pub fn x_range_max(&self) -> f64 {
            self.x_range_max
        }
        pub fn node_count(&self) -> usize{
            self.c_values.len()
        }

        pub fn calc_u(&self, x: f64) -> f64 {
            let a = self.x_range_min;
            let b = self.x_range_max;
            ((2.0 * x) - a - b) / (b - a)
        }
        pub fn calc_x(&self, u: f64) -> f64 {
            let a = self.x_range_min;
            let b = self.x_range_max;
            ((b - a) / 2.0) * u + ((a + b) / 2.0)
        }
        pub fn approximate(&self, x: f64) -> f64{
            let u = self.calc_u(x);
            let mut t_prev = 1.0;
            let mut t = u;
            let mut y = self.c_values[0];
            for i in 1..self.node_count() as usize {
                y += t * self.c_values[i];
                let t_next = 2.0 * u * t - t_prev;
                t_prev = t;
                t = t_next;
            }
            y
        }
    }
}

fn calc_u(x: f64) -> f64 {
    return (x - 1.0) / 2.0;
}

fn calc_x(u: f64) -> f64 {
    return 2.0 * u + 1.0;
}

fn chebyshev_node(i: u64, n: u64) -> f64 {
    return (((2.0 * i as f64 - 1.0) * f64::consts::PI) / (2.0 * n as f64)).cos();
}

fn fx(x: f64) -> f64 {
    return (1.0 / 3.0) * x.powi(3) + 2.0 * x.powi(2) + x - 10.0;
}

fn erf_0(x: f32) -> f32 {
    let x = x.abs();
    let a1 = 0.278393_f32;
    let a2 = 0.230389_f32;
    let a3 = 0.000972_f32;
    let a4 = 0.078108_f32;
    let x2 = x * x;
    let x3 = x2 * x;
    let x4 = x3 * x;

    let denominator = 1_f32 + a1 * x + a2 * x2 + a3 * x3 + a4 * x4;
    let denominator = denominator.powi(4);

    return 1_f32 - (1_f32 / denominator);
}


fn main() {
    let mut nodes = [0.0; 5];
    (1..=5_u64).for_each(|i| nodes[(i - 1) as usize] = chebyshev_node(i, 5));
    println!("nodes");
    let nodes = nodes;
    (0..5_usize).for_each(|i| println!("{}", nodes[i]));

    println!("converted");
    let mut xs = [0.0; 5];
    (0..5_usize).for_each(|i| xs[i] = calc_x(nodes[i]));
    (0..5_usize).for_each(|i| println!("{}", xs[i]));

    println!("pure");
    let mut ys0 = [0.0; 5];
    (0..5_usize).for_each(|i| ys0[i] = fx(nodes[i]));
    (0..5_usize).for_each(|i| println!("{}", ys0[i]));

    println!("converted");
    let mut ys1 = [0.0; 5];
    (0..5_usize).for_each(|i| ys1[i] = fx(xs[i]));
    (0..5_usize).for_each(|i| println!("{}", ys1[i]));

    println!("{}", erf_0(3.0));
    println!("{}", erf_0(3.0));
    println!("Hello, world!");

    // let test = cheb::Chebyshev::new(-1.0,3.0,4_u64,(|x| (1.0 / 3.0) * x.powi(3) + 2.0 * x.powi(2) + x - 10.0));
    // println!("{:#?}", test);
    // println!("{}", test.approximate(0.0));
    let test = cheb::gen_chebyshev_approximator(-1.0,3.0,4_u64,(|x| (1.0 / 3.0) * x.powi(3) + 2.0 * x.powi(2) + x - 10.0));
    // println!("{:#?}", test);
    println!("{}", test(0.0));

}
