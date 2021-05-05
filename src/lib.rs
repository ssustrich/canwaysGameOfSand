use log::{debug, error, info, warn};
use pixels::{Error, Pixels, SurfaceTexture};
use rand::Rng;
use winit::dpi::{LogicalPosition, LogicalSize, PhysicalSize};
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit_input_helper::WinitInputHelper;
use env_logger::Env;

const PARTICLETYPES: [&str; 3] = ["NONE", "SAND", "WATER"];

pub struct SandGrid {
        particles: Vec<Particle>,
        width: usize,
        height: usize,
        active_type: usize,
        // Should always be the same size as `cells`. When updating, we read from
        // `cells` and write to `scratch_cells`, then swap. Otherwise it's not in
        // use, and `cells` should be updated directly.
        scratch_particles: Vec<Particle>,
    }
    
    impl SandGrid {
        pub fn new_empty(width: usize, height: usize) -> Self {
            assert!(width != 0 && height != 0);
            let size = width.checked_mul(height).expect("too big");
            Self {
                particles: vec![Particle::default(); size],
                scratch_particles: vec![Particle::default(); size],
                active_type: 1,
                width,
                height,
            }
        }
        pub fn clear(&mut self){
            for x in 0..self.particles.len(){
                self.particles[x] = Particle::default();
            }
        }
        pub   fn set_brush_type(&mut self, brush_type: usize){
            self.active_type = brush_type;
        }
    
        pub   fn new_random(width: usize, height: usize) -> Self {
            let mut result = Self::new_empty(width, height);
            result.randomize();
            result
        }
        
    
        pub   fn randomize(&mut self) {
            let mut rng: randomize::PCG32 = generate_seed().into();
            for c in self.particles.iter_mut() {
                let mut rng1 = rand::thread_rng();
                //println!("Integer: {}", rng1.gen_range(0..PARTICLETYPES.len());
                *c = Particle::new(rng1.gen_range(0..PARTICLETYPES.len()), true);
            }
            // run a few simulation iterations for aesthetics (If we don't, the
            // noise is ugly)
            for _ in 0..3 {
                self.update();
            }
        }
    
        pub    fn update_water(&mut self, idx: usize){
            log::debug!("{:?}", self.particles[idx]);
            //check to see if we can move down
            let mut v: Vec<isize> = self.getEightNeighbors(idx);
            let mut ui = v[6];
            let mut ul = v[5];
            let mut ur = v[7];
            
            //we hit the bottom
            if ui == -1 {
                self.particles[idx].active = true; 
                self.particles[idx].already_updated = true;
            }
            else if ui > -1 && self.particles[ui as usize].active {
                if rand::random(){
                    ul ^= ur;
                    ur ^= ul;
                    ul ^= ur;
                }
           
                //check bl  (which may be swapped)
                if ul > -1 && !self.particles[ui as usize].active {
                    self.particles[idx].already_updated = true;
                    self.particles[idx].active = false;
                    self.particles[ul as usize].active = true;
                    self.particles[ul as usize].p_type = self.particles[idx].p_type;
                    self.particles[idx].p_type = 0;
                }
                else if ur > -1 && !self.particles[ur as usize].active {
                    self.particles[idx].already_updated = true;
                    self.particles[idx].active = false;
                    self.particles[ur as usize].active = true;
                    self.particles[ur as usize].p_type = self.particles[idx].p_type;
                    self.particles[idx].p_type = 0;
    
    
                } else{
                self.particles[idx].already_updated = true;
                self.particles[idx].active = true; 
                }
            }
            else {
                self.particles[idx].already_updated = true;
                self.particles[idx].active = false;
    
                self.particles[ui as usize].active = true;
                self.particles[ui as usize].p_type = self.particles[idx].p_type;
                self.particles[idx].p_type = 0;
            }
    
    
        }
        pub fn update_sand(&mut self) {
            let mut counter = 0;
            for idx in (0..self.particles.len()).rev() {
                //let neibs = self.count_neibs(x, y);
                //println!("Checking for alive cell at index {:?}", self.particles[idx]);
                //self.printCrazy8(v, idx);
    
                if !self.particles[idx].already_updated && self.particles[idx].active && self.particles[idx].p_type ==1  {
                    
                    //if at base level set to inactive go to next particle
                    //if there is not a particle below it then move it down one akd keep it active
                    //if there a dir particle below it can it go left or right
                    //if so then go bl or br and set to inactince
                    //else stop and set to inactive
    
                    
                    counter+=1;
                    //println!("executing loop counter {}",counter);
                    log::debug!("{:?}", self.particles[idx]);
                    //check to see if we can move down
                    let mut v: Vec<isize> = self.getEightNeighbors(idx);
                    let mut bi = v[2] as isize;
                    //we hit the bottom
                    if bi == -1 {
                        self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                    } else if  self.particles[bi as usize].p_type == 0{
                        self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                        self.particles[idx].p_type =0;
                        self.particles[bi as usize].active = true;
                        self.particles[bi as usize].already_updated = false;
                        self.particles[bi as usize].p_type =1;
                    } else if self.particles[bi as usize].p_type == 1 && self.particles[bi as usize].active == false {
                        let mut bl = v[3] as isize;
                        let mut br = v[1] as isize;
                        if rand::random() {
                                bl ^= br;
                                br ^= bl;
                                bl ^= br;
                        }
                        if  bl > -1 && self.particles[bl as usize].p_type == 0  {
                            self.particles[idx].active = false;
                            self.particles[idx].already_updated = true;
                            self.particles[idx].p_type =0;
                            self.particles[bl as usize].active = true;
                            self.particles[bl as usize].already_updated = false;
                            self.particles[bl as usize].p_type =1;
                        }
                        else if  br > -1 && self.particles[br as usize].p_type == 0{
                            self.particles[idx].active = false;
                            self.particles[idx].already_updated = true;
                            self.particles[idx].p_type =0;
                            self.particles[br as usize].active = true;
                            self.particles[br as usize].already_updated = false;
                            self.particles[br as usize].p_type =1;
                        }
                         else  {
                             self.particles[idx].active = false;
                            self.particles[idx].already_updated = true;
                        }
                    }
                    else {
                        println!("come here often?");
                        self.particles[idx].active = false;
                        self.particles[idx].already_updated = true;
                        self.particles[idx].p_type =0;
                        self.particles[bi as usize].active = true;
                        self.particles[bi as usize].already_updated = false;
                        self.particles[bi as usize].p_type =1;
                    }
                    
                }    
            }
        }
    
           pub  fn update(&mut self) {
            for idx in (0..self.particles.len()).rev() {
                   if !self.particles[idx].already_updated && self.particles[idx].active {
    
                    if self.particles[idx].p_type ==1 {
                        self.update_sand();
                    }
                   }
            }
   

        }
    
      pub  fn toggle(&mut self, x: isize, y: isize) -> bool {
            if let Some(i) = self.grid_idx(x, y) {
                let was_alive = self.particles[i].active;
                self.particles[i].set_active(!was_alive);
                !was_alive
            } else {
                false
            }
        }

       pub fn draw(&self, screen: &mut [u8]) {
            debug_assert_eq!(screen.len(), 4 * self.particles.len());
            for (c, pix) in self.particles.iter().zip(screen.chunks_exact_mut(4)) {
                let color = if c.p_type == 1 {
                    [0, 0xff, 0xff, 0xff]
                } else {
                    [0, 0, 0x00, 0xff]
                };
                pix.copy_from_slice(&color);
            }
        }






        // pub fn toggle(&mut self, x: isize, y: isize) -> bool {
        //     if let Some(i) = self.grid_idx(x, y) {
        //         let was_active = self.particles[i].active;
        //         self.particles[i].set_active(!was_active);
        //         self.particles[i].already_updated = false;
        //         self.particles[i].p_type = self.active_type;
        //         !was_active
        //     } else {
        //         false
        //     }
        // }
    
        // pub fn draw(&self, screen: &mut [u8]) {
        //     debug_assert_eq!(screen.len(), 4 * self.particles.len());
        //     for (c, pix) in self.particles.iter().zip(screen.chunks_exact_mut(4)) {
        //         let mut color = [0, 0x00, 0x00, 0x00];
    
        //         if c.p_type ==1{
        //             color = [0x96, 0x4b, 0x00, 0xff];    
        //         }
        //         else if c.p_type == 2{
        //             color = [0, 0xff, 0xff, 0xff];
        //         }
    
        //         pix.copy_from_slice(&color);
        //     }
        // }
    
      pub   fn set_line(&mut self, x0: isize, y0: isize, x1: isize, y1: isize, active: bool) {
            // probably should do sutherland-hodgeman if this were more serious.
            // instead just clamp the start pos, and draw until moving towards the
            // end pos takes us out of bounds.
            let x0 = x0.max(0).min(self.width as isize);
            let y0 = y0.max(0).min(self.height as isize);
            for (x, y) in line_drawing::Bresenham::new((x0, y0), (x1, y1)) {
                if let Some(i) = self.grid_idx(x, y) {
                    self.particles[i].set_active(active);
                    self.particles[i].already_updated = false;
                    self.particles[i].p_type = self.active_type;
                } else {
                    break;
                }
            }
        }
    
        pub  fn getXYfromInx(&self, idx: usize) -> (usize, usize) {
            let row: usize = idx % self.width;
            let column: usize = idx / self.width;
            (row, column)
        }
    
    
    
            /* Given an index in the array of X this function will return the index on the 8
            neighbors in an array of len 8 where array[0] is the cell to the immediate right of ?
            and the continue in a clowise fashion. If the cell is touching an edge of the game board the value 
            for neighbors that are off the board is -1
    
            [5]  [6]  [7]
    
            [4]   X   [0]
    
            [3]  [2]  [1]
    
            Example 1) Cell X is in the middle of a 3x3 game board
    
            idx = 4 and the board would look like this
            0  1  2
            3  4  5
            6  7  8
            The returned vector is <5,8,7,6,3,0,1,2>
    
            Example 2) Cell X is in the upper left corner of a 3x3 board
            idx = 0 and the board would look like this, remember -1 indicataes a wall or edge
            -1 -1 -1
            -1  0  1
            -1  3  4
            The returned vector is <1,4,3,-1,-1,-1,-1,-1>
            */
            pub     fn getEightNeighbors(&self, idx: usize) -> Vec<isize> {
            let coord = self.getXYfromInx(idx);
            let mut v: Vec<isize> = vec![-1; 8];
            log::info!("[getEightNeighbors] Called with indx:{}, which maps to x,y:{:?}",idx, coord);
            //println!("{}{:?}",idx, coord );
    
            v[0] = if coord.0 == self.width-1 { -1 } else { (idx+1) as isize };
            v[1] = if coord.0 == self.width-1 || coord.1 == self.height -1 { -1 } else {(idx + 1 + self.width) as isize};
            v[2] = if coord.1 == self.height -1 { -1 } else { (idx+ self.width) as isize };
            v[3] = if coord.0 == 0 || coord.1 == self.height -1 { -1 } else { (idx+ self.width -1 ) as isize };
            v[4] = if coord.0 == 0 {-1} else {(idx-1) as isize};
            v[5] = if coord.1 == 0 || coord.0 == 0 {-1} else {(idx - 1 - self.width) as isize};
            v[6] = if coord.1 == 0 {-1} else {(idx - self.width) as isize};
            v[7] = if coord.0 == self.width-1 || coord.1 == 0 {-1} else {(idx + 1 - self.width) as isize};
            log::debug!("[getEightNeighbors] {:?}",v);
            v
        }
    
        pub   fn printCrazy8(&self, x: Vec<isize>, idx: usize){
    
            println!();
            println!("{} {} {}",x[5], x[6], x[7] );
            println!("{} {} {}",x[4], idx, x[0]);
            println!("{} {} {}",x[3],x[2], x[1]);
            println!();
        }
    
        pub  fn grid_idx<I: std::convert::TryInto<usize>>(&self, x: I, y: I) -> Option<usize> {
            if let (Ok(x), Ok(y)) = (x.try_into(), y.try_into()) {
                if x < self.width && y < self.height {
                    Some(x + y * self.width)
                } else {
                    None
                }
            } else {
                None
            }
        }
    }
    /// Generate a pseudorandom seed for the game's PRNG.
fn generate_seed() -> (u64, u64) {
    use byteorder::{ByteOrder, NativeEndian};
    use getrandom::getrandom;

    let mut seed = [0_u8; 16];

    getrandom(&mut seed).expect("failed to getrandom");

    (
        NativeEndian::read_u64(&seed[0..8]),
        NativeEndian::read_u64(&seed[8..16]),
    )
}

     #[derive(Clone, Copy, Debug, Default)]
        /// The most basic element in teh game 
        pub struct Particle {
            p_type: usize,
            active: bool,
            already_updated: bool,
            velocity: f32,
        }
        
        impl Particle {
            pub    fn new(p_type: usize, active: bool) -> Self {
                Self {
                    p_type,
                    active,
                    already_updated: false,
                    velocity: 1.0,
                }
            }
        
            #[must_use]
            pub   fn next_state(mut self, active: bool) -> Self {
                self.active = active;
                self
            }
        
            pub    fn set_active(&mut self, active: bool) {
                *self = self.next_state(active);
            }
        }